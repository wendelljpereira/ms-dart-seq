require(biomaRt)
require(dplyr)
require(clusterProfiler)

############################################
# Annotation using BioMart - Configuration #
############################################

annot_file <- snakemake@params[["annotation_file_rda"]]

if (file.exists(annot_file) == TRUE){
    load(annot_file)
}else{
    # parameters to set BioMart
    its_in_list_marts <- snakemake@params[["its_in_listMarts"]] # yes or no
    biomart_name <- snakemake@params[["biomart_name"]]
    biomart_dataset <- snakemake@params[["biomart_dataset"]]
    
    # if Mart is not in listMarts, also set:
    biomart_host <- snakemake@params[["biomart_host"]]
    biomart_vschema <- snakemake@params[["biomart_vschema"]]
    
    # Set what division in gff3 file column 9 containing the gene names
    c_name <- as.numeric(snakemake@params[["c_name"]])
    
    # Set the Mart to use.
    if (its_in_list_marts == "yes"){
        
        myusemart <- useMart(biomart_name, dataset = biomart_dataset)
        
    }else if (its_in_list_marts == "no"){
        
        ## For use BioMarts from Phytozome
        Mart <- new("Mart",
                    biomart = biomart_name,
                    vschema = biomart_vschema, 
                    host = biomart_host)
        
        ## Select the dataset to use
        mysets <- listDatasets(Mart)
        mydataset <- mysets$dataset[mysets$dataset == biomart_dataset] # Select the dataset "Phytozome 12 Genomes"
        myusemart <- useDataset(as.character(mydataset), mart = Mart)
    }
    
    ######################################
    # Annotation of all E. grandis genes #
    ######################################
    
    # List all atributes available in selected Mart
    attributes <- listAttributes(mart = myusemart)
    
    # Select the attributes useful to annotation
    my_attributes <- attributes[c(1:35), 1] 
    
    # List all filters available in selected Mart
    filters <- listFilters(myusemart)
    
    # Select the information which will be use as input (filter)
    ## As our list containing the Gene id, we selected the "gene_name_filter", filter 11 in the list of filters
    my_filter <- filters[11, 1]
    
    # Reaf the list of genes from gff3 file
    genes_IDs <- read.table(snakemake@input[[1]], sep = "\t")
    
    ## Extract the names from the Tags column
    genes_IDs <- genes_IDs[genes_IDs[, 3] == "gene", ]
    genes_names <- data.frame(do.call("rbind",
                                      strsplit(as.character(genes_IDs[,9]),
                                               ";",
                                               fixed = TRUE)))
    
    genes_names <- data.frame(do.call("rbind", strsplit(as.character(genes_names[, c_name]), "=", fixed = TRUE)))
    genes_names <- as.character(genes_names$X2)
    
    # Execute the annotation for the genes in "genes_list" vector
    biomart_annot <- getBM(attributes = my_attributes,
                           filters = my_filter,
                           mart = myusemart,
                           values = genes_names)
    
    save(file = annot_file, biomart_annot)
}

# Summarise to one row per gene. Multiples annotations will be separate by comma
annotation_file_res <- biomart_annot %>%
    group_by(gene_name1) %>%
    summarise_all(funs(paste(unique(.), collapse = ", ")))

# Write annotation table
write.table(annotation_file_res,
            snakemake@output[[1]],
            col.names = T,
            row.names = F,
            quote = F,
            sep = "\t")