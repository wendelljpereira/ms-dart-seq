require(AnnotationForge)
require(tidyr)
require(dplyr)
require(foreach)
require(doMC)
require(gdata)
registerDoMC(7)  # Allow to use seven processors in parallel.

###############################################################
###### Combines the annotations from Blast2GO and BioMart ######
###############################################################

# Blast2GO

##Reads the Blast2GO file with the annotation of all genes
blast2GO_complete <- read.table(snakemake@input[[1]],
                                sep = "\t",
                                header = T,
                                quote = "\"")

## Changes the Enzyme codes (E.C) for the same pattern used by BioMart
blast2GO_complete$Enzyme.Code <- gsub("EC:", "", blast2GO_complete$Enzyme.Code)

##Filters the annotations columns which are not informative.
blast2GO <- blast2GO_complete[, -c(3:15, 17, 18)]

# BioMart

## Reads the BioMart file with the annotation of all genes
biomart_complete <- read.table(snakemake@input[[2]],
                               sep = "\t",
                               header = T,
                               quote = "\"",
                               fill = T)

## Filters the annotations columns which are not informative.
biomart_data <- biomart_complete[, -c(9:13, 16, 17, 29:33)]

## Remove Nas
biomart_data <- sapply(biomart_data, as.character)
biomart_data[is.na(biomart_data)] <- ""
biomart_data <- as.data.frame(biomart_data)

# Selects the compatible annotation between the both sources, BioMart and Blast2GO and merge the datas in an unique data.frame based in the genes names.

blast2GO <- blast2GO[, c("Sequence.Name",
                         "Sequence.Description",
                         "Annotation.GO.ID",
                         "Enzyme.Code",
                         "Enzyme.Name")]

colnames(blast2GO) <- c("gene_name",
                        "b2g_gene_description",
                        "b2g_go_id",
                        "b2g_kegg_enzyme_id",
                        "b2g_kegg_enzyme_desc")

biomart_data <- biomart_data[, c("gene_name1",
                                 "gene_description",
                                 "go_id",
                                 "kegg_enzyme_id",
                                 "kegg_enzyme_desc")]

colnames(biomart_data) <- c("gene_name",
                            "bioM_gene_description",
                            "bioM_go_id",
                            "bioM_kegg_enzyme_id",
                            "bioM_kegg_enzyme_desc")

all_annot <- merge(biomart_data, blast2GO, by = "gene_name")

## Compares the two annotations. The final results are the union between the annotations. Additionaly, the terms which was exclusive for each source are insert in new columns, for manual inspection.

# vector with the columns for comparison.
feature <- c("gene_description",
             "go_id",
             "kegg_enzyme_id",
             "kegg_enzyme_desc")

annot_final <- data.frame(gene_name = all_annot[, "gene_name"])

count <- as.numeric(0)
for (f in feature){
  combine_data_all <- data.frame()

  combine_data_all_2 <- foreach(i = 1:nrow(all_annot), .combine = "rbind") %dopar% {
    # Extract the columns to comparison using the feature vector. In each loop, one annotation is compared.
    teste_set <- all_annot[, c(1, grep(f, colnames(all_annot)))]

    # Determines the terms for each source, Blast2GO and BioMart.
    gene_name <- paste(teste_set[i, 1])
    bioM_set_names <- strsplit(as.character(teste_set[i, grep("bioM", colnames(teste_set))]),
                               ";",
                               fixed = TRUE)[[1]]

    b2g_set_names <- strsplit(as.character(teste_set[i, grep("b2g", colnames(teste_set))]),
                              ";",
                              fixed = TRUE)[[1]]

    # Compares the terms and build the data.frame.
    comb_merge <- paste(unique(c(bioM_set_names, b2g_set_names)),
                        collapse = " & ")
    bioM_uniq <- paste(unique(bioM_set_names[!bioM_set_names %in% b2g_set_names]),
                       collapse = " & ")
    b2g_uniq <- paste(unique(b2g_set_names[!b2g_set_names %in% bioM_set_names]),
                      collapse = " & ")
    common_terms <- paste(unique(b2g_set_names[b2g_set_names %in% bioM_set_names]),
                          collapse = " & ")

    combine_data <- data.frame(gene_name,
                               comb_merge,
                               bioM_uniq,
                               b2g_uniq,
                               common_terms)

    colnames(combine_data) <- c("gene_name",
                                paste(f, "merged", sep = "_"),
                                paste(f, "unique_in_BioMart", sep = "_"),
                                paste(f, "unique_in_blastGO", sep = "_"),
                                paste(f, "common_to_both", sep = "_"))

    combine_data_all <- rbind(combine_data_all, combine_data)

  }
  annot_final <- merge(annot_final, combine_data_all_2, by = "gene_name")
}

save(file = snakemake@output[[1]], annot_final)

annot_final <- data.frame(lapply(annot_final, function(x) {gsub("NA &", "", x)} ))

# Inserts the annotation unique generated in the blast2GO.
blast2GO_unique <- blast2GO_complete[, c("Sequence.Name",
                                         "InterPro.Accession",
                                         "InterPro.Type",
                                         "InterPro.Name",
                                         "InterPro.Signatures",
                                         "InterPro.GO.ID",
                                         "InterPro.GO.Term",
                                         "InterPro.GO.Category")]

colnames(blast2GO_unique) <- c("gene_name",
                               "InterPro.Accession",
                               "InterPro.Type",
                               "InterPro.Name",
                               "InterPro.Signatures",
                               "InterPro.GO.ID",
                               "InterPro.GO.Term",
                               "InterPro.GO.Category")

annot_final <- merge(annot_final, blast2GO_unique, by = "gene_name")

# reads the file with the distances between the tested sites and the closest gene
distance <- read.table(snakemake@input[[3]], sep = "\t")

# keeps only the closest element of each methylation site
distance_filtered <- data.frame()
for( i in unique(distance[, 4]) ){

  dist_sub <- distance[distance[, 4] == i, ]
  min_dist <- as.numeric(min(abs(dist_sub[, 16])))

  dist_sub <- dist_sub[abs(dist_sub[, 16]) == min_dist, ]

  distance_filtered <- rbind(distance_filtered, dist_sub)
}

distance <- distance_filtered[, c(4, 2, 3, 6, 16, 15, 7, 10, 11, 13)]

# Select the columns of interest and split the name of genes.
genes_names <- data.frame(do.call("rbind",
                                  strsplit(as.character(distance[, 6]), "=", fixed = TRUE)))

genes_names <- as.character(genes_names$X3)
distance[, 6] <- genes_names

colnames(distance) <- c("locus",
                        "locus_start",
                        "locus_end",
                        "locus_strand",
                        "distance_to_gene(kb)",
                        "gene_name",
                        "chromosome",
                        "gene_start",
                        "gene_end",
                        "gene_strand" )

# For each loci, merge the distance to the closest gene with the annotation.
annot_table <- merge(distance, annot_final, by = "gene_name")
annot_table <- annot_table[, c(2:6, 1, 7:33)]

# Writes the table with the combined annotation.
write.table(annot_table,
            snakemake@output[[2]],
            col.names = T,
            row.names = F,
            sep = "\t",
            quote = F)
