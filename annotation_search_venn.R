require(gdata)
require(dplyr)
require(clusterProfiler)

# Loads the OrgDb package containing the GO terms.
pkg <- snakemake@params[["OrgDb_package_name"]]
require(pkg, character.only = TRUE)

########################
###### Annotation ######
########################

# Loads the function that will make the annotation
anotation_search <- function(query = query,
                             annotation_file = annotation_file,
                             list_by = list_by,
                             output_name = output_name,
                             query_is_file = "FALSE",
                             query_is = query_is,
                             isoform_name = "TRUE"){

    # Reads the MSD-Methylated sites. Check if it is a file or a vector with the MSD Sites
    if (query_is_file == "TRUE"){
        lista <- read.table(query, colClasses = "character")
        lista <- as.character(lista$V1)
    }else if (query_is_file == "FALSE"){
        lista <- as.character(query[!is.na(query)])
    }else{
        print("Invalid 'query_is_file' argument!")
    }

    # Reads the annotation table of the genome
    anot <- read.table(annotation_file,
                       header = T,
                       sep = "\t",
                       quote = "\"",
                       colClasses = "character",
                       fill = T)

    colnames(anot) <- c("locus",
                        "locus_start",
                        "locus_end",
                        "locus_strand",
                        "distance_to_gene_kb",
                        "gene_name",
                        "chromosome",
                        "gene_start",
                        "gene_end",
                        "gene_strand",
                        "gene_description",
                        "gene_desc_unique_BioMart",
                        "gene_desc_unique_B2GO",
                        "gene_descriptiom_intersect",
                        "annotation_go_id",
                        "go_id_unique_BioMart",
                        "go_id_unique_B2GO",
                        "go_id_intersect",
                        "enzyme_code",
                        "enzyme_code_unique_BioMart",
                        "enzyme_code_unique_B2GO",
                        "enzyme_code_intersect",
                        "enzyme_description",
                        "enzyme_description_unique_BioMart",
                        "enzyme_description_unique_B2GO",
                        "enzyme_description_intersect",
                        "interpro_accession",
                        "interpro_type",
                        "interpro_name",
                        "interpro_signatures",
                        "interpro_go_id",
                        "interpro_go_term",
                        "interpro_go_category")

    # Verify how the output should be generated
    if (query_is == "mark"){
        if (list_by == "mark"){

            anotation <- anot[anot$locus %in% lista, ]
            without_anot <- lista[lista %in% anot$locus == "FALSE"]

            write.table(anotation,
                        paste(output_name, "by_marks.txt", sep = "_"),
                        col.names = T,
                        row.names = F,
                        quote = F,
                        sep = "\t")

            # Warning about the sites without methylation and write a file with their names
            if (length(without_anot) > 0){
                print(paste("Warning!: There are",
                            length(without_anot),
                            "marks without annotations for the sample",
                            output_name,
                            "!"))

                write.table(without_anot,
                            paste(output_name, "without_annotations.txt", sep = "_"),
                            col.names = F,
                            row.names = F,
                            quote = F,
                            sep = "\t")
            }
        }
        else if (list_by == "gene"){

            anotation <- anot[anot$locus %in% lista, ]
            without_anot <- lista[lista %in% anot$locus == "FALSE"]

            # Join the MSD-Methylated sites of the same genes.
            gene_anot_complete <- data.frame()
            for (i in unique(anotation$gene_name)){
                gene <- anotation[anotation$gene_name == i, ]
                locus_gene <- paste(gene$locus, collapse = ",")
                locus_strand <- paste(gene$locus_strand, collapse = ",")
                locus_distance <- paste(gene$distance_to_gene_kb, collapse = ",")
                gene_anot <- gene[1, ]
                gene_anot$locus <- locus_gene
                gene_anot$locus_strand <- locus_strand
                gene_anot$distance_to_gene_kb <- locus_distance
                gene_anot_complete <- rbind(gene_anot_complete, gene_anot)
            }

            # Formats the file and remove unecessary columns
            genes_names <- as.data.frame(gene_anot_complete$gene_name)
            colnames(genes_names) <- "gene_name"
            genes_annotations <- cbind(genes_names, gene_anot_complete)
            genes_annotations <- genes_annotations[, -c(13:15, 17:19, 21:23, 25:27)]

            if (isoform_name == "TRUE"){
                isoform <- as.data.frame(paste(genes_annotations$gene_name, ".1", sep = ""))
                colnames(isoform) <- "isoform_name"
                genes_annotations <- cbind(isoform, genes_annotations)
            }

            print(paste(output_name, "by_genes.txt", sep = "_"))

            write.table(genes_annotations,
                        paste(output_name, "by_genes.txt", sep = "_"),
                        col.names = T,
                        row.names = F,
                        quote = F,
                        sep = "\t")

            if (length(without_anot) > 0){
                print(paste("Warning!: There are",
                            length(without_anot),
                            "marks without annotations for the sample",
                            output_name,
                            "!"))

                write.table(without_anot,
                            paste(output_name, "without_annotations.txt", sep = "_"),
                            col.names = F,
                            row.names = F,
                            quote = F,
                            sep = "\t")
            }
        }else{
            print("Invalid 'list_by' argument!")
            print("'list_by' accepts only 'gene' or 'mark'")
        }
    }else if (query_is == "gene"){
        if (list_by == "mark"){

            anotation <- anot[anot$gene_name %in% lista, ]
            without_anot <- lista[lista %in% anot$gene_name == "FALSE"]

            write.table(anotation,
                        paste(output_name, "by_marks.txt", sep = "_"),
                        col.names = T,
                        row.names = F,
                        quote = F,
                        sep = "\t")

            if (length(without_anot) > 0){
                print(paste("Warning!: There are",
                            length(without_anot),
                            "marks without annotations for the sample",
                            output_name,
                            "!"))

                write.table(without_anot,
                            paste(output_name, "without_annotations.txt", sep = "_"),
                            col.names = F,
                            row.names = F,
                            quote = F,
                            sep = "\t")
            }
        }
        else if (list_by == "gene"){

            anotation <- anot[anot$gene_name %in% lista, ]
            without_anot <- lista[lista %in% anot$gene_name == "FALSE"]

            gene_anot_complete <- data.frame()
            for (i in unique(anotation$gene_name)){
                gene <- anotation[anotation$gene_name == i, ]
                locus_gene <- paste(gene$locus, collapse = ",")
                locus_strand <- paste(gene$locus_strand, collapse = ",")
                locus_distance <- paste(gene$distance_to_gene_kb, collapse = ",")
                gene_anot <- gene[1, ]
                gene_anot$locus <- locus_gene
                gene_anot$locus_strand <- locus_strand
                gene_anot$distance_to_gene_kb <- locus_distance
                gene_anot_complete <- rbind(gene_anot_complete, gene_anot)
            }

            genes_names <- as.data.frame(gene_anot_complete$gene_name)
            colnames(genes_names) <- "gene_name"
            genes_annotations <- cbind(genes_names, gene_anot_complete)
            genes_annotations <- genes_annotations[, -c(13:15, 17:19, 21:23, 25:27)]

            if (isoform_name == "TRUE"){
                isoform <- as.data.frame(paste(genes_annotations$gene_name, ".1", sep = ""))
                colnames(isoform) <- "isoform_name"
                genes_annotations <- cbind(isoform, genes_annotations)
            }

            write.table(genes_annotations,
                        paste(output_name, "by_genes.txt", sep = "_"),
                        col.names = T,
                        row.names = F,
                        quote = F,
                        sep = "\t")

            if (length(without_anot) > 0){
                print(paste("Warning!: There are",
                            length(without_anot),
                            "marks without annotations for the sample",
                            output_name,
                            "!"))

                write.table(without_anot,
                            paste(output_name, "without_annotations.txt", sep = "_"),
                            col.names = F,
                            row.names = F,
                            quote = F,
                            sep = "\t")
            }
        }
        else{
            print("Invalid 'list_by' argument!")
            print("'list_by' accepts only 'gene' or 'mark'")
        }
    }else if (query_is == "transcript"){

        lista2 <- character()
        for (i in lista){
            new_id <- unlist(strsplit(i, "[.]"))
            new_id <- paste(new_id[1], new_id[2], sep = ".")
            lista2 <- rbind(lista2, new_id)
        }
        lista <- as.character(lista2)

        if (list_by == "mark"){

            anotation <- anot[anot$gene_name %in% lista, ]
            without_anot <- lista[lista %in% anot$gene_name == "FALSE"]

            write.table(anotation,
                        paste(output_name, "by_marks.txt", sep = "_"),
                        col.names = T,
                        row.names = F,
                        quote = F,
                        sep = "\t")

            if (length(without_anot) > 0){
                print(paste("Warning!: There are",
                            length(without_anot),
                            "marks without annotations for the sample",
                            output_name,
                            "!"))

                write.table(without_anot,
                            paste(output_name, "without_annotations.txt", sep = "_"),
                            col.names = F,
                            row.names = F,
                            quote = F,
                            sep = "\t")
            }
        }
        else if (list_by == "gene"){

            anotation <- anot[anot$gene_name %in% lista, ]
            without_anot <- lista[lista %in% anot$gene_name == "FALSE"]

            gene_anot_complete <- data.frame()
            for (i in unique(anotation$gene_name)){
                gene <- anotation[anotation$gene_name == i, ]
                locus_gene <- paste(gene$locus, collapse = ",")
                locus_strand <- paste(gene$locus_strand, collapse = ",")
                locus_distance <- paste(gene$distance_to_gene_kb, collapse = ",")
                gene_anot <- gene[1, ]
                gene_anot$locus <- locus_gene
                gene_anot$locus_strand <- locus_strand
                gene_anot$distance_to_gene_kb <- locus_distance
                gene_anot_complete <- rbind(gene_anot_complete, gene_anot)
            }

            genes_names <- as.data.frame(gene_anot_complete$gene_name)
            colnames(genes_names) <- "gene_name"
            genes_annotations <- cbind(genes_names, gene_anot_complete)
            genes_annotations <- genes_annotations[, -c(13:15, 17:19, 21:23, 25:27)]

            if (isoform_name == "TRUE"){
                isoform <- as.data.frame(paste(genes_annotations$gene_name, ".1", sep = ""))
                colnames(isoform) <- "isoform_name"
                genes_annotations <- cbind(isoform, genes_annotations)
            }

            write.table(genes_annotations,
                        paste(output_name, "by_genes.txt", sep = "_"),
                        col.names = T,
                        row.names = F,
                        quote = F,
                        sep = "\t")

            if (length(without_anot) > 0){
                print(paste("Warning!: There are",
                            length(without_anot),
                            "marks without annotations for the sample",
                            output_name,
                            "!"))

                write.table(without_anot,
                            paste(output_name, "without_annotations.txt", sep = "_"),
                            col.names = F,
                            row.names = F,
                            quote = F,
                            sep = "\t")
            }
        }
        else{
            print("Invalid 'list_by' argument!")
            print("'list_by' accepts only 'gene' or 'mark'")
        }
    }
}

for (i in 2:(length(snakemake@input) - 1)){

    anotation_search(query = paste(snakemake@input[[i]]),
                     annotation_file = snakemake@input[[1]],
                     list_by = "gene",
                     output_name = paste(unlist(strsplit(
                         snakemake@input[[i]],
                         split = "_genes",
                         fixed = TRUE))[1],
                         "annotation",
                         sep = "_"),
                     query_is_file = "TRUE",
                     query_is = "transcript",
                     isoform_name = "TRUE")

    print(
        paste(
            unlist(
                strsplit(snakemake@input[[i]],
                         split = "_genes",
                         fixed = TRUE))[1],
            "annotation",
            sep = "_"))
}

#################################
###### Go terms enrichment ######
#################################

# Reads the list of genes from gff3 file
genes_IDs <- read.table(snakemake@input[[length(snakemake@input)]],
                        sep = "\t")

col_gene_name <- snakemake@params[["col_gene_name"]]

## Extracts the names from the Tags column
genes_IDs <- genes_IDs[genes_IDs[, 3] == "gene", ]
genes_names <- data.frame(
    do.call("rbind",
            strsplit(as.character(genes_IDs[, 9]), ";", fixed = TRUE)))

genes_names <- data.frame(
    do.call("rbind",
            strsplit(as.character(genes_names[, col_gene_name]), "=", fixed = TRUE)))
genes_names <- as.character(genes_names$X2)

# Verifies if all genes or a subset of then should be used as universe in enrichment analysis.
if (snakemake@params[["change_universe"]] == TRUE){

    universe_subset <- read.table(snakemake@params[["universe_subset"]],
                                  sep = "\t")
    universe <- genes_names[genes_names %in% universe_subset$V1]

}else if (snakemake@params[["change_universe"]] == FALSE){
    universe <- genes_names
}

# Enrichment parameters
padjustmethod <- snakemake@params[["pAdjustMethod"]]
pvaluecutoff <- snakemake@params[["pvalueCutoff"]]
qvalue <- snakemake@params[["qvalue"]]
sim_cutoff <- snakemake@params[["sim_cutoff"]] # cutoff to removes similarity

# plot parameters
number_of_terms_by_plot <- snakemake@params[["number_of_terms_by_plot"]]
plot_by <- snakemake@params[["plot_by"]]

# Makes a list to use in "compareCluster" function. Each vector in the lista is a sample.
file_genes_all <- data.frame()
for (i in 2:(length(snakemake@input) - 1)){
    genes_list <- read.table(snakemake@input[[i]], sep = "\t")
    colnames(genes_list) <- strsplit(paste(snakemake@input[[i]]), ".txt")
    if (i == 2){
        file_genes_all <- genes_list
    }else{
        file_genes_all <- cbindX(file_genes_all, genes_list)
    }
}

compare_list <- as.list(file_genes_all)

# Remove NAs
compare_list <- lapply(compare_list, function(x) x[!is.na(x)])

#
ont <- c("BP", "MF", "CC")
all_rich_ont_terms <- data.frame()
for (o in ont){
    # Executes the enrichment analysis and compare the most significant terms for each sample.
    comparison <- compareCluster(compare_list,
                                 fun = "enrichGO",
                                 OrgDb = "org.Egrandis.eg.db",
                                 keyType = "GID", # Nome da coluna com os genes no OrgDb
                                 ont = paste(o),
                                 pAdjustMethod = padjustmethod,
                                 universe = universe,
                                 pvalueCutoff = pvaluecutoff,
                                 qvalue = qvalue)

    # Removes the redundance of terms
    comparison_simp <- simplify(comparison,
                                cutoff = sim_cutoff,
                                by = "p.adjust",
                                select_fun = min)

    terms_plot <- dotplot(comparison_simp,
                          showCategory = number_of_terms_by_plot,
                          font.size = 14,
                          by = plot_by)

    svg(filename = paste0("images/cluster_profiler",
                          paste("GO_terms_venn_sub_set",
                                o,
                                ".svg",
                                sep = "")),
        width = 16,
        height = 10,
        pointsize = 12)

    print(terms_plot)
    dev.off()

    rich_ont_terms <- as.data.frame(comparison_simp)
    rich_ont_terms$Ontology <- rep(o, nrow(rich_ont_terms))
    rich_ont_terms <- rich_ont_terms[, c(1:3, 11, 4:10)]
    all_rich_ont_terms <- rbind(all_rich_ont_terms, rich_ont_terms)
}

for (i in unique(all_rich_ont_terms[, 1])){
    enriched_terms <- all_rich_ont_terms[all_rich_ont_terms[, 1] == i, ]
    write.table(enriched_terms,
                paste(i, "enriched_GO_terms.txt", sep = "_"),
                row.names = F,
                col.names = T,
                sep = "\t",
                quote = F)

    print(paste(i, "enriched_GO_terms.txt", sep = "_"))
}
