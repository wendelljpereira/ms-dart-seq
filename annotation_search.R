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
                        paste(output_name, "by_marks.tst", sep = "_"),
                        col.names = T,
                        row.names = F,
                        quote = F,
                        sep = "\t")

            # Warning about the sites without methylation and write a file with their names
            if ( length(without_anot) > 0){
                print(paste("Warning!: There are",
                            length(without_anot),
                            "marks without annotations for the sample",
                            output_name,
                            "!"))

                write.table(without_anot,
                            paste(output_name, "without_annotations.tst", sep = "_"),
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

            print(paste(output_name, "by_genes.tst", sep = "_"))

            write.table(genes_annotations,
                        paste(output_name, "by_genes.tst", sep = "_"),
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
                            paste(output_name, "without_annotations.tst", sep = "_"),
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
                        paste(output_name, "by_marks.tst", sep = "_"),
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
                            paste(output_name, "without_annotations.tst", sep = "_"),
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
                        paste(output_name, "by_genes.tst", sep = "_"),
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
                            paste(output_name, "without_annotations.tst", sep = "_"),
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
                        paste(output_name, "by_marks.tst", sep = "_"),
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
                            paste(output_name, "without_annotations.tst", sep = "_"),
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
                        paste(output_name, "by_genes.tst", sep = "_"),
                        col.names = T,
                        row.names = F,
                        quote = F,
                        sep = "\t")

            if (length(without_anot) > 0){
                print(paste("Warning!: There are",
                            length(without_anot),
                            "marks without annotations for the sample",
                            output_name, "!"))

                write.table(without_anot,
                            paste(output_name, "without_annotations.tst", sep = "_"),
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

# Reads the MDS-Methylated sites
G_int <- read.table(snakemake@input[[1]],
                    header = T,
                    na.strings = "NA",
                    colClasses = "character")

## MDS-Methylated sites in genes
in_g <- read.table(snakemake@input[[2]], colClasses = "character")[4]
in_g_noncod <- read.table(snakemake@input[[3]], colClasses = "character")[4]
marks_in_genes <- unique(c(in_g$V4, in_g_noncod$V4))

for (i in 1:length(names(G_int))){

    # Filter the sites on genes
    query <- G_int[G_int[, i] %in% marks_in_genes, i]

    # Executes the analysis
    anotation_search(query = query,
                     annotation_file = snakemake@input[[4]],
                     list_by = snakemake@params[["list_by"]],
                     output_name = paste0("annotation/sample_annot/",
                                          paste(names(G_int)[i],
                                                "group4",
                                                sep = "_")),
                     query_is_file = snakemake@params[["query_is_file"]],
                     query_is = snakemake@params[["query_is"]],
                     isoform_name = snakemake@params[["isoform_name"]])

    # Creates a file with the names of the methylated genes in each tissue
    file_genes <- read.table(paste0("annotation/sample_annot/",
                                    paste(names(G_int)[i],
                                          snakemake@params[["groups"]],
                                          sep = "_", "by_genes.tst")),
                             header = T,
                             sep = "\t",
                             quote = "\"")[1]

    colnames(file_genes) <- paste(names(G_int)[i])
    if (i == 1){
        file_genes_all <- file_genes
    }else if (i > 1){
        file_genes_all <- cbindX(file_genes_all, file_genes)
    }
}

write.table(file_genes_all,
            snakemake@output[[1]],
            sep = "\t",
            col.names = T,
            row.names = F,
            quote = F)

#################################
###### Go terms enrichment ######
#################################

# Reads the list of genes from gff3 file
genes_IDs <- read.table(snakemake@input[[5]], sep = "\t")
col_gene_name <- snakemake@params[["col_gene_name"]]

## Extracts the names from the Tags column
genes_IDs <- genes_IDs[genes_IDs[, 3] == "gene", ]
genes_names <- data.frame(do.call("rbind",
                                  strsplit(as.character(genes_IDs[, 9]),
                                           ";",
                                           fixed = TRUE)
)
)

genes_names <- data.frame(do.call("rbind",
                                  strsplit(as.character(genes_names[, col_gene_name]),
                                           "=",
                                           fixed = TRUE)
)
)

genes_names <- as.character(genes_names$X2)

# Verifies if all genes or a subset of then should be used as universe in enrichment analysis.
if (snakemake@params[["change_universe"]] == TRUE){
    universe_subset <- read.table(snakemake@params[["universe_subset"]], sep = "\t")
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

# Makes a list to use in "compareCluster" function. Each vector in the list is a sample.
compare_list <- as.list(file_genes_all)

# Remove NAs
compare_list <- lapply(compare_list, function(x) x[!is.na(x)])

ont <- c("BP", "MF", "CC")
all_rich_ont_terms <- data.frame()

for (o in ont){
    # Executes the enrichment analysis and compare the most significant terms for each sample.
    comparison <- compareCluster(compare_list,
                                 fun = "enrichGO",
                                 OrgDb = "org.Egrandis.eg.db",
                                 keyType = "GID",
                                 ont = paste(o),
                                 pAdjustMethod = padjustmethod,
                                 universe = universe,
                                 pvalueCutoff = pvaluecutoff,
                                 qvalue = qvalue)

    # Removes the redundance of terms
    comparison_simp <- simplify(x = comparison,
                                cutoff = sim_cutoff,
                                by = "p.adjust",
                                select_fun = min)

    terms_plot <- dotplot(comparison_simp,
                          showCategory = number_of_terms_by_plot,
                          font.size = 14,
                          by = plot_by)

    svg(filename = paste("images/cluster_profiler/GO_terms_", o, ".svg", sep = ""),
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
                paste0("annotation/sample_annot/", paste(i, "enriched_GO_terms.tst", sep = "_")),
                row.names = F,
                col.names = T,
                sep = "\t",
                quote = F)

    print(paste0("annotation/sample_annot/", paste(i, "enriched_GO_terms.tst", sep = "_")))
}
