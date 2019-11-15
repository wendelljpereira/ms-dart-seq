require(ggplot2)
require(DESeq2)
require(plyr)
require(gdata)
require(tidyr)

# Help function to do the differential expression analysis
deseq2_with_counts <- function(file = file,
                               prefix = prefix,
                               clone_name = clone_name,
                               tissue = tissue,
                               fdr = 0.001,
                               log_fc = 1,
                               sep_into = sep_into,
                               subset_model = "normal",
                               no_bio_rep = "FALSE",
                               dispersion = NULL,
                               filter = "posit",
                               min_msp = min_msp,
                               intersect_file = intersect_file){

    data <- read.table(as.character(file),
                       header = T,
                       sep = ",",
                       check.names = F)

    marcas <- as.character(data[, 1])

    prefixo <- as.character(prefix)
    clone_name <- as.character(clone_name)
    tecido <- as.character(tissue)
    fdr <- as.numeric(fdr)
    log_fc <- as.numeric(log_fc)
    sep_into <- as.numeric(sep_into)
    subset_model <- as.character(subset_model)
    filter <- as.character(filter)
    min_msp <- as.numeric(min_msp)

    # Selects the name of the sample and tissue that should be used
    clone <- data[, grep(pattern = clone_name, colnames(data))]
    clone <- clone[, grep(pattern = tecido, colnames(clone))]

    clone <- as.data.frame(clone, row.names = marcas)

    # Selects the model of the analysis
    if (subset_model == "nonzero"){

        clone <- clone[rowSums(clone) > 0, ]
        print("subset_model = nonzero")

    } else if (subset_model == "intersect"){

        intersect <- read.table(intersect_file,
                                header = T,
                                quote = "\"",
                                sep = "\t",
                                check.names = F,
                                colClasses = "character")[1]

        clone <- clone[rownames(clone) %in% intersect$BRASUZ1, ]

        print("subset_model = intersect")

    } else {

        print("subset_model = normal")

    }

    # Writes a file with the names of the sites that were texted for methylation
    tested_sites <- rownames(clone)

    write.table(tested_sites,
                "tested_sites_to_methylation.txt",
                col.names = F,
                row.names = F,
                quote = F,
                sep = "\t")

    # Filters by the minimum counts to be considerated as a true site
    if (filter == "posit"){
        if (no_bio_rep == "FALSE"){
            filtred_clone <- data.frame()
            for (i in 1:nrow(clone)){
                count_average <- clone[i, ]
                count_average <- count_average[, grep(pattern = "ms", colnames(clone))]
                if ( (sum(count_average[1, ]) / length(count_average)) >= min_msp){
                    filtred_clone <- rbind(filtred_clone, clone[i, ])
                }
            }
        }else if (no_bio_rep == "TRUE"){
            filtred_clone <- data.frame()
            for (i in 1:nrow(clone)){
                count_average <- clone[i, ]
                count_average <- as.data.frame(count_average[, grep(pattern = "ms", colnames(clone))])
                if ( (sum(count_average[1, ]) / length(count_average)) >= min_msp){
                    filtred_clone <- rbind(filtred_clone, clone[i, ])
                }
            }
        }
    }else if (filter == "negat"){
        if (no_bio_rep == "FALSE"){
            filtred_clone <- data.frame()
            for (i in 1:nrow(clone)){
                count_average <- clone[i, ]
                count_average <- count_average[, grep(pattern = "hp", colnames(clone))]
                if ( (sum(count_average[1, ]) / length(count_average)) >= min_msp){
                    filtred_clone <- rbind(filtred_clone, clone[i, ])
                }
            }
        }else if (no_bio_rep == "TRUE"){
            filtred_clone <- data.frame()
            for (i in 1:nrow(clone)){
                count_average <- clone[i, ]
                count_average <- as.data.frame(count_average[, grep(pattern = "hp", colnames(clone))])
                if ( (sum(count_average[1, ]) / length(count_average)) >= min_msp){
                    filtred_clone <- rbind(filtred_clone, clone[i, ])
                }
            }
        }
    }

    print(paste("For sample",
                paste("'", clone_name, "_", tecido, "'", sep = ""),
                nrow(filtred_clone),
                "marks passed for all filters!"))

    msp_10_sites <- as.data.frame(rownames(filtred_clone))

    write.table(msp_10_sites,
                paste(clone_name, "_", tecido, "_", "msp_bigger_than_10.txt", sep = ""),
                row.names = F,
                col.names = F,
                quote = F)

    clone <- lapply(filtred_clone, as.integer)
    clone <- as.data.frame(clone, check.names = FALSE)

    clone <- as.data.frame(clone, row.names = rownames(filtred_clone))

    # Help function to generate some descriptive stats
    descritivas <- function(y) {
        data.frame(n = length(y),
                   mean = mean(y),
                   IC95inf = t.test(y)$conf.int[1],
                   IC95sup = t.test(y)$conf.int[2],
                   var = var(y),
                   std.dev. = sd(y),
                   std.err. = sd(y) / sqrt(length(y)),
                   CV = 100 * sd(y) / mean(y),
                   min = min(y),
                   max = max(y),
                   Q1 = quantile(y, 0.25),
                   median = median(y),
                   Q3 = quantile(y, 0.75),
                   sum(y))
    }

    # Determines the groups of samples using the first section of the name
    names <- as.data.frame(names(clone))
    colnames(names) <- "nome"

    sep <- paste0("x", 1:(sep_into + 1))

    names_group <- separate(names,
                            nome,
                            into = c(paste0("null", 1:(sep_into)), "Enzime"),
                            sep = "_",
                            remove = T,
                            extra = "merge")

    names_group <- as.character(names_group$Enzime)

    x <- as.character()
    for (i in 1:length(names_group)){
        x <- rbind(x, substr(names_group[i], 0, 2))
    }
    groups <- as.character(x)
    groups <- gsub(".", "", groups, fixed = T)

    # Load the processed data to DEseq2 object
    coldata <- data.frame(condition = groups,
                          type = c("single-read",
                                   "single-read",
                                   "single-read",
                                   "single-read",
                                   "single-read",
                                   "single-read"),
                          row.names = paste(colnames(clone)))

    data_deseq <- DESeqDataSetFromMatrix(countData = clone,
                                         colData = coldata,
                                         design = ~ condition)

    # Determine the reference group
    data_deseq$condition <- relevel(data_deseq$condition, ref = "hp")

    # Executes the analysis
    data_deseq <- DESeq(data_deseq)

    # Extract the results
    raw_results <- as.data.frame(results(data_deseq))

    # Filter by FoldChange and FDR
    if (filter == "posit"){

        results_sig <- subset(raw_results, raw_results$padj <= fdr)
        results_sig <- subset(results_sig, results_sig$log2FoldChange >= log_fc)

    }else if (filter == "negat"){

        results_sig <- subset(raw_results, raw_results$padj <= fdr)
        results_sig <- subset(results_sig, results_sig$log2FoldChange <= 0)
        results_sig <- subset(results_sig, abs(results_sig$log2FoldChange) >= log_fc)

    }else if (filter == "full"){

        results_sig <- subset(raw_results, raw_results$padj <= fdr)
        results_sig <- subset(results_sig, abs(results_sig$log2FoldChange) >= log_fc)

    }else{

        print("error: filter value is incorrectly defined.")

    }

    if (nrow(results_sig) == 0){

        print(paste("There is no differentialy expressed sites for",
                    prefix,
                    clone_name,
                    "tissue:",
                    tissue,
                    "!",
                    sep = " "))

    }else{

        marcks_sig <- rownames(results_sig)
        Data_DE <- as.data.frame(marcks_sig)
        colnames(Data_DE) <- paste(clone_name, tecido, sep = "_")

        if (file.exists(paste(prefixo, "DE_marks.txt", sep = "_")) == "FALSE"){

            write.table(Data_DE,
                        file = paste(prefixo, "DE_marks.txt", sep = "_"),
                        sep = "\t",
                        quote = F,
                        row.names = F,
                        col.names = T)

        }else if (file.exists(paste(prefixo, "DE_marks.txt", sep = "_")) == "TRUE"){

            Data_DE_g <- read.table(paste(prefixo, "DE_marks.txt", sep = "_"),
                                    header = T,
                                    sep = "\t")

            Data_DE_g <- cbindX(Data_DE_g, Data_DE)

            write.table(Data_DE_g,
                        file = paste(prefixo, "DE_marks.txt", sep = "_"),
                        sep = "\t",
                        quote = F,
                        row.names = F,
                        col.names = T)
        }

        factor <- rownames(results_sig)
        counts <- data.frame()
        for (i in 1:length(factor)){
            counts <- rbind(counts, subset(clone, row.names(clone) == factor[i]))
        }

        data_stats <- apply(counts, 2, descritivas)
        data_stats <- do.call("rbind", data_stats)

        if (file.exists(paste(prefixo, "DE_stats.txt", sep = "_")) == "FALSE"){

            write.table(data_stats,
                        file = paste(prefixo, "DE_stats.txt", sep = "_"),
                        sep = "\t",
                        quote = F,
                        row.names = T,
                        col.names = T)

        }else if (file.exists(paste(prefixo, "DE_stats.txt", sep = "_")) == "TRUE"){

            data_stats_g <- read.table(paste(prefixo, "DE_stats.txt", sep = "_"),
                                       header = T,
                                       sep = "\t")

            data_stats_g <- rbind(data_stats_g, data_stats)

            write.table(data_stats_g,
                        file = paste(prefixo, "DE_stats.txt", sep = "_"),
                        sep = "\t",
                        quote = F,
                        row.names = T,
                        col.names = T)
        }
    }
}

with_tec_reps <- snakemake@params[["samples_with_tec_reps"]]
number_of_tec_rep <- as.numeric(snakemake@params[["number_of_tec_rep"]])
samples_without_rep <- snakemake@params[["samples_without_rep"]]
g <- as.numeric(snakemake@params[["groups"]])

# Executes the function for the samples of BRASUZ1
for (t in 1:length(snakemake@params[["tissues"]])){
    deseq2_with_counts(file = snakemake@input[2],
                       prefix = paste(snakemake@params[["prefix"]], g, sep = ""),
                       clone_name = snakemake@params[["genotypes"]],
                       tissue = snakemake@params[["tissues"]][t],
                       fdr = snakemake@params[["fdr"]],
                       log_fc = snakemake@params[["log_fold_change"]],
                       sep_into = snakemake@params[["sep_into"]],
                       subset_model = snakemake@params[["subset_model"]],
                       no_bio_rep = snakemake@params[["no_bio_rep"]],
                       dispersion = snakemake@params[["dispersion"]],
                       filter = snakemake@params[["filtration_mode"]],
                       min_msp = number_of_tec_rep * as.numeric(snakemake@params[["min_msp"]]),
                       intersect_file = paste(snakemake@input[1]))
}
