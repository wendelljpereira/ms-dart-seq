require(ggplot2)
require(edgeR)
require(plyr)
require(gdata)
require(tidyr)

# Help function to do the differential expression analysis
edgeR_with_counts <- function(file = file,
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
                                colClasses = "character")

        intersect <- as.character(intersect[, 1])

        clone <- clone[rownames(clone) %in% intersect, ]
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
                paste(clone_name, "_", tecido, "_", "msp_bigger_than_threshold.txt", sep = ""),
                row.names = F,
                col.names = F,
                quote = F)
    clone <- filtred_clone

    # Help function to generate some descriptive stats
    descritivas <- function(y) {
        data.frame(n = length(y),
                   mean = mean(y),
                   IC95inf = t.test(y)$conf.int[1],
                   IC95sup = t.test(y)$conf.int[2],
                   var = var(y), std.dev. = sd(y),
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

    # Loads the processed data to edgeR object
    data_edgeR <- DGEList(counts = clone, group = groups)

    if (no_bio_rep == "FALSE"){

        # Determines the dispersion and save the values in a new file
        data_edgeR <- estimateDisp(data_edgeR)

        if (file.exists(paste(prefixo, "dispersions.txt", sep = "_")) == "FALSE"){

            disp_value <- as.data.frame(rbind(data_edgeR$common.dispersion))

            colnames(disp_value) <- "Common_dispersion"
            rownames(disp_value) <- paste(clone_name, tecido, sep = "_")

            write.table(disp_value,
                        file = paste(prefixo, "dispersions.txt", sep = "_"),
                        sep = "\t",
                        quote = F,
                        row.names = T,
                        col.names = T)

        } else if (file.exists(paste(prefixo, "dispersions.txt", sep = "_")) == "TRUE"){

            disp_g <- read.table(paste(prefixo, "dispersions.txt", sep = "_"),
                                 sep = "\t",
                                 header = T,
                                 row.names = 1)

            disp_value <- rbind(data_edgeR$common.dispersion)

            rownames(disp_value) <- paste(clone_name, tecido, sep = "_")
            colnames(disp_value) <- "Common_dispersion"

            disp_g <- rbind(disp_g, disp_value)

            write.table(disp_g,
                        file = paste(prefixo, "dispersions.txt", sep = "_"),
                        sep = "\t",
                        quote = F,
                        row.names = T,
                        col.names = T)
        }

        # Determines the significantly different counts between MspI and HpaII
        tags <- exactTest(data_edgeR, pair = c("hp", "ms"))

    } else if (no_bio_rep == "TRUE"){

        if (file.exists(paste(prefixo, "dispersions.txt", sep = "_")) == "FALSE"){

            disp_value <- as.data.frame(rbind(as.data.frame(dispersion)))

            colnames(disp_value) <- "Common_dispersion"
            rownames(disp_value) <- paste(clone_name, tecido, sep = "_")

            write.table(disp_value,
                        file = paste(prefixo, "dispersions.txt", sep = "_"),
                        sep = "\t",
                        quote = F,
                        row.names = T,
                        col.names = T)

        } else if (file.exists(paste(prefixo, "dispersions.txt", sep = "_")) == "TRUE"){

            disp_g <- read.table(paste(prefixo, "dispersions.txt", sep = "_"),
                                 sep = "\t",
                                 header = T,
                                 row.names = 1)

            disp_value <- rbind(as.data.frame(dispersion))

            rownames(disp_value) <- paste(clone_name, tecido, sep = "_")
            colnames(disp_value) <- "Common_dispersion"

            disp_g <- rbind(disp_g, disp_value)

            write.table(disp_g,
                        file = paste(prefixo, "dispersions.txt", sep = "_"),
                        sep = "\t",
                        quote = F,
                        row.names = T,
                        col.names = T)

        }

        tags <- exactTest(data_edgeR,
                          pair = c("hp", "ms"),
                          dispersion = as.numeric(dispersion))

    }else{
        print("wrong value of the no_bio_rep argument value!")
    }

    tTags <- topTags(tags, n = NULL)
    results <- as.data.frame(tTags)

    # Filter by FoldChange and FDR
    if (filter == "posit"){

        results_sig <- subset(results, results$FDR <= fdr)
        results_sig <- subset(results_sig, results_sig$logFC >= log_fc)

    }else if (filter == "negat"){

        results_sig <- subset(results, results$FDR <= fdr)
        results_sig <- subset(results_sig, results_sig$logFC <= 0)
        results_sig <- subset(results_sig, abs(results_sig$logFC) >= log_fc)

    } else if (filter == "full"){

        results_sig <- subset(results, results$FDR <= fdr)
        results_sig <- subset(results_sig, abs(results_sig$logFC) >= log_fc)

    } else {

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

    } else {

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

        # Descriptive stats
        factor <- rownames(results_sig)
        counts <- data.frame()

        for (i in 1:length(factor)){

            counts <- rbind(counts,
                            subset(clone, row.names(clone) == factor[i]))

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
    edgeR_with_counts(file = snakemake@input[2],
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
