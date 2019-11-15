"Usage: cpgs_validation.R [--cpg_file <cpg>] (--out1 <O1>) (--out2 <O2>) <input1> <input2>
-h --help    show this
--out1   name1    specify the name for the first output file
cpgs_validation.R -h | --help  show this message
" -> doc

# load the docopt library
require(docopt)

# retrieve the command-line arguments
opts <- docopt(doc)
save.image(file = "docopt.Rda")

require(plyr)
require(berryFunctions)
require(alluvial)
require(reshape2)
require(googleVis)
require(data.table)

cpg_info_file <- paste(opts$`<cpg>`) #cpgs_validation_methy_info_CpG_raw.rda

# Help function
correct_sites <- function(sites_bed_file = sites_bed_file, data = data){

    all_sites <- read.table(sites_bed_file, sep = "\t")
    all_sites <- unique(as.character(all_sites[, 4]))

    dup_sites <- all_sites[grep(pattern = "|", all_sites, fixed = T)]

    dup_sites_names <- strsplit(x = dup_sites, split = "|", fixed = T)
    dup_sites_names <- unique(c(do.call("rbind", dup_sites_names)))

    new_data_corrected <- data.frame()
    for (i in 1:ncol(data)){

        new_data_elements <- data[i]
        new_data_renamed_all <- data.frame()
        for (j in 1:nrow(new_data_elements)){

            if (as.character(new_data_elements[j, 1]) %in% dup_sites_names){

                x <- grep(pattern = paste0("\\<", new_data_elements[j, 1], "\\|"), x = dup_sites)
                y <- grep(pattern = paste0("\\|", new_data_elements[j, 1], "\\>"), x = dup_sites)

                if (length(x) > 0){
                    z <- x
                }else if (length(y) > 0){
                    z <- y
                }else{
                    print("Error!")
                }

                new_data_renamed <- as.data.frame(dup_sites[z])
                colnames(new_data_renamed) <- "new_data_renamed"
                new_data_renamed_all <- rbind(new_data_renamed_all, new_data_renamed)

            }else{

                new_data_renamed <- as.data.frame(new_data_elements[j, 1])
                colnames(new_data_renamed) <- "new_data_renamed"
                new_data_renamed_all <- rbind(new_data_renamed_all,
                                              new_data_renamed)

            }

        }
        if (i == 1){
            colnames(new_data_renamed_all) <- colnames(data)[i]
            new_data_corrected <- as.data.frame(unique(new_data_renamed_all))

        }else{
            colnames(new_data_renamed_all) <- colnames(data)[i]
            new_data_corrected <- cbindX(new_data_corrected,
                                         as.data.frame(unique(new_data_renamed_all)))
        }
    }

    return(new_data_corrected)
}

# Check if the analysis was already realized
if (file.exists(cpg_info_file) == TRUE){

    load("cpgs_validation_methy_info_CpG_raw.rda")

}else{

    clusters_CpG <- as.data.frame(fread(opts$`<input1>`, sep = "\t"))

    true_clusters_CpG_names <- unique(clusters_CpG[duplicated(clusters_CpG$V7) == TRUE, 7])
    true_clusters_CpG <- clusters_CpG[clusters_CpG$V7 %in% true_clusters_CpG_names, ]

    brasuz_sites_CpG <- true_clusters_CpG

    min_fdr <- 0.01

    brasuz_sites_CpG$classification <- ifelse(brasuz_sites_CpG$V5 < min_fdr,
                                              "Methylated",
                                              "Unmethylated")

    methy_info_CpG <- data.frame()

    for (i in 1:length(unique(brasuz_sites_CpG$V7))){

        cluster <- brasuz_sites_CpG[brasuz_sites_CpG$V7 == unique(brasuz_sites_CpG$V7)[i], ]

        # Selects only the MSD sites where both cytosines were sampled by BS-seq
        if (length(grep(cluster$V4, pattern = "%")) == 2){
            dart_site <- paste(as.character(cluster[-grep(cluster$V4, pattern = "%"), 4]),
                               sep = "-",
                               collapse = "-")

            CpGs <- cluster[grep(cluster$V4, pattern = "%"), ]
            CpG_plus_level <- as.character(CpGs[CpGs$V6 == "+", 4])
            CpG_plus_status <- as.character(CpGs[CpGs$V6 == "+", 8])
            CpG_minus_level <- as.character(CpGs[CpGs$V6 == "-", 4])
            CpG_minus_status <- as.character(CpGs[CpGs$V6 == "-", 8])

            methy_info_CpG <- rbind(methy_info_CpG,
                                    data.frame(msdartseq = dart_site,
                                               CpG_plus_level = CpG_plus_level,
                                               CpG_plus_status = CpG_plus_status,
                                               CpG_minus_level = CpG_minus_level,
                                               CpG_minus_status = CpG_minus_status))
        }
    }

    save(file = "cpgs_validation_methy_info_CpG_raw.rda", methy_info_CpG)
}

## Methylated sites ##
methy_info <- methy_info_CpG

CpG_plus <- colsplit(string = methy_info[, 2],
                     pattern = "[(]|[%]|[/]",
                     names = c("a",
                               "n_of_reads_CpG_plus",
                               "meth_level_CpG1",
                               "d"))[c(3, 2)]

CpG_minus <- colsplit(string = methy_info[, 4],
                      pattern = "[(]|[%]|[/]",
                      names = c("a",
                                "n_of_reads_CpG_minus",
                                "meth_level_CpG2",
                                "d"))[c(3, 2)]

methy_info_detailed <- data.frame(msdartseq = methy_info[, 1],
                                  CpG_plus = CpG_plus,
                                  CpG_plus_met_status = methy_info[, 3],
                                  CpG_minus = CpG_minus,
                                  CpG_minus_met_status = methy_info[, 5])

# Select only the cytosines supported for at least 3 reads
methy_info_detailed <- methy_info_detailed[methy_info_detailed[, 3] >= 3 & methy_info_detailed[, 6] >= 3, ]

# Makes the classification of the MSD-sites
restrict_site_status <- character()
for (i in 1:nrow(methy_info_detailed)){

    methy_info_site <- methy_info_detailed[i, ]

    if (methy_info_site[, 4] == "Methylated" & methy_info_site[, 7] == "Methylated"){

        restrict_site_status <- c(restrict_site_status, "Full-methylated")

    } else if ((methy_info_site[, 4] == "Unmethylated" &
               methy_info_site[, 7] == "Methylated") |
               (methy_info_site[, 4] == "Methylated" &
               methy_info_site[, 7] == "Unmethylated")) {

        restrict_site_status <- c(restrict_site_status, "Hemi-methylated")

    } else if (methy_info_site[, 4] == "Unmethylated" &
               methy_info_site[, 7] == "Unmethylated"){

        restrict_site_status <- c(restrict_site_status, "Unmethylated")

    }
}

methy_info_detailed$restrict_site_status <- restrict_site_status

most_methylated_CpG <- character()
for (c in 1:nrow(methy_info_detailed)){

    methy_info_site <- methy_info_detailed[c, ]

    if (methy_info_site[, 8] == "Hemi-methylated"){

        meth_level <- ifelse(methy_info_site[, 4] == "Methylated",
                             methy_info_site[, 2],
                             methy_info_site[, 5])

        most_methylated_CpG <- c(most_methylated_CpG, meth_level)

    } else {

        meth_level <- ifelse(methy_info_site[, 2] > methy_info_site[, 5],
                             methy_info_site[, 2],
                             methy_info_site[, 5])

        most_methylated_CpG <- c(most_methylated_CpG, meth_level)

    }
}

methy_info_detailed$most_methylated_CpG <- most_methylated_CpG

write.table(methy_info_detailed,
    paste(opts$O1),
    col.names = T,
    row.names = F,
    quote = F,
    sep = "\t")

##########################################################################################
##                            Evaluation of the Adult leaves                            ##
##########################################################################################

brasuz_methylated_sites <- read.table(opts$`<input2>`, header = T, sep = "\t")
brasuz_methylated_sites <- brasuz_methylated_sites[, 1] # adult leaves

dart_meth_sites_validation <- methy_info_detailed[methy_info_detailed$msdartseq %in% brasuz_methylated_sites, ]

##########################################
##         Methylated Sites             ##
##########################################

classification_meth_CpG <- classify(dart_meth_sites_validation$most_methylated_CpG,
                                    "c",
                                    breaks = seq(0, 100, 10))

dart_meth_sites_validation$classification_meth_CpG <- classification_meth_CpG$index

plot_data_p1 <- count(dart_meth_sites_validation,
                      c("restrict_site_status", "classification_meth_CpG"))

write.table(plot_data_p1,
            paste(opts$O2),
            col.names = T,
            row.names = F,
            quote = F,
            sep = "\t")

###############################################################################
##                            Unmethylated sites                             ##
###############################################################################

brasuz_sites <- read.table("BRASUZ1_leaf.ad_msp_bigger_than_10.txt",
                                        header = F,
                                        sep = "\t")

sites_bed_file <- "position_of_the_sampled_sites/msdartseq_methylation_sites_of_sequenced_fragments_merged.bed"

merged_frags <- read.table(sites_bed_file, sep = "\t")
merged_frags <- merged_frags[merged_frags$V6 == "+", ]

write.table(merged_frags,
            "position_of_the_sampled_sites/msdartseq_methylation_sites_of_sequenced_fragments_merged_plus.bed",
            col.names = F,
            row.names = F,
            sep = "\t",
            quote = F)

sites_bed_file <- "position_of_the_sampled_sites/msdartseq_methylation_sites_of_sequenced_fragments_merged_plus.bed"

brasuz_unmethylated_sites <- correct_sites(sites_bed_file = sites_bed_file,
                                         data = brasuz_sites)

teste <- brasuz_unmethylated_sites[brasuz_unmethylated_sites$V1 %in% brasuz_methylated_sites, ]

brasuz_unmethylated_sites <- brasuz_unmethylated_sites[!brasuz_unmethylated_sites$V1 %in% brasuz_methylated_sites, ]

unmetylated_bs <- methy_info_detailed[methy_info_detailed$msdartseq %in% brasuz_unmethylated_sites, ]

count(unmetylated_bs, vars = "restrict_site_status")

write.table(plot_data_p1,
            paste(opts$O2),
            col.names = T,
            row.names = F,
            quote = F,
            sep = "\t")
