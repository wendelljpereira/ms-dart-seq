"Usage: marks_closest_restriction_site_search.R (--out1 <O1>) (--out2 <O2>) (--out3 <O3>) <input1> <input2> <input3>
-h --help    show this help
--out1  name1   bed file with the position of all sequenced fragments
--out2  name2   bed file with the position of the methylation sited of the sequenced fragments
--out3  name3   csv with the counts
input1  input1  bed file with the position of all MS-DArT tags
input2  input2  bed file with the position of all restriction sites in the genome
input3  input3  tst file with the counts of each Tag for each sample
marks_closest_restriction_site_search.R -h | --help  show this message
" -> doc

# Loads the docopt library
require(docopt)
require(tidyverse)

# Retrieves the command-line arguments
opts <- docopt(doc)

# Reads the bed file with mapping position for each mark.
# This will be used to search the closest restriction site upstream and
# dowstream of each mark.
marks_pos <- read_tsv(opts$`<input1>`,
                      col_names = FALSE,
                      col_types = "ciiccc")

# Reads a bed file containing the restriction sites positions.
sites_bed <- read_tsv(opts$`<input2>`,
                      col_names = FALSE,
                      col_types = "ciiccc",
                      progress = F)

sites <- sites_bed[sites_bed[, 6] == "+", ]
names(sites) <- c("chr", "start", "end", "name", "score", "strand")

# Searches for the closest restriction sites of each MSD-Tag
closest_sites_search <- function(i, marks_pos, sites){

    # update progress bar
    Sys.sleep(0.1)
    setTxtProgressBar(pb, i)

    mark_target <- as.character(marks_pos[i, 4])
    start_position <- as.numeric(marks_pos[i, 2])
    end_position <- as.numeric(marks_pos[i, 3])
    feature <- as.character(marks_pos[i, 1])
    strand_target <- as.character(marks_pos[i, 6])

    sites_feature <- sites[sites[, 1] == feature, ]

    if (strand_target == "+"){

        if (as.numeric(min(sites_feature$end[sites_feature$end >= end_position])) == "Inf"){

            enzime_site_downstream <- as.data.frame(rbind(rep(NA, 6)))

        }else{

            closest_downstream <- as.numeric(min(sites_feature$end[sites_feature$end >= end_position]))
            enzime_site_downstream <- sites_feature[sites_feature$end == closest_downstream, ]

        }

        # Organizes to bed format
        enzime_site_downstream <- cbind(enzime_site_downstream[, c(1, 2, 3, 4)],
                                        as.data.frame(0),
                                        as.data.frame(paste(enzime_site_downstream[, 6])))

        target_fragment <- cbind(as.data.frame(feature),
                                 start_position,
                                 enzime_site_downstream[, 3],
                                 as.data.frame(mark_target),
                                 as.data.frame(0),
                                 as.data.frame(strand_target))

        enzimes_frag_position <- cbind(target_fragment, enzime_site_downstream)

        colnames(enzimes_frag_position) <- c("mark_feature",
                                             "mark_frag_start",
                                             "mark_frag_end",
                                             "mark ID",
                                             "score",
                                             "mark strand",
                                             "complement_site_feature",
                                             "complement_site_start",
                                             "complement_site_end",
                                             "complement_site_enzime",
                                             "score",
                                             "complement_site_strand")


    } else if (strand_target == "-") {

        if (as.numeric(max(sites_feature$start[sites_feature$start <= start_position])) == "-Inf"){

            enzime_site_upstream <- as.data.frame(rbind(rep(NA, 6)))

        } else {

            closest_upstream <- as.numeric(max(sites_feature$start[sites_feature$start <= start_position]))
            enzime_site_upstream <- as.data.frame(sites_feature[sites_feature$start == closest_upstream, ])

        }

        enzime_site_upstream <- cbind(enzime_site_upstream[, c(1, 2, 3, 4)],
                                      as.data.frame(0),
                                      as.data.frame(paste(enzime_site_upstream[, 6])))

        target_fragment <- cbind(as.data.frame(feature),
                                 enzime_site_upstream[, 3],
                                 end_position,
                                 as.data.frame(mark_target),
                                 as.data.frame(0),
                                 as.data.frame(strand_target))

        enzimes_frag_position <- cbind(target_fragment, enzime_site_upstream)

        colnames(enzimes_frag_position) <- c("mark_feature",
                                             "mark_frag_start",
                                             "mark_frag_end",
                                             "mark ID",
                                             "score",
                                             "mark strand",
                                             "complement_site_feature",
                                             "complement_site_start",
                                             "complement_site_end",
                                             "complement_site_enzime",
                                             "score",
                                             "complement_site_strand")
    }

    enzimes_frag_position_complete <- as.data.frame(enzimes_frag_position)

}

# Executes the function
# Set a progress bar
pb <- txtProgressBar(min = 0, max = as.numeric(nrow(marks_pos)), style = 3)
closest_sites <- lapply(1:as.numeric(nrow(marks_pos)),
                        closest_sites_search,
                        marks_pos = marks_pos,
                        sites = sites)
close(pb)

closest_sites_df <- do.call("rbind", closest_sites)

save.image(file = "fragments.Rda")

# Removes the sites where it is not possible to identify restriction sites at the ends of the fragment
closest_sites_df <- as.data.frame(closest_sites_df[complete.cases(closest_sites_df), ])

closest_sites_df$mark_frag_start <- as.numeric(closest_sites_df$mark_frag_start)
closest_sites_df$mark_frag_end <- as.numeric(closest_sites_df$mark_frag_end)

# Removes the sites that do not have a MspI/HpaII sites
closest_sites_df_pst_msp <- closest_sites_df[closest_sites_df$complement_site_enzime == "MspI", ]

# sort the files by genomic position
closest_sites_df_pst_msp <- closest_sites_df_pst_msp[with(closest_sites_df_pst_msp, order(mark_feature, mark_frag_start)), ]

write.table(closest_sites_df_pst_msp,
            paste(opts$O1),
            col.names = F,
            row.names = F,
            quote = F,
            sep = "\t")

# Determines the methylation position

## frags in plus strand
pst_mspI_plus <- closest_sites_df_pst_msp[closest_sites_df_pst_msp$`mark strand` == "+", ]

methyl_plus_cpg_plus <- pst_mspI_plus[pst_mspI_plus$complement_site_enzime == "MspI", ]
methyl_plus_cpg_plus$mark_frag_start <- methyl_plus_cpg_plus$complement_site_start + 1
methyl_plus_cpg_plus$mark_frag_end <- methyl_plus_cpg_plus$complement_site_start + 2

methyl_plus_cpg_minus <- pst_mspI_plus[pst_mspI_plus$complement_site_enzime == "MspI", ]
methyl_plus_cpg_minus$mark_frag_start <- methyl_plus_cpg_minus$complement_site_start + 2
methyl_plus_cpg_minus$mark_frag_end <- methyl_plus_cpg_minus$complement_site_start + 3
methyl_plus_cpg_minus$`mark strand` <- "-"

methyl_plus_df <- rbind(methyl_plus_cpg_plus[, 1:6], methyl_plus_cpg_minus[, 1:6])

## frags in minus strand
pst_mspI_minus <- closest_sites_df_pst_msp[closest_sites_df_pst_msp$`mark strand` == "-", ]

methyl_minus_cpg_plus <- pst_mspI_minus[pst_mspI_minus$complement_site_enzime == "MspI", ]
methyl_minus_cpg_plus$mark_frag_start <- methyl_minus_cpg_plus$complement_site_start + 1
methyl_minus_cpg_plus$mark_frag_end <- methyl_minus_cpg_plus$complement_site_start + 2
methyl_minus_cpg_plus$`mark strand` <- "+"

methyl_minus_cpg_minus <- pst_mspI_minus[pst_mspI_minus$complement_site_enzime == "MspI", ]
methyl_minus_cpg_minus$mark_frag_start <- methyl_minus_cpg_minus$complement_site_start + 2
methyl_minus_cpg_minus$mark_frag_end <- methyl_minus_cpg_minus$complement_site_start + 3

methyl_minus_df <- rbind(methyl_minus_cpg_plus[, 1:6], methyl_minus_cpg_minus[, 1:6])

sites_tested_to_methylation <- rbind(methyl_plus_df, methyl_minus_df)

# sort by genomic position
sites_tested_to_methylation <- sites_tested_to_methylation[with(sites_tested_to_methylation, order(mark_feature, mark_frag_start)), ]

write.table(sites_tested_to_methylation,
            paste(opts$O2),
            col.names = F,
            row.names = F,
            quote = F,
            sep = "\t")

# Removes sites for what was not possible to define de restrictions sites from the counts file
counts <- read.table(opts$`<input3>`, sep = "\t", header = T, check.names = F)
counts <- counts[counts$Geneid %in% closest_sites_df_pst_msp[, 4], ]

write.table(counts, paste(opts$O3), col.names = T, row.names = F, sep = ",")
