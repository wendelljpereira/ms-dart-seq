"Usage: counts_correction.R (--out1 <O1>) (--out2 <O2>) <input1> <input2>
-h --help    show this
--out1   bed_file    bed file with MS-DArT-seq tags with corrected positions
--out2   tst_file    table with the counts of the MS-DArT-seq tags with corrected positions
<input1>   bed_file    bed file with MS-DArT-seq tags positions
<input2>   tst_file    table with the counts of the MS-DArT-seq tags
counts_correction.R -h | --help  show this message
" -> doc

require(docopt)
require(tidyverse)

# retrieve the command-line arguments
opts <- docopt(doc)

all_positions <- read.table(opts$`<input1>`,
                            sep = "\t")

all_counts <- read.table(opts$`<input2>`,
                         sep = "\t",
                         header = T,
                         check.names = F)

clusters_2_or_more_seqs_names <- unique(all_positions[duplicated(all_positions$V7) == TRUE, 7])

temp_bed <- all_positions[all_positions$V7 %in% clusters_2_or_more_seqs_names, ]
temp_counts <- all_counts[all_counts$Geneid %in% temp_bed$V4, ]

## Corrects counts using the MSD-Tags as a reference
new_features <- data.frame()
for (i in unique(temp_bed$V7)){

    temp_bed_subset <- temp_bed[temp_bed$V7 == i, ]
    counts_subset <- temp_counts[temp_counts$Geneid %in% temp_bed_subset$V4, ]

    # Gives a way to separate the MSD-Tags from the Help-Tags
    counts_subset$reg_class <- temp_bed_subset$V5

    if(unique(counts_subset$Strand) == "+"){

        counts_subset <- counts_subset[with(counts_subset, order(reg_class, desc(End))), ]

    }else if(unique(counts_subset$Strand) == "-"){

        counts_subset <- counts_subset[with(counts_subset, order(reg_class, Start)), ]

    }

    counts_subset <- counts_subset[, -ncol(counts_subset)]

    ## 2 overlaping Tags
    if(nrow(counts_subset) == 3){

        biggest_feature <- cbind(counts_subset[1, c(1:6)],
                                 counts_subset[3, -c(1:6)])

        colnames(biggest_feature) <- names(counts_subset)

        smallest_feature <- cbind(counts_subset[2, c(1:6)],
                                  counts_subset[1, -c(1:6)] - counts_subset[3, -c(1:6)])
        colnames(smallest_feature) <- names(counts_subset)

        features_normalized <- rbind(biggest_feature, smallest_feature)

        new_features <- rbind(new_features, features_normalized)

    }else if(nrow(counts_subset) == 5){

        biggest_feature <- cbind(counts_subset[1, c(1:6)], counts_subset[4, -c(1:6)])
        colnames(biggest_feature) <- names(counts_subset)

        middle_feature <- cbind(counts_subset[2, c(1:6)], counts_subset[5, -c(1:6)] - counts_subset[4, -c(1:6)])
        colnames(middle_feature) <- names(counts_subset)

        smallest_feature <- cbind(counts_subset[3, c(1:6)], counts_subset[3, -c(1:6)] - counts_subset[5, -c(1:6)])
        colnames(smallest_feature) <- names(counts_subset)

        features_normalized <- rbind(biggest_feature, middle_feature, smallest_feature)

        new_features <- rbind(new_features, features_normalized)

    }else if(nrow(counts_subset) == 7){

        biggest_feature <- cbind(counts_subset[1, c(1:6)], counts_subset[5, -c(1:6)])
        colnames(biggest_feature) <- names(counts_subset)

        middle_feature <- cbind(counts_subset[2, c(1:6)], counts_subset[6, -c(1:6)] - counts_subset[5, -c(1:6)])
        colnames(middle_feature) <- names(counts_subset)

        middle_feature_2 <- cbind(counts_subset[3, c(1:6)], counts_subset[7, -c(1:6)] - counts_subset[6, -c(1:6)])
        colnames(middle_feature_2) <- names(counts_subset)

        smallest_feature <- cbind(counts_subset[4, c(1:6)], counts_subset[4, -c(1:6)] - counts_subset[7, -c(1:6)])
        colnames(smallest_feature) <- names(counts_subset)

        features_normalized <- rbind(biggest_feature, middle_feature, middle_feature_2, smallest_feature)

        new_features <- rbind(new_features, features_normalized)

    }else{
        print(paste("Warning: More than 5 features in cluster ", i, "!", sep = ""))
    }
}

# Removes the original counts of the sites that were corrected
all_counts <- all_counts[!all_counts$Geneid %in% temp_counts$Geneid, ]

# Adds the new counts (corrected)
colnames(new_features) <- colnames(all_counts)
all_counts <- rbind(all_counts, new_features)

all_counts <- all_counts[order(all_counts$Chr, all_counts$Start), ]

all_counts$Geneid <- paste("MS-DArT_site", seq(1, nrow(all_counts)), sep = "_")

all_counts <- all_counts[, -c(2,3,4,5,6)]

all_positions <- all_positions[!all_positions$V4 %in% temp_bed$V4, ]
new_features_bed <- new_features[, c(2,3,4,1,6,5)]
new_features_bed$Length <- 0

all_positions <- all_positions[,1:6]
all_positions$V5 <- 0

colnames(new_features_bed) <- names(all_positions)
all_positions <- rbind(all_positions, new_features_bed)

all_positions <- arrange(all_positions, V1, V2)

all_positions$V4 <- paste("MS-DArT_site", seq(1, nrow(all_positions)), sep = "_")

write.table(all_positions,
            paste(opts$O1),
            col.names = F,
            row.names = F,
            sep = "\t",
            quote = F)

write.table(all_counts,
            paste(opts$O2),
            col.names = T,
            row.names = F,
            sep = "\t",
            quote = F)
