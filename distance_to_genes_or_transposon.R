require(ggplot2)
require(gdata)
options(scipen=999)

# Closest gene of each mark
all_genes <- read.table(snakemake@input[[1]], colClasses = "character")
all_genes$V16 <- as.numeric(all_genes$V16)

# Closest TE of each mark
all_transposons <- read.table(snakemake@input[[2]],
                              sep = "\t",
                              colClasses = "character")

all_transposons$V13 <- as.numeric(all_transposons$V13)

# Cleaning
transposons <- cbind(as.data.frame(rep("transposon", nrow(all_transposons))),
                     all_transposons)
transposons <- transposons[!transposons$V11 == ".", ]
transposons <- transposons[, c(1, 5, 14)]
colnames(transposons) <- c("feature", "mark", "distance")

genes <- cbind(as.data.frame(rep("gene", nrow(all_genes))), all_genes)
genes <- genes[genes$V9 == "gene", ]
genes <- genes[, c(1, 5, 17)]
colnames(genes) <- c("feature", "mark", "distance")

features <- rbind(genes, transposons)

features$distance <- abs(features$distance)
features <- cbind(as.data.frame(rep("Mapped reads", nrow(features))), features)
colnames(features) <- c("type", "feature", "mark", "distance")

features_10kb <- features[features$distance <= 10000, ]

# To use in the combined plot
methy_sites_10kb <- features_10kb

bin <- c(1, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10001)
features_10kb <- cbind(features_10kb, findInterval(features_10kb$distance, bin))
colnames(features_10kb) <- c("type", "feature", "mark", "distance", "groups")
features_10kb$groups <- as.factor(features_10kb$groups)

#### >>>>>>>>>>> Integrar no snakemake
write.table(features_10kb,
            "distance_to_genes_and_tes_methylated_sites.tst",
            col.names = T,
            row.names = F,
            sep = "\t")

interval_labels <- c("0" = "0",
                     "1" = "0.001-1",
                     "2" = "1-2",
                     "3" = "2-3",
                     "4" = "3-4",
                     "5" = "4-5",
                     "6" = "5-6",
                     "7" = "6-7",
                     "8" = "7-8",
                     "9" = "8-9",
                     "10" = "9-10")

svg(filename = snakemake@output[[1]], width = 8, height = 6, pointsize = 12)
ggplot(features_10kb, aes(x = groups)) +
    geom_bar(aes(fill = feature), position = "dodge", alpha = 0.7) +
    xlab("Distance (Kbp)") +
    ylab("Number of methylated sites") +
    scale_y_continuous(breaks = seq(0, 800, 100), limits = c(0, 800)) +
    scale_x_discrete(labels = interval_labels) +
    theme_bw() +
    theme(legend.text = element_text(size = 15),
          axis.title.y = element_text(size = 16, vjust = 2),
          axis.title.x = element_text(size = 16, vjust = 0),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          legend.title = element_blank(),
          strip.text = element_text(size = 15),
          plot.title = element_text(size = 17, face = "bold")) +
    scale_fill_manual(
        labels = c("Genes", "Transposons", "Transposons in genes"),
        values = c("Darkblue", "gold3", "red")) +
    ggtitle(paste("Full set - Distance to genes and transposons"))
dev.off()

# Reads the methylated sites
metil <- read.table(snakemake@input[[3]], header = T, colClasses = "character")

colnames(metil) <- c("BRASUZ1 Adult Leaves",
                     "BRASUZ1 Juvenile Leaves",
                     "BRASUZ1 Xylem")

for (i in 1:length(metil)){
    tissue <- colnames(metil)[i]
    marks <- unique(metil[, i])
    bar_plot_sub <- features_10kb[features_10kb$mark %in% as.character(marks), ]

    write.table(bar_plot_sub,
                paste0(tissue, "_distance_to_genes_and_tes_methylated_sites.tst"),
                col.names = T,
                row.names = F,
                sep = "\t")

    bar_plot <- ggplot(bar_plot_sub, aes(x = groups)) +
        geom_bar(aes(fill = feature), position = "dodge", alpha = 0.7) +
        xlab("Distance (Kbp)") +
        ylab("Number of methylated sites") +
        scale_y_continuous(breaks = seq(0, 800, 100), limits = c(0, 800)) +
        scale_x_discrete(labels = interval_labels) +
        theme_bw() +
        theme(legend.text = element_text(size = 18),
              axis.title.y = element_text(size = 16, vjust = 2),
              axis.title.x = element_text(size = 16, vjust = 0),
              axis.text.x = element_text(size = 14),
              axis.text.y = element_text(size = 14),
              legend.title = element_blank(),
              strip.text = element_text(size = 15),
              plot.title = element_text(size = 17, face = "bold")) +
        scale_fill_manual(
            labels = c("Genes", "Transposons", "Transposons in genes"),
            values = c("#E69F00", "#0072B2")) +
        ggtitle(paste(tissue, " - Distance to genes and transposons"))

    # Exports the plots
    svg(paste("images/distance_to_genes_and_TEs/",
              tissue,
              "Distance to genes and transposons.svg"),
        width = 8,
        height = 6,
        pointsize = 12)
    print(bar_plot)
    dev.off()
}

####################################
#### Distance of the mspI sites ####
####################################

# Loads the information of the closest genes of each MSD sites
all_genes <- read.table(snakemake@input[[4]],
                        colClasses = "character",
                        sep = "\t")
all_genes$V16 <- as.numeric(all_genes$V16)

# Loads the information of the closest TE of each MSD sites
all_transposons <- read.table(snakemake@input[[5]],
                              sep = "\t",
                              colClasses = "character")
all_transposons$V13 <- as.numeric(all_transposons$V13)

# Cleaning
transposons <- cbind(as.data.frame(rep("transposon", nrow(all_transposons))),
                     all_transposons)
transposons <- transposons[!transposons$V11 == ".", ]
transposons <- transposons[, c(1, 5, 14)]
colnames(transposons) <- c("feature", "mark", "distance")

genes <- cbind(as.data.frame(rep("gene", nrow(all_genes))), all_genes)
genes <- genes[genes$V9 == "gene", ]
genes <- genes[, c(1, 5, 17)]
colnames(genes) <- c("feature", "mark", "distance")

# Join genes and TEs information
features <- rbind(genes, transposons)

# calculates the absolute value of the distance
features$distance <- abs(features$distance)
features <- cbind(as.data.frame(rep("Mapped reads", nrow(features))), features)
colnames(features) <- c("type", "feature", "mark", "distance")

# all_marks_10Kb
features_10kb <- features[features$distance <= 10000, ]

# Removes the marks that are in genes or TEs
features_10kb <- features_10kb[!features_10kb$distance == 0, ]
all_sites_10Kb <- features_10kb

bin <- c(1, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10001)
features_10kb <- cbind(features_10kb, findInterval(features_10kb$distance, bin))
colnames(features_10kb) <- c("type", "feature", "mark", "distance", "groups")
features_10kb$groups <- as.factor(features_10kb$groups)

#### >>>>>>>>>>> Integrar no snakemake
write.table(features_10kb,
            "distance_to_genes_and_tes_MspI_sites.tst",
            col.names = T,
            row.names = F,
            sep = "\t")

svg(filename = snakemake@output[[2]], width = 8, height = 6, pointsize = 12)
ggplot(features_10kb, aes(x = groups)) +
    geom_bar(aes(fill = feature), position = "dodge", alpha = 0.8) +
    xlab("Distance (Kbp)") +
    ylab("Number of MspI sites") +
    scale_y_continuous(breaks = seq(0, 300000, 50000), limits = c(0, 300000)) +
    scale_x_discrete(labels = interval_labels) +
    theme_bw() +
    theme(legend.text = element_text(size = 15),
          axis.title.y = element_text(size = 16, vjust = 2),
          axis.title.x = element_text(size = 16,vjust = 0),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          legend.title = element_blank(),
          strip.text = element_text(size = 15),
          plot.title = element_text(size = 17, face = "bold")) +
    scale_fill_manual(
        labels = c("Genes", "Transposons", "Transposons in genes"),
        values = c("#E69F00", "#0072B2")) +
    ggtitle(paste("MspI - Distance to genes and transposons"))
dev.off()

### Combining both plots in one ###

## ====> Insert the code below in snakemake <==== ##

methy_sites_10kb$site_class <- "methylated_sites"
all_sites_10Kb$site_class <- "all_mspI"

# To use in the combined plot
all_features_comb <- rbind(methy_sites_10kb,
                           all_sites_10Kb)


# Defines the interval of the samples
bin <- c(1, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10001)

all_features_comb <- cbind(all_features_comb,
                           findInterval(all_features_comb$distance, bin))

colnames(all_features_comb) <- c("type",
                                 "feature",
                                 "mark",
                                 "distance",
                                 "site_class",
                                 "groups")

all_features_comb$groups <- as.factor(all_features_comb$groups)

labels <- c(methylated_sites = "Methylated sites",
            all_mspI = "All MspI sites")

## Plot
svg(filename = "images/distance_to_genes_and_TEs/distance_to_genes_and_TEs_combined_plot.svg",
    width = 12,
    height = 6,
    pointsize = 12)

ggplot(all_features_comb, aes(x = groups)) +
    geom_bar(aes(fill = feature), position = "dodge", alpha = 0.7) +
    xlab("Distance (kbp)") +
    ylab("Number of sites") +
    facet_wrap(~ site_class, scales = "free_y", labeller = labeller(site_class = labels)) +
    scale_y_continuous(breaks = scales::pretty_breaks(7), limits = c(0, NA)) +
    scale_x_discrete(labels = interval_labels) +
    theme_bw() +
    theme(legend.text = element_text(size = 15),
          axis.title.y = element_text(size = 16, vjust = 2),
          axis.title.x = element_text(size = 16,vjust = 0),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          legend.title = element_blank(),
          strip.text = element_text(size = 15),
          plot.title = element_text(size = 17, face = "bold")) +
    scale_fill_manual(
        labels = c("Genes", "Transposons", "Transposons in genes"),
        values = c("#E69F00", "#0072B2"))
dev.off()
