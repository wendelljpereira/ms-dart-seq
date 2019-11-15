require(edgeR)
require(stringr)
require(ggfortify)
require(gridExtra)
require(plyr)
require(VennDiagram)
require(ggthemes)

# All samples and marks with true sites
counts_data <- read.table(snakemake@input[[1]],
                          header = T, sep = ",",
                          check.names = F,
                          colClasses = c("character", rep("numeric", 18)))

# File with information about marks with MspI counts larger than 0 in all samples.
true_sites <- read.table(snakemake@input[[2]], header = T, sep = "\t")[, 1]
counts_data <- counts_data[counts_data[, 1] %in% true_sites, ][, -1]

counts_data_t <- as.data.frame(t(counts_data))

clones_names <- as.data.frame(row.names(counts_data_t))
colnames(clones_names) <- "clone"
clones_names <- as.data.frame(str_split_fixed(clones_names$clone, "_", 4))
colnames(clones_names) <- c("group", "clone", "tissue", "enzime")
clones_names$enzime <- substr(clones_names$enzime, 1, 2)

clones_names$clone_tissue <- paste(clones_names$clone,
                                   clones_names$tissue,
                                   sep = "_")

clones_names$clone_tissue_enzime <- paste(clones_names$clone,
                                          clones_names$tissue,
                                          clones_names$enzime,
                                          sep = "_")

# Removes sites with counts equal to zero
subset_index <- as.data.frame(lapply(counts_data_t, function(x) sum(x) > 0))
counts <- counts_data_t[, subset_index == TRUE]

## Executes the PCA
res.pca <- prcomp(counts, scale = TRUE, center = TRUE)
x <- summary(res.pca)

A <- autoplot(res.pca,
              data = clones_names,
              colour = "clone_tissue",
              shape = "enzime", size = 3,
              frame = TRUE,
              frame.type = "norm") +
    theme(legend.text = element_text(size = 14),
          axis.title.y = element_text(size = 15),
          axis.title.x = element_text(size = 15),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          legend.title = element_blank(),
          strip.text = element_text(size = 15))+
    scale_colour_colorblind()+
    scale_fill_colorblind()+
    theme_bw()

svg(filename = snakemake@output[[1]], width = 10, height = 7, pointsize = 12)
A
dev.off()

# All samples, but only DM marks
counts_data <- read.table(snakemake@input[[1]],
                          header = T,
                          sep = ",",
                          check.names = F,
                          colClasses = c("character", rep("numeric", 18)))

# File with information about DM marks in the samples.
DM_sites <- read.table(snakemake@input[[3]],
                       header = T,
                       sep = "\t",
                       colClasses = c("character"))

DM_sites <- unique(c(DM_sites[, 1], DM_sites[, 2], DM_sites[, 3]))

counts_data <- counts_data[counts_data[, 1] %in% DM_sites, ][, -1]

counts_data_t <- as.data.frame(t(counts_data))

clones_names <- as.data.frame(row.names(counts_data_t))
colnames(clones_names) <- "clone"
clones_names <- as.data.frame(str_split_fixed(clones_names$clone, "_", 4))
colnames(clones_names) <- c("group", "clone", "tissue", "enzime")
clones_names$enzime <- substr(clones_names$enzime, 1, 2)
clones_names$clone_tissue <- paste(clones_names$clone,
                                   clones_names$tissue,
                                   sep = "_")

clones_names$clone_tissue_enzime <- paste(clones_names$clone,
                                          clones_names$tissue,
                                          clones_names$enzime,
                                          sep = "_")

subset_index <- as.data.frame(lapply(counts_data_t, function(x) sum(x) > 0))
counts <- counts_data_t[, subset_index == TRUE]

res.pca <- prcomp(counts, scale = TRUE, center = TRUE)
x <- summary(res.pca)

B <- autoplot(res.pca,
              data = clones_names,
              colour = "clone_tissue",
              shape = "enzime",
              size = 3,
              frame = TRUE,
              frame.type = "norm") +
    theme(legend.text = element_text(size = 14),
          axis.title.y = element_text(size = 15),
          axis.title.x = element_text(size = 15),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          legend.title = element_blank(),
          strip.text = element_text(size = 15))+
    scale_colour_colorblind()+
    scale_fill_colorblind()+
    theme_bw()

svg(filename = snakemake@output[[2]], width = 10, height = 7, pointsize = 12)
B
dev.off()

## PCA using only MspI counts of the methylated sites.

counts_msp <- counts[-c(1:9), ]

res.pca <- prcomp(counts_msp, scale = TRUE, center = TRUE)
x <- summary(res.pca)

clones_names_ms <- clones_names[clones_names$enzime == "ms", ]

C <- autoplot(res.pca,
              data = clones_names_ms,
              colour = "clone_tissue",
              shape = "enzime",
              size = 3,
              frame = TRUE,
              frame.type = "norm") +
    theme(legend.text = element_text(size = 14),
          axis.title.y = element_text(size = 15),
          axis.title.x = element_text(size = 15),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          legend.title = element_blank(),
          strip.text = element_text(size = 15))+
    scale_colour_colorblind()+
    scale_fill_colorblind()+
    theme_bw()

svg(filename = "images/PCAs/pca_group4_only_DM_sites_only_MS.svg", width = 10, height = 7, pointsize = 12)
C
dev.off()
