"Usage:
    reproducibility.R (-o <out> | --output=<out>) (-co <out2> | --correlation_output=<out2>) (<input1>) (<input2>) (<input3>)
    reproducibility.R -h | --help
Options:
    -h --help  Show this options panel.
    -o, --output <out>  Specify the name for the output file. This output is a image showing the distributions curves for the three replicates from each tissue.
    -co, --correlation_output <out2>  Specify the prefiz for the correlation plots.
Files:
    input1  Marks without overlap file
    input2  Counts file
    input3  File with all methylated sites
" -> doc

# loads the docopt library
require(docopt)
# retrieves the command-line arguments
opts <- docopt(doc)

# Loads the necessary packages
require(reshape2)
require(ggplot2)
require(ggsci)
require(gridExtra)
require(corrplot)
require(RColorBrewer)
require(gridGraphics)

# Loads marks without overlaps.
marks <- read.table(opts$input1)
marks <- unique(marks)

# Loads counts file.
counts <- read.table(opts$input2,
                     header = T,
                     sep = "\t",
                     quote = "\"",
                     check.names = F)

counts_file <- counts[counts[, 1] %in% marks[, 1], ]
counts_file[] <- lapply(counts_file, function(x) as.numeric(x))

counts_file_ms <- counts_file[, grep("_ms", colnames(counts_file))]
counts_file_hp <- counts_file[, grep("_hp", colnames(counts_file))]

ad.leaves <- counts_file_ms[, grep("leaf.ad", colnames(counts_file_ms))]
juv.leaves <- counts_file_ms[, grep("leaf.juv", colnames(counts_file_ms))]
xylem <- counts_file_ms[, grep("wood", colnames(counts_file_ms))]

ad.leaves_plot_rep1 <- data.frame(sample = rep("rep1", nrow(ad.leaves)), counts = ad.leaves[, 1])
ad.leaves_plot_rep2 <- data.frame(sample = rep("rep2", nrow(ad.leaves)), counts = ad.leaves[, 2])
ad.leaves_plot_rep3 <- data.frame(sample = rep("rep3", nrow(ad.leaves)), counts = ad.leaves[, 3])

ad.leaves_plot <- rbind(ad.leaves_plot_rep1, ad.leaves_plot_rep2, ad.leaves_plot_rep3)

juv.leaves_plot_rep1 <- data.frame(sample = rep("rep1", nrow(juv.leaves)), counts = juv.leaves[, 1])
juv.leaves_plot_rep2 <- data.frame(sample = rep("rep2", nrow(juv.leaves)), counts = juv.leaves[, 2])
juv.leaves_plot_rep3 <- data.frame(sample = rep("rep3", nrow(juv.leaves)), counts = juv.leaves[, 3])

juv.leaves_plot <- rbind(juv.leaves_plot_rep1, juv.leaves_plot_rep2, juv.leaves_plot_rep3)

xylem_plot_rep1 <- data.frame(sample = rep("rep1", nrow(xylem)), counts = xylem[, 1])
xylem_plot_rep2 <- data.frame(sample = rep("rep2", nrow(xylem)), counts = xylem[, 2])
xylem_plot_rep3 <- data.frame(sample = rep("rep3", nrow(xylem)), counts = xylem[, 3])

xylem_plot <- rbind(xylem_plot_rep1, xylem_plot_rep2, xylem_plot_rep3)

A <- ggplot(ad.leaves_plot, aes(x = counts, colour = sample)) +
    geom_freqpoly(alpha = 0.5) +
    xlim(0, 300) +
    ylim(0, 5000) +
    ylab("Number of marks") +
    labs(title = "Adult Leaves") +
    geom_vline(aes(xintercept = 10), color = "red", linetype = "dashed", size = 0.4) +
    scale_fill_manual(values = c("red", "blue", "green")) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          legend.title = element_blank(),
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 14),
          legend.text = element_text(size = 14))

B <- ggplot(juv.leaves_plot, aes(x = counts, colour = sample)) +
    geom_freqpoly(alpha = 0.5) +
    xlim(0, 300) +
    ylim(0, 5000) +
    ylab("Number of marks") +
    labs(title = "Juvenile Leaves") +
    geom_vline(aes(xintercept = 10), color = "red", linetype = "dashed", size = 0.4) +
    scale_fill_manual(values = c("red", "blue", "green")) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          legend.title = element_blank(),
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 14),
          legend.text = element_text(size = 14))

C <- ggplot(xylem_plot, aes(x = counts, colour = sample)) +
    geom_freqpoly(alpha = 0.5) +
    xlim(0, 300) +
    ylim(0, 5000) +
    ylab("Number of marks") +
    labs(title = "Xylem") +
    geom_vline(aes(xintercept = 10), color = "red", linetype = "dashed", size = 0.4) +
    scale_fill_manual(values = c("red", "blue", "green")) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          legend.title = element_blank(),
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 14),
          legend.text = element_text(size = 14))

svg(filename = paste(opts$output), width = 16, height = 6, pointsize = 12)
grid.arrange(A, B, C, ncol = 3)
dev.off()

# Correlation plot (pearson)

# mat : is a matrix of data
# ... : further arguments to pass to the native R cor.test function
cor.mtest <- function(mat, ...) {
    mat <- as.matrix(mat)
    n <- ncol(mat)
    p.mat <- matrix(NA, n, n)
    diag(p.mat) <- 0
    for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
            tmp <- cor.test(mat[, i], mat[, j], ...)
            p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
        }
    }
    colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
    p.mat
}

# Load evaluated sites
# Load marks without overlaps.
# Load marks without overlaps.
marks <- read.table(opts$input1)
marks <- unique(marks)

# Load counts file.
counts <- read.table(opts$input2, header = T, sep = "\t", quote = "\"", check.names = F)

# Extracts the counts of the evaluated sites
counts_file <- counts[counts[, 1] %in% marks[, 1], ]
counts_file[] <- lapply(counts_file, function(x) as.numeric(x))

counts_file_ms <- counts_file[, grep("_ms", colnames(counts_file))]
counts_file_hp <- counts_file[, grep("_hp", colnames(counts_file))]

# Loads methylated sites
methylated_sites <- read.table(opts$input3, header = T, sep = "\t", quote = "\"", check.names = F)

# Removes NAs
methylated_sites <- lapply(methylated_sites, function(x) x[!is.na(x)])

# Identify all methylated sites
methylated_sites_names <- as.character()
for (i in 1:length(methylated_sites)){
    methylated_sites_names <- c(methylated_sites_names, as.character(unlist(methylated_sites[i])))
}

methylated_sites_names <- unique(methylated_sites_names)

# Extracts the counts of the methylated sites
meth_counts_file <- counts[counts[, 1] %in% methylated_sites_names, ]
counts_file[] <- lapply(counts_file, function(x) as.numeric(x))

meth_counts_file_ms <- meth_counts_file[, grep("_ms", colnames(meth_counts_file))]
meth_counts_file_hp <- meth_counts_file[, grep("_hp", colnames(meth_counts_file))]

colnames(counts_file_ms) <- c("Adult leaves - rep1",
                              "Juvenile leaves - rep1",
                              "Adult leaves - rep2",
                              "Juvenile leaves - rep2",
                              "Adult leaves - rep3",
                              "Juvenile leaves - rep3",
                              "Xylem - rep1",
                              "Xylem - rep2",
                              "Xylem - rep3")

colnames(counts_file_hp) <- c("Adult leaves - rep1",
                              "Juvenile leaves - rep1",
                              "Adult leaves - rep2",
                              "Juvenile leaves - rep2",
                              "Xylem - rep1",
                              "Xylem - rep2",
                              "Adult leaves - rep3",
                              "Juvenile leaves - rep3",
                              "Xylem - rep3")

colnames(meth_counts_file_ms) <- c("Adult leaves - rep1",
                                   "Juvenile leaves - rep1",
                                   "Adult leaves - rep2",
                                   "Juvenile leaves - rep2",
                                   "Adult leaves - rep3",
                                   "Juvenile leaves - rep3",
                                   "Xylem - rep1",
                                   "Xylem - rep2",
                                   "Xylem - rep3")

colnames(meth_counts_file_hp) <- c("Adult leaves - rep1",
                                   "Juvenile leaves - rep1",
                                   "Adult leaves - rep2",
                                   "Juvenile leaves - rep2",
                                   "Xylem - rep1",
                                   "Xylem - rep2",
                                   "Adult leaves - rep3",
                                   "Juvenile leaves - rep3",
                                   "Xylem - rep3")

# Correlation matrix
correlation_mspI <- cor(counts_file_ms)
p.mat_mspI <- cor.mtest(counts_file_ms)
correlation_hpaII <- cor(counts_file_hp)
p.mat_hpaII <- cor.mtest(counts_file_hp)

meth_correlation_mspI <- cor(meth_counts_file_ms)
meth_p.mat_mspI <- cor.mtest(meth_counts_file_ms)
meth_correlation_hpaII <- cor(meth_counts_file_hp)
meth_p.mat_hpaII <- cor.mtest(meth_counts_file_hp)

# Function to extract the names of the samples (second e third columns in the ID)
extract_names <- function(x){
    new_names <- as.data.frame(matrix(unlist(strsplit(x = x, "_")), nrow = length(x), byrow = T))
    new_names <- paste(new_names[, 2], new_names[, 3])
}

colnames(correlation_mspI) <- colnames(correlation_mspI)
colnames(correlation_hpaII) <- colnames(correlation_hpaII)
colnames(meth_correlation_mspI) <- colnames(meth_correlation_mspI)
colnames(meth_correlation_hpaII) <- colnames(meth_correlation_hpaII)

rownames(correlation_mspI) <- rownames(correlation_mspI)
rownames(correlation_hpaII) <- rownames(correlation_hpaII)
rownames(meth_correlation_mspI) <- rownames(meth_correlation_mspI)
rownames(meth_correlation_hpaII) <- rownames(meth_correlation_hpaII)

# Correlation plots
col <- colorRampPalette(c("#EE3B3B", "#FF6EB4", "#FFFFFF", "#00FFFF", "#6A5ACD"))

svg(filename = paste(opts$correlation_output, "methylated_sites.svg", sep = "_"), width = 16, height = 10, pointsize = 12)
par(mfrow=c(1, 2))
corrplot(meth_correlation_mspI, method = "square", col = col(200),
         type = "upper", order = "hclust",
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 45, #Text label color and rotation
         # Combine with significance
         p.mat = meth_p.mat_mspI, sig.level = 0.01, insig = "blank",
         # hide correlation coefficient on the principal diagonal
         diag = T, title = "a) Methylated sites - Counts of the PstI-MspI library", mar = c(0, 0, 2, 0))

corrplot(meth_correlation_hpaII, method = "square", col = col(200),
         type = "upper", order = "hclust", addrect = 3,
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 45, #Text label color and rotation
         # Combine with significance
         p.mat = meth_p.mat_mspI, sig.level = 0.01, insig = "blank",
         # hide correlation coefficient on the principal diagonal
         diag = T, title = "b) Methylated sites - Counts of the PstI-HpaII library", mar = c(0, 0, 2, 0))
dev.off()

svg(filename = paste(opts$correlation_output, "sampled_sites.svg", sep = "_"), width = 16, height = 10, pointsize = 12)

par(mfrow=c(1, 2))

corrplot(correlation_mspI, method = "square", col = col(200),
         type = "upper", order = "hclust",
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 45, #Text label color and rotation
         # Combine with significance
         p.mat = p.mat_mspI, sig.level = 0.01, insig = "blank",
         # hide correlation coefficient on the principal diagonal
         diag = T, title = "a) Sampled sites - Counts of the PstI-MspI library", mar = c(0, 0, 2, 0))

corrplot(correlation_hpaII, method = "square", col = col(200),
         type = "upper", order = "hclust",
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 45, #Text label color and rotation
         # Combine with significance
         p.mat = p.mat_mspI, sig.level = 0.01, insig = "blank",
         # hide correlation coefficient on the principal diagonal
         diag = T, title = "b) Sampled sites - Counts of the PstI-HpaII library", mar = c(0, 0, 2, 0))
dev.off()

svg(filename = paste(opts$correlation_output, ".svg", sep = ""), width = 12, height = 12, pointsize = 12)
par(mfrow=c(2, 2))
corrplot(correlation_mspI, method = "square", col = col(200),
         type = "upper", order = "hclust",
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 45, #Text label color and rotation
         # Combine with significance
         p.mat = p.mat_mspI, sig.level = 0.01, insig = "blank",
         # hide correlation coefficient on the principal diagonal
         diag = T, title = "a) Sampled sites - Counts of the PstI-MspI library", mar = c(0, 0, 2, 0))

corrplot(correlation_hpaII, method = "square", col = col(200),
         type = "upper", order = "hclust",
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 45, #Text label color and rotation
         # Combine with significance
         p.mat = p.mat_mspI, sig.level = 0.01, insig = "blank",
         # hide correlation coefficient on the principal diagonal
         diag = T, title = "b) Sampled sites - Counts of the PstI-HpaII library", mar = c(0, 0, 2, 0))

corrplot(meth_correlation_mspI, method = "square", col = col(200),
         type = "upper", order = "hclust",
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 45, #Text label color and rotation
         # Combine with significance
         p.mat = meth_p.mat_mspI, sig.level = 0.01, insig = "blank",
         # hide correlation coefficient on the principal diagonal
         diag = T, title = "c) Methylated sites - Counts of the PstI-MspI library", mar = c(0, 0, 2, 0))

corrplot(meth_correlation_hpaII, method = "square", col = col(200),
         type = "upper", order = "hclust", addrect = 3,
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 45, #Text label color and rotation
         # Combine with significance
         p.mat = meth_p.mat_mspI, sig.level = 0.01, insig = "blank",
         # hide correlation coefficient on the principal diagonal
         diag = T, title = "d) Methylated sites - Counts of the PstI-HpaII library", mar = c(0, 0, 2, 0))
dev.off()

### Correlation plot (pearson) ###

# mat : is a matrix of data
# ... : further arguments to pass to the native R cor.test function
cor.mtest <- function(mat, ...) {
    mat <- as.matrix(mat)
    n <- ncol(mat)
    p.mat <- matrix(NA, n, n)
    diag(p.mat) <- 0
    for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
            tmp <- cor.test(mat[, i], mat[, j], ...)
            p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
        }
    }
    colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
    p.mat
}

# Loads evaluated sites
# Loads marks without overlaps.
# Loads marks without overlaps.
marks <- read.table(opts$input1)
marks <- unique(marks)

# Loads counts file.
counts <- read.table(opts$input2, header = T, sep = "\t", quote = "\"", check.names = F)

# Extracts the counts of the evaluated sites
counts_file <- counts[counts[, 1] %in% marks[, 1], ]
counts_file[] <- lapply(counts_file, function(x) as.numeric(x))

counts_file_ms <- counts_file[, grep("_ms", colnames(counts_file))]
counts_file_hp <- counts_file[, grep("_hp", colnames(counts_file))]

# Loads methylated sites
methylated_sites <- read.table(opts$input3, header = T, sep = "\t", quote = "\"", check.names = F)

# Removes NAs
methylated_sites <- lapply(methylated_sites, function(x) x[!is.na(x)])

# Identify all methylated sites
methylated_sites_names <- as.character()
for (i in 1:length(methylated_sites)){
    methylated_sites_names <- c(methylated_sites_names, as.character(unlist(methylated_sites[i])))
}

methylated_sites_names <- unique(methylated_sites_names)

# Extracts the counts of the methylated sites
meth_counts_file <- counts[counts[, 1] %in% methylated_sites_names, ]
counts_file[] <- lapply(counts_file, function(x) as.numeric(x))

meth_counts_file_ms <- meth_counts_file[, grep("_ms", colnames(meth_counts_file))]
meth_counts_file_hp <- meth_counts_file[, grep("_hp", colnames(meth_counts_file))]

colnames(counts_file_ms) <- c("Al - rep1",
                              "Jl - rep1",
                              "Al - rep2",
                              "Jl - rep2",
                              "Al - rep3",
                              "Jl - rep3",
                              "Xy - rep1",
                              "Xy - rep2",
                              "Xy - rep3")

colnames(counts_file_hp) <- c("Al - rep1",
                              "Jl - rep1",
                              "Al - rep2",
                              "Jl - rep2",
                              "Xy - rep1",
                              "Xy - rep2",
                              "Al - rep3",
                              "Jl - rep3",
                              "Xy - rep3")

colnames(meth_counts_file_ms) <- c("Al - rep1",
                                   "Jl - rep1",
                                   "Al - rep2",
                                   "Jl - rep2",
                                   "Al - rep3",
                                   "Jl - rep3",
                                   "Xy - rep1",
                                   "Xy - rep2",
                                   "Xy - rep3")

colnames(meth_counts_file_hp) <- c("Al - rep1",
                                   "Jl - rep1",
                                   "Al - rep2",
                                   "Jl - rep2",
                                   "Xy - rep1",
                                   "Xy - rep2",
                                   "Al - rep3",
                                   "Jl - rep3",
                                   "Xy - rep3")

# Correlation matrix
correlation_mspI <- cor(counts_file_ms)
p.mat_mspI <- cor.mtest(counts_file_ms)
correlation_hpaII <- cor(counts_file_hp)
p.mat_hpaII <- cor.mtest(counts_file_hp)

meth_correlation_mspI <- cor(meth_counts_file_ms)
meth_p.mat_mspI <- cor.mtest(meth_counts_file_ms)
meth_correlation_hpaII <- cor(meth_counts_file_hp)
meth_p.mat_hpaII <- cor.mtest(meth_counts_file_hp)

# Function to extract the names of the samples (second e third columns in the ID)
extract_names <- function(x){
    new_names <- as.data.frame(matrix(unlist(strsplit(x = x, "_")), nrow = length(x), byrow = T))
    new_names <- paste(new_names[, 2], new_names[, 3])
}

colnames(correlation_mspI) <- colnames(correlation_mspI)
colnames(correlation_hpaII) <- colnames(correlation_hpaII)
colnames(meth_correlation_mspI) <- colnames(meth_correlation_mspI)
colnames(meth_correlation_hpaII) <- colnames(meth_correlation_hpaII)

rownames(correlation_mspI) <- rownames(correlation_mspI)
rownames(correlation_hpaII) <- rownames(correlation_hpaII)
rownames(meth_correlation_mspI) <- rownames(meth_correlation_mspI)
rownames(meth_correlation_hpaII) <- rownames(meth_correlation_hpaII)

# Correlation plots
col <- colorRampPalette(c("#EE3B3B", "#FF6EB4", "#FFFFFF", "#00FFFF", "#6A5ACD"))

svg(filename = paste(opts$correlation_output, "methylated_sites_short_name.svg", sep = "_"), width = 16, height = 10, pointsize = 12)
par(mfrow=c(1, 2))
corrplot(meth_correlation_mspI, method = "square", col = col(200),
         type = "upper", order = "hclust",
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 45, #Text label color and rotation
         # Combine with significance
         p.mat = meth_p.mat_mspI, sig.level = 0.01, insig = "blank",
         # hide correlation coefficient on the principal diagonal
         diag = T, title = "a) Methylated sites - Counts of the PstI-MspI library", mar = c(0, 0, 2, 0))

corrplot(meth_correlation_hpaII, method = "square", col = col(200),
         type = "upper", order = "hclust", addrect = 3,
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 45, #Text label color and rotation
         # Combine with significance
         p.mat = meth_p.mat_mspI, sig.level = 0.01, insig = "blank",
         # hide correlation coefficient on the principal diagonal
         diag = T, title = "b) Methylated sites - Counts of the PstI-HpaII library", mar = c(0, 0, 2, 0))
dev.off()

svg(filename = paste(opts$correlation_output, "sampled_sites_short_name.svg", sep = "_"), width = 16, height = 10, pointsize = 12)

par(mfrow=c(1, 2))

corrplot(correlation_mspI, method = "square", col = col(200),
         type = "upper", order = "hclust",
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 45, #Text label color and rotation
         # Combine with significance
         p.mat = p.mat_mspI, sig.level = 0.01, insig = "blank",
         # hide correlation coefficient on the principal diagonal
         diag = T, title = "a) Sampled sites - Counts of the PstI-MspI library", mar = c(0, 0, 2, 0))

corrplot(correlation_hpaII, method = "square", col = col(200),
         type = "upper", order = "hclust",
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 45, #Text label color and rotation
         # Combine with significance
         p.mat = p.mat_mspI, sig.level = 0.01, insig = "blank",
         # hide correlation coefficient on the principal diagonal
         diag = T, title = "b) Sampled sites - Counts of the PstI-HpaII library", mar = c(0, 0, 2, 0))
dev.off()

svg(filename = paste(opts$correlation_output, "_short_name.svg", sep = ""), width = 12, height = 12, pointsize = 12)
par(mfrow=c(2, 2))
corrplot(correlation_mspI, method = "square", col = col(200),
         type = "upper", order = "hclust",
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 45, #Text label color and rotation
         # Combine with significance
         p.mat = p.mat_mspI, sig.level = 0.01, insig = "blank",
         # hide correlation coefficient on the principal diagonal
         diag = T, title = "a) Sampled sites - Counts of the PstI-MspI library", mar = c(0, 0, 2, 0))

corrplot(correlation_hpaII, method = "square", col = col(200),
         type = "upper", order = "hclust",
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 45, #Text label color and rotation
         # Combine with significance
         p.mat = p.mat_mspI, sig.level = 0.01, insig = "blank",
         # hide correlation coefficient on the principal diagonal
         diag = T, title = "b) Sampled sites - Counts of the PstI-HpaII library", mar = c(0, 0, 2, 0))

corrplot(meth_correlation_mspI, method = "square", col = col(200),
         type = "upper", order = "hclust",
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 45, #Text label color and rotation
         # Combine with significance
         p.mat = meth_p.mat_mspI, sig.level = 0.01, insig = "blank",
         # hide correlation coefficient on the principal diagonal
         diag = T, title = "c) Methylated sites - Counts of the PstI-MspI library", mar = c(0, 0, 2, 0))

corrplot(meth_correlation_hpaII, method = "square", col = col(200),
         type = "upper", order = "hclust", addrect = 3,
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 45, #Text label color and rotation
         # Combine with significance
         p.mat = meth_p.mat_mspI, sig.level = 0.01, insig = "blank",
         # hide correlation coefficient on the principal diagonal
         diag = T, title = "d) Methylated sites - Counts of the PstI-HpaII library", mar = c(0, 0, 2, 0))
dev.off()
