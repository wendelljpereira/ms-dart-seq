"Usage: id_part1.R (--out1 <O1>) <input1> <input2> <input3> <input4> <input5> <input6> <input7> <input8>
-h --help    show this
--out1   name1    specify the name for the first output file
id_part1.R -h | --help  show this message
" -> doc

# load the docopt library
require(docopt)
require(tidyverse)
require(plyr)
require(VennDiagram)
require(gridExtra)
require(gdata)
# retrieve the command-line arguments
opts <- docopt(doc)

outside_genes_or_tes <- unique(read.table(opts$`<input1>`, colClasses = "character")[4])
transposons <- unique(read.table(opts$`<input2>`, colClasses = "character")[4])
exons <- unique(read.table(opts$`<input3>`, colClasses = "character")[4])
exons_containing_transposons <- unique(read.table(opts$`<input4>`, colClasses = "character")[4])
exon_overlaping_genes <- unique(read.table(opts$`<input5>`, colClasses = "character")[4])
exon_overlaping_genes <- unique(unique(exon_overlaping_genes))
te_within_intron_utr <- unique(read.table(opts$`<input6>`, colClasses = "character")[4])
intron_or_utr <- unique(read.table(opts$`<input7>`, colClasses = "character")[4])

# make the data.frame to make the layers of the plot
dat <- data.frame(count = c(nrow(exons) + nrow(intron_or_utr),
                            nrow(transposons) + nrow(te_within_intron_utr),
                            nrow(outside_genes_or_tes),
                            nrow(exons_containing_transposons) + nrow(exon_overlaping_genes),
                            nrow(intron_or_utr),
                            nrow(exons),
                            nrow(transposons),
                            nrow(te_within_intron_utr),
                            nrow(outside_genes_or_tes),
                            nrow(exon_overlaping_genes),
                            nrow(exons_containing_transposons)),
                  ring = c("x", "x", "x", "x", "y", "y", "y", "y", "y", "y", "y"),
                  category = c("Gene", "TEs", "Intergenic", "Unknow", "Intron or UTR", "Exon", "TE - intergenic", "TE inside a intron or UTR", "Intergenic and outside of the TEs", "Overlaps of genes", "Large overlaps of exons with TEs"))

# compute fractions
dat %<>% group_by(ring) %>% mutate(fraction = count / sum(count),
                                   ymax = cumsum(fraction),
                                   ymin = c(0, ymax[1:length(ymax) - 1]))

# Add x limits
baseNum <- 3
dat$xmax <- as.numeric(dat$ring) + baseNum
dat$xmin <- dat$xmax - 1

# write a table with the counts in each category
write.table(dat,
            "Mehylated_sites_distribution_graph_table.tst",
            col.names = T,
            row.names = F,
            sep = "\t",
            quote = F)

color_values <- c("#a672c9",
                  "#724193",
                  "#FFF68F",
                  "#FFF68F",
                  "#c255fc",
                  "#75dd81",
                  "#00ff1d",
                  "#1cb9ee",
                  "#7db9f2",
                  "#399cf9",
                  "#01aa15")

break_labels <- c("Gene",
                  "TEs",
                  "Intergenic",
                  "Unknow",
                  "Intron or UTR",
                  "Exon",
                  "TE - intergenic",
                  "TE inside a intron or UTR",
                  "Intergenic and outside of the TEs",
                  "Overlaps of genes",
                  "Large overlaps of exons with TEs")

# plot
chart <- ggplot(dat, aes(fill = category,
                         ymax = ymax,
                         ymin = ymin,
                         xmax = xmax,
                         xmin = xmin)) +
    geom_rect(colour = "grey30") +
    coord_polar(theta = "y") +
    xlim(c(0, 6)) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    theme(axis.text = element_blank(), axis.title = element_blank()) +
    theme(axis.ticks = element_blank(), panel.border = element_blank()) +
    theme(legend.text = element_text(size = 12), legend.title = element_blank()) +
    labs(title = "Genomic context of methylated sites") +
    scale_fill_manual(values = color_values,
                      breaks = break_labels) +
    geom_label(aes(label = paste(count, "\n", paste0("(", round(fraction * 100, 1), "%)")),
                   x = xmax,
                   y = ymax,
                   size = 11),
               inherit.aes = T,
               show.legend = F)

svg(filename = paste(opts$O1), width = 8, height = 8, pointsize = 12)
chart
dev.off()

# Reads the files with samples data
metil <- read.table(opts$`<input8>`, colClasses = "character", header = T)
colnames(metil) <- c("BRASUZ1 Adult Leaf", "BRASUZ1 Juvenile Leaf", "BRASUZ1 Xylem")

for (i in 1:length(metil)){

    tissue <- colnames(metil)[i]
    marks <- unique(metil[, i])

    fora_de_genes_ou_transposons_sub <- outside_genes_or_tes[outside_genes_or_tes[, 1] %in% as.character(marks), ]
    length(unique(fora_de_genes_ou_transposons_sub))

    transposons_sub <- transposons[transposons[, 1] %in% as.character(marks), ]
    length(unique(transposons_sub))

    exons_sub <- exons[exons[, 1] %in% as.character(marks), ]
    length(unique(exons_sub))

    exons_com_transposons_sub <- exons_containing_transposons[exons_containing_transposons[, 1] %in% as.character(marks), ]
    length(unique(exons_com_transposons_sub))

    exons_genes_duplos_sub <- exon_overlaping_genes[exon_overlaping_genes[, 1] %in% as.character(marks), ]
    length(unique(exons_genes_duplos_sub))

    nao_codantes_com_transposons_sub <- te_within_intron_utr[te_within_intron_utr[, 1] %in% as.character(marks), ]
    length(unique(nao_codantes_com_transposons_sub))

    nao_codantes_sem_transposons_sub <- intron_or_utr[intron_or_utr[, 1] %in% as.character(marks), ]
    length(unique(nao_codantes_sem_transposons_sub))

    category_labels <- c("Gene",
                         "TEs",
                         "Intergenic",
                         "Unknow",
                         "Intron or UTR",
                         "Exon",
                         "TE - intergenic",
                         "TE inside a intron or UTR",
                         "Intergenic and outside of the TEs",
                         "Overlaps of genes",
                         "Large overlaps of exons with TEs")

    ring <- ring_labels <- c("x", "x", "x", "x", "y", "y", "y", "y", "y", "y", "y")

    # make the data.frame to make the layers of the plot
    dat_sub <- data.frame(count = c(length(exons_sub) + length(nao_codantes_sem_transposons_sub),
                                    length(transposons_sub) + length(nao_codantes_com_transposons_sub),
                                    length(fora_de_genes_ou_transposons_sub),
                                    length(exons_com_transposons_sub) + length(exons_genes_duplos_sub),
                                    length(nao_codantes_sem_transposons_sub),
                                    length(exons_sub),
                                    length(transposons_sub),
                                    length(nao_codantes_com_transposons_sub),
                                    length(fora_de_genes_ou_transposons_sub),
                                    length(exons_genes_duplos_sub),
                                    length(exons_com_transposons_sub)),
                          ring = ring_labels,
                          category = category_labels)

    # compute fractions
    dat_sub %<>% group_by(ring) %>% mutate(fraction = count / sum(count),
                                           ymax = cumsum(fraction),
                                           ymin = c(0, ymax[1:length(ymax) - 1]))

    # Add x limits
    baseNum <- 3
    dat_sub$xmax <- as.numeric(dat_sub$ring) + baseNum
    dat_sub$xmin <- dat_sub$xmax - 1

    # write a table with the counts in each category
    write.table(dat_sub,
                paste(tissue, "Mehylated_sites_distribution_graph_table.tst", sep = "_"),
                col.names = T,
                row.names = F,
                sep = "\t",
                quote = F)

    # plot
    chart_sub <- ggplot(dat_sub, aes(fill = category,
                                     ymax = ymax,
                                     ymin = ymin,
                                     xmax = xmax,
                                     xmin = xmin)) +
        geom_rect(colour = "grey30") +
        coord_polar(theta = "y") +
        xlim(c(0, 6)) +
        theme_bw() +
        theme(panel.grid = element_blank()) +
        theme(axis.text = element_blank(), axis.title = element_blank()) +
        theme(axis.ticks = element_blank(), panel.border = element_blank()) +
        theme(legend.text = element_text(size = 12), legend.title = element_blank()) +
        labs(title="Genomic context of methylated sites") +
        scale_fill_manual(values = color_values,
                          breaks = break_labels) +
        geom_label(aes(label = paste(count, "\n", paste0("(", round(fraction*100, 1), "%)")),
                       x = xmax,
                       y = ymax,
                       size = 11),
                   inherit.aes = T,
                   show.legend = F)

    svg(paste(paste0("images/genomic_context/", tissue), "Mehylated_sites_distribution.svg", sep = "_"),
        width = 8,
        height = 8,
        pointsize = 12)

    print(chart_sub)

    dev.off()
}

########################################################
### Genomic context of the marks in the intersection ###
########################################################

tissue_methylation <- read.table(opts$`<input8>`,
                                 colClasses = "character",
                                 header = T)

sample1 <- as.character(tissue_methylation[, 1])
sample2 <- as.character(tissue_methylation[, 2])
sample3 <- as.character(tissue_methylation[, 3])

sample1 <- unique(sample1[complete.cases(sample1)])
sample2 <- unique(sample2[complete.cases(sample2)])
sample3 <- unique(sample3[complete.cases(sample3)])

inter_12_full <- unique(intersect(sample1, sample2))
inter_13_full <- unique(intersect(sample1, sample3))
inter_23_full <- unique(intersect(sample2, sample3))

inter_123_full <- unique(intersect(inter_12_full, sample3))

marks <- unique(inter_123_full)

fora_de_genes_ou_transposons_sub <- outside_genes_or_tes[outside_genes_or_tes[, 1] %in% as.character(marks), ]
length(unique(fora_de_genes_ou_transposons_sub))

transposons_sub <- transposons[transposons[, 1] %in% as.character(marks), ]
length(unique(transposons_sub))

exons_sub <- exons[exons[, 1] %in% as.character(marks), ]
length(unique(exons_sub))

exons_com_transposons_sub <- exons_containing_transposons[exons_containing_transposons[, 1] %in% as.character(marks), ]
length(unique(exons_com_transposons_sub))

exons_genes_duplos_sub <- exon_overlaping_genes[exon_overlaping_genes[, 1] %in% as.character(marks), ]
length(unique(exons_genes_duplos_sub))

nao_codantes_com_transposons_sub <- te_within_intron_utr[te_within_intron_utr[, 1] %in% as.character(marks), ]
length(unique(nao_codantes_com_transposons_sub))

nao_codantes_sem_transposons_sub <- intron_or_utr[intron_or_utr[, 1] %in% as.character(marks), ]
length(unique(nao_codantes_sem_transposons_sub))

category_labels <- c("Gene",
                     "TEs",
                     "Intergenic",
                     "Unknow",
                     "Intron or UTR",
                     "Exon",
                     "TE - intergenic",
                     "TE inside a intron or UTR",
                     "Intergenic and outside of the TEs",
                     "Overlaps of genes",
                     "Large overlaps of exons with TEs")

ring <- ring_labels <- c("x", "x", "x", "x", "y", "y", "y", "y", "y", "y", "y")

# make the data.frame to make the layers of the plot
dat_sub <- data.frame(count = c(length(exons_sub) + length(nao_codantes_sem_transposons_sub),
                                length(transposons_sub) + length(nao_codantes_com_transposons_sub),
                                length(fora_de_genes_ou_transposons_sub),
                                length(exons_com_transposons_sub) + length(exons_genes_duplos_sub),
                                length(nao_codantes_sem_transposons_sub),
                                length(exons_sub),
                                length(transposons_sub),
                                length(nao_codantes_com_transposons_sub),
                                length(fora_de_genes_ou_transposons_sub),
                                length(exons_genes_duplos_sub),
                                length(exons_com_transposons_sub)),
                      ring = ring_labels,
                      category = category_labels)

# compute fractions
dat_sub %<>% group_by(ring) %>% mutate(fraction = count / sum(count),
                                       ymax = cumsum(fraction),
                                       ymin = c(0, ymax[1:length(ymax) - 1]))

# Add x limits
baseNum <- 3
dat_sub$xmax <- as.numeric(dat_sub$ring) + baseNum
dat_sub$xmin <- dat_sub$xmax - 1

# write a table with the counts in each category
write.table(dat_sub,
            paste("intersection", "Mehylated_sites_distribution_graph_table.tst", sep = "_"),
            col.names = T,
            row.names = F,
            sep = "\t",
            quote = F)

# plot
chart_sub <- ggplot(dat_sub, aes(fill = category,
                                 ymax = ymax,
                                 ymin = ymin,
                                 xmax = xmax,
                                 xmin = xmin)) +
    geom_rect(colour = "grey30") +
    coord_polar(theta = "y") +
    xlim(c(0, 6)) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    theme(axis.text = element_blank(), axis.title = element_blank()) +
    theme(axis.ticks = element_blank(), panel.border = element_blank()) +
    theme(legend.text = element_text(size = 12), legend.title = element_blank()) +
    labs(title="Genomic context of methylated sites") +
    scale_fill_manual(values = color_values,
                      breaks = break_labels) +
    geom_label(aes(label = paste(count, "\n", paste0("(", round(fraction*100, 1), "%)")),
                   x = xmax,
                   y = ymax,
                   size = 11),
               inherit.aes = T,
               show.legend = F)

svg(paste(paste0("images/genomic_context/", "intersection"), "Mehylated_sites_distribution.svg", sep = "_"),
    width = 8,
    height = 8,
    pointsize = 12)

print(chart_sub)

dev.off()

###############################################################
### Comparing the subsets of DNA methylation among tissues. ###
###############################################################

## venn plot function
venn_3_samples <- function(sample1 = sample1, sample2 = sample2, sample3 = sample3, name1 = name1, name2 = name2, name3 = name3, clone_name = clone_name, save_ids = "FALSE"){

    sample1 <- unique(sample1[complete.cases(sample1)])
    sample2 <- unique(sample2[complete.cases(sample2)])
    sample3 <- unique(sample3[complete.cases(sample3)])

    inter_12_full <- unique(intersect(sample1, sample2))
    inter_13_full <- unique(intersect(sample1, sample3))
    inter_23_full <- unique(intersect(sample2, sample3))

    inter_123_full <- unique(intersect(inter_12_full, sample3))
    inter_12 <- unique(inter_12_full[!inter_12_full %in% inter_123_full])
    inter_13 <- unique(inter_13_full[!inter_13_full %in% inter_123_full])
    inter_23 <- unique(inter_23_full[!inter_23_full %in% inter_123_full])
    unic_s1 <- unique(sample1[!sample1 %in% c(inter_12, inter_13, inter_123_full)])
    unic_s2 <- unique(sample2[!sample2 %in% c(inter_12, inter_23, inter_123_full)])
    unic_s3 <- unique(sample3[!sample3 %in% c(inter_13, inter_23, inter_123_full)])

    if(save_ids == "TRUE"){
        write.table(inter_123_full, paste(clone_name, name1, "vs", name2, "vs", name3, "intersection_marks.txt", sep = "_"), row.names = F, col.names = F, quote = F)
        write.table(inter_12, paste(clone_name, "only", name1, "vs", name2, "intersection_marks.txt", sep = "_"), row.names = F, col.names = F, quote = F)
        write.table(inter_13, paste(clone_name, "only", name1, "vs", name3, "intersection_marks.txt", sep = "_"), row.names = F, col.names = F, quote = F)
        write.table(inter_23, paste(clone_name, "only", name2, "vs", name3, "intersection_marks.txt", sep = "_"), row.names = F, col.names = F, quote = F)
        write.table(inter_12_full, paste(clone_name, "all", name1, "vs", name2, "intersection_marks.txt", sep = "_"), row.names = F, col.names = F, quote = F)
        write.table(inter_13_full, paste(clone_name, "all", name1, "vs", name3, "intersection_marks.txt", sep = "_"), row.names = F, col.names = F, quote = F)
        write.table(inter_23_full, paste(clone_name, "all", name2, "vs", name3, "intersection_marks.txt", sep = "_"), row.names = F, col.names = F, quote = F)
        write.table(unic_s1, paste(clone_name, name1, "exclusive_marks.txt", sep = "_"), row.names = F, col.names = F, quote = F)
        write.table(unic_s2, paste(clone_name, name2, "exclusive_marks.txt", sep = "_"), row.names = F, col.names = F, quote = F)
        write.table(unic_s3, paste(clone_name, name3, "exclusive_marks.txt", sep = "_"), row.names = F, col.names = F, quote = F)
    }else if(save_ids == "FALSE"){
    }else{
        print("Invalid save_ids option!")
    }

    grid.newpage();
    venn.plot <- draw.triple.venn(
        area1 = length(sample1),
        area2 = length(sample2),
        area3 = length(sample3),
        n12 = length(inter_12_full),
        n13 = length(inter_13_full),
        n23 = length(inter_23_full),
        n123 = length(inter_123_full),
        cross.area = length(inter),
        alpha = 0.75,
        category = c(deparse(name1), deparse(name2), deparse(name3)),
        fill = c("darkgreen", "green", "#8B4513"),
        lty = "blank",
        cex = 4,
        cat.cex = 2.5,
        cat.dist = 0.055,
        ext.pos = 0,
        ext.dist = -0.05,
        ext.length = 0.85,
        ext.line.lwd = 2,
        ext.line.lty = "dashed",
        scaled = T,
        print.mode = c("raw", "percent"),
        rotation.degree = 0)
}

# Reads the methylation in each tissue
tissue_methylation <- read.table(opts$`<input8>`, colClasses = "character", header = T)

names <- c("BRASUZ1 Adult Leaf", "BRASUZ1 Juvenile Leaf", "BRASUZ1 Xylem")
colnames(tissue_methylation) <- names

dfList <- list(outside_genes_or_tes,
               exons,
               transposons,
               te_within_intron_utr,
               intron_or_utr)

dfnames <- c("intergenic",
             "exons",
             "transposons_within_intergenic_regions",
             "transposons_within_intron_or_utr",
             "within_intron_utrs_of_genes")

for(i in 1:length(dfList)){

    AL <- tissue_methylation[tissue_methylation[, 1] %in% dfList[[i]][, 1], ][1]
    JL <- tissue_methylation[tissue_methylation[, 2] %in% dfList[[i]][, 1], ][2]
    XY <- tissue_methylation[tissue_methylation[, 3] %in% dfList[[i]][, 1], ][3]

    A <- venn_3_samples(sample1 = AL[,1],
                        sample2 = JL[,1],
                        sample3 = XY[,1],
                        name1 = names[1],
                        name2 = names[2],
                        name3 = names[3],
                        clone_name = "")

    svg(filename = paste("images/genomic_context/venn_plot",
                         paste(dfnames[i]), ".svg", sep = "_"),
        width = 12,
        height = 12,
        pointsize = 12)

    grid.arrange(grobTree(A),
                 ncol = 1,
                 top = textGrob(dfnames[i], gp = gpar(fontsize = 30, font = 8)))
    dev.off()
}
