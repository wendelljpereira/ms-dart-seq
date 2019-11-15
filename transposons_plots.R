require(ggplot2)
require(plyr)
require(VennDiagram)
require(gridExtra)
require(gdata)
require(scales)
require(splitstackshape)

# Help function to create a venn plot
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

    # write files with the Ids of the sites in each subset of the venn plots
    if(save_ids == "TRUE"){

        write.table(inter_123_full,
             paste(clone_name, name1, "vs", name2, "vs", name3, "intersection_transposons.txt", sep = "_"), row.names = F, col.names = F, quote = F)
        write.table(inter_12,
             paste(clone_name, "only", name1, "vs", name2, "intersection_transposons.txt", sep = "_"),
             row.names = F, col.names = F, quote = F)
        write.table(inter_13,
             paste(clone_name, "only", name1, "vs", name3, "intersection_transposons.txt", sep = "_"),
             row.names = F, col.names = F, quote = F)
        write.table(inter_23,
             paste(clone_name, "only", name2, "vs", name3, "intersection_transposons.txt", sep = "_"),
              row.names = F, col.names = F, quote = F)
        write.table(inter_12_full,
             paste(clone_name, "all", name1, "vs", name2, "intersection_transposons.txt", sep = "_"),
              row.names = F, col.names = F, quote = F)
        write.table(inter_13_full,
             paste(clone_name, "all", name1, "vs", name3, "intersection_transposons.txt", sep = "_"),
              row.names = F, col.names = F, quote = F)
        write.table(inter_23_full,
             paste(clone_name, "all", name2, "vs", name3, "intersection_transposons.txt", sep = "_"),
              row.names = F, col.names = F, quote = F)
        write.table(unic_s1,
             paste(clone_name, name1, "exclusive_transposons.txt", sep = "_"),
              row.names = F, col.names = F, quote = F)
        write.table(unic_s2,
             paste(clone_name, name2, "exclusive_transposons.txt", sep = "_"),
              row.names = F, col.names = F, quote = F)
        write.table(unic_s3,
             paste(clone_name, name3, "exclusive_transposons.txt", sep = "_"),
              row.names = F, col.names = F, quote = F)

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
        alpha = 0.5,
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

transposons <- read.table(snakemake@input[[1]], sep = "\t")
transposon_in_genes <- read.table(snakemake@input[[2]], sep = "\t")

transposons <- transposons[, c(4, 10)]
transposon_in_genes <- transposon_in_genes[, c(4, 10)]

marks_in_transp <- c(unique(as.character(transposons$V4)),
 unique(as.character(transposon_in_genes$V4)))

# Reads the MSD-Methylated sites
DM_marks <- read.table(snakemake@input[[3]],
                       header = T,
                       na.strings = "NA",
                       colClasses = "character")


# Defines the Methylated TEs of each sample
transposons_all <- data.frame()
for(i in 1:length(names(DM_marks))){
    if(i == 1){

        query <- DM_marks[DM_marks[, i] %in% marks_in_transp, i]

        transposons_group <- transposons[transposons$V4 %in% query,][2]
        colnames(transposons_group) <- paste(names(DM_marks)[i])
        transposons_all <- cbind(transposons_group)

    }else{

        query <- DM_marks[DM_marks[,i] %in% marks_in_transp, i]

        transposons_group <- transposons[transposons$V4 %in% query,][2]
        colnames(transposons_group) <- paste(names(DM_marks)[i])
        transposons_all <- cbindX(transposons_all, transposons_group)
    }
}

# Writes a file with the TEs of each sample
write.table(transposons_all,
            snakemake@output[[1]],
            sep = "\t",
            col.names = T,
            row.names = F,
            quote = F)

# Venn plots
A <- venn_3_samples(sample1 = transposons_all[,1],
                    sample2 = transposons_all[,2],
                    sample3 = transposons_all[,3],
                    name1 = colnames(transposons_all)[1],
                    name2 = colnames(transposons_all)[2],
                    name3 = colnames(transposons_all)[3],
                    clone_name = "BRASUZ1",
                    save_ids = "TRUE")
dev.off()

## Export as a svg file
svg(filename = snakemake@output[[2]], width = 12, height = 12, pointsize = 12)
grid.arrange(grobTree(A), ncol = 1, top = textGrob("BRASUZ1", gp = gpar(fontsize = 30, font = 8)))
dev.off()

# Select the subset of TEs of the venn plots to annotation
sample1 = transposons_all[, 1]
sample2 = transposons_all[, 2]
sample3 = transposons_all[, 3]

sample1 <- unique(as.character(sample1[complete.cases(sample1)]))
sample2 <- unique(as.character(sample2[complete.cases(sample2)]))
sample3 <- unique(as.character(sample3[complete.cases(sample3)]))

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

venn_subset_transp <- cbindX(as.data.frame(inter_123_full),
                             as.data.frame(unic_s1),
                             as.data.frame(unic_s2),
                             as.data.frame(unic_s3))

transp_genome <- read.table(snakemake@input[[4]], sep = "\t")[,4]

transposons_all <- cbindX(as.data.frame(sample1),
                          as.data.frame(sample2),
                          as.data.frame(sample3),
                          venn_subset_transp)

transposons_all <- cbindX(transposons_all, as.data.frame(transp_genome))

colnames(transposons_all) <- c("Adult leaves",
                               "Juvenile leaves",
                               "Xylem",
                               "Intersection",
                               "Unique in Adult",
                               "Unique in Juvenile",
                               "Unique in Xylem",
                               "Genome")

# Verify the TEs that are not possible to classify correctly.
TE_not_class <- read.table(snakemake@input[[4]], sep = "\t")
TE_not_class <- cSplit(indt = TE_not_class, splitCols = "V4", sep = ",", drop = F)

TE_not_class <- TE_not_class[is.na(TE_not_class$V4_02) == FALSE, ]
TE_not_class <- as.character(TE_not_class$V4)

# Defines the groups to be used in the plot
use_intersections <- snakemake@params[["use_intersections"]]

if(use_intersections == TRUE){

    transp_class_all <- data.frame()
    for(i in 1:length(names(transposons_all))){

        transp_subset <- unique(as.character(transposons_all[, i]))

        # check TEs with overlaps
        transp_unclass <- unique(as.character(transp_subset[transp_subset %in% TE_not_class]))

        # check if in the regions with more than one TE, the TEs are of differently classified
        transp_unclass_final <- character()
        for( t in 1:length(transp_unclass)){

            sum_of_class <- (length(grep("DMX", transp_unclass[t])) +
                                 length(grep("DTX", transp_unclass[t])) +
                                 length(grep("DHX", transp_unclass[t])) +
                                 length(grep("DXX-MtTE", transp_unclass[t])) +
                                 length(grep("DXX_Blc", transp_unclass[t])) +
                                 length(grep("RYX", transp_unclass[t])) +
                                 length(grep("RtX", transp_unclass[t])) +
                                 length(grep("RLX", transp_unclass[t])) +
                                 length(grep("RXX-LARD", transp_unclass[t])) +
                                 length(grep("RXX-TRtM", transp_unclass[t])) +
                                 length(grep("RSX", transp_unclass[t])) +
                                 length(grep("RXX_Blc", transp_unclass[t])))

            if(sum_of_class > 1){

                transp_unclass_final <- c(transp_unclass_final,  transp_unclass[t])

            }
        }

        unknow_TE <- length(transp_unclass_final)

        transp_subset <- transp_subset[!transp_subset %in% transp_unclass_final]

        #MAVERICK | DMX
        maverick <- length(grep("DMX", transp_subset))
        #CACTA | DTX
        cacta <- length(grep("DTX", transp_subset))
        #HELITRON | DHX
        helitron <- length(grep("DHX", transp_subset))
        #MITE | DXX-MITE
        mite <- length(grep("DXX-MITE", transp_subset))
        #DNA_general | DXX_Blc
        dna_general <- length(grep("DXX_Blc", transp_subset))
        #DIRS/VIPER | RYX
        dirs_viper <- length(grep("RYX", transp_subset))
        #LINE | RIX
        line <- length(grep("RIX", transp_subset))
        #LTR | RLX
        ltr <- length(grep("RLX", transp_subset))
        #LARD | RXX-LARD
        lard <- length(grep("RXX-LARD", transp_subset))
        #TRIM | RXX-TRIM
        trim <- length(grep("RXX-TRIM", transp_subset))
        #SINE | RSX
        sine <- length(grep("RSX", transp_subset))
        #general | RXX_Blc
        rna_general <- length(grep("RXX_Blc", transp_subset))

        category <- c("MAVERICK",
                      "CACTA",
                      "HELITRON",
                      "MITE",
                      "DNA_general",
                      "DIRS/VIPER",
                      "LINE",
                      "LTR",
                      "LARD",
                      "TRIM",
                      "SINE",
                      "RNA_general",
                      "Nested TEs")

        len <- c(maverick,
                 cacta,
                 helitron,
                 mite,
                 dna_general,
                 dirs_viper,
                 line,ltr,
                 lard,
                 trim,
                 sine,
                 rna_general,
                 unknow_TE)

        transp_class <- data.frame(rep(names(transposons_all)[i], length(category)), category, len)
        transp_class_all <- rbind(transp_class_all ,transp_class)
    }

    colnames(transp_class_all) <- c("samples", "classif", "quantif")

    transp_class_all <- transp_class_all[!transp_class_all$quantif == "0",]

    group_plot <-character()
    for(i in 1:nrow(transp_class_all)){
        if(paste(transp_class_all$samples[i]) == "Genome"){
            group_plot <- append(group_plot,"E. grandis")
        }else{
            group_plot <- append(group_plot,"Samples")
        }
    }

    transp_class_all <- cbind(group_plot, transp_class_all)

    write.table(transp_class_all, snakemake@output[[3]], quote = F, row.names = F, col.names = F, sep = "\t")

    ## Defines the colors of the bars.
    fill <- c("#000000",
              "#E69F00",
              "#56B4E9",
              "#009E73",
              "#F0E442",
              "#0072B2",
              "#D55E00",
              "#CC79A7",
              "#2C7417",
              "#CBE75F",
              "#A877EA",
              "#6BF4ED",
              "#AC6AC9")

    transp_plot <- ggplot(transp_class_all, aes(x = samples, y = quantif, fill = classif)) +
        geom_bar(stat = "identity", alpha = 0.9, width = 0.5) +
        scale_x_discrete(name="")+
        scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
        facet_wrap( ~ group_plot, ncol=2, scales="free")+
        ylab("Number of transposons")+
        theme_bw()+
        theme(legend.text=element_text(size=14), axis.title.y=element_text(size = 20, vjust=2), axis.title.x=element_text(size=14,vjust=0), axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15), legend.title = element_blank(), strip.text.x = element_text(size = 15))+
        scale_fill_manual(values=fill)

    transp_plot_perc <- ggplot(transp_class_all, aes(x = samples, y = quantif, fill = classif)) +
        geom_bar(position = "fill", stat = "identity", alpha = 0.9, width = 0.5) +
        scale_x_discrete(name="")+
        scale_y_continuous(labels = percent_format())+
        facet_wrap( ~ group_plot, ncol=2, scales="free")+
        ylab("Number of transposons")+
        theme_bw()+
        theme(legend.text=element_text(size=14), axis.title.y=element_text(size = 20, vjust=2), axis.title.x=element_text(size=14,vjust=0), axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15), legend.title = element_blank(), strip.text.x = element_text(size = 15))+
        scale_fill_manual(values=fill)

    # Get the ggplot grob
    transp_plot_grob = ggplotGrob(transp_plot)
    transp_plot_perc_grob = ggplotGrob(transp_plot_perc)

    # Builds and save the bar plot.
    svg(snakemake@output[[4]], width=10, height=7)
    grid.newpage()
    grid.draw(transp_plot_grob)
    dev.off()

    # Builds and save the bar plot.
    svg(snakemake@output[[5]], width=10, height=7)
    grid.newpage()
    grid.draw(transp_plot_perc_grob)
    dev.off()

}else if(use_intersections == FALSE){

    transposons_all <- transposons_all[, colnames(transposons_all) %in% c("Adult leaves", "Juvenile leaves", "Xylem", "Genome") ]

    transp_class_all <- data.frame()
    for(i in 1:length(names(transposons_all))){

        transp_subset <- unique(as.character(transposons_all[, i]))

        # check TEs with overlaps
        transp_unclass <- unique(as.character(transp_subset[transp_subset %in% TE_not_class]))

        transp_unclass_final <- character()
        for( t in 1:length(transp_unclass)){

            sum_of_class <- (length(grep("DMX", transp_unclass[t])) +
                                 length(grep("DTX", transp_unclass[t])) +
                                 length(grep("DHX", transp_unclass[t])) +
                                 length(grep("DXX-MtTE", transp_unclass[t])) +
                                 length(grep("DXX_Blc", transp_unclass[t])) +
                                 length(grep("RYX", transp_unclass[t])) +
                                 length(grep("RtX", transp_unclass[t])) +
                                 length(grep("RLX", transp_unclass[t])) +
                                 length(grep("RXX-LARD", transp_unclass[t])) +
                                 length(grep("RXX-TRtM", transp_unclass[t])) +
                                 length(grep("RSX", transp_unclass[t])) +
                                 length(grep("RXX_Blc", transp_unclass[t])))

            if(sum_of_class > 1){
                transp_unclass_final <- c(transp_unclass_final,  transp_unclass[t])
            }
        }

        unknow_TE <- length(transp_unclass_final)

        transp_subset <- transp_subset[!transp_subset %in% transp_unclass_final]

        #MAVERICK | DMX
        maverick <- length(grep("DMX", transp_subset))
        #CACTA | DTX
        cacta <- length(grep("DTX", transp_subset))
        #HELITRON | DHX
        helitron <- length(grep("DHX", transp_subset))
        #MITE | DXX-MITE
        mite <- length(grep("DXX-MITE", transp_subset))
        #DNA_general | DXX_Blc
        dna_general <- length(grep("DXX_Blc", transp_subset))
        #DIRS/VIPER | RYX
        dirs_viper <- length(grep("RYX", transp_subset))
        #LINE | RIX
        line <- length(grep("RIX", transp_subset))
        #LTR | RLX
        ltr <- length(grep("RLX", transp_subset))
        #LARD | RXX-LARD
        lard <- length(grep("RXX-LARD", transp_subset))
        #TRIM | RXX-TRIM
        trim <- length(grep("RXX-TRIM", transp_subset))
        #SINE | RSX
        sine <- length(grep("RSX", transp_subset))
        #general | RXX_Blc
        rna_general <- length(grep("RXX_Blc", transp_subset))

        category <- c("MAVERICK",
                      "CACTA",
                      "HELITRON",
                      "MITE",
                      "DNA_general",
                      "DIRS/VIPER",
                      "LINE",
                      "LTR",
                      "LARD",
                      "TRIM",
                      "SINE",
                      "RNA_general",
                      "Nested TEs")

        len <- c(maverick,
                 cacta,
                 helitron,
                 mite,
                 dna_general,
                 dirs_viper,
                 line,ltr,
                 lard,
                 trim,
                 sine,
                 rna_general,
                 unknow_TE)

        transp_class <- data.frame(rep(names(transposons_all)[i], length(category)), category, len)
        transp_class_all <- rbind(transp_class_all ,transp_class)
    }

    colnames(transp_class_all) <- c("samples", "classif", "quantif")

    transp_class_all <- transp_class_all[!transp_class_all$quantif == "0",]

    group_plot <-character()
    for(i in 1:nrow(transp_class_all)){
        if(paste(transp_class_all$samples[i]) == "Genome"){
            group_plot <- append(group_plot, "E. grandis")
        }else{
            group_plot <- append(group_plot, "Samples")
        }
    }

    transp_class_all <- cbind(group_plot, transp_class_all)

    write.table(transp_class_all, snakemake@output[[3]], quote = F, row.names = F, col.names = F, sep = "\t")

    ## Defines the colors of the bars.
    fill <- c("#000000",
              "#E69F00",
              "#56B4E9",
              "#009E73",
              "#F0E442",
              "#0072B2",
              "#D55E00",
              "#CC79A7",
              "#2C7417",
              "#CBE75F",
              "#A877EA",
              "#6BF4ED",
              "#AC6AC9")

    # labels of the grids
    levels(transp_class_all$group_plot) <- c("E. grandis" = "E. grandis - all TEs", "Samples" = "Samples - Methylated TEs")

    transp_plot <- ggplot(transp_class_all, aes(x = samples, y = quantif, fill = classif)) +
        geom_bar(stat = "identity", alpha = 0.9, width = 0.4) +
        scale_x_discrete(name="")+
        scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
        facet_wrap( ~ group_plot, ncol = 2, scales = "free") +
        ylab("Number of transposons")+
        theme_bw()+
        theme(legend.text=element_text(size=14),
              axis.title.y = element_text(size = 20, vjust = 2),
              axis.title.x=element_text(size = 14, vjust = 0),
              axis.text.x = element_text(size = 15),
              axis.text.y = element_text(size = 15),
              legend.title = element_blank(),
              strip.text.x = element_text(size = 15))+
        scale_fill_manual(values = fill)

    # Get the ggplot grob
    transp_plot_grob = ggplotGrob(transp_plot)

    # changes the dimension of the the second box
    transp_plot_grob$widths[9] = 2*transp_plot_grob$widths[9]

    # Builds and save the bar plot.
    svg(snakemake@output[[4]], width=10, height=7)
    # Draw the plot
    grid.newpage()
    grid.draw(transp_plot_grob)
    dev.off()

    transp_plot_perc <- ggplot(transp_class_all, aes(x = samples, y = quantif, fill = classif)) +
        geom_bar(position = "fill", stat = "identity", alpha = 0.9, width = 0.5) +
        scale_x_discrete(name="")+
        scale_y_continuous(labels = percent_format())+
        facet_wrap( ~ group_plot, ncol=2, scales="free")+
        ylab("Relative frequency")+
        theme_bw()+
        theme(legend.text=element_text(size=14), axis.title.y=element_text(size = 20, vjust=2), axis.title.x=element_text(size=14,vjust=0), axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15), legend.title = element_blank(), strip.text.x = element_text(size = 15))+
        scale_fill_manual(values=fill)

    # Get the ggplot grob
    transp_plot_perc_grob = ggplotGrob(transp_plot_perc)

    # changes the dimension of the the second box
    transp_plot_perc_grob$widths[9] = 2*transp_plot_perc_grob$widths[9]

    # Builds and save the bar plot.
    svg(snakemake@output[[5]], width=10, height=7)
    # Draw the plot
    grid.newpage()
    grid.draw(transp_plot_perc_grob)
    dev.off()
}
