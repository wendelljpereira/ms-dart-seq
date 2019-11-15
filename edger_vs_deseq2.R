require(gdata)
require(VennDiagram)
require(gridExtra)

# Function to remove redundance between MSD-Tags
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
                    print("Erro!")
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

# Help function to create a venn plot
venn_2_samples <- function(sample1 = sample1,
                           sample2 = sample2,
                           name1 = name1,
                           name2 = name2,
                           clone_name = clone_name,
                           save_ids = "FALSE"){

    inter <- unique(intersect(sample1[complete.cases(sample1)],
                              sample2[complete.cases(sample2)]))

    unic_s1 <- unique(setdiff(sample1[complete.cases(sample1)],
                              sample2[complete.cases(sample2)]))

    unic_s2 <- unique(setdiff(sample2[complete.cases(sample2)],
                              sample1[complete.cases(sample1)]))

    if (length(unic_s1) > length(unic_s2)){
        grid.newpage();
        venn.plot <- draw.pairwise.venn(
            area1 = length(unic_s1) + length(inter),
            area2 = length(unic_s2) + length(inter),
            cross.area = length(inter),
            alpha = 0.65,
            category = c(deparse(name1), deparse(name2)),
            fill = c("#8470FF", "#FF7F24"),
            lty = "blank",
            cex = 2,
            cat.cex = 2,
            cat.pos = c(0, 0),
            cat.dist = 0.055,
            cat.just = list(c(1, 0), c(0, 0)),
            ext.pos = 0,
            ext.dist = -0.05,
            ext.length = 0.85,
            ext.line.lwd = 2,
            ext.line.lty = "dashed",
            scaled = T,
            print.mode = c("raw", "percent"),
            rotation.degree = 0)
    }else{
        grid.newpage();
        venn.plot <- draw.pairwise.venn(
            area1 = length(unic_s1) + length(inter),
            area2 = length(unic_s2) + length(inter),
            cross.area = length(inter),
            alpha = 0.65,
            category = c(deparse(name1), deparse(name2)),
            fill = c("#8470FF", "#FF7F24"),
            lty = "blank",
            cex = 2,
            cat.cex = 2,
            cat.pos = c(0, 0),
            cat.dist = 0.055,
            cat.just = list(c(1, 0), c(0, 0)),
            ext.pos = 0,
            ext.dist = -0.05,
            ext.length = 0.85,
            ext.line.lwd = 2,
            ext.line.lty = "dashed",
            scaled = T,
            print.mode = c("raw", "percent"),
            rotation.degree = 180)
    }
}

# Reads the methylated sites
edger_dm_marks <- read.table(snakemake@input[[1]],
                             header = T,
                             colClasses = "character")

deseq_dm_marks <- read.table(snakemake@input[[2]],
                             header = T,
                             colClasses = "character")

# Determines the intersections
new_data_edger_deseq_intersect <- data.frame()
for (i in 1:ncol(edger_dm_marks)){

    intersect_df <- as.data.frame(intersect(edger_dm_marks[, i],
                                            deseq_dm_marks[, i]))

    colnames(intersect_df) <- paste(colnames(edger_dm_marks)[i])

    if (i == 1){

        new_data_edger_deseq_intersect <- intersect_df

    }else{

        new_data_edger_deseq_intersect <- cbindX(new_data_edger_deseq_intersect,
                                                 intersect_df)

    }
}

write.table(new_data_edger_deseq_intersect,
            snakemake@output[[1]],
            row.names = F,
            col.names = T,
            sep = "\t",
            quote = F)

sites_bed_file <- snakemake@input[[3]]

edger_deseq_intersect_corrected <- correct_sites(sites_bed_file = sites_bed_file,
                                                 data = new_data_edger_deseq_intersect)

write.table(edger_deseq_intersect_corrected,
            snakemake@output[[2]],
            row.names = F,
            col.names = T,
            sep = "\t",
            quote = F)

###############
## Venn plot ##
###############

edger_dm_marks_corrected <- correct_sites(sites_bed_file = sites_bed_file,
                                          data = edger_dm_marks)

deseq_dm_marks_corrected <- correct_sites(sites_bed_file = sites_bed_file,
                                          data = deseq_dm_marks)

## Plots
plot1 <- venn_2_samples(sample1 = edger_dm_marks_corrected[, 1],
                        sample2 = deseq_dm_marks_corrected[, 1],
                        name1 = "BRASUZ1 adult leaf - edgeR",
                        name2 = "BRASUZ1 adulte leaf - DEseq2")

plot2 <- venn_2_samples(sample1 = edger_dm_marks_corrected[, 2],
                        sample2 = deseq_dm_marks_corrected[, 2],
                        name1 = "BRASUZ1 juvenile leaf - edgeR",
                        name2 = "BRASUZ1 juvenile leaf - DEseq2")

plot3 <- venn_2_samples(sample1 = edger_dm_marks_corrected[, 3],
                        sample2 = deseq_dm_marks_corrected[, 3],
                        name1 = "BRASUZ1 xylem - edgeR",
                        name2 = "BRASUZ1 xylem - DEseq2")

svg(filename = snakemake@output[[3]],
    width = 18,
    height = 6,
    pointsize = 12)

grid.arrange(grobTree(plot1),
             grobTree(plot2),
             grobTree(plot3),
             ncol = 3,
             top = textGrob("edgeR vs DEseq2 methylated sites",
                          gp = gpar(fontsize = 40, font = 8)))
dev.off()
