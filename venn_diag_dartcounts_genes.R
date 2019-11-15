require(plyr)
require(VennDiagram)
require(gridExtra)

G4_int <- read.table(snakemake@input[[1]], header = T, na.strings = "NA", colClasses = "character")

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

    if(save_ids == "TRUE"){
        write.table(inter_123_full, paste(clone_name, name1, "vs", name2, "vs", name3, "intersection_genes.tst", sep = "_"), row.names = F, col.names = F, quote = F)
        write.table(inter_12, paste(clone_name, "only", name1, "vs", name2, "intersection_genes.tst", sep = "_"), row.names = F, col.names = F, quote = F)
        write.table(inter_13, paste(clone_name, "only", name1, "vs", name3, "intersection_genes.tst", sep = "_"), row.names = F, col.names = F, quote = F)
        write.table(inter_23, paste(clone_name, "only", name2, "vs", name3, "intersection_genes.tst", sep = "_"), row.names = F, col.names = F, quote = F)
        write.table(inter_12_full, paste(clone_name, "all", name1, "vs", name2, "intersection_genes.tst", sep = "_"), row.names = F, col.names = F, quote = F)
        write.table(inter_13_full, paste(clone_name, "all", name1, "vs", name3, "intersection_genes.tst", sep = "_"), row.names = F, col.names = F, quote = F)
        write.table(inter_23_full, paste(clone_name, "all", name2, "vs", name3, "intersection_genes.tst", sep = "_"), row.names = F, col.names = F, quote = F)
        write.table(unic_s1, paste(clone_name, name1, "exclusive_genes.tst", sep = "_"), row.names = F, col.names = F, quote = F)
        write.table(unic_s2, paste(clone_name, name2, "exclusive_genes.tst", sep = "_"), row.names = F, col.names = F, quote = F)
        write.table(unic_s3, paste(clone_name, name3, "exclusive_genes.tst", sep = "_"), row.names = F, col.names = F, quote = F)
    }else if(save_ids == "FALSE"){
    }else{
        print("Invalid save_ids option!")
    }

    # venn plot
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

# Venn plot of the methylated genes
A <- venn_3_samples(sample1 = G4_int[,1],
                    sample2 = G4_int[,2],
                    sample3 = G4_int[,3],
                    name1 = snakemake@params[["tissues"]][1],
                    name2 = snakemake@params[["tissues"]][2],
                    name3 = snakemake@params[["tissues"]][3],
                    clone_name = snakemake@params[["clone_name"]],
                    save_ids = snakemake@params[["save_ids"]])

svg(filename = snakemake@output[[1]], width = 12, height = 12, pointsize = 12)
grid.arrange(grobTree(A), ncol=1, top=textGrob(snakemake@params[["clone_name"]], gp=gpar(fontsize=30,font=8)))
dev.off()
