# Fisher Test - Contingency table #

############################
## venn plot of all marks ##
############################

print("Fisher Test - all MSD-Sites")
# Reads the data
msd_marks <- read.table(snakemake@input[[1]],
                       sep = "\t",
                       header = T,
                       colClasses = "character")

# Removes NAs
msd_marks <- lapply(msd_marks, function(x) x[!is.na(x)])

total_tested_marks <- 19526
dm_brasuz.ad <- length(msd_marks[[1]])
dm_brasuz.juv <- length(msd_marks[[2]])
dm_brasuz.wood <- length(msd_marks[[3]])
dm_int_ad_juv <- length(intersect(msd_marks[[1]], msd_marks[[2]]))
dm_int_ad_wood <- length(intersect(msd_marks[[1]], msd_marks[[3]]))
dm_int_juv_wood <- length(intersect(msd_marks[[2]], msd_marks[[3]]))

#BRASUZ1 Leaf.juv vs Leaf.ad
A <- dm_int_ad_juv
B <- dm_brasuz.juv - dm_int_ad_juv
C <- dm_brasuz.ad - dm_int_ad_juv
D <- total_tested_marks - (length(unique(as.character(c(msd_marks[[1]], msd_marks[[2]])))))

cont_table <- matrix(c(A, B, C, D),
                     nrow = 2,
                     dimnames = list(Leaf.juv = c("Methylated", "ñ Methylated"),
                                     Leaf.ad = c("Methylated", "ñ Methylated")))

print("BRASUZ1")
print(cont_table)
fisher.test(cont_table, alternative = "greater")

#BRASUZ1 Leaf.juv vs Leaf.wood
A <- dm_int_juv_wood
B <- dm_brasuz.juv - dm_int_juv_wood
C <- dm_brasuz.wood - dm_int_juv_wood
D <- total_tested_marks - (length(unique(as.character(c(msd_marks[[2]], msd_marks[[3]])))))

cont_table <- matrix(c(A, B, C, D),
                     nrow = 2,
                     dimnames = list(Leaf.juv = c("Methylated", "ñ Methylated"),
                                     Leaf.wood = c("Methylated", "ñ Methylated")))

print("BRASUZ1")
print(cont_table)
fisher.test(cont_table, alternative = "greater")

#BRASUZ1 Leaf.ad vs Leaf.wood
A <- dm_int_ad_wood
B <- dm_brasuz.ad - dm_int_ad_wood
C <- dm_brasuz.wood - dm_int_ad_wood
D <- total_tested_marks - (length(unique(as.character(c(msd_marks[[1]], msd_marks[[3]])))))

cont_table <- matrix(c(A, B, C, D),
                     nrow = 2,
                     dimnames = list(Leaf.juv = c("Methylated", "ñ Methylated"),
                                     Leaf.wood = c("Methylated", "ñ Methylated")))

print("BRASUZ1")
print(cont_table)
fisher.test(cont_table, alternative = "greater")

###################################################
## venn plot marks without a gene or transposon  ##
###################################################

print("Fisher Test - MSD-Sites outside genes or tranposons")
# Reads the data
marks_out_gt <- read.table(snakemake@input[[2]],
                           sep = "\t",
                           header = T,
                           colClasses = "character")

# Removes NAs
marks_out_gt <- lapply(marks_out_gt, function(x) x[!is.na(x)])

#
total_tested_marks <- 20386
dm_brasuz.ad <- length(marks_out_gt[[1]])
dm_brasuz.juv <- length(marks_out_gt[[2]])
dm_brasuz.wood <- length(marks_out_gt[[3]])
dm_int_ad_juv <- length(intersect(marks_out_gt[[1]], marks_out_gt[[2]]))
dm_int_ad_wood <- length(intersect(marks_out_gt[[1]], marks_out_gt[[3]]))
dm_int_juv_wood <- length(intersect(marks_out_gt[[2]], marks_out_gt[[3]]))

#BRASUZ1 Leaf.juv vs Leaf.ad
A <- dm_int_ad_juv
B <- dm_brasuz.juv - dm_int_ad_juv
C <- dm_brasuz.ad - dm_int_ad_juv
D <- total_tested_marks - (length(unique(as.character(c(marks_out_gt[[1]], marks_out_gt[[2]])))))

cont_table <- matrix(c(A, B, C, D),
                     nrow = 2,
                     dimnames = list(Leaf.juv = c("Methylated", "ñ Methylated"),
                                     Leaf.ad = c("Methylated", "ñ Methylated")))

print("BRASUZ1")
print(cont_table)
fisher.test(cont_table, alternative = "greater")

#BRASUZ1 Leaf.juv vs Leaf.wood
A <- dm_int_juv_wood
B <- dm_brasuz.juv - dm_int_juv_wood
C <- dm_brasuz.wood - dm_int_juv_wood
D <- total_tested_marks - (length(unique(as.character(c(marks_out_gt[[2]], marks_out_gt[[3]])))))

cont_table <- matrix(c(A, B, C, D),
                     nrow = 2,
                     dimnames = list(Leaf.juv = c("Methylated", "ñ Methylated"),
                                     Leaf.wood = c("Methylated", "ñ Methylated")))

print("BRASUZ1")
print(cont_table)
fisher.test(cont_table, alternative = "greater")

#BRASUZ1 Leaf.ad vs Leaf.wood
A <- dm_int_ad_wood
B <- dm_brasuz.ad - dm_int_ad_wood
C <- dm_brasuz.wood - dm_int_ad_wood
D <- total_tested_marks - (length(unique(as.character(c(marks_out_gt[[1]], marks_out_gt[[3]])))))

cont_table <- matrix(c(A, B, C, D),
                     nrow = 2,
                     dimnames = list(Leaf.juv = c("Methylated", "ñ Methylated"),
                                     Leaf.wood = c("Methylated", "ñ Methylated")))

print("BRASUZ1")
print(cont_table)
fisher.test(cont_table, alternative = "greater")

######################
## Methylated genes ##
######################

print("Fisher Test - Methylated genes")

# Reads the data
dm_genes <- read.table(snakemake@input[[3]],
                       sep = "\t",
                       header = T,
                       colClasses = "character")

# Removes NAs
dm_genes <- lapply(dm_genes, function(x) x[!is.na(x)])

#
total_tested_genes <- 36349
dm_brasuz.ad <- length(dm_genes[[1]])
dm_brasuz.juv <- length(dm_genes[[2]])
dm_brasuz.wood <- length(dm_genes[[3]])
dm_int_ad_juv <- length(intersect(dm_genes[[1]], dm_genes[[2]]))
dm_int_ad_wood <- length(intersect(dm_genes[[1]], dm_genes[[3]]))
dm_int_juv_wood <- length(intersect(dm_genes[[2]], dm_genes[[3]]))

#BRASUZ1 Leaf.juv vs Leaf.ad
A <- dm_int_ad_juv
B <- dm_brasuz.juv - dm_int_ad_juv
C <- dm_brasuz.ad - dm_int_ad_juv
D <- total_tested_genes - (length(unique(as.character(c(dm_genes[[1]], dm_genes[[2]])))))

cont_table <- matrix(c(A, B, C, D),
                     nrow = 2,
                     dimnames = list(Leaf.juv = c("Methylated", "ñ Methylated"),
                                     Leaf.ad = c("Methylated", "ñ Methylated")))

print("BRASUZ1")
print(cont_table)
fisher.test(cont_table, alternative = "two.sided")

#BRASUZ1 Leaf.juv vs Leaf.wood
A <- dm_int_juv_wood
B <- dm_brasuz.juv - dm_int_juv_wood
C <- dm_brasuz.wood - dm_int_juv_wood
D <- total_tested_genes - (length(unique(as.character(c(dm_genes[[2]], dm_genes[[3]])))))

cont_table <- matrix(c(A, B, C, D),
                     nrow = 2,
                     dimnames = list(Leaf.juv = c("Methylated", "ñ Methylated"),
                                     Leaf.wood = c("Methylated", "ñ Methylated")))

print("BRASUZ1")
print(cont_table)
fisher.test(cont_table, alternative = "two.sided")

#BRASUZ1 Leaf.ad vs Leaf.wood
A <- dm_int_ad_wood
B <- dm_brasuz.ad - dm_int_ad_wood
C <- dm_brasuz.wood - dm_int_ad_wood
D <- total_tested_genes - (length(unique(as.character(c(dm_genes[[1]], dm_genes[[3]])))))

cont_table <- matrix(c(A, B, C, D),
                     nrow = 2,
                     dimnames = list(Leaf.juv = c("Methylated", "ñ Methylated"),
                                     Leaf.wood = c("Methylated", "ñ Methylated")))

print("BRASUZ1")
print(cont_table)
fisher.test(cont_table, alternative = "two.sided")

####################
## Methylated TEs ##
####################

print("Fisher Test - Methylated transposons")

# Reads the data
dm_transposons <- read.table(snakemake@input[[4]],
                             sep = "\t",
                             header = T,
                             colClasses = "character")

# Removes NAs
dm_transposons <- lapply(dm_transposons, function(x) x[!is.na(x)])

#
total_tested_transposons <- 223479
dm_brasuz.ad <- length(dm_transposons[[1]])
dm_brasuz.juv <- length(dm_transposons[[2]])
dm_brasuz.wood <- length(dm_transposons[[3]])
dm_int_ad_juv <- length(intersect(dm_transposons[[1]], dm_transposons[[2]]))
dm_int_ad_wood <- length(intersect(dm_transposons[[1]], dm_transposons[[3]]))
dm_int_juv_wood <- length(intersect(dm_transposons[[2]], dm_transposons[[3]]))

#BRASUZ1 Leaf.juv vs Leaf.ad
A <- dm_int_ad_juv
B <- dm_brasuz.juv - dm_int_ad_juv
C <- dm_brasuz.ad - dm_int_ad_juv
D <- total_tested_transposons - (length(unique(as.character(c(dm_transposons[[1]], dm_transposons[[2]])))))

cont_table <- matrix(c(A, B, C, D),
                     nrow = 2,
                     dimnames = list(Leaf.juv = c("Methylated", "ñ Methylated"),
                                     Leaf.ad = c("Methylated", "ñ Methylated")))

print("BRASUZ1")
print(cont_table)
fisher.test(cont_table, alternative = "two.sided")

#BRASUZ1 Leaf.juv vs Leaf.wood
A <- dm_int_juv_wood
B <- dm_brasuz.juv - dm_int_juv_wood
C <- dm_brasuz.wood - dm_int_juv_wood
D <- total_tested_transposons - (length(unique(as.character(c(dm_transposons[[2]], dm_transposons[[3]])))))

cont_table <- matrix(c(A, B, C, D),
                     nrow = 2,
                     dimnames = list(Leaf.juv = c("Methylated", "ñ Methylated"),
                                     Leaf.wood = c("Methylated", "ñ Methylated")))

print("BRASUZ1")
print(cont_table)
fisher.test(cont_table, alternative = "two.sided")

#BRASUZ1 Leaf.ad vs Leaf.wood
A <- dm_int_ad_wood
B <- dm_brasuz.ad - dm_int_ad_wood
C <- dm_brasuz.wood - dm_int_ad_wood
D <- total_tested_transposons - (length(unique(as.character(c(dm_transposons[[1]], dm_transposons[[3]])))))

cont_table <- matrix(c(A, B, C, D),
                     nrow = 2,
                     dimnames = list(Leaf.juv = c("Methylated", "ñ Methylated"),
                                     Leaf.wood = c("Methylated", "ñ Methylated")))

print("BRASUZ1")
print(cont_table)
fisher.test(cont_table, alternative = "two.sided")

#####################
## edgeR_vs_deseq2 ##
#####################

#BRASUZ1 Adult Leaf
A <- 3166
B <- 98
C <- 276
D <- 13322

cont_table <- matrix(c(A, B, C, D),
                     nrow = 2,
                     dimnames = list(edgeR = c("Methylated", "ñ Methylated"),
                                     DEseq2 = c("Methylated", "ñ Methylated")))

print("BRASUZ1")
print(cont_table)
fisher.test(cont_table, alternative = "two.sided")

#BRASUZ1 Juvenile Leaf
A <- 3257
B <- 100
C <- 227
D <- 13211

cont_table <- matrix(c(A, B, C, D),
                     nrow = 2,
                     dimnames = list(edgeR = c("Methylated", "ñ Methylated"),
                                     DEseq2 = c("Methylated", "ñ Methylated")))

print("BRASUZ1")
print(cont_table)
fisher.test(cont_table, alternative = "two.sided")

#BRASUZ1  Xylem
A <- 2966
B <- 92
C <- 295
D <- 12887

cont_table <- matrix(c(A, B, C, D),
                     nrow = 2,
                     dimnames = list(edgeR = c("Methylated", "ñ Methylated"),
                                     DEseq2 = c("Methylated", "ñ Methylated")))

print("BRASUZ1")
print(cont_table)
fisher.test(cont_table, alternative = "two.sided")
