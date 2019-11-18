"Usage: statistical_tests.R (--out1 <O1>) (--out2 <O2>) (--out3 <O3>) (--out4 <O4>) <input1> <input2> <input3> <input4> <input5> <input6> <input7> <input8> <input9> <input10> <input11> <input12> <input13> <input14> <input15> <input16> <input17> <input18>
-h --help show this error message.
<input1>  name1 name1
<input2>  name2 name2
<input3>  name3 name3
<input4>  name4 name4
<input5>  name5 name5
<input6>  name6 name6
<input7>  name7 name7
<input8>  name8 name8
<input9>  name9 name9
<input10> name10  name10
<input11> name11  name11
<input12> name12  name12
<input13> name13  name13
<input14> name14  name14
<input15> name15  name15
<input16> name16  name16
<input17> name17  name17
<input18> name18  name18
--out1  name1 specify the name for the first output file
--out2  name2 specify the name for the first output file
--out3  name3 specify the name for the first output file
--out4  name4 specify the name for the first output file
statistical_tests.R -h | --help  show this message
" -> doc

require(docopt)
require(vcd)
require(ggplot2)
require(scales)
require(RVAideMemoire)
require(data.table)

# retrieve the command-line arguments
opts <- docopt(doc)
save.image("statistical_tests.Rda")

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

###############################################################
#### Tests the homogenety of the sampled sites distribution ####
###############################################################

counts_1Mb <- read.table(opts$`<input1>`, sep = "\t")
counts_500Kb <- read.table(opts$`<input2>`, sep = "\t")
counts_250Kb <- read.table(opts$`<input3>`, sep = "\t")
counts_100Kb <- read.table(opts$`<input4>`, sep = "\t")

# tests if the counts follow a uniform distribution
for(i in unique(counts_1Mb$V1)){
    subset_counts <- counts_1Mb[counts_1Mb$V1 == i, 4]
    print(i)
    print(chisq.test(subset_counts))
}

for(i in unique(counts_500Kb$V1)){
    subset_counts <- counts_500Kb[counts_500Kb$V1 == i, 4]
    print(i)
    print(chisq.test(subset_counts))
}

for(i in unique(counts_250Kb$V1)){
    subset_counts <- counts_250Kb[counts_250Kb$V1 == i, 4]
    print(i)
    print(chisq.test(subset_counts))
}

for(i in unique(counts_100Kb$V1)){
    subset_counts <- counts_100Kb[counts_100Kb$V1 == i, 4]
    print(i)
    print(chisq.test(subset_counts))
}

#####################################################################################
####       Test of the methylated sites distribution among BRASUZ1 tissues       ####
#####################################################################################

# This section applies the Cochran's Q test to test if there is differences between the number of methylated sites between the tissues.

# Read the file with all sampled sites in the intersection
inters_sites <- read.table(opts$`<input5>`,
                           header = T)[1]

sites_bed_file <- opts$`<input6>`

# Applies the function to remove redundance of the tags representing the same MS site
corrected_int_sites <- correct_sites(sites_bed_file = sites_bed_file,
                                                 data = inters_sites)

all_sites <- read.table(opts$`<input6>`)

selected_sites <- all_sites[all_sites$V4 %in% corrected_int_sites[, 1], ]

# Reads the methylated sites
methyl_sites <- read.table(opts$`<input7>`,
                           header = T)

# Organize the table to apply the Cochran's Q test. 1 if the site was methylated, 0 if not.
ad_leaves <- data.frame(tissue = "adult leaves",
                        MS_site = corrected_int_sites[, 1],
                        methylation_status = ifelse(corrected_int_sites[, 1] %in% methyl_sites[, 1], 1, 0))


juv_leaves <- data.frame(tissue = "juvenile leaves",
                         MS_site = corrected_int_sites[, 1],
                         methylation_status = ifelse(corrected_int_sites[, 1] %in% methyl_sites[, 2], 1, 0))

xylem <- data.frame(tissue = "xylem",
                    MS_site = corrected_int_sites[, 1],
                    methylation_status = ifelse(corrected_int_sites[, 1] %in% methyl_sites[, 3], 1, 0))

# Join the data in an new data.frame
ms_table <- rbind(ad_leaves, juv_leaves, xylem)

# Show the number of methylated sites (1) and unmethylated sites (0)
(plot_data <- xtabs( ~ tissue + methylation_status, data = ms_table))

cochran.qtest(methylation_status ~ tissue | MS_site, data = ms_table)

#######################################################################
####       Test of the methylated genes and methylated TEs         ####
#######################################################################

### genes ###

all_genes <- read.table(opts$`<input8>`)[9]
all_genes <- strsplit(as.character(all_genes[, 1]), "=")

all_genes <- data.frame(matrix(unlist(all_genes),
                               ncol = 3,
                               byrow = T))[3]

methyl_genes <- read.table(opts$`<input9>`,
                           header = T)

ad_leaves <- data.frame(tissue = "adult leaves",
                        genes = all_genes[, 1],
                        methy_status = ifelse(all_genes[, 1] %in% methyl_genes[, 1], 1, 0))

juv_leaves <- data.frame(tissue = "juvenile leaves",
                         genes = all_genes[, 1],
                         methy_status = ifelse(all_genes[, 1] %in% methyl_genes[, 2], 1, 0))

xylem  <- data.frame(tissue = "xylem",
                     genes = all_genes[, 1],
                     methy_status = ifelse(all_genes[, 1] %in% methyl_genes[, 3], 1, 0))

genes_table <- rbind(ad_leaves,
                     juv_leaves,
                     xylem)

# Show the number of methylated sites (1) and unmethylated sites (0)
(plot_data <- xtabs( ~ tissue + methy_status, data = genes_table))

cochran.qtest(methy_status ~ tissue | genes, data = genes_table)

methy_TEs <- read.table(opts$`<input10>`,
                        header = T)

all_TEs <- read.table(opts$`<input11>`)[4]
all_TEs <- unique(all_TEs)

ad_leaves <- data.frame(tissue = "adult leaves",
                        TEs = all_TEs[, 1],
                        methy_status = ifelse(all_TEs[, 1] %in% methy_TEs[, 1], 1, 0))

juv_leaves <- data.frame(tissue = "juvenile leaves",
                         TEs = all_TEs[, 1],
                         methy_status = ifelse(all_TEs[, 1] %in% methy_TEs[, 2], 1, 0))

xylem <- data.frame(tissue = "xylem",
                    TEs = all_TEs[, 1],
                    methy_status = ifelse(all_TEs[, 1] %in% methy_TEs[, 3], 1, 0))

TEs_table <- rbind(ad_leaves,
                   juv_leaves,
                   xylem)

# Shows the number of methylated sites (1) and unmethylated sites (0)
(plot_data <- xtabs( ~ tissue + methy_status, data = TEs_table))

cochran.qtest(methy_status ~ tissue | TEs, data = TEs_table)

###########################################
#### Test of proportions - TEs classes ####
###########################################

# reads the file with the TEs classification
table <- read.table(opts$`<input12>`,
                    sep = "\t",
                    header = F,
                    stringsAsFactors = FALSE)

table <- table[, c(2:4)]

colnames(table) <- c("Tissues", "TE_class", "counts")

# makes the matrix to chi-square test

adult_leaves <- table[table[, 1] == "Adult leaves", -1]
colnames(adult_leaves) <- c("TE_class", "Adult leaves")

juvenile_leaves <- table[table[, 1] == "Juvenile leaves", -1]
colnames(juvenile_leaves) <- c("TE_class", "Juvenile leaves")

xylem <- table[table[, 1] == "Xylem", -1]
colnames(xylem) <- c("TE_class", "Xylem")

genome <- table[table[, 1] == "Genome", -1]
colnames(genome) <- c("TE_class", "Genome")

table_class <- Reduce(function(x, y) merge(x, y, all = TRUE),
                      list(adult_leaves, juvenile_leaves, xylem, genome))

# changes the data to integer, so is possible change the NAs to 0
table_class_2 <- data.frame(as.character(table_class[, 1]),
                            as.integer(table_class[, 2]),
                            as.integer(table_class[, 3]),
                            as.integer(table_class[, 4]),
                            as.integer(table_class[, 5]),
                            stringsAsFactors = FALSE)

table_class_2[is.na(table_class_2)] <- 0

matrix <- matrix(unlist(table_class_2), 13)

matrix_int <- apply(matrix[, 2:5], 2, as.integer)
row.names(matrix_int) <- matrix[, 1]
colnames(matrix_int) <- colnames(table_class[2:5])

#Since some of the classes has an expected value lower than 5
## Uses all classes but simulating the p-value by MC
(result_int_mc <- chisq.test(matrix_int, simulate.p.value =  T, B = 10000))

#Uses the residuals to show the classes which differ among the samples
svg(paste(opts$O1),
    width = 8,
    height = 8,
    pointsize = 12)

assoc(matrix_int,
      main = "BRASUZ tissues",
      shade = TRUE,
      labeling = labeling_border(rot_labels = c(90,0,0,0),
                                 just_labels = c("left", "center", "center", "right")),
      labeling_args=list(gp_labels=gpar(fontsize=16),
                         gp_varnames=gpar(fontsize=16)),
      legend_args=list(fontsize=14),
      margins = unit(6, "lines"),
      legend_width = unit(7, "lines"),
      spacing = spacing_equal(unit(0.4, "lines")))

dev.off()

# Removes the classes with expected values < 5
matrix_int_sub <- matrix_int[-c(2,3,8,11), ]

(result_int_sub <- chisq.test(matrix_int_sub, correct = F))

####################################
###### Using only the tissues ######
####################################

# Excludes the date of the genome
tissues_matrix_int <- matrix_int[, !colnames(matrix_int) == "Genome"]
tissues_matrix_int_sub <- matrix_int_sub[, !colnames(matrix_int_sub) == "Genome"]

(result_int <- chisq.test(tissues_matrix_int_sub, correct = F))

# Removes the class with less than 5 counts
tissues_matrix_int_sub <- tissues_matrix_int_sub[-6, ]

(result_int <- chisq.test(tissues_matrix_int_sub, correct = F))

###############################################
#### Tests of proportions - Genomic context ####
###############################################

adult_leaves <- read.table(opts$`<input13>`,
                           header = T, sep = "\t")

juvenile_leaves <- read.table(opts$`<input14>`,
                              header = T, sep = "\t")

xylem <- read.table(opts$`<input15>`,
                    header = T, sep = "\t")

adult_leaves <- adult_leaves[, c(1,3)]
juvenile_leaves <- juvenile_leaves[, c(1,3)]
xylem <- xylem[, c(1,3)]

gene_class <- Reduce(function(x, y) merge(x, y, by = "category"),
                     list(adult_leaves, juvenile_leaves, xylem))

gene_class_matrix <- matrix(unlist(gene_class[, -1]), 11)
colnames(gene_class_matrix) <- c("Adult leaves", "Juvenile leaves", "Xylem")
rownames(gene_class_matrix) <- gene_class[, 1]

# Selects the specific classes
gene_class_matrix <- gene_class_matrix[-c(2, 3, 10, 11), ]

# Chi.test
(results_gene_class <- chisq.test(gene_class_matrix, correct = F))

#Uses the residuals to show the classes which differ among the samples
svg(paste(opts$O2),
    width = 8,
    height = 6,
    pointsize = 12)

assoc(gene_class_matrix,
      main = "BRASUZ tissues",
      shade = TRUE,
      labeling = labeling_border(rot_labels = c(90,0,0,0),
                                 just_labels = c("left", "center", "center", "right")),
      labeling_args=list(gp_labels=gpar(fontsize=16),
                         gp_varnames=gpar(fontsize=16)),
      legend_args=list(fontsize=14),
      margins = unit(6, "lines"),
      legend_width = unit(7, "lines"),
      spacing = spacing_equal(unit(0.4, "lines")))

dev.off()

######################################################################
### Tests of the distribution of the marks in the vicinity of genes ###
######################################################################

# Closests gene of each mark
all_genes_meth <- read.table("distance_to_genes_and_TEs/distance_to_genes_features.txt", colClasses = "character")
all_genes_meth$V16 <- as.numeric(all_genes_meth$V16)

# Closests TE of each mark
all_transposons_meth <- read.table("distance_to_genes_and_TEs/distance_to_transposons_features.txt",
                              sep = "\t",
                              colClasses = "character")

all_transposons_meth$V13 <- as.numeric(all_transposons_meth$V13)

# Cleaning
transposons_meth <- cbind(as.data.frame(rep("transposon", nrow(all_transposons_meth))),
                     all_transposons_meth)

transposons_meth <- transposons_meth[!transposons_meth$V11 == ".", ]
transposons_meth <- transposons_meth[, c(1, 5, 14)]
colnames(transposons_meth) <- c("feature", "mark", "distance")

genes_meth <- cbind(as.data.frame(rep("gene", nrow(all_genes_meth))), all_genes_meth)
genes_meth <- genes_meth[genes_meth$V9 == "gene", ]
genes_meth <- genes_meth[, c(1, 5, 17)]
colnames(genes_meth) <- c("feature", "mark", "distance")

# combines the genes and Tes
features_meth <- rbind(genes_meth, transposons_meth)

# Filters the regions of interest (10Kb)
features_meth$distance <- abs(features_meth$distance)
features_meth <- cbind(as.data.frame(rep("Mapped reads", nrow(features_meth))), features_meth)
colnames(features_meth) <- c("type", "feature", "mark", "distance")

features_10kb_meth <- features_meth[features_meth$distance <= 10000, ]

bin <- c(1, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10001)
features_10kb_meth <- cbind(features_10kb_meth, findInterval(features_10kb_meth$distance, bin))
colnames(features_10kb_meth) <- c("type", "feature", "mark", "distance", "groups")
features_10kb_meth$groups <- as.factor(features_10kb_meth$groups)

# Closest gene of each marks
all_genes <- read.table("distance_to_genes_and_TEs/mspI_distance_to_genes_features.tst",
                        colClasses = "character",
                        sep = "\t")
all_genes$V16 <- as.numeric(all_genes$V16)

# Closest TEs of each marks
all_transposons <- read.table("distance_to_genes_and_TEs/mspI_distance_to_transposons_features.tst",
                              sep = "\t",
                              colClasses = "character")
all_transposons$V13 <- as.numeric(all_transposons$V13)

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
features_10kb <- features_10kb[!features_10kb$distance == 0, ]

all_sites_10Kb <- features_10kb

bin <- c(1, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10001)
features_10kb <- cbind(features_10kb, findInterval(features_10kb$distance, bin))
colnames(features_10kb) <- c("type", "feature", "mark", "distance", "groups")
features_10kb$groups <- as.factor(features_10kb$groups)

features_10kb_meth$site_class <- "methylated_sites"
features_10kb$site_class <- "all_mspI"

all_features_comb <- rbind(features_10kb_meth,
                           features_10kb)

bin <- c(1, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10001)

all_features_comb <- cbind(all_features_comb,
                           findInterval(all_features_comb$distance, bin))

all_features_comb <- all_features_comb[, -c(1,5)]

colnames(all_features_comb) <- c("feature",
                                 "mark",
                                 "distance",
                                 "site_class",
                                 "groups")

all_features_comb$groups <- as.factor(all_features_comb$groups)

###########
## Genes ##
###########
all_features_genes <- all_features_comb[all_features_comb$feature == "gene", ]

gene_df <- as.data.frame(all_features_genes %>%
  group_by(site_class, groups) %>%
  summarise(count = n()))

gene_class_matrix <- matrix(unlist(gene_df[, -1]), 10 )
gene_class_matrix <- gene_class_matrix[, -c(1, 2)]

colnames(gene_class_matrix) <- c("MspI sites", "All methylated sites")
rownames(gene_class_matrix) <- c("0.001-1 kb",
                                 "1-2 kb",
                                 "2-3 kb",
                                 "3-4 kb",
                                 "4-5 kb",
                                 "5-6 kb",
                                 "6-7 kb",
                                 "7-8 kb",
                                 "8-9 kb",
                                 "9-10 kb")

# Chi.test
(results_gene_class <- chisq.test(gene_class_matrix, correct = F))

# Plots the residuals
svg("images/association_plots/association_plot_mspI_vicinity_genes.svg",
    width = 8,
    height = 8,
    pointsize = 12)

assoc(gene_class_matrix,
      main = "BRASUZ tissues",
      shade = TRUE,
      labeling = labeling_border(rot_labels = c(90,0,0,0),
                                 just_labels = c("left", "center", "center", "right")),
      labeling_args=list(gp_labels=gpar(fontsize=16),
                         gp_varnames=gpar(fontsize=16)),
      legend_args=list(fontsize=14),
      margins = unit(6, "lines"),
      legend_width = unit(7, "lines"),
      spacing = spacing_equal(unit(0.4, "lines")))

dev.off()

###########
##  TEs  ##
###########
all_features_tes <- all_features_comb[all_features_comb$feature == "transposon", ]

tes_df <- as.data.frame(all_features_tes %>%
                           group_by(site_class, groups) %>%
                           summarise(count = n()))

tes_class_matrix <- matrix(unlist(tes_df[, -1]), 10 )
tes_class_matrix <- tes_class_matrix[, -c(1, 2)]

colnames(tes_class_matrix) <- c("MspI sites", "All methylated sites")
rownames(tes_class_matrix) <- c("0.001-1 kb",
                                 "1-2 kb",
                                 "2-3 kb",
                                 "3-4 kb",
                                 "4-5 kb",
                                 "5-6 kb",
                                 "6-7 kb",
                                 "7-8 kb",
                                 "8-9 kb",
                                 "9-10 kb")

# Chi.test
(results_tes_class <- chisq.test(tes_class_matrix, correct = F))

# Plots the residuals
svg("images/association_plots/association_plot_mspI_vicinity_tes.svg",
    width = 8,
    height = 8,
    pointsize = 12)

assoc(tes_class_matrix,
      main = "BRASUZ tissues",
      shade = TRUE,
      labeling = labeling_border(rot_labels = c(90,0,0,0),
                                 just_labels = c("left", "center", "center", "right")),
      labeling_args=list(gp_labels=gpar(fontsize=16),
                         gp_varnames=gpar(fontsize=16)),
      legend_args=list(fontsize=14),
      margins = unit(6, "lines"),
      legend_width = unit(7, "lines"),
      spacing = spacing_equal(unit(0.4, "lines")))

dev.off()
