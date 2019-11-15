"Usage: genes_with_GBM.R <input1> <input2> <input3> <input4> <input5> <input6> <input7> <input8>
-h --help    show this
genes_with_GBM.R -h | --help  show this message
" -> doc

# load the necessary packages
require(tidyverse)
require(gdata)
require(docopt)

# retrieve the command-line arguments
opts <- docopt(doc)
save.image(file="docopt.rda")

correct_sites <- function(sites_bed_file = sites_bed_file, data = data){

    all_sites <- read.table(sites_bed_file, sep = "\t")
    all_sites <- all_sites[all_sites$V6 == "+", ]

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

############################################
#### % of genes sampled and that are GbM ###
############################################

# read the list of genes sampled for at least one MS-DArT-seq Tag.
all_sampled_genes <- read.table(opts$`<input1>`,
                                sep = "\t")

# Determines how many of the sampled genes were considered as GbM.
bisulfito <- read.table(opts$`<input2>`,
                        header = T,
                        sep = ",")

bisulfito_gbm <- bisulfito[bisulfito$class == "CG-gbm", ]

sampled_genes_gbm <- all_sampled_genes[all_sampled_genes$V1 %in% bisulfito_gbm$gene, ]

print(paste("Of the ", nrow(all_sampled_genes), " sampled genes, ", length(sampled_genes_gbm), " are genes with GbM.", sep = ""))

############################################
#### % of sampled site in genes with GbM ###
############################################

# reads the file with the sampled sites in genes
sampled_sites_in_genes <- read.table(opts$`<input3>`,
                                     sep = "\t")

sites_intersection <- read.table(opts$`<input4>`,
                                 header = T,
                                 sep = "\t")

intersect_corrected <- correct_sites(sites_bed_file = opts$`<input5>`,
                                     data = sites_intersection)

sampled_sites_in_genes <- sampled_sites_in_genes[sampled_sites_in_genes$V4 %in% intersect_corrected$BRASUZ1_leaf.ad, ]


print(paste("Of the ", nrow(all_sampled_genes), " sampled genes, ", length(sampled_genes_gbm), " are genes with GbM.", sep = ""))

# number of samples sites in genes (after remove the marks outside the set of marks in the intersection of the three tissues.)
length(unique(sampled_sites_in_genes$V4))

# Generates a file with all genes that have at least one of the sequenced mspI sites.
genes_names <- strsplit(paste(sampled_sites_in_genes$V15), "=")
genes_names_df <- as.data.frame(do.call("rbind", genes_names))

sampled_sites_in_genes$gene_names <- genes_names_df$V3

# Removes the duplicated marks that are in genes which overlaps other gene.
sampled_sites_in_genes <- sampled_sites_in_genes[sampled_sites_in_genes$gene_names %in% all_sampled_genes$V1, ]

sampled_sites_in_genes_gbm <- sampled_sites_in_genes[sampled_sites_in_genes$gene_names %in% bisulfito_gbm$gene, ]

names_sampled_sites_in_genes_gbm <- unique(sampled_sites_in_genes_gbm$V4)

# Number of sampled sites within genes with gbm
print(paste("Number of sampled sites within genes with gbm: ", length(names_sampled_sites_in_genes_gbm), sep = ""))

# Number of genes with GBM that contains at least one of the sampled sites
names_genes_with_sampled_sites <- unique(sampled_sites_in_genes_gbm$gene_names)
print(paste("Number of gbm genes that contains sequenced sites: ", length(names_genes_with_sampled_sites), sep = ""))


## Reads the files with the information of the methylated sites within genes for each tissue.
folha_adulta <- read.table(opts$`<input6>`,
                           header = T,
                           sep = "\t",
                           quote = "\"")

folha_juvenil <- read.table(opts$`<input7>`,
                            header = T,
                            sep = "\t",
                            quote = "\"")

xilema <- read.table(opts$`<input8>`,
                     header = T,
                     sep = "\t",
                     quote = "\"")

all_tissues <- rbind(folha_adulta, folha_juvenil, xilema)

# Filter to keep the methylated sites in genes with GbM
all_meth_sites_in_gbm <- all_tissues[all_tissues$gene_name %in% bisulfito_gbm$gene, ]

all_meth_sites_in_gbm <- separate(data = all_meth_sites_in_gbm,
                                  col = "locus",
                                  into = c("col1",
                                           "col2",
                                           "col3",
                                           "col4",
                                           "col5",
                                           "col6",
                                           "col7",
                                           "col8",
                                           "col9",
                                           "col10"),
                                  sep = ",",
                                  fill = "right")

all_meth_sites_in_gbm <- c(all_meth_sites_in_gbm$col1,
                           all_meth_sites_in_gbm$col2,
                           all_meth_sites_in_gbm$col3,
                           all_meth_sites_in_gbm$col4,
                           all_meth_sites_in_gbm$col5,
                           all_meth_sites_in_gbm$col6,
                           all_meth_sites_in_gbm$col7,
                           all_meth_sites_in_gbm$col8,
                           all_meth_sites_in_gbm$col9,
                           all_meth_sites_in_gbm$col10)

all_meth_sites_in_gbm <- unique(all_meth_sites_in_gbm)

# Number of methylated sites within genes with gbm
print(paste("Number of methylated sites within genes with gbm: ", length(all_meth_sites_in_gbm), sep = ""))
