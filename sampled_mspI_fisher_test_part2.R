"Usage: sampled_mspI_fisher_test_part2.R (--out1 <O1>) <input1> <input2> <input3> <input4> <input5> <input6>
-h --help    show this
--out1   tab_file    List of genes with sequenced MspI sites
<input1>   bed_file    position_of_the_sampled_sites/msdartseq_methylation_sites_of_sequenced_fragments_merged.bed
<input2>   tab_file    tested_sites_to_methylation.txt
<input3>   tab_file    sampled_mspI_fisher_test/mspI_genes_intersect.txt
<input4>   tab_file    sampled_mspI_fisher_test/tested_mspI_genes_intersect.txt
<input5>   bed_file    methylated_sites/Brasuz1_methylatied_sites.bed
<input6>   bed_file    restriction_sites/msp_pst_sites_positions.bed
sampled_mspI_fisher_test_part2.R -h | --help  show this message
" -> doc

require(docopt)
require(data.table)

# retrieve the command-line arguments
opts <- docopt(doc)
save.image(file = "docopt.rda")

# Help function to remove the redunction of the MS tags
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
                x <- grep(pattern = paste0("\\<", new_data_elements[j, 1], "\\|"),
                          x = dup_sites)
                y <- grep(pattern = paste0("\\|", new_data_elements[j, 1], "\\>"),
                          x = dup_sites)

                if (length(x) > 0){
                    z <- x

                }else if (length(y) > 0){
                    z <- y

                }else{
                    print("Erro!")
                }

                new_data_renamed <- as.data.frame(dup_sites[z])
                colnames(new_data_renamed) <- "new_data_renamed"
                new_data_renamed_all <- rbind(new_data_renamed_all,
                                              new_data_renamed)
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

# Reads the set of sampled sites.
tested_sites_mspI <- read.table(opts$`<input1>`)
tested_sites_mspI <- tested_sites_mspI[tested_sites_mspI$V6 == "+", ]

# Reads the set of marks that were used by edgeR and Deseq2
tested_sites_mspI_intersection <- read.table(opts$`<input2>`,
                                             header = F)

# Removes the redundance between tags
tested_sites_mspI_intersection <- correct_sites(sites_bed_file = opts$`<input1>`,
                                                data = tested_sites_mspI_intersection)

# Reads the MspI sites in genes
all_mspI_in_genes <- read.table(opts$`<input3>`, sep = "\t")

all_mspI_in_genes <- all_mspI_in_genes[!duplicated(paste(all_mspI_in_genes$V1,
                                                         all_mspI_in_genes$V2,
                                                         all_mspI_in_genes$V3,
                                                         all_mspI_in_genes$V4)), ]

# Sequenced mspI sites in genes
tested_mspI_in_genes <- read.table(opts$`<input4>`, sep = "\t")
tested_mspI_in_genes <- tested_mspI_in_genes[!duplicated(tested_mspI_in_genes$V4), ]

# List of genes with at least one sequenced MspI site
genes_names <- strsplit(paste(tested_mspI_in_genes$V15), "=")
genes_names_df <- as.data.frame(do.call("rbind", genes_names))
names_genes_with_sequenced_mspI <- unique(genes_names_df$V3)

write.table(names_genes_with_sequenced_mspI,
            paste(opts$O1),
            col.names = F,
            row.names = F,
            quote = F,
            sep = "\t")

# Reads all restriction sites detected in the genome
genome_sites <- as.data.frame(fread(opts$`<input5>`))
genome_sites_mspI <- genome_sites[genome_sites$V4 == "MspI" & genome_sites$V6 == "+", ]

#################################
# Fisher exact test - all sites #
#################################

A <- length(unique(tested_mspI_in_genes$V4))
B <- nrow(tested_sites_mspI) - A
C <- nrow(all_mspI_in_genes) - A
D <- (nrow(genome_sites_mspI) - nrow(all_mspI_in_genes)) - B

cont_table <- matrix(c(A, B, C, D),
                     nrow = 2,
                     dimnames = list(c("Genes", "Intergenic"),
                                     c("Sampled", "Ñ Sampled")))

print("Sampled MspI sites")
print(cont_table)
fisher.test(cont_table, alternative = "two.side")
fisher.test(cont_table, alternative = "greater")

############################################
#  Fisher exact test - intersection sites  #
############################################
tested_mspI_in_genes_intersect <- tested_mspI_in_genes[tested_mspI_in_genes$V4 %in% tested_sites_mspI_intersection$V1, ]

A <- length(unique(tested_mspI_in_genes_intersect$V4))
B <- nrow(tested_sites_mspI_intersection) - A
C <- nrow(all_mspI_in_genes)- A
D <- (nrow(genome_sites_mspI)-nrow(all_mspI_in_genes)) - B

cont_table <- matrix(c(A, B, C, D),
                     nrow = 2,
                     dimnames = list(c("Genes", "Intergenic"),
                                     c("Sampled", "Ñ Sampled")))

print("Sampled MspI sites - Intersection of samples")
print(cont_table)
fisher.test(cont_table, alternative = "two.side")
fisher.test(cont_table, alternative = "greater")

#####################################
# MDS Methylated sites vs MSD Sites #
#####################################

tested_sites_mspI <- tested_sites_mspI[tested_sites_mspI$V4 %in% tested_sites_mspI_intersection$V1, ]
tested_mspI_in_genes <- tested_mspI_in_genes[tested_mspI_in_genes$V4 %in% tested_sites_mspI_intersection$V1, ]

testes_mspI_sites_intergenic <- tested_sites_mspI[!tested_sites_mspI$V4 %in% tested_mspI_in_genes$V4, ]

methylated_sites <- read.table(opts$`<input6>`)

E <- nrow(tested_mspI_in_genes[tested_mspI_in_genes$V4 %in% methylated_sites$V4, ])
F. <- nrow(testes_mspI_sites_intergenic[testes_mspI_sites_intergenic$V4 %in% methylated_sites$V4, ])
G <- nrow(tested_mspI_in_genes) - E
H <- nrow(testes_mspI_sites_intergenic) - F.

# Remove os sítios da classe unknow
E <- E - 307

cont_table <- matrix(c(E, F., G, H),
                     nrow = 2,
                     dimnames = list(c("Genes", "Intergenic"),
                                     c("Methylated", "Ñ Methylated")))

print("Methylated sites")
print(cont_table)
fisher.test(cont_table, alternative = "two.side")
fisher.test(cont_table, alternative = "greater")
