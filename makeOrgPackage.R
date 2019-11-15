require(AnnotationForge)
require(tidyverse)

## Prepares the data.frame to build the annotation package.
load(snakemake@input[[1]])

go_annot <- annot_final[, c(1,6)]
colnames(go_annot) <- c("gene_name", "GO")

go_annot_per_row <- go_annot %>%
    mutate(V3 = strsplit(as.character(GO), "&")) %>%
    unnest(V3)

go_annot_per_row <- go_annot_per_row[, c(1, 3)]
colnames(go_annot_per_row)[2] <- "GO"

go_annot_per_row <- go_annot_per_row %>%
    mutate(V3 = strsplit(as.character(GO), ",")) %>%
    unnest(V3)

go_annot_per_row <- go_annot_per_row[, c(1, 3)]

go_annot_per_row$V3 <- gsub(pattern = " ", replacement = "", go_annot_per_row$V3)

go_annot_per_row <- cbind(go_annot_per_row,
                          rep("Blast2GO_BioMart", nrow(go_annot_per_row)))

colnames(go_annot_per_row) <- c("GID", "GO", "EVIDENCE")

go_annot_per_row$GO <- as.factor(go_annot_per_row$GO)
go_annot_per_row$GID <- as.character(noquote(go_annot_per_row$GID))

# Removes NAs
go_annot_final <- go_annot_per_row[!go_annot_per_row$GO == "NA",]

# Removes duplicate rows
go_annot_final$comb <- paste(go_annot_final$GID,
                             go_annot_final$GO,
                             go_annot_final$EVIDENCE,
                             sep = "_")

go_annot_final <- go_annot_final[duplicated(go_annot_final$comb) == FALSE, c(1:3)]

go_annot_final <- go_annot_final[go_annot_final[, 2] != "", ]
go_annot_final <- go_annot_final[go_annot_final[, 3] != "", ]

## Calls the function to build the annotation package
makeOrgPackage(go=go_annot_final,
               version="0.1",
               maintainer="wendell pereira <wendelljpereira@gmail.com>",
               author="wendell pereira",
               outputDir = ".",
               tax_id="71139",
               genus="Eucalyptus",
               species="grandis",
               goTable="go")

# Installs the generated package
install.packages("./org.Egrandis.eg.db", repos=NULL)
