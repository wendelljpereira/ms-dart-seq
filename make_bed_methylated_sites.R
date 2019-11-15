sites <- read.table(snakemake@input[[1]], sep = "\t")

# Reads the methylated sites
marks <- read.table(snakemake@input[[2]], header = T, sep = "\t")

marks_unique <- as.data.frame(unique(as.character(as.matrix(marks))))
marks_unique <- na.omit(marks_unique)

# Extract the position of the methylated sites
marks_pos <- sites[sites$V4 %in% marks_unique[, 1], ]

write.table(marks_pos,
            snakemake@output[[1]],
            col.names = F,
            row.names = F,
            quote = F,
            sep = "\t")

# Selects the positions in plus strand as representative of the methylation position
marks_pos_plus <- marks_pos[marks_pos$V6 == "+", ]

write.table(marks_pos_plus,
            snakemake@output[[2]],
            col.names = F,
            row.names = F,
            quote = F,
            sep = "\t")
