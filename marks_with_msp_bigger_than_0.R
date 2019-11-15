"Usage: marks_with_msp_bigger_than_0.R (--out1 <O1>) <input1>
-h --help    show this
--out1  name1   bed file with all the corrected positions of the MS-DArT Tags
input1  input1  bed file with all MS-DArT Tags
marks_with_msp_bigger_than_0.R -h | --help  show this message
" -> doc

# load the docopt library
require(docopt)
# retrieve the command-line arguments
opts <- docopt(doc)

require(gdata)

dados <- read.table(opts$`<input1>`,
                    header = T,
                    sep = "\t",
                    check.names = F)

# Select all sites with counts bigger than one in each column.

for (i in 2:ncol(dados)){
    if (i == 2){
        data_without_zero_all_samples <- dados[dados[, i] > 0, ]
        all_marks <- as.data.frame(data_without_zero_all_samples[, 1])
        colnames(all_marks) <- paste(names(dados)[i])
     }else{
         data_without_zero_sub <- dados[dados[, i] > 0, ]
         all_marks_sub <- as.data.frame(data_without_zero_sub[, 1])
         colnames(all_marks_sub) <- paste(names(dados)[i])
         all_marks <- cbindX(all_marks, all_marks_sub)
     }
}

# Separates the columns (samples) accordingly with the experimental groups (first name before the separator "_").
## Writes a file for each group.
    y <- data.frame()
for (i in 1:length(names(all_marks))){
    x <- as.data.frame(paste(strsplit(names(all_marks)[i], "_")[[1]][1]))
    colnames(x) <- "group"
    y <- rbind(y, x)
}
y$group <- as.factor(y$group)

for (i in 1:length(levels(y$group))){
    group <- all_marks[y$group == sort(levels(y$group))[i]]
    write.table(group,
                paste(opts$O1),
                row.names = F,
                quote = F,
                sep = "\t")
}