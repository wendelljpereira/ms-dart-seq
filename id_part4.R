"Usage: id_part4.R (--out1 <O1>) (--out2 <O2>) (--out3 <O3>) (--out4 <O4>) <input1> <input2> <input3> <input4>
-h --help    show this
--out1   name1    specify the name for the first output file
--out2   name2    specify the name for the second output file
--out3   name3    specify the name for the third output file
--out4   name4    specify the name for the third output file
id_part4.R -h | --help  show this message
" -> doc

# load the docopt library
require(docopt)
# retrieve the command-line arguments
opts <- docopt(doc)

marks_exons <- read.table(opts$`<input1>`, colClasses = "character")

marks_exons_with_transposons <- read.table(opts$`<input2>`,
                                           colClasses = "character")

marks_exons_with_transposons_names <- unique(marks_exons_with_transposons$V4)
marks_exons_with_transposons <- marks_exons[marks_exons$V4 %in% marks_exons_with_transposons_names, ]

marks_exons_without_transposons <- marks_exons[!marks_exons$V4 %in% marks_exons_with_transposons_names, ]

write.table(marks_exons_without_transposons,
            paste(opts$O1),
            col.names = F,
            row.names = F,
            quote = F,
            sep = "\t")

write.table(marks_exons_with_transposons,
            paste(opts$O2),
            col.names = F,
            row.names = F,
            quote = F,
            sep = "\t")

marks_outside_exons <- read.table(opts$`<input3>`, colClasses = "character")

marks_outside_exons_with_transposons <- read.table(opts$`<input4>`, colClasses = "character")
marks_outside_exons_with_transposons_names <- unique(marks_outside_exons_with_transposons$V4)

non_coding_marks_without_transposons <- marks_outside_exons[!marks_outside_exons$V4 %in% marks_outside_exons_with_transposons_names, ]

non_coding_marks_with_transposons <- marks_outside_exons[marks_outside_exons$V4 %in% marks_outside_exons_with_transposons_names, ]

write.table(non_coding_marks_with_transposons,
            paste(opts$O3),
            col.names = F,
            row.names = F,
            quote = F,
            sep = "\t")

write.table(non_coding_marks_without_transposons,
            paste(opts$O4),
            col.names = F,
            row.names = F,
            quote = F,
            sep = "\t")
