"Usage: id_part2.R (--out1 <O1>) (--out2 <O2>) <input1> <input2>
-h --help    show this
--out1   name1    specify the name for the first output file
--out2   name2    specify the name for the second output file
id_part2.R -h | --help  show this message
" -> doc

# load the docopt library
require(docopt)
# retrieve the command-line arguments
opts <- docopt(doc)
save.image()

MSD <- read.table(opts$`<input1>`, colClasses = "character")

MSD_transposon_int <- read.table(opts$`<input2>`,
                                    sep = "\t",
                                    colClasses = "character")

MSD_TEs_int_names <- unique(MSD_transposon_int$V4)

# Generates the liste of MSD-Methylated sites in TEs
MSD_in_TEs <- MSD[MSD$V4 %in% MSD_TEs_int_names, ]

write.table(MSD_in_TEs,
            paste(opts$O1),
            sep = "\t",
            col.names = F,
            row.names = F,
            quote = F)

MSD_outside_genes_and_TEs <- MSD[!MSD$V4 %in% MSD_TEs_int_names, ]

write.table(MSD_outside_genes_and_TEs,
            paste(opts$O2),
            sep = "\t",
            col.names = F,
            row.names = F,
            quote = F)
