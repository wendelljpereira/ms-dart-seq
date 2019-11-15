"Usage: id_part1.R (--out1 <O1>) (--out2 <O2>) <input1> <input2>
-h --help    show this
--out1   name1    specify the name for the first output file
--out2   name2    specify the name for the second output file
id_part1.R -h | --help  show this message
" -> doc

# loads the docopt library
require(docopt)
# retrieves the command-line arguments
opts <- docopt(doc)

# Reads the file with the genomic features information of the E. grandis genome
marcas <- read.table(opts$`<input1>`, colClasses = "character")

# Reads the sites that overlaps genes
marcas_genes_int <- read.table(opts$`<input2>`,
                               sep = "\t",
                               colClasses = "character")

marcas_genes_int_names <- unique(marcas_genes_int$V4)
marcas_em_genes <- marcas[marcas$V4 %in% marcas_genes_int_names, ]

write.table(marcas_em_genes,
            paste(opts$O1),
            sep = "\t",
            col.names = F,
            row.names = F,
            quote = F)

marcas_fora_dos_genes <- marcas[!marcas$V4 %in% marcas_genes_int_names, ]

write.table(marcas_fora_dos_genes,
            paste(opts$O2),
            sep = "\t",
            col.names = F,
            row.names = F,
            quote = F)
