"Usage: gff_to_gtf_biotype.R (--out1 <O1>) (--out2 <O2>) <input1> 
-h --help    show this
--out1   gft_file    GTF file of the GFF3
--out1   gft_file    GTF file with the biotype
<input1>   name1    GFF3 file
gff_to_gtf_biotype.R -h | --help  show this message
" -> doc

require(docopt)
require(rtracklayer)

# retrieve the command-line arguments
opts <- docopt(doc)
save.image()

## Convert gff to gtf
gff <- import(opts$`<input1>`)

export(gff, paste(opts$O1), "gtf")


## Add a biotype description for the genes of the Eucalyptus grandis genome. It is necessary to run AlFA(). Also, it is necessary to convert from GFF3 to GTF format.
gtf <- read.table(paste(opts$O1),
                  sep = "\t",
                  header = F)

# In the case of the E. grandis genome all feature in the gff3 are genes coding for proteins. So, it is possible to add the same biotype for all of them 
gtf$V9 <- paste(gtf$V9, " gene_biotype \"protein_coding\";", sep = "")

write.table(gtf,
            paste(opts$O2),
            sep = "\t",
            quote = F,
            col.names = F,
            row.names = F)

