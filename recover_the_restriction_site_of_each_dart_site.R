"Usage: recover_the_restriction_site_of_each_dart_site.R (--out1 <O1>) <input1>
-h --help    show this
--out1   name1    specify the name for the first output file
recover_the_restriction_site_of_each_dart_site.R -h | --help  show this message
" -> doc

# load the docopt library
require(docopt)
# retrieve the command-line arguments
opts <- docopt(doc)

brasuz_methylated_sites <- read.table(opts$`<input1>`, header = F, sep = "\t")
brasuz_methylated_sites <- brasuz_methylated_sites[brasuz_methylated_sites$V6 == "+", ]

brasuz_methylated_sites_mspI <- brasuz_methylated_sites
brasuz_methylated_sites_mspI$V2 <- brasuz_methylated_sites_mspI$V2 - 1
brasuz_methylated_sites_mspI$V3 <- brasuz_methylated_sites_mspI$V3 + 2

write.table(brasuz_methylated_sites_mspI,
            paste(opts$O1),
            col.names = F,
            row.names = F,
            sep = "\t",
            quote = F)
