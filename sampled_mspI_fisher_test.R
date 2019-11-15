"Usage: sampled_mspI_fisher_test.R (--out1 <O1>) <input1>
-h --help    show this
--out1   bed_file    Position of the MspI sites in the plus strand
<input1>   bed_file    Position of the restriction sites of MspI and PstI
sampled_mspI_fisher_test.R -h | --help  show this message
" -> doc

require(docopt)
require(data.table)

# retrieves the command-line arguments
opts <- docopt(doc)

# Reads the set of restriction sites
genome_sites <- as.data.frame(fread(opts$`<input1>`))

# Keeps only the MspI sites in the "+" strand
genome_sites_mspI <- genome_sites[genome_sites$V4 == "MspI" & genome_sites$V6 == "+", ]

# writes a file with the MspI sites in the plus strand
write.table(genome_sites_mspI,
            paste(opts$O1),
            col.names = F,
            row.names = F,
            quote = F,
            sep = "\t")
