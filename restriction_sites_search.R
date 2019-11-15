"Usage: restriction_sites_search.R (--genome <genome>) (--out1 <O1>) (--out2 <O2>) <input1>
-h --help    show this help
--genome    genome  Name of the BSgenome package with the reference genome
--out1  name1   fasta file with the restriction site sequence of the enzymes
--out2  name2   bed file with the position of all restriction sites in the genome
input1  input1  bed file with the position of all mspI restriction sites in the genome
restriction_sites_search.R -h | --help  show this message
" -> doc

# load the docopt library
require(docopt)
require(tidyverse)

# retrieve the command-line arguments
opts <- docopt(doc)

# load the genome
require(opts$`<genome>`,
        character.only = T)

genome <- get(opts$`<genome>`)
genome

# Loads the functions
find_restriction_sites <- function(i,
                                   restriction_site_list=restriction_site_list,
                                   target_sequence=target_sequence,
                                   target_sequence_name = target_sequence_name){
  enzime_name <- names(restriction_site_list)[i]
  restriction_site <- restriction_site_list[[i]]

  # search in plus strand
  plus_matches <- matchPattern(restriction_site, target_sequence)
  rep.int(enzime_name, length(plus_matches))

  data.frame(rep(target_sequence_name, length(plus_matches)),
             start(plus_matches) - 1,
             end(plus_matches),
             rep(enzime_name, length(plus_matches)),
             rep(0, length(plus_matches)),
             rep("+", length(plus_matches)))
}

restrict_site_search <- function(j,
                                 seqnames = seqnames,
                                 restriction_site_list = restriction_site_list,
                                 genome = genome){
  target_sequence <- genome[[seqnames[j]]]
  target_sequence_name <- seqnames[j]
  cat("> Finding all restriction sites in", seqnames[j], "\n")

  # Execute the function find_restriction_sites for each enzime.
  sites_plus <- lapply(1:length(restriction_site_list),
                       find_restriction_sites,
                       restriction_site_list,
                       target_sequence,
                       target_sequence_name)

  sites_plus_df <- as.data.frame(do.call(rbind, sites_plus))
}

# fasta file with the sequence of the restrictions sites
restriction_site_list <- readDNAStringSet(opts$`<input1>`, "fasta")
restriction_site_list

# Reads the names of chromosomes and scafoolds and put on a vector
seqnames <- seqnames(genome)

# Executes the function to find all restriction sites
sites <- lapply(1:length(seqnames),
                restrict_site_search,
                seqnames,
                restriction_site_list,
                genome)

sites <- as.data.frame(do.call(rbind, sites))
colnames(sites) <- c("V1", "V2", "V3", "V4", "V5", "V6")

# Writes the bed file containing the restriction sites.
## sort
sites <- sites[with(sites, order(V1, V2)), ]

write.table(sites,
            paste(opts$O1),
            col.names = F,
            row.names = F,
            quote = F,
            sep = "\t")

# Save a bed file with the MspI positions
mspI_sites <- sites[sites$V4 == "MspI", ]

write.table(mspI_sites,
            paste(opts$O2),
            col.names = F,
            row.names = F,
            quote = F,
            sep = "\t")
