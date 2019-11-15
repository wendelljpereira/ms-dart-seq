"Usage: split_clusters.R (--out1 <O1>) <input1>
-h --help    show this help information
<input1>    input1  file with the clusters of reads and PstI or MspI sites
--out1      name1   File with the clusters that intersect at least one of the restriction sites
cns_atac_processing.R -h | --help  show this message
" -> doc

# loads the docopt library
require(docopt)
require(data.table)

# retrieves the command-line arguments
opts <- docopt(doc)

clusters <- as.data.frame(fread(opts$`<input1>`,
                                sep="\t",
                                header = F))

# Extracts clusters with one element and without any restriction site
clusters_2_or_more_seqs_names <- unique(clusters[duplicated(clusters$V7) == TRUE, 7])

## Extracts the useful clusters ##
clusters_2_or_more_seqs <- clusters[clusters$V7 %in% clusters_2_or_more_seqs_names, ]

write.table(clusters_2_or_more_seqs,
            opts$O1,
            sep = "\t",
            col.names = F,
            row.names = F,
            quote = F)
