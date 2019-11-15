require(GenomicRanges)
require(ggbio)
require(Gviz)
require(rtracklayer)
require(gdata)

correct_sites <- function(sites_bed_file = sites_bed_file, data = data){

    all_sites <- read.table(sites_bed_file, sep = "\t")
    all_sites <- unique(as.character(all_sites[, 4]))

    dup_sites <- all_sites[grep(pattern = "|", all_sites, fixed = T)]

    dup_sites_names <- strsplit(x = dup_sites, split = "|", fixed = T)
    dup_sites_names <- unique(c(do.call("rbind", dup_sites_names)))

    new_data_corrected <- data.frame()
    for (i in 1:ncol(data)){
        new_data_elements <- data[i]
        new_data_renamed_all <- data.frame()
        for (j in 1:nrow(new_data_elements)){

            if (as.character(new_data_elements[j, 1]) %in% dup_sites_names){

                x <- grep(pattern = paste0("\\<", new_data_elements[j, 1], "\\|"),
                          x = dup_sites)
                y <- grep(pattern = paste0("\\|", new_data_elements[j, 1], "\\>"),
                          x = dup_sites)

                if (length(x) > 0){
                    z <- x
                }else if (length(y) > 0){
                    z <- y
                }else{
                    print("Erro!")
                }

                new_data_renamed <- as.data.frame(dup_sites[z])
                colnames(new_data_renamed) <- "new_data_renamed"
                new_data_renamed_all <- rbind(new_data_renamed_all,
                                              new_data_renamed)
            }else{
                new_data_renamed <- as.data.frame(new_data_elements[j, 1])
                colnames(new_data_renamed) <- "new_data_renamed"
                new_data_renamed_all <- rbind(new_data_renamed_all,
                                              new_data_renamed)
            }

        }
        if (i == 1){
            colnames(new_data_renamed_all) <- colnames(data)[i]
            new_data_corrected <- as.data.frame(unique(new_data_renamed_all))

        }else{
            colnames(new_data_renamed_all) <- colnames(data)[i]
            new_data_corrected <- cbindX(new_data_corrected,
                                         as.data.frame(unique(new_data_renamed_all)))
        }
    }

    return(new_data_corrected)
}

if (file.exists("file_to_plot")) file.remove("file_to_plot")
if (file.exists("euk.txt")) file.remove("euk.txt")

#####################################
# Creates the file in Grange format #
#####################################

data <- read.table(snakemake@input[[1]], header = F)
colnames(data) <- c("chr", "start", "end", "id", "score", "strand")

data <- data[data$strand == "+", ]

data_intersect <- read.table(snakemake@input[[2]], header = T)[1]

intersect_corrected <- correct_sites(sites_bed_file = snakemake@input[[1]], data = data_intersect)

data <- data[data$id %in% intersect_corrected[, 1], ]

bed <- with(data, GRanges(chr, IRanges(start + 1, end), strand, id = id))
bed_merged <- reduce(bed)

## configure the chr sizes
vezes <- seq(1, 360)
ini <- -249999
fim <- 0
cromo <- sprintf("%02d", 1:11)
for (i in vezes) {
    ini <- ini + 250000
    fim <- fim + 250000
    result <- restrict(bed_merged, start = ini, end = fim)
    for (j in cromo) {
        atual <- paste("Chr", j, sep = "")
        scor <- length(grep(atual, result))
        line <- paste(atual, ini, fim, "*", scor, sep = "\t")
        write(line, file = "file_to_plot", append = TRUE)
    }
}

###############################
# Creates the plot using Gviz #
###############################
data2 <- read.table("file_to_plot", header = F)
colnames(data2) <- c("chr", "start", "end", "strand", "score")
bed2 <- with(data2, GRanges(chr, IRanges(start, end), strand, score))

####################
# Sets X axis size #
####################
gtrack <- GenomeAxisTrack()
displayPars(gtrack)$size <- 1.2
displayPars(gtrack)$cex <- 6
displayPars(gtrack)$col.line <- "black"
displayPars(gtrack)$col <- "black"

pos <- rep("0", 11)
chr <- c("Chr01",
         "Chr02",
         "Chr03",
         "Chr04",
         "Chr05",
         "Chr06",
         "Chr07",
         "Chr08",
         "Chr09",
         "Chr10",
         "Chr11")

len <- c(40297282,
         64237462,
         80088348,
         41978404,
         74731017,
         53893726,
         52447651,
         74330457,
         39019482,
         39359118,
         45510589)

chr_file <- data.frame(chr, pos, len)

write.table(chr_file,
            "euk.txt",
            col.names = F,
            row.names = F,
            quote = F,
            sep = "\t")

chr <- import.bed("euk.txt")
ncols <- 3
nrows <- 4
grid.newpage()
png(snakemake@output[[1]], width = 5000, height = 3000)

# Join the plots
pushViewport(viewport(layout = grid.layout(nrows, ncols)))
i <- 0
ii <- 0
k <- 1
for (j in cromo) {
    i <- i + 1
    ii <- ii + 1
    atual <- paste("Chr", j, sep = "")
    nome <- paste("Bands_distribution_in_", atual, sep = "")
    chr_idx <- levels(seqnames(chr)) == atual
    f <- 1
    t <- end(chr[chr_idx])
    dTrack <- DataTrack(bed2,
                        name = atual,
                        chromosome = atual,
                        ylim = c(0, 50),
                        fill.histogram = c("#EE7600", "#556B2F"),
                        cex.legend = 10)

    displayPars(dTrack)$cex.axis <- 4
    displayPars(dTrack)$fontsize.legend <- 4
    displayPars(dTrack)$cex <- 4
    displayPars(dTrack)$fontcolor.legend <- "black"
    pushViewport(viewport(layout.pos.col = ii, layout.pos.row = k))

    plotTracks(list(dTrack, gtrack),
               littleTicks = FALSE,
               add = TRUE,
               from = f,
               to = t,
               chromosome = atual,
               type = "histogram",
               col.axis = "black",
               background.panel = "transparent",
               showTitle = TRUE,
               fontcolor = "black",
               background.title = "transparent",
               main = atual,
               cex.main = 8,
               fontcolor.legend = "black",
               legend = FALSE)
    popViewport(1)
    if (i %% 3 == 0){
        k <- k + 1
        ii <- 0
    }
}
dev.off()
