require(docopt)
require(data.table)

# Create the function to return the mode of a vector
getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
}

files <- list.files(pattern="*.bed", full.names=T, recursive=FALSE)

system("mkdir good_clusters")
system("mkdir clusters_without_pstI")

for(i in 1:length(files)){

    # Read the file with all clusters formed by the combination of reads and restriction sites positions
    cluster_with_overlaps <- as.data.frame(fread(paste(files[i]), sep="\t", header = F))

    # Verity if there is a PstI site
    if("PstI" %in% unique(cluster_with_overlaps$V4)){

        if(unique(cluster_with_overlaps$V6) == "+"){

            pos_pstI <- max(cluster_with_overlaps[cluster_with_overlaps$V4 == "PstI", 2])

            new_positions <- cluster_with_overlaps[cluster_with_overlaps$V2 >= pos_pstI, ]
            biggest_frag <- new_positions[new_positions$V3 - new_positions$V2 == max(new_positions$V3 - new_positions$V2), ]
            biggest_frag_unique <- biggest_frag[1, ]

            if("PstI" %in% unique(new_positions$V4) & "MspI" %in% unique(new_positions$V4)){

                new_positions <- new_positions[new_positions$V4 %in% c("PstI", "MspI"), ]

                new_positions$V2 <- max(new_positions[new_positions$V4 == "PstI", 2])
                new_positions_frags <- new_positions[new_positions$V4 == "MspI", ]

                # Is the intervals PstI-MspI suported for at least one read?
                ## Test if each of the PstI-MspI intervals containg reads that start on the pstI site and end in the mspI site.
                new_positions_frags_checked <- data.frame()
                for(r in 1:nrow(new_positions_frags)){
                    if(any(cluster_with_overlaps$V2 %in% seq(new_positions_frags[r, 2], new_positions_frags[r, 2] + 4) & cluster_with_overlaps$V3 %in% seq(new_positions_frags[r, 3] -3, new_positions_frags[r, 3])) == TRUE){
                        new_positions_frags_checked <- rbind(new_positions_frags_checked, new_positions_frags[r, ])
                    }
                }

                # If any PstI-MspI frag is supported by reads, test if there is any reads bigger than this frag and, If true, add this read as a new feature.
                if(nrow(new_positions_frags_checked) >= 1){

                    # Test if the there is any frags bigger than the biggest PstI-MspI position.
                    if(biggest_frag_unique$V2 >= min(new_positions_frags$V2) && biggest_frag_unique$V3 <= max(new_positions_frags$V3)){

                    }else{
                        new_positions_frags_checked <- rbind(new_positions_frags_checked, biggest_frag_unique)
                    }
                }else{
                    new_positions_frags_checked <- biggest_frag_unique
                }
            }else{
                new_positions_frags_checked <- biggest_frag_unique
            }
            # Creates a positions that represents the difference between the overlap features. Will be use to normalize the counts.
            new_positions_frags_checked$V5 <- strsplit(files[i], "./")[[1]][[2]]

            if(nrow(new_positions_frags_checked) > 1){
                end_positions <- sort(new_positions_frags_checked$V3, decreasing=F)
                diff_pos <- data.frame()
                for(n in 1:(nrow(new_positions_frags_checked)-1)){
                    diff_pos_V1 = unique(new_positions_frags_checked$V1)
                    diff_pos_V2 <- end_positions[n]
                    diff_pos_V3 <- end_positions[n+1]
                    diff_pos_V4 <- unique(new_positions_frags_checked$V5)
                    diff_pos_V5 <- paste(strsplit(files[i], "./")[[1]][[2]], "diff_frag", sep = "_")
                    diff_pos_V6 <- unique(new_positions_frags_checked$V6)
                    diff_pos_V7 <- unique(new_positions_frags_checked$V7)
                    diff_pos <- rbind(diff_pos, data.frame(V1 = diff_pos_V1,
                                                           V2 = diff_pos_V2,
                                                           V3 = diff_pos_V3,
                                                           V4 = diff_pos_V4,
                                                           V5 = diff_pos_V5,
                                                           V6 = diff_pos_V6,
                                                           V7 = diff_pos_V7))
                }
                new_positions_frags_checked <- rbind(new_positions_frags_checked, diff_pos)
            }

            write.table(new_positions_frags_checked, file = paste('./',"good_clusters",'/', strsplit(files[i], "./")[[1]][[2]], sep = ""), col.names = F, row.names = F, quote = F, sep = "\t")

        }else if(unique(cluster_with_overlaps$V6) == "-"){
            pos_pstI <- min(cluster_with_overlaps[cluster_with_overlaps$V4 == "PstI", 3])

            new_positions <- cluster_with_overlaps[cluster_with_overlaps$V2 <= pos_pstI, ]
            biggest_frag <- new_positions[abs(new_positions$V3 - new_positions$V2) == max(abs(new_positions$V3 - new_positions$V2)), ]
            biggest_frag_unique <- biggest_frag[1, ]

            if("PstI" %in% unique(new_positions$V4) & "MspI" %in% unique(new_positions$V4)){

                new_positions <- new_positions[new_positions$V4 %in% c("PstI", "MspI"), ]

                new_positions$V3 <- min(new_positions[new_positions$V4 == "PstI", 3])
                new_positions_frags <- new_positions[new_positions$V4 == "MspI", ]

                # Is the intervals PstI-MspI suported for at least one read?
                ## Test if each of the PstI-MspI intervals containg reads that start on the pstI site and end in the mspI site.
                new_positions_frags_checked <- data.frame()
                for(r in 1:nrow(new_positions_frags)){
                    if(any(cluster_with_overlaps$V2 %in% seq(new_positions_frags[r, 2], new_positions_frags[r, 2] + 3) & cluster_with_overlaps$V3 %in% seq(new_positions_frags[r, 3] -4, new_positions_frags[r, 3])) == TRUE){
                        new_positions_frags_checked <- rbind(new_positions_frags_checked, new_positions_frags[r, ])
                    }
                }

                # If any PstI-MspI frag is supported by reads, test if there is any reads bigger than this frag and, If true, add this read as a new feature.
                if(nrow(new_positions_frags_checked) >= 1){

                    # Tests if the there is any frags bigger than the biggest PstI-MspI position.
                    if(biggest_frag_unique$V2 >= min(new_positions_frags$V2) && biggest_frag_unique$V3 <= max(new_positions_frags$V3)){

                    }else{
                        new_positions_frags_checked <- rbind(new_positions_frags_checked, biggest_frag_unique)
                    }
                }else{
                    new_positions_frags_checked <- biggest_frag_unique
                }
            }else{
                new_positions_frags_checked <- biggest_frag_unique
            }

            # Creates a positions that represents the difference between the overlap features. Will be use to normalize the counts.
            new_positions_frags_checked$V5 <- strsplit(files[i], "./")[[1]][[2]]

            if(nrow(new_positions_frags_checked) > 1){
                end_positions <- sort(new_positions_frags_checked$V2, decreasing=F)
                diff_pos <- data.frame()
                for(n in 1:(nrow(new_positions_frags_checked)-1)){
                    diff_pos_V1 = unique(new_positions_frags_checked$V1)
                    diff_pos_V2 <- end_positions[n]
                    diff_pos_V3 <- end_positions[n+1]
                    diff_pos_V4 <- unique(new_positions_frags_checked$V5)
                    diff_pos_V5 <- paste(strsplit(files[i], "./")[[1]][[2]], "diff_frag", sep = "_")
                    diff_pos_V6 <- unique(new_positions_frags_checked$V6)
                    diff_pos_V7 <- unique(new_positions_frags_checked$V7)
                    diff_pos <- rbind(diff_pos, data.frame(V1 = diff_pos_V1,
                                                           V2 = diff_pos_V2,
                                                           V3 = diff_pos_V3,
                                                           V4 = diff_pos_V4,
                                                           V5 = diff_pos_V5,
                                                           V6 = diff_pos_V6,
                                                           V7 = diff_pos_V7))
                }
                new_positions_frags_checked <- rbind(new_positions_frags_checked, diff_pos)
            }

            write.table(new_positions_frags_checked, file = paste('./',"good_clusters",'/', strsplit(files[i], "./")[[1]][[2]], sep = ""), col.names = F, row.names = F, quote = F, sep = "\t")
        }
    }else if(!"PstI" %in% unique(cluster_with_overlaps$V4)){

        # Is the position without PstI sites supported by at least 10 reads?
        if(nrow(cluster_with_overlaps) < 10){

        }else{
            if(unique(cluster_with_overlaps$V6) == "+"){
                pos_pstI <- max(getmode(cluster_with_overlaps[cluster_with_overlaps$V4 != "MspI", 2]))

                new_positions <- cluster_with_overlaps[cluster_with_overlaps$V2 >= pos_pstI, ]
                biggest_frag <- new_positions[new_positions$V3 - new_positions$V2 == max(new_positions$V3 - new_positions$V2), ]
                biggest_frag_unique <- biggest_frag[1, ]

                # There is a MspI site cover by the reads?
                if("MspI" %in% unique(new_positions$V4)){
                    new_positions <- new_positions[new_positions$V4 %in% "MspI", ]

                    # Necessary to deal with reads where the start position containg a MspI in the genome (caused by a 'SNP')
                    ## Select the MspI closest to the begin of the feature, test it the 'pstI position" is the same and ignore the mspI frag if true.
                    new_positions_mspI <- new_positions[new_positions$V2 == min(new_positions$V2), ]

                    if(nrow(new_positions) == 1){
                        if(pos_pstI %in% seq(new_positions_mspI[, 2], new_positions_mspI[, 3])){
                            new_positions_frags <- biggest_frag_unique
                        }else{
                            new_positions$V2 <- pos_pstI
                            new_positions_frags <- new_positions
                        }
                    }else if(nrow(new_positions) > 1){
                        if(pos_pstI %in% seq(new_positions_mspI[, 2], new_positions_mspI[, 3])){
                            new_positions <- new_positions[!new_positions$V2 == new_positions_mspI$V2, ]
                            new_positions$V2 <- pos_pstI
                            new_positions_frags <- new_positions
                        }
                    }

                    # Is the intervals PstI-MspI suported for at least one read?
                    ## Test if each of the PstI-MspI intervals containg reads that start on the pstI site and end in the mspI site.
                    new_positions_frags_checked <- data.frame()
                    for(r in 1:nrow(new_positions_frags)){
                        if(any(cluster_with_overlaps$V2 %in% seq(new_positions_frags[r, 2], new_positions_frags[r, 2] + 4) & cluster_with_overlaps$V3 %in% seq(new_positions_frags[r, 3] -3, new_positions_frags[r, 3])) == TRUE){
                            new_positions_frags_checked <- rbind(new_positions_frags_checked, new_positions_frags[r, ])
                        }
                    }

                    # If any PstI-MspI frag is supported by reads, test if there is any reads bigger than this frag and, If true, add this read as a new feature.
                    if(nrow(new_positions_frags_checked) >= 1){

                        # Test if the there is any frags bigger than the biggest PstI-MspI position.
                        if(biggest_frag_unique$V2 >= min(new_positions_frags$V2) && biggest_frag_unique$V3 <= max(new_positions_frags$V3)){

                        }else{
                            new_positions_frags_checked <- rbind(new_positions_frags_checked, biggest_frag_unique)
                        }
                    }else{
                        new_positions_frags_checked <- biggest_frag_unique
                    }
                }else{
                    new_positions_frags_checked <- biggest_frag_unique
                }

                # Create a positions that represents the difference between the overlap features. Will be use to normalize the counts.
                new_positions_frags_checked$V5 <- strsplit(files[i], "./")[[1]][[2]]

                if(nrow(new_positions_frags_checked) > 1){
                    end_positions <- sort(new_positions_frags_checked$V3, decreasing=F)
                    diff_pos <- data.frame()
                    for(n in 1:(nrow(new_positions_frags_checked)-1)){
                        diff_pos_V1 = unique(new_positions_frags_checked$V1)
                        diff_pos_V2 <- end_positions[n]
                        diff_pos_V3 <- end_positions[n+1]
                        diff_pos_V4 <- unique(new_positions_frags_checked$V5)
                        diff_pos_V5 <- paste(strsplit(files[i], "./")[[1]][[2]], "diff_frag", sep = "_")
                        diff_pos_V6 <- unique(new_positions_frags_checked$V6)
                        diff_pos_V7 <- unique(new_positions_frags_checked$V7)
                        diff_pos <- rbind(diff_pos, data.frame(V1 = diff_pos_V1,
                                                               V2 = diff_pos_V2,
                                                               V3 = diff_pos_V3,
                                                               V4 = diff_pos_V4,
                                                               V5 = diff_pos_V5,
                                                               V6 = diff_pos_V6,
                                                               V7 = diff_pos_V7))
                    }
                    new_positions_frags_checked <- rbind(new_positions_frags_checked, diff_pos)
                }

                write.table(new_positions_frags_checked, file = paste('./',"clusters_without_pstI",'/', strsplit(files[i], "./")[[1]][[2]], sep = ""), col.names = F, row.names = F, quote = F, sep = "\t")

            }else if(unique(cluster_with_overlaps$V6) == "-"){
                pos_pstI <- min(getmode(cluster_with_overlaps[cluster_with_overlaps$V4 != "MspI", 3]))

                new_positions <- cluster_with_overlaps[cluster_with_overlaps$V2 <= pos_pstI, ]
                biggest_frag <- new_positions[abs(new_positions$V3 - new_positions$V2) == max(abs(new_positions$V3 - new_positions$V2)), ]
                biggest_frag_unique <- biggest_frag[1, ]

                # There is a MspI site cover by the reads?
                if("MspI" %in% unique(new_positions$V4)){
                    new_positions <- new_positions[new_positions$V4 %in% "MspI", ]

                    # Necessary to deal with reads where the start position containg a MspI in the genome (caused by a 'SNP')
                    ## Select the MspI closest to the end (+) of the feature, test it the 'pstI position" is the same and ignore the mspI frag if true.
                    new_positions_mspI <- new_positions[new_positions$V3 == max(new_positions$V3), ]

                    if(nrow(new_positions) == 1){
                        if(pos_pstI %in% seq(new_positions_mspI[, 2], new_positions_mspI[, 3])){
                            new_positions_frags <- biggest_frag_unique
                        }else{
                            new_positions$V3 <- pos_pstI
                            new_positions_frags <- new_positions
                        }
                    }else if(nrow(new_positions) > 1){
                        if(pos_pstI %in% seq(new_positions_mspI[, 2], new_positions_mspI[, 3])){
                            new_positions <- new_positions[!new_positions$V3 == new_positions_mspI$V3, ]
                            new_positions$V3 <- pos_pstI
                            new_positions_frags <- new_positions
                        }
                    }

                    # Is the intervals PstI-MspI suported for at least one read?
                    ## Test if each of the PstI-MspI intervals containg reads that start on the pstI site and end in the mspI site.
                    new_positions_frags_checked <- data.frame()
                    for(r in 1:nrow(new_positions_frags)){
                        if(any(cluster_with_overlaps$V2 %in% seq(new_positions_frags[r, 2], new_positions_frags[r, 2] + 3) & cluster_with_overlaps$V3 %in% seq(new_positions_frags[r, 3] -4, new_positions_frags[r, 3])) == TRUE){
                            new_positions_frags_checked <- rbind(new_positions_frags_checked, new_positions_frags[r, ])
                        }
                    }

                    # If any PstI-MspI frag is supported by reads, test if there is any reads bigger than this frag and, If true, add this read as a new feature.
                    if(nrow(new_positions_frags_checked) >= 1){

                        # Test if the there is any frags bigger than the biggest PstI-MspI position.
                        if(biggest_frag_unique$V2 >= min(new_positions_frags$V2) && biggest_frag_unique$V3 <= max(new_positions_frags$V3)){

                        }else{
                            new_positions_frags_checked <- rbind(new_positions_frags_checked, biggest_frag_unique)
                        }
                    }else{
                        new_positions_frags_checked <- biggest_frag_unique
                    }
                }else{
                    new_positions_frags_checked <- biggest_frag_unique
                }

                # Create a positions that represents the difference between the overlap features. Will be use to normalize the counts.
                new_positions_frags_checked$V5 <- strsplit(files[i], "./")[[1]][[2]]

                if(nrow(new_positions_frags_checked) > 1){
                    end_positions <- sort(new_positions_frags_checked$V2, decreasing=F)
                    diff_pos <- data.frame()
                    for(n in 1:(nrow(new_positions_frags_checked)-1)){
                        diff_pos_V1 = unique(new_positions_frags_checked$V1)
                        diff_pos_V2 <- end_positions[n]
                        diff_pos_V3 <- end_positions[n+1]
                        diff_pos_V4 <- unique(new_positions_frags_checked$V5)
                        diff_pos_V5 <- paste(strsplit(files[i], "./")[[1]][[2]], "diff_frag", sep = "_")
                        diff_pos_V6 <- unique(new_positions_frags_checked$V6)
                        diff_pos_V7 <- unique(new_positions_frags_checked$V7)
                        diff_pos <- rbind(diff_pos, data.frame(V1 = diff_pos_V1,
                                                               V2 = diff_pos_V2,
                                                               V3 = diff_pos_V3,
                                                               V4 = diff_pos_V4,
                                                               V5 = diff_pos_V5,
                                                               V6 = diff_pos_V6,
                                                               V7 = diff_pos_V7))
                    }
                    new_positions_frags_checked <- rbind(new_positions_frags_checked, diff_pos)
                }

                write.table(new_positions_frags_checked, file = paste('./',"clusters_without_pstI",'/', strsplit(files[i], "./")[[1]][[2]], sep = ""), col.names = F, row.names = F, quote = F, sep = "\t")
            }
        }
    }
}
