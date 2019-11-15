require(docopt)
require(gdata)

intersect_marks <- function(data_intersect = data_intersect,
                            clone_name = "name",
                            tissue = "name",
                            group = "name",
                            prefix = "file",
                            enzyme = "ms",
                            int_name = "name",
                            non_overlap=non_overlap,
                            non_overlap_is_file=non_overlap_is_file
                            ){

    clone_name <- as.character(c(clone_name))
    tissue <- as.character(c(tissue))
    group <- as.character(c(group))
    prefix <- as.character(prefix)
    enzime <- as.character(c(enzyme))
    int_name <- as.character(int_name)

    if (non_overlap_is_file == "TRUE"){
        non_overlap <- read.table(non_overlap)[, 5]
    }else if (non_overlap_is_file == "FALSE"){
        non_overlap <- non_overlap
    }

    clone <- data_intersect

    # Checks if is necessary to filter the marks of a group os samples.
    if (clone_name == "name"){
        print("all clones will be use!")
    }else{
        for (i in 1:length(clone_name)){
            y <- clone[, grep(pattern = clone_name[i], colnames(clone))]
            if (i == 1){
                x <- as.data.frame(y)
            }else{
                x <- cbindX(x, y)
            }
        }
        clone <- x
    }

    # Checks if is necessary to filter the marks of a tissue.
    if (tissue == "name"){
        print("all tissues will be use!")
    }else{
        for (i in 1:length(tissue)){
            y <- clone[, grep(pattern = tissue[i], colnames(clone))]
            if (i == 1){
                x <- as.data.frame(y)
            }else{
                x <- cbindX(x, y)
            }
        }
        clone <- x
    }

    if (group == "name"){
        print("all groups will be use!")
    }else{
        for (i in 1:length(group)){
            y <- clone[, substr(names(clone), 1, 2) == paste(group[i], "_", sep = "")]
            if (i == 1){
                x <- as.data.frame(y)
            }else{
                x <- cbindX(x, y)
            }
        }
        clone <- x
    }

    # Checks if is necessary to filter the marks by the enzymes.
    if (enzyme == "names"){
        print("all enzymes data will be use!")
    }else{
        for (i in 1:length(enzyme)){
            y <- clone[, grep(pattern = enzyme[i], colnames(clone))]
            if (i == 1){
                x <- as.data.frame(y)
            }else{
                x <- cbindX(x, y)
            }
        }
        clone <- x
    }

    x <- as.list(clone)

    # Removes NAs
    x <- lapply(x, function(x) x[!is.na(x)])

    # Determine the intersection of MSD-Sites among the filtered samples
    Int <- Reduce(intersect, x)
    I <- as.data.frame(Int)
    colnames(I) <- paste(int_name)

    if (file.exists(paste(prefix, "intersect_marks.txt", sep = "_")) == "FALSE"){

        write.table(I,
                    file = paste(prefix, "intersect_marks.txt", sep = "_"),
                    sep = "\t",
                    quote = F,
                    row.names = F,
                    col.names = T)

    } else if (file.exists(paste(prefix, "intersect_marks.txt", sep = "_")) == "TRUE"){

        Data_int_g <- read.table(paste(prefix, "intersect_marks.txt", sep = "_"), header = T, sep = "\t")
        Data_int_g <- cbindX(Data_int_g, as.data.frame(I))

        write.table(Data_int_g, file = paste(prefix, "intersect_marks.txt", sep = "_"), sep = "\t", quote = F, row.names = F, col.names = T)
    }
}

intersect_group_files <- snakemake@params[["grupos_intersect"]]

all_sites <- data.frame()
for (i in 1:length(intersect_group_files)){
    if (i == 1){

        sites <- read.table(intersect_group_files[i],
                                header = T,
                                sep = "\t",
                                quote = "\"",
                                check.names = F,
                                na.strings = "NA",
                                colClasses = "character")
        all_sites <- sites

    } else{

        sites <- read.table(intersect_group_files[i],
                                header = T,
                                sep = "\t",
                                quote = "\"",
                                check.names = F,
                                na.strings = "NA",
                                colClasses = "character")
        all_sites <- cbindX(all_sites, sites)

    }
}

for (i in 1:length(snakemake@params[["names"]])){

    if (i == 1){

        x <- as.data.frame(read.table(paste(snakemake@input[[i]]))[, 4])
        colnames(x) <- snakemake@params[["names"]][i]
        without_overlaps <- x

    }else{

        x <- as.data.frame(read.table(paste(snakemake@input[[i]]))[, 4])
        colnames(x) <- snakemake@params[["names"]][i]
        without_overlaps <- cbindX(without_overlaps, x)

    }
}

x <- as.list(without_overlaps)

## Removes NAs
x <- lapply(x, function(x) x[!is.na(x)])

without_overlaps <- unique(Reduce(intersect, x))

for (i in 1:length(snakemake@params[["names"]])){
    for (j in 1:length(snakemake@params[["tissue"]])){
        intersect_marks(data_intersect = all_sites,
                        clone_name = snakemake@params[["names"]],
                        tissue = "name",
                        prefix = snakemake@params[["prefix"]][1],
                        enzyme = snakemake@params[["enzyme"]][1],
                        int_name = paste(snakemake@params[["names"]][i],
                                         snakemake@params[["tissue"]][j],
                                         sep = "_"),
                        non_overlap = without_overlaps,
                        non_overlap_is_file = "FALSE")
    }
}
