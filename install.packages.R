## Promotes the installation of all R packages necessary to BRASUZ1 analysis ##

# Help function to install the R packages available in the CRAN
use_package <- function(p){
    if (!is.element(p, installed.packages()[, 1])){
        install.packages(p, dep = TRUE)
        require(p, character.only = TRUE)
    }else{
        require(p, character.only = TRUE)
    }
}

############
### CRAN ###
############

# List of necessary R packages from CRAN
packages <- c("alluvial",
              "berryFunctions",
              "corrplot",
              "data.table",
              "docopt",
              "doMC",
              "dplyr",
              "foreach",
              "gdata",
              "ggfortify",
              "ggplot2",
              "ggsci",
              "ggthemes",
              "googleVis",
              "gridExtra",
              "gridGraphics",
              "Matching",
              "plyr",
              "RColorBrewer",
              "reshape2",
              "sjstats",
              "splitstackshape",
              "stringr",
              "tidyr",
              "tidyverse",
              "VennDiagram")

# install and/or load the packages from CRAN
for (i in packages){
    use_package(i)
}

####################
### Bioconductor ###
####################

# List of necessary packages from bioconductor
bioc_packages <- c("AnnotationForge",
                   "biomaRt",
                   "clusterProfiler",
                   "DESeq2",
                   "edgeR",
                   "GenomicRanges",
                   "ggbio",
                   "Gviz",
                   "rtracklayer")

#Install and/or load the Bioconductor packages
for (p in bioc_packages){
    if (is.element(p, installed.packages()[, 1]) == FALSE){
        source("https://bioconductor.org/biocLite.R")
        biocLite(paste(p))
        require(p, character.only = TRUE)
    }else{
        require(p, character.only = TRUE)
    }
}
