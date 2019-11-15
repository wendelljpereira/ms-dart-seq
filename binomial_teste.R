require(data.table)
require(reshape2)

meth_cpg <- as.data.frame(fread(snakemake@input[[1]], sep = "\t"))

# BS-seq data processing
cpg_info <- colsplit(string = meth_cpg[, 4],
                     pattern = "[(]|[%]|[/]|[']",
                     names = c("NA",
                               "n_of_meth",
                               "n_of_reads",
                               "meth_level_cpg1",
                               "d"))[c(2, 3)]

meth_cpg <- cbind(meth_cpg, cpg_info)

# Executes the binomial tests
meth_cpg_pvalues <- numeric()
for (i in 1:nrow(meth_cpg)) {
    meth_cpg_pvalues[i] <- binom.test(meth_cpg[i, 10],
                                      meth_cpg[i, 11],
                                      p = 0.06,
                                      alternative = "greater",
                                      conf.level = 0.95)$p.value
}

meth_cpg$pvalues <- meth_cpg_pvalues

# Corrects the FDR (cpg + CHG)
meth_cpg$fdr <- p.adjust(meth_cpg_pvalues, "fdr", n = nrow(meth_cpg))

meth_cpg <- meth_cpg[, c(1, 2, 3, 4, 13, 6)]

write.table(meth_cpg,
            snakemake@output[[1]],
            col.names = F,
            row.names = F,
            sep = "\t",
            quote = F)
