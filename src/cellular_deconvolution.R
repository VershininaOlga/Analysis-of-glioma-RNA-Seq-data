rm(list=ls())

library(EPIC)

args <- commandArgs(trailingOnly = TRUE)
path <- args[1]
setwd(path)
data_fn <- args[2]


datalog2TPMnorm <- as.matrix(read.csv(data_fn, sep = ";", 
                                      check.names = FALSE, row.names = 'Gene'))

dataTPM <- (2 ^ datalog2TPMnorm) - 1

epic_dec <- EPIC(dataTPM)
cell_frac <- epic_dec$cellFractions
res <- data.frame(data.frame("sample" = rownames(cell_frac)), cell_frac, check.names = FALSE)

write.table(res, file="data/TCGA-LGG-epic-deconvolution.csv",
            sep = ";", quote = F, row.names = FALSE)
