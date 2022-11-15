rm(list=ls())

library(EPIC)

args <- commandArgs(trailingOnly = TRUE)
database <- args[1]
path <- args[2]
setwd(path)
tpmData_fn <- args[3]


datalog2TPMnorm <- read.csv(tpmData_fn, sep = ";", check.names = FALSE, row.names = "gene_name")
datalog2TPMnorm <- as.matrix(datalog2TPMnorm)
dataTPM <- (2 ^ datalog2TPMnorm) - 1


epic_dec <- EPIC(dataTPM)
cell_frac <- epic_dec$cellFractions
res <- data.frame(data.frame("sample" = rownames(cell_frac)), cell_frac, check.names = FALSE)

write.table(res, file=paste("data/", database, "-epic-deconvolution.csv", sep=""),
            sep = ";", quote = F, row.names = FALSE)
