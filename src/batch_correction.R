rm(list=ls())

library(dplyr)
library(sva)

args <- commandArgs(trailingOnly = TRUE)
path <- args[1]
setwd(path)
database <- args[2]
countData_fn <- args[3]

countData <- read.csv(countData_fn, sep = ";", check.names = FALSE)
countData <- distinct(countData, gene_name, .keep_all = TRUE)
countData <- na.omit(countData)
rownames(countData) <- countData$gene_name
countData["gene_name"] <- NULL
countData <- as.matrix(countData)

if (database == "GEO")
{
  batch <- c(rep(1, 18), rep(2, 24))
}
if (database == "CGGA")
{
  batch <- c(rep(1, 271), rep(2, 137))
}

adjusted <- ComBat_seq(countData, batch = batch, group = NULL)

write.table(data.frame(data.frame("gene_name" = rownames(adjusted)),
                       adjusted, check.names = FALSE),
            file = paste(unlist(strsplit(countData_fn, split = "-counts.csv", 
                                         fixed = TRUE))[1],
                         "-ComBatSeq-counts.csv", sep = ""),
            sep = ";", quote = FALSE, row.names = FALSE)
