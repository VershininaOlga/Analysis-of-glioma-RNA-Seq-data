rm(list=ls())

library(dplyr)
library(sva)
library(DESeq2)

args <- commandArgs(trailingOnly = TRUE)
path <- args[1]
setwd(path)
dataset <- args[2]
countData_fn <- args[3]
geneLengths_fn <- args[4]


if ((dataset == "CGGA-LGG") | (dataset == "CGGA-GBM")){
  
  countDataComBatSeq_fn <- paste("data/", dataset, "-ComBatSeq-counts.csv", sep="")
  
  if (!file.exists(countDataComBatSeq_fn)){
    countData <- read.csv(countData_fn, sep = ";", check.names = FALSE)
    rownames(countData) <- countData$gene_name
    countData["gene_name"] <- NULL
    countData <- as.matrix(countData)
    if (dataset == "CGGA-LGG"){
      batch <- c(rep(1, 271), rep(2, 137))
    }
    if (dataset == "CGGA-GBM"){
      batch <- c(rep(1, 133), rep(2, 85))
    }
    adjusted <- ComBat_seq(countData, batch = batch, group = NULL)
    countData <- data.frame(data.frame("gene_name" = rownames(adjusted)),
                            adjusted, check.names = FALSE)
    write.table(countData, file = countDataComBatSeq_fn, sep = ";", quote = FALSE, row.names = FALSE)
  }
  
  countData_fn <- countDataComBatSeq_fn
}

print(countData_fn)

geneLengths <- read.csv(geneLengths_fn, sep = ";", check.names = FALSE)
countData <- read.csv(countData_fn, sep = ";", check.names = FALSE)
if (dataset == "TCGA-GBM"){
  countData["gene_id"] <- NULL
  countData["gene_type"] <- NULL
}
countData <- distinct(countData, gene_name, .keep_all = TRUE)
countData <- merge(countData, geneLengths, by.x = "gene_name", by.y = "gene_name")
geneLengths <- countData[c("gene_name", "feature_interval_length")]
countData["feature_interval_length"] <- NULL


# Normalization
rownames(countData) <- countData$gene_name
countData[1] <- NULL

if (all(rownames(countData) == geneLengths$gene_name))
{
  RPK <- sweep(countData, 1, geneLengths$feature_interval_length, "/") * 10^3
  TPM <- sweep(RPK, 2, colSums(RPK), "/") * 10^6
  log2TPM <- log2(TPM + 1)
  res <- data.frame(data.frame("gene_name" = rownames(countData)), log2TPM, 
                    check.names = FALSE)
  write.table(res, paste("data/", dataset, "-log2TPMnorm-data.csv", sep = ""),
              sep = ";", quote = FALSE, row.names = FALSE)
} else{
  stop("Gene names must be in the same order")
}
