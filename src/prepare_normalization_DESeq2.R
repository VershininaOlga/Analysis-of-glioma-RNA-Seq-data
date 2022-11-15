rm(list=ls())

library(dplyr)
library(tidyr)
library(DESeq2)
library(apeglm)
library(sva)
library(GenomicFeatures)


# Samples in the data file and the attribute file must be in the same order

# Example of data (counts matrix) file:
# gene_name;C1;C2;E1;E2;...    
# A1BG;12;25;57;77;...
# A1CF;176;201;15;28;...
# ...

# Example of sample attributes file:
# sample;condition;[other columns]
# C1;control
# C2;control
# E1;experiment
# E2;experiment
# ...


read_data <- function(countData_fn, attributes_fn){
  countData <- read.csv(countData_fn, sep = ";", check.names = FALSE)
  countData <- countData[order(countData$gene_name),]
  countData <- na.omit(countData)
  
  attributes <- read.csv(attributes_fn,  sep = ";", check.names = FALSE)
  attributes$condition <- factor(attributes$condition)
  
  return(list(countData = countData, attributes = attributes))
}


save_data <- function(data, data_fn){
  write.table(data, file = data_fn, sep = ";", quote = FALSE, row.names = FALSE)
}


get_protein_coding_genes_1 <- function(countData, attributes, annotation_fn){
  annotation <- read.csv(annotation_fn, row.names = NULL, sep = "\t")
  colnames(annotation)[1] <- "feature"
  annotation <- annotation[annotation$feature == "gene",]
  annotation <- annotation[annotation$class == "protein_coding",]
  
  #leaving the last record
  while (any(duplicated(annotation$symbol))){
    annotation <- annotation %>% group_by(symbol) %>% filter(duplicated(symbol) | n()==1)
  }
  
  n_subjects <- ncol(countData)
  countData <- merge(countData, annotation, by.x = "gene_name", by.y = "symbol")
  geneLengths <- countData[c("gene_name", "feature_interval_length")]
  countData <- countData[,0:n_subjects]
  
  return(list(countData = countData, attributes = attributes, lengths = geneLengths))
}


get_protein_coding_genes_2 <- function(countData, attributes, annotation_fn){
  txdb <- makeTxDbFromGFF(annotation_fn, format = "gtf")
  exonsListPerGene <- exonsBy(txdb, by = "gene")
  geneLengths <- as.data.frame(sum(width(reduce(exonsListPerGene))))
  colnames(geneLengths)[1] <- "feature_interval_length"
  geneLengths["gene_id"] <- rownames(geneLengths)
  
  countData <- countData[countData$gene_type == "protein_coding",]
  countData <- distinct(countData, gene_name, .keep_all = TRUE) # delete dupl genes
  countData <- merge(countData, geneLengths, by.x = "gene_id", by.y = "gene_id")
  countData <- countData[order(countData$gene_name),]
  countData["gene_id"] <- NULL
  countData["gene_type"] <- NULL
  
  geneLengths <- countData[c("gene_name", "feature_interval_length")]
  countData["feature_interval_length"] <- NULL
  
  save_data(geneLengths, "data/gencodev36-protein-coding-gene-lengths.csv")
  
  return(list(countData = countData, attributes = attributes, lengths = geneLengths))
}


prepare_data_by_groups <- function(countData, attributes, groups){
  rownames(countData) <- countData$gene_name
  countData[1] <- NULL
  if (all(colnames(countData) == attributes$sample)){
    group1 <- attributes[attributes$condition == groups[1],]
    group2 <- attributes[attributes$condition == groups[2],]
    attributesGroups <- rbind(group1, group2)
    countDataGroups <- cbind(countData[, group1$sample], countData[, group2$sample])
    return(list(countData = countDataGroups, attributes = attributesGroups))
  }
  else{
    stop("Sample names must be in the same order")
  }
}


run_diff_expression <- function(countData, attributes, groups, refGroup, paired){
  
  data <- prepare_data_by_groups(countData, attributes, groups)
  countDataGroups <- data$countData
  attributesGroups <- data$attributes
  
  if (all(colnames(countDataGroups) == attributesGroups$sample)){
    if (paired == TRUE){
      attributesGroups$subject <- factor(attributesGroups$subject)
      dds <- DESeqDataSetFromMatrix(countData = countDataGroups,
                                    colData = attributesGroups,
                                    design = ~ subject + condition)
    }
    else{
      dds <- DESeqDataSetFromMatrix(countData = countDataGroups,
                                    colData = attributesGroups,
                                    design = ~ condition)
    }
    
    dds$condition <- relevel(dds$condition, ref = refGroup)
    dds <- DESeq(dds)
    print(paste("condition_", groups[1], "_vs_", groups[2], sep=""))
    resDEseq2 <- lfcShrink(dds, coef = paste("condition_", groups[1], "_vs_",
                                             groups[2], sep = ""), type = "apeglm")
    resDEseq2 <- as.data.frame(resDEseq2)
    resDEseq2 <- data.frame(data.frame("gene_name" = rownames(resDEseq2)), resDEseq2)
    return(resDEseq2) 
  }
  else{
    stop("Sample names must be in the same order")
  }
}


normalization <- function(countData, attributes, geneLengths, normalizationType){
  rownames(countData) <- countData$gene_name
  countData[1] <- NULL
  
  if (normalizationType == "deseq2"){
    if (all(colnames(countData) == attributes$sample)){
      dds <- DESeqDataSetFromMatrix(countData = countData, 
                                    colData = attributes, 
                                    design = ~ 1)
      dds <- estimateSizeFactors(dds)
      normalizedCountData <- counts(dds, normalized = TRUE)
      log2NormalizedCountData <- log2(normalizedCountData + 1)
      res <- data.frame(data.frame("gene_name" = rownames(countData)), 
                        log2NormalizedCountData, check.names = FALSE)
      return(res)
    }
    else{
      stop("Sample names must be in the same order")
    }
  }
  
  if (normalizationType == "tpm"){
    if (all(rownames(countData) == geneLengths$gene_name)){
      RPK <- sweep(countData, 1, geneLengths$feature_interval_length, "/") * 10^3
      TPM <- sweep(RPK, 2, colSums(RPK), "/") * 10^6
      log2TPM <- log2(TPM + 1)
      res <- data.frame(data.frame("gene_name" = rownames(countData)), log2TPM, 
                        check.names = FALSE)
      return(res)
    }
    else{
      stop("Gene names must be in the same order")
    }
  }
}


pipeline1 <- function(countData_fn, attributes_fn, annotation_fn){
  # 1 - Read data
  data <- read_data(countData_fn, attributes_fn)
  
  # 2 - Isolation of protein coding genes
  dataProt <- get_protein_coding_genes_1(data$countData, data$attributes, annotation_fn)
  save_data(dataProt$countData, paste(unlist(strsplit(countData_fn, split = ".csv", 
                                                      fixed = TRUE))[1],
                                      "-protein-coding-genes.csv", sep = ""))
  rm(data)
  
  # 3 - Normalization
  #countDataNorm <- normalization(dataProt$countData, dataProt$attributes, 
                                 #dataProt$lengths, "deseq2")
  #save_data(countDataNorm, paste(unlist(strsplit(countData_fn, split = "-counts.csv", 
                                                 #fixed = TRUE))[1],
                                 #"-log2DESeq2norm-data.csv", sep = ""))

  countDataNorm <- normalization(dataProt$countData, dataProt$attributes, 
                                 dataProt$lengths, "tpm")
  save_data(countDataNorm, paste(unlist(strsplit(countData_fn, split = "-counts.csv", 
                                                 fixed = TRUE))[1],
                                 "-log2TPMnorm-data.csv", sep = ""))
  
  # 4 - Differential Expression Analysis
  DCPSvsDC <- run_diff_expression(dataProt$countData, dataProt$attributes, 
                                  c("DCPS", "DC"), "DC", TRUE)
  save_data(DCPSvsDC, "data/DESeq2res-DCPS-vs-DC.csv")
  DCMTXvsDC <- run_diff_expression(dataProt$countData, dataProt$attributes, 
                                   c("DCMTX", "DC"), "DC", TRUE)
  save_data(DCMTXvsDC, "data/DESeq2res-DCMTX-vs-DC.csv")
}


pipeline2 <- function(countData_fn, attributes_fn, annotation_fn,
                      countDataControls_fn, attributesControls_fn){
  # 1 - Read data
  data <- read_data(countData_fn, attributes_fn)
  
  # 2 - Isolation of protein coding genes
  dataProt <- get_protein_coding_genes_2(data$countData, data$attributes, annotation_fn)
  save_data(dataProt$countData, paste(unlist(strsplit(countData_fn, split = ".csv", 
                                                      fixed = TRUE))[1],
                                      "-protein-coding-genes.csv", sep = ""))
  rm(data)
  
  # 3 - Normalization
  #countDataNorm <- normalization(dataProt$countData, dataProt$attributes, 
                                 #dataProt$lengths, "deseq2")
  #save_data(countDataNorm, paste(unlist(strsplit(countData_fn, split = "-counts.csv", 
                                                 #fixed = TRUE))[1],
                                 #"-log2DESeq2norm-data.csv", sep = ""))

  countDataNorm <- normalization(dataProt$countData, dataProt$attributes, 
                                 dataProt$lengths, "tpm")
  save_data(countDataNorm, paste(unlist(strsplit(countData_fn, split = "-counts.csv", 
                                                 fixed = TRUE))[1],
                                 "-log2TPMnorm-data.csv", sep = ""))
  
  # 4 - Differential Expression Analysis
  dataControls <- read_data(countDataControls_fn, attributesControls_fn)
  countDataComb <- merge(dataProt$countData, dataControls$countData,
                         by.x = "gene_name", by.y = "gene_name")
  attrComb <- rbind(dataProt$attributes, dataControls$attributes)
  rm(dataControls)
  
  DeadvsAlive <- run_diff_expression(countDataComb, attrComb, 
                                     c("Dead", "Alive"), "Alive", FALSE)
  save_data(DeadvsAlive, "data/DESeq2res-Dead-vs-Alive.csv")
  AlivevsControl <- run_diff_expression(countDataComb, attrComb, 
                                        c("Alive", "Control"), "Control", FALSE)
  save_data(AlivevsControl, "data/DESeq2res-Alive-vs-Control.csv")
  DeadvsControl <- run_diff_expression(countDataComb, attrComb, 
                                       c("Dead", "Control"), "Control", FALSE)
  save_data(DeadvsControl, "data/DESeq2res-Dead-vs-Control.csv")
  
}


args <- commandArgs(trailingOnly = TRUE)
path <- args[1]
setwd(path)
countData_fn <- args[2]
attributes_fn <- args[3]
annotation_fn  <- args[4]
type_exp <- as.numeric((args[5]))
countDataControls_fn  <- args[6]
attributesControls_fn  <- args[7]

if (type_exp == 1){
  pipeline1(countData_fn, attributes_fn, annotation_fn)
}
if (type_exp == 2){
  pipeline2(countData_fn, attributes_fn, annotation_fn, 
            countDataControls_fn, attributesControls_fn)
}
