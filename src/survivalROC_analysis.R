rm(list=ls())

library(survivalROC)

args <- commandArgs(trailingOnly = TRUE)
path <- args[1]
setwd(path)
data_fn <- args[2]

data <- read.csv(data_fn, sep = ";", check.names = FALSE)

ROC1 = survivalROC(Stime = data$duration,
                   status = data$observed,
                   marker = data$risk_score,
                   predict.time = 365,
                   method = "KM")
ROC1df <- data.frame(specificity = ROC1$FP, sensitivity = ROC1$TP)
write.table(ROC1df, file="data/TCGA-LGG-ROC-1-year.csv", sep = ';', quote = F, row.names = FALSE)

ROC3 = survivalROC(Stime = data$duration,
                   status = data$observed,
                   marker = data$risk_score,
                   predict.time = 1095,
                   method = "KM")
ROC3df <- data.frame(specificity = ROC3$FP, sensitivity = ROC3$TP)
write.table(ROC3df, file="data/TCGA-LGG-ROC-3-year.csv", sep = ';', quote = F, row.names = FALSE)

ROC5 = survivalROC(Stime = data$duration,
                   status = data$observed,
                   marker = data$risk_score,
                   predict.time = 1825,
                   method = "KM")
ROC5df <- data.frame(specificity = ROC5$FP, sensitivity = ROC5$TP)
write.table(ROC5df, file="data/TCGA-LGG-ROC-5-year.csv", sep = ';', quote = F, row.names = FALSE)

AUCdf <- data.frame(year = c(1, 3, 5), AUC = c(ROC1$AUC, ROC3$AUC, ROC5$AUC))
write.table(AUCdf, file="data/TCGA-LGG-AUC-1-3-5-years.csv", sep = ';', quote = F, row.names = FALSE)
