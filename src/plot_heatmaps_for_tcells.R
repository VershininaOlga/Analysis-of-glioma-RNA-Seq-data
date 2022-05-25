rm(list=ls())

library("gplots")
library("ggplot2")
library(ComplexHeatmap)
library(latex2exp)
library(circlize)

args <- commandArgs(trailingOnly = TRUE)
path <- args[1]
setwd(path)


tcells = c('Th1', 'Th17', 'Treg', 'CTL')

for (i in 1:length(tcells))
{
  if (tcells[i] == 'Th1')
    center_genes = c('IFNG')
  if (tcells[i] == 'Th17')
    center_genes = c('IL17A', 'RORC')
  if (tcells[i] == 'CTL')
    center_genes = c('CD8B')
  if (tcells[i] == 'Treg')
    center_genes = c('CD2', 'CD3D', 'CD3G')
  
  corr_matrix = read.table(paste('data/TCGA-LGG-', tcells[i], '-corr-matrix.csv', sep=''),
                           sep=';', header=TRUE, row.names='X')
  genes = colnames(corr_matrix)
  
  sizes <- c()
  faces <- c()
  for (j in (1:length(genes)))
  {
    if (genes[j] %in% center_genes)
    {
      sizes[j] = 19
      faces[j] = 'bold'
    }
    else
    {
      sizes[j] = 16
      faces[j] = 'plain'
    }
  }
  
  col_fun = colorRamp2(c(-1.0, 0.0, 1.0), c("green", "black", "red"))
  lgd = Legend(col_fun = col_fun, title = "Pearson's coefficient", direction = "horizontal",
               title_position = "topcenter", title_gp = gpar(fontsize = 18, fontface = 'bold'),
               labels_gp = gpar(fontsize = 17),
               legend_width = unit(15, "cm"), grid_height = unit(1, "cm"))
  
  HM = Heatmap(corr_matrix, 
               col = circlize::colorRamp2(c(-1, -0.5, 0, 0.5, 1), c("green", "green3", "black", "red4", "red")),
               row_labels = genes, column_labels = genes,
               row_names_gp = gpar(fontsize = sizes, fontface = faces),
               column_names_gp = gpar(fontsize = sizes, fontface = faces),
               column_title = paste(tcells[i],'associated metagene', sep='-'), 
               column_title_gp = gpar(fontsize = 22, fontface = 'bold'),
               cluster_rows = T, cluster_columns = T,
               clustering_method_rows = "average", clustering_method_columns = "average",
               column_dend_height = unit(30, "mm"),
               row_dend_width = unit(30, "mm"),
               show_heatmap_legend = FALSE,
  ) 
  
  png(file=paste("pictures/corr_heatmap_for_", tcells[i], ".png", sep = ""), 
      width=2500, height=2500, res=300)
  draw(HM, auto_adjust = FALSE, 
       #heatmap_legend_side = "bottom", 
       #annotation_legend_list = lgd, 
       #annotation_legend_side = "bottom"
       )
  dev.off()
}
