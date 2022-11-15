rm(list=ls())

library("gplots")
library("ggplot2")
library(ComplexHeatmap)
library(latex2exp)
library(circlize)

args <- commandArgs(trailingOnly = TRUE)
path <- args[1]
setwd(path)
corr_PS_fn <- args[2]
corr_MTX_fn <- args[3]


corr_PS = read.table(corr_PS_fn, sep = ';', header = TRUE, row.names = 'X')
corr_MTX = read.table(corr_MTX_fn, sep = ';', header = TRUE, row.names = 'X')
genes = colnames(corr_PS)


col_fun = colorRamp2(c(-1.0, 0.0, 1.0), c("green", "black", "red"))
lgd = Legend(col_fun = col_fun, title = "Pearson's coefficient", direction = "horizontal",
             title_position = "topcenter", title_gp = gpar(fontsize = 17),
             labels_gp = gpar(fontsize = 17),
             legend_width = unit(10, "cm"), grid_height = unit(0.7, "cm"))

HM_PS = Heatmap(corr_PS, 
                col = greenred(100),
                row_labels = genes, column_labels = genes,
                row_names_gp = gpar(fontsize = 10, fontface = "italic"),
                column_names_gp = gpar(fontsize = 10, fontface = "italic"),
                cluster_rows = F, cluster_columns = F,
                show_heatmap_legend = FALSE
) 

HM_MTX = Heatmap(corr_MTX, 
                 col = greenred(100),
                 row_labels = genes, column_labels = genes,
                 row_names_gp = gpar(fontsize = 10, fontface = "italic"),
                 column_names_gp = gpar(fontsize = 10, fontface = "italic"),
                 cluster_rows = F, cluster_columns = F,
                 show_heatmap_legend = FALSE
)

ht_list = HM_PS + HM_MTX
png(file = "pictures/corr_heatmap_PS_MTX_61_marker_genes.png",
    width = 5000, height = 2700, res = 300)
draw(ht_list, annotation_legend_list = lgd, 
     heatmap_legend_side = "bottom", annotation_legend_side = "bottom", 
     auto_adjust = FALSE)
dev.off()
