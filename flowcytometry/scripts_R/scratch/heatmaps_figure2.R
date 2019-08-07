library(gtable) ; library(grid)

heatmaps.tphe <- list(CD4.TPHE = pheatmap(scale(heat$CD4.TPHE[, sort(colnames(heat$CD4.TPHE))]), 
                                          color = viridis(color.cut),
                                          annotation_row = data.frame(row.names = rownames(heat$CD4.TPHE), Dataset = c(rep("CD4.TPHE", nrow(heat$CD4.TPHE)))),
                                          annotation_colors = list(Dataset = c(CD4.TPHE = "orange")),
                                          annotation_names_row = F,
                                          scale = "row", 
                                          cluster_rows = F, 
                                          cluster_cols = F, 
                                          show_colnames = F),
                      CD8.TPHE = pheatmap(scale(heat$CD8.TPHE[, sort(colnames(heat$CD4.TPHE))]),
                                          color = viridis(color.cut),
                                          annotation_row = data.frame(row.names = rownames(heat$CD8.TPHE), Dataset = c(rep("CD8.TPHE", nrow(heat$CD8.TPHE)))),
                                          annotation_colors = list(Dataset = c(CD8.TPHE = "blue")),
                                          annotation_names_row = F,
                                          scale = "row", 
                                          cluster_rows = F, 
                                          cluster_cols = F,
                                          fontsize_col = 14))

plot.grob <- p$gtable$grob[[1]]
xlab.grob <- p$gtable$grob[[2]]  
ylab.grob <- p$gtable$grob[[3]]  
legend.grob <- p$gtable$grob[[4]]  