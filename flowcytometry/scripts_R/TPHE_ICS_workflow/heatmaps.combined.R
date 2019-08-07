library(gridExtra) ; library(pheatmap) ; library(viridis)

#### unify heatmaps for annotation key

## read in flowSOM cluster MFI .csvs
heat.data <- list.files("./results_R/fSOMs", full.names = TRUE, recursive = TRUE, pattern = "MetaCluster_MFI.csv")

heat <- list(CD4.TPHE = data.frame(read.csv(grep("TPHE_fSOM.CD4", heat.data, value = T), row.names = 1)),
             CD4.ICS = data.frame(read.csv(grep("ICS_fSOM.CD4", heat.data, value = T), row.names = 1)),
             CD8.TPHE = data.frame(read.csv(grep("TPHE_fSOM.CD8", heat.data, value = T), row.names = 1)),
             CD8.ICS = data.frame(read.csv(grep("ICS_fSOM.CD8", heat.data, value = T), row.names = 1))
)

## for legibility, just use meta cluster number
heat <- sapply(names(heat), function(i) {
  rownames(heat[[i]]) <- sub("Meta.Cluster_", "", rownames(heat[[i]]))
  return(heat[[i]])
})

#OR

## generate unique rownames
heat <- sapply(names(heat), function(i) {
  rownames(heat[[i]]) <- sub("Meta.Cluster_", "mc.", paste0(i, "_", rownames(heat[[i]])))
  return(heat[[i]])
})



## match columns; per panel -- TPHE and ICS
heat$CD4.TPHE[, colnames(heat$CD8.TPHE)[!(colnames(heat$CD8.TPHE) %in% colnames(heat$CD4.TPHE))]] <- NA
heat$CD8.TPHE[, colnames(heat$CD4.TPHE)[!(colnames(heat$CD4.TPHE) %in% colnames(heat$CD8.TPHE))]] <- NA

heat$CD4.ICS[, colnames(heat$CD8.ICS)[!(colnames(heat$CD8.ICS) %in% colnames(heat$CD4.ICS))]] <- NA
heat$CD8.ICS[, colnames(heat$CD4.ICS)[!(colnames(heat$CD4.ICS) %in% colnames(heat$CD8.ICS))]] <- NA

## TPHE markers as cols, clusters as rows; vertical stack using grid.arrange
color.cut <- 4 # parametize color scale in final heatmap

heatmaps.tphe <- list(CD4.TPHE = pheatmap(scale(heat$CD4.TPHE[, sort(colnames(heat$CD4.TPHE))]), 
                                     color = viridis(color.cut),
                                     annotation_row = data.frame(row.names = rownames(heat$CD4.TPHE), Dataset = c(rep("CD4.TPHE", nrow(heat$CD4.TPHE)))),
                                     annotation_colors = list(Dataset = c(CD4.TPHE = "blue")),
                                     annotation_names_row = F,
                                     scale = "row", 
                                     cluster_rows = F, 
                                     cluster_cols = F, 
                                     show_colnames = F)[[4]],
                      CD8.TPHE = pheatmap(scale(heat$CD8.TPHE[, sort(colnames(heat$CD4.TPHE))]),
                                     color = viridis(color.cut),
                                     annotation_row = data.frame(row.names = rownames(heat$CD8.TPHE), Dataset = c(rep("CD8.TPHE", nrow(heat$CD8.TPHE)))),
                                     annotation_colors = list(Dataset = c(CD8.TPHE = "red")),
                                     annotation_names_row = F,
                                     scale = "row", 
                                     cluster_rows = F, 
                                     cluster_cols = F,
                                     fontsize_col = 14)[[4]])

tphe.maps <- do.call(grid.arrange, heatmaps.tphe)

pdf("./results_R/plots/TPHE_ICS/Heatmaps_TPHE_CD4.CD8.combined.pdf", width = 8.5, height = 11)
do.call(grid.arrange, heatmaps.tphe)
dev.off()

## ICS markers as cols, clusters as rows; vertical stack using grid.arrange; ics.col_order <- c(8,5,2,6,3,1,4,7) to avoid 'NA' in middle
ics.col_order.CD4 <- c(8,5,2,6,3,1,4,7)
ics.col_order.CD8 <- c(3,6,2,7,4,1,5,8)
color.cut <- 4 # parametize color scale in final heatmap

heatmaps.ics <- list(CD4.ICS = pheatmap(scale(heat$CD4.ICS[, ics.col_order.CD4]), 
                                          color = viridis(color.cut),
                                          annotation_row = data.frame(row.names = rownames(heat$CD4.ICS), Dataset = c(rep("CD4.ICS", nrow(heat$CD4.ICS)))),
                                          annotation_colors = list(Dataset = c(CD4.ICS = "orange")),
                                          annotation_names_row = F,
                                          scale = "row", 
                                          cluster_rows = F, 
                                          cluster_cols = F, 
                                          show_colnames = F)[[4]],
                      CD8.ICS = pheatmap(scale(heat$CD8.ICS[, ics.col_order.CD8]),
                                          color = viridis(color.cut),
                                          annotation_row = data.frame(row.names = rownames(heat$CD8.ICS), Dataset = c(rep("CD8.ICS", nrow(heat$CD8.ICS)))),
                                          annotation_colors = list(Dataset = c(CD8.ICS = "purple")),
                                          annotation_names_row = F,
                                          scale = "row", 
                                          cluster_rows = F, 
                                          cluster_cols = F,
                                          fontsize_col = 14)[[4]])

ics.maps <- do.call(grid.arrange, heatmaps.ics)

pdf("./results_R/plots/TPHE_ICS/Heatmaps_ICS_CD4.CD8.combined.pdf", width = 8.5, height = 11)
do.call(grid.arrange, heatmaps.ics)
dev.off()

##
pdf("./results_R/plots/TPHE_ICS/Heatmaps_TPHE.ICS_CD4.CD8.combined.pdf", width = 11, height = 8.5)
grid.arrange(heatmaps.tphe$CD4.TPHE, heatmaps.ics$CD4.ICS, heatmaps.tphe$CD8.TPHE, heatmaps.ics$CD8.ICS)
dev.off()
