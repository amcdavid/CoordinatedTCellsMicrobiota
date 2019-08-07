library(gridExtra) ; library(pheatmap) ; library(viridis)

#### unify heatmaps for annotation key

## read in flowSOM cluster MFI .csvs
input.mfi <- list.files("./results_R/fSOMs/TPHE_ICS/", full.names = T, pattern = "MFI")

fsom.mfi <- list(tphe.cd4 = read.csv(input.mfi[grep("tphe.cd4", input.mfi)], row.names = 1),
                 tphe.cd8 = read.csv(input.mfi[grep("tphe.cd8", input.mfi)], row.names = 1),
                 ics.cd4 = read.csv(input.mfi[grep("ics.cd4", input.mfi)], row.names = 1),
                 ics.cd8 = read.csv(input.mfi[grep("ics.cd8", input.mfi)], row.names = 1)
)

## for legibility, just use meta cluster number
fsom.mfi <- sapply(names(fsom.mfi), function(i) {
  rownames(fsom.mfi[[i]]) <- sub("_meta.cluster", "", rownames(fsom.mfi[[i]]))
  return(fsom.mfi[[i]])
})

#OR

# ## generate unique rownames
# heat <- sapply(names(heat), function(i) {
#   rownames(heat[[i]]) <- sub("_meta.cluster", "mc.", paste0(i, "_", rownames(heat[[i]])))
#   return(heat[[i]])
# })

## match columns; per panel -- TPHE and ICS
fsom.mfi$tphe.cd4[, colnames(fsom.mfi$tphe.cd8)[!(colnames(fsom.mfi$tphe.cd8) %in% colnames(fsom.mfi$tphe.cd4))]] <- NA
fsom.mfi$tphe.cd8[, colnames(fsom.mfi$tphe.cd4)[!(colnames(fsom.mfi$tphe.cd4) %in% colnames(fsom.mfi$tphe.cd8))]] <- NA

fsom.mfi$ics.cd4[, colnames(fsom.mfi$ics.cd8)[!(colnames(fsom.mfi$ics.cd8) %in% colnames(fsom.mfi$ics.cd4))]] <- NA
fsom.mfi$ics.cd8[, colnames(fsom.mfi$ics.cd4)[!(colnames(fsom.mfi$ics.cd4) %in% colnames(fsom.mfi$ics.cd8))]] <- NA

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