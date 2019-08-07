library(pheatmap); library(viridis); library(grid); library(gridExtra)

fsom.heat <- function(fsom.object){
  heat <- aggregate(fsom.object$map$codes, list(meta.cluster = fsom.object$metaclustering), median)[, -1]
  rownames(heat) <- paste0(rownames(heat), "_meta.cluster")
  colnames(heat) <- fsom.object$markers[fsom.object$map$colsUsed]
  heat
}

##
tphe.cd4.col_order <- c("CD127", "FOXP3", "CD31", "CD28", "CD45RO", "CD197", "CD185", "CD57" )
tphe.cd8.col_order <- c("CD122", "CD127", "FOXP3", "CD31", "CD28", "CD45RO", "CD197", "CD185", "CD57", "GRZB", "PERFORIN")
ics.cd4.col_order <- c("CD45RA", "IL8", "IL2", "TNFA", "IL4", "IL17", "IFNG")
ics.cd8.col_order <- c("CD45RA", "IL8", "IL2", "TNFA", "IL4", "IL17", "IFNG", "CD107A")

fsom.mfi <- list(tphe.cd4 = fsom.heat(readRDS("./results_R/fSOMs/TPHE/2018-11-14.fSOM_PASS3_CD4p_rerun.rds"))[-19, ][tphe.cd4.col_order],
                       tphe.cd8 = fsom.heat(readRDS("./results_R/fSOMs/TPHE/2019-01-18.fSOM_PASS3_CD8p.rds"))[tphe.cd8.col_order],
                       ics.cd4 = fsom.heat(readRDS("./analyses_R/FSOM_analysis/03_CD4_fSOM.CD4Cytokines.rds"))[ics.cd4.col_order],
                       ics.cd8 = fsom.heat(readRDS("./analyses_R/FSOM_analysis/03B_CD8_fSOM.CD8Cytokines.RDS"))[ics.cd8.col_order]
                       )
fsom.mfi$ics.cd8 <- fsom.mfi$ics.cd8[-grep("IL17", colnames(fsom.mfi$ics.cd8))]
fsom.mfi$ics.cd8[12, ] <- (fsom.mfi$ics.cd8[12, ] + fsom.mfi$ics.cd8[13, ])/2
fsom.mfi$ics.cd8 <- fsom.mfi$ics.cd8[-13, ]

sapply(seq(fsom.mfi), function(i) write.csv(fsom.mfi[[i]], file.path("./results_R/fSOMs/TPHE_ICS/", paste0("fsom_MetaCluster.MFI_", names(fsom.mfi)[i], ".csv"))))
##

input.mfi <- list.files("./results_R/fSOMs/TPHE_ICS/", full.names = T, pattern = "MFI")

fsom.mfi <- list(tphe.cd4 = read.csv(input.mfi[grep("tphe.cd4", input.mfi)], row.names = 1),
                 tphe.cd8 = read.csv(input.mfi[grep("tphe.cd8", input.mfi)], row.names = 1),
                 ics.cd4 = read.csv(input.mfi[grep("ics.cd4", input.mfi)], row.names = 1),
                 ics.cd8 = read.csv(input.mfi[grep("ics.cd8", input.mfi)], row.names = 1)
                 )

fsom.heatmaps <- list(tphe.cd4.heat = pheatmap(scale(fsom.mfi$tphe.cd4), 
                                               scale = "row", 
                                               clustering_method = "mcquitty",
                                               cluster_cols = T,
                                               cutree_rows = 4,
                                               color = viridis(20),
                                               main = "CD4 T-cell Phenotyping Panel FlowSOM Meta-cluster MFI"),
                      tphe.cd8.heat = pheatmap(scale(fsom.mfi$tphe.cd8), 
                                               scale = "row", 
                                               clustering_method = "mcquitty",
                                               cluster_cols = T,
                                               cutree_rows = 5,
                                               color = viridis(5),
                                               main = "CD8 T-cell Phenotyping Panel FlowSOM Meta-cluster MFI"),
                      ics.cd4.heat = pheatmap(scale(fsom.mfi$ics.cd4), 
                                              scale = "row", 
                                              clustering_method = "mcquitty",
                                              cluster_cols = T,
                                              cutree_rows = 4,
                                              color = viridis(5),
                                              main = "CD4 T-cell Intracellular Cytokine Panel FlowSOM Meta-cluster MFI"),
                      ics.cd8.heat = pheatmap(scale(fsom.mfi$ics.cd8), 
                                              scale = "row", 
                                              clustering_method = "mcquitty",
                                              cluster_cols = T,
                                              cutree_rows = 7,
                                              color = viridis(20),
                                              main = "CD8 T-cell Intracellular Cytokine Panel FlowSOM Meta-cluster MFI")
)

##
meta.cluster.row_order.tphe.cd4 <- data.frame(tphe.cd4 = fsom.heatmaps$tphe.cd4.heat$tree_row$order)
meta.cluster.row_order.tphe.cd4[meta.cluster.row_order.tphe.cd4$tphe.cd4 >= 19, ] <- meta.cluster.row_order.tphe.cd4[meta.cluster.row_order.tphe.cd4$tphe.cd4 >= 19, ] + 1
write.csv(meta.cluster.row_order.tphe.cd4, "./results_R/fSOMs/TPHE_ICS/meta_cluster.row_order.tphe_CD4.csv", row.names = F)

meta.cluster.row_order.tphe.cd8 <- data.frame(fsom.heatmaps$tphe.cd8 = tphe.cd8.heat$tree_row$order)
write.csv(meta.cluster.row_order.tphe.cd8, "./results_R/fSOMs/TPHE_ICS/meta_cluster.row_order.tphe_CD8.csv", row.names = F)

meta.cluster.row_order.ics.cd4 <- data.frame(fsom.heatmaps$ics.cd4 = ics.cd4.heat$tree_row$order)
write.csv(meta.cluster.row_order.ics.cd4, "./results_R/fSOMs/TPHE_ICS/meta_cluster.row_order.ics_CD4.csv", row.names = F)

meta.cluster.row_order.ics.cd8 <- data.frame(fsom.heatmaps$ics.cd8 = ics.cd8.heat$tree_row$order)
write.csv(meta.cluster.row_order.ics.cd8, "./results_R/fSOMs/TPHE_ICS/meta_cluster.row_order.ics_CD8.csv", row.names = F)
##

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

## match columns; per panel -- TPHE and ICS
fsom.mfi$tphe.cd4[, colnames(fsom.mfi$tphe.cd8)[!(colnames(fsom.mfi$tphe.cd8) %in% colnames(fsom.mfi$tphe.cd4))]] <- NA
fsom.mfi$tphe.cd8[, colnames(fsom.mfi$tphe.cd4)[!(colnames(fsom.mfi$tphe.cd4) %in% colnames(fsom.mfi$tphe.cd8))]] <- NA

fsom.mfi$ics.cd4[, colnames(fsom.mfi$ics.cd8)[!(colnames(fsom.mfi$ics.cd8) %in% colnames(fsom.mfi$ics.cd4))]] <- NA
fsom.mfi$ics.cd8[, colnames(fsom.mfi$ics.cd4)[!(colnames(fsom.mfi$ics.cd4) %in% colnames(fsom.mfi$ics.cd8))]] <- NA

tphe.col_order <- c("CD122", "CD127", "FOXP3", "CD31", "CD28", "CD45RO", "CD197", "CD185", "CD57", "GRZB", "PERFORIN")
ics.col_order <- c("CD45RA", "IL8", "IL2", "TNFA", "IL4", "IL17", "IFNG", "CD107A")

## generate heatmaps; use pheatmap list element '4' for gtable
fsom.heatmaps <- list(tphe.cd4.heat = pheatmap(scale(fsom.mfi$tphe.cd4[, tphe.col_order]), 
                                               scale = "row", 
                                               clustering_method = "mcquitty",
                                               cluster_rows = T,
                                               cluster_cols = F,
                                               cutree_rows = 4,
                                               color = viridis(20),
                                               main = "CD4 T-cell Phenotyping Panel FlowSOM Meta-cluster MFI",
                                               show_colnames = F)[[4]],
                      tphe.cd8.heat = pheatmap(scale(fsom.mfi$tphe.cd8[, tphe.col_order]), 
                                               scale = "row", 
                                               clustering_method = "mcquitty",
                                               cluster_rows = T,
                                               cluster_cols = F,
                                               cutree_rows = 5,
                                               color = viridis(5),
                                               main = "CD8 T-cell Phenotyping Panel FlowSOM Meta-cluster MFI",
                                               fontsize_col = 12)[[4]],
                      ics.cd4.heat = pheatmap(scale(fsom.mfi$ics.cd4[ics.col_order]), 
                                              scale = "row", 
                                              clustering_method = "mcquitty",
                                              cluster_cols = F,
                                              cutree_rows = 4,
                                              color = viridis(5),
                                              main = "CD4 T-cell Intracellular Cytokine Panel FlowSOM Meta-cluster MFI",
                                              show_colnames = F)[[4]],
                      ics.cd8.heat = pheatmap(scale(fsom.mfi$ics.cd8[ics.col_order]), 
                                              scale = "row", 
                                              clustering_method = "mcquitty",
                                              cluster_cols = F,
                                              cutree_rows = 7,
                                              color = viridis(20),
                                              main = "CD8 T-cell Intracellular Cytokine Panel FlowSOM Meta-cluster MFI",
                                              fontsize_col = 12)[[4]]
)

do.call(grid.arrange, fsom.heatmaps[c(1,2)])
do.call(grid.arrange, fsom.heatmaps[c(3,4)])

pdf(file.path("./results_R/plots/TPHE_ICS", paste0(Sys.Date(), "_", "Heatmaps_TPHE_CD4.CD8.combined.pdf")), width = 8.5, height = 11)
do.call(grid.arrange, fsom.heatmaps[c(1,2)])
dev.off()

pdf(file.path("./results_R/plots/TPHE_ICS", paste0(Sys.Date(), "_", "Heatmaps_ICS_CD4.CD8.combined.pdf")), width = 8.5, height = 11)
do.call(grid.arrange, fsom.heatmaps[c(3,4)])
dev.off()

pdf(file.path("./results_R/plots/TPHE_ICS", paste0(Sys.Date(), "_", "Heatmaps_TPHE.ICS_CD4.CD8.combined.pdf")), width = 8.5, height = 11)
grid.arrange(fsom.heatmaps$tphe.cd4.heat, fsom.heatmaps$ics.cd4.heat, fsom.heatmaps$tphe.cd8.heat, fsom.heatmaps$ics.cd8.heat)
dev.off()
