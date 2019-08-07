library(data.table)
library(pheatmap)

dat <- readRDS("./analyses_R/ICS_tSNE.Explore/SEB.CD4CD69pos.RDS")

heat <- dat$fsom.mclust.agg[, dat$dims.used]

rownames(heat) <- paste0("Meta.Cluster_", rownames(heat))

pheatmap(scale(heat), scale = "row", cluster_rows = FALSE, treeheight_col = 0, fontsize_col = 18)

pdf("ICS_fSOM_CD4pCD69p.metaclusters.HEATMAP.pdf")
pheatmap(scale(heat), scale = "row", cluster_rows = FALSE, treeheight_col = 0, fontsize_col = 18)
dev.off()


dat <- readRDS("./analyses_R/ICS_tSNE.Explore/SEB.CD8CD69pos.RDS")

heat <- dat$fsom.mclust.agg[, dat$dims.used]

rownames(heat) <- paste0("Meta.Cluster_", rownames(heat))

pheatmap(scale(heat), scale = "row", cluster_rows = FALSE, treeheight_col = 0, fontsize_col = 18)

pdf("ICS_fSOM_CD8pCD69p.metaclusters.HEATMAP.pdf")
pheatmap(scale(heat), scale = "row", cluster_rows = FALSE, treeheight_col = 0, fontsize_col = 18)
dev.off()


heat <- read.csv("./results_R/fSOMs/ICS/fSOM_CD8cytokines.MetaClustersMFI.csv", row.names = 1)
