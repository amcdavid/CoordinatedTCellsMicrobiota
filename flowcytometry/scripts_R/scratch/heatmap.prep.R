TPHE.CD4.heat <- read.csv("./results_R/fSOMs/TPHE/TPHE_fSOM.CD4_MetaCluster_MFI.csv", row.names = 1)
TPHE.CD4.heat <- TPHE.CD4.heat[-19, ] # bad cluster

TPHE.CD8.heat <- read.csv("./results_R/fSOMs/TPHE/TPHE_fSOM.CD8_MetaCluster_MFI.csv", row.names = 1)

ICS.CD4.heat <- read.csv("./results_R/fSOMs/ICS/ICS_fSOM.CD4cytokines_MetaClusters_MFI.csv", row.names = 1)

ICS.CD8.heat <- read.csv("./results_R/fSOMs/ICS/ICS_fSOM.CD8cytokines_MetaClusters_MFI.csv", row.names = 1)
ICS.CD8.heat[12, ] <- (ICS.CD8.heat[12, ] + ICS.CD8.heat[13, ])/2
ICS.CD8.heat <- ICS.CD8.heat[-13, ]

write.csv(TPHE.CD4.heat, "./results_R/fSOMs/TPHE/FINAL/TPHE_fSOM.CD4_MetaCluster_MFI.csv", row.names = TRUE)
write.csv(TPHE.CD8.heat, "./results_R/fSOMs/TPHE/FINAL/TPHE_fSOM.CD8_MetaCluster_MFI.csv", row.names = TRUE)
write.csv(ICS.CD4.heat, "./results_R/fSOMs/ICS/FINAL/ICS_fSOM.CD4cytokines_MetaClusters_MFI.csv", row.names = TRUE)
write.csv(ICS.CD8.heat, "./results_R/fSOMs/ICS/FINAL/ICS_fSOM.CD8cytokines_MetaClusters_MFI.csv", row.names = TRUE)

library(pheatmap)

heat.data <- list.files("./results_R/fSOMs", full.names = TRUE, recursive = TRUE, pattern = "MetaCluster_MFI.csv")

heat <- list(CD4.TPHE = data.frame(read.csv(grep("TPHE_fSOM.CD4", heat.data, value = T), row.names = 1)),
             CD4.ICS = data.frame(read.csv(grep("ICS_fSOM.CD4", heat.data, value = T), row.names = 1)),
             CD8.TPHE = data.frame(read.csv(grep("TPHE_fSOM.CD8", heat.data, value = T), row.names = 1)),
             CD8.ICS = data.frame(read.csv(grep("ICS_fSOM.CD8", heat.data, value = T), row.names = 1))
             )

pheatmap(scale(heat$CD4.TPHE), scale = "row")
pheatmap(scale(heat$CD8.TPHE), scale = "row")
pheatmap(scale(heat$CD4.ICS), scale = "row")
pheatmap(scale(heat$CD8.ICS), scale = "row")

pheatmap(scale(heat$CD4.TPHE), scale = "row", cluster_rows = FALSE, cluster_cols = FALSE)
pheatmap(scale(heat$CD8.TPHE), scale = "row", cluster_rows = FALSE, cluster_cols = FALSE)
pheatmap(scale(heat$CD4.ICS), scale = "row", cluster_rows = FALSE, cluster_cols = FALSE)
pheatmap(scale(heat$CD8.ICS), scale = "row", cluster_rows = FALSE, cluster_cols = FALSE)


pdf("./temp/Heatmaps_FlowSOM.MetaClusterMFI_TPHE.ICS.pdf", width = 11, height = 8.5)
for(i in seq(heat)){
    pheatmap(scale(heat[[i]]), scale = "row", main = paste0(names(heat)[i], " ", "FlowSOM Meta-Cluster MFI"))
}
for(i in seq(heat)){
  pheatmap(scale(heat[[i]]), scale = "row", cluster_rows = FALSE, cluster_cols = FALSE, main = paste0(names(heat)[i], " ", "FlowSOM Meta-Cluster MFI"))
}
dev.off()

##
library(ComplexHeatmap)

library(gridExtra); library(pheatmap)
m <- matrix(c(1:4), ncol=2)
n <- matrix(c(1,1,1,2), ncol=2)
a <- list(pheatmap(m)[[4]])
a[[2]] <- pheatmap(n)[[4]]
z <- do.call(grid.arrange,a)
plot(z)

heatmaps <- list(CD4.TPHE = pheatmap(scale(heat$CD4.TPHE), scale = "row", cluster_rows = FALSE, cluster_cols = FALSE)[[4]],
                 CD8.TPHE = pheatmap(scale(heat$CD8.TPHE), scale = "row", cluster_rows = FALSE, cluster_cols = FALSE)[[4]])

heatmaps.tpose <- list(CD4.TPHE = pheatmap(t(scale(heat$CD4.TPHE)), scale = "column", cluster_rows = FALSE, cluster_cols = FALSE)[[4]],
                       CD8.TPHE = pheatmap(t(scale(heat$CD8.TPHE)), scale = "column", cluster_rows = FALSE, cluster_cols = FALSE)[[4]])

z <- do.call(grid.arrange, heatmaps)
