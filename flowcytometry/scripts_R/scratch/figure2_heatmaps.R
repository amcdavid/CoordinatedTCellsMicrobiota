library(pheatmap); library(viridis); library(grid)

input.mfi <- list.files("./results_R/fSOMs/TPHE_ICS/", full.names = T, pattern = "MFI")


# tphe.cd4.col_order <- c("CD127", "FOXP3", "CD31", "CD28", "CD45RO", "CD197", "CD185", "CD57" )
# tphe.cd8.col_order <- c("CD122", "CD127", "FOXP3", "CD31", "CD28", "CD45RO", "CD197", "CD185", "CD57", "GRZB", "PERFORIN")
# ics.cd4.col_order <- c("CD45RA", "IL8", "IL2", "TNFA", "IL4", "IL17", "IFNG")
# ics.cd8.col_order <- c("CD45RA", "IL8", "IL2", "TNFA", "IL4", "IL17", "IFNG", "CD107A")

fsom.mfi <- list(tphe.cd4 = read.csv(input.mfi[grep("tphe.cd4", input.mfi)], row.names = 1),
                 tphe.cd8 = read.csv(input.mfi[grep("tphe.cd8", input.mfi)], row.names = 1),
                 ics.cd4 = read.csv(input.mfi[grep("ics.cd4", input.mfi)], row.names = 1),
                 ics.cd8 = read.csv(input.mfi[grep("ics.cd8", input.mfi)], row.names = 1)
)

pdf("./results_R/fSOMs/TPHE_ICS/meta.clusters_heatmaps.pdf")
fsom.heatmaps <- list(tphe.cd4.heat = pheatmap(scale(fsom.mfi$tphe.cd4), 
                                               scale = "row", 
                                               clustering_method = "mcquitty",
                                               cluster_cols = F,
                                               cutree_rows = 4,
                                               color = viridis(20),
                                               main = "CD4 T-cell Phenotyping Panel FlowSOM Meta-cluster MFI"),
                      tphe.cd8.heat = pheatmap(scale(fsom.mfi$tphe.cd8), 
                                               scale = "row", 
                                               clustering_method = "mcquitty",
                                               cluster_cols = F,
                                               cutree_rows = 5,
                                               color = viridis(5),
                                               main = "CD8 T-cell Phenotyping Panel FlowSOM Meta-cluster MFI"),
                      ics.cd4.heat = pheatmap(scale(fsom.mfi$ics.cd4), 
                                              scale = "row", 
                                              clustering_method = "mcquitty",
                                              cluster_cols = F,
                                              cutree_rows = 5,
                                              color = viridis(5),
                                              main = "CD4 T-cell Intracellular Cytokine Panel FlowSOM Meta-cluster MFI"),
                      ics.cd8.heat = pheatmap(scale(fsom.mfi$ics.cd8), 
                                              scale = "row", 
                                              clustering_method = "mcquitty",
                                              cluster_cols = F,
                                              cutree_rows = 7,
                                              color = viridis(20),
                                              main = "CD8 T-cell Intracellular Cytokine Panel FlowSOM Meta-cluster MFI")
)
dev.off()
