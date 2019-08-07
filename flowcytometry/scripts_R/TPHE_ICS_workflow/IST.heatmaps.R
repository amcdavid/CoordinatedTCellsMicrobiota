##IST composition heatmaps

library("pheatmap")
library("viridis")

ist.md <- list(tphe = read.csv("./ISTs/tphe_md.txt", sep = "\t", row.names = 1),
               ics = read.csv("./ISTs/ics_md.txt", sep = "\t", row.names = 1))

ist.md <- lapply(ist.md, function(i) i[order(i), , drop = F])

ist <- list(tphe.cd4 = read.csv("./ISTs/tphe_cd4_comp.txt", sep = "\t", row.names = 1),
            tphe.cd8 = read.csv("./ISTs/tphe_cd8_comp.txt", sep = "\t", row.names = 1),
            ics.cd4 = read.csv("./ISTs/ics_cd4_comp.txt", sep = "\t", row.names = 1),
            ics.cd8 = read.csv("./ISTs/ics_cd8_comp.txt", sep = "\t", row.names = 1))

ist[grep("tphe", names(ist))] <- lapply(ist[grep("tphe", names(ist))], function(i) i[, rownames(ist.md[[grep("tphe", names(ist.md))]])])
ist[grep("ics", names(ist))] <- lapply(ist[grep("ics", names(ist))], function(i) i[, rownames(ist.md[[grep("ics", names(ist.md))]])])

rownames(ist$tphe.cd4) <- sub("TPHE CD4 C", "", rownames(ist$tphe.cd4))
rownames(ist$tphe.cd8) <- sub("TPHE CD8 C", "", rownames(ist$tphe.cd8))
rownames(ist$ics.cd4) <- sub("ICS CD4 C", "", rownames(ist$ics.cd4))
rownames(ist$ics.cd8) <- sub("ICS CD8 C", "", rownames(ist$ics.cd8))

ist.scale <- lapply(ist, function(i) t(apply(i, 1L, scales::rescale)))

#pheatmap(ist.scale$tphe.cd4, color = colors, annotation_col = ist.md$tphe, cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = F, annotation_legend = F, gaps_col = unname(cumsum(table(ist.md$tphe))))
#BETTER:
#pheatmap(ist.scale$tphe.cd4, color = viridis(20), cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE, gaps_col = unname(cumsum(table(ist.md$tphe))))


ist.heatmaps <- list(ics.cd4 = pheatmap(ist.scale$ics.cd4, 
                                        color = colorRampPalette(rev(RColorBrewer::brewer.pal(n=7, name="RdYlBu")), bias=3)(100), 
                                        #annotation_col = ist.md$ics,
                                        border_color = "black",
                                        cluster_rows = FALSE, 
                                        cluster_cols = FALSE, 
                                        show_colnames = F, 
                                        annotation_legend = F, 
                                        gaps_col = unname(cumsum(table(ist.md$ics))),
                                        silent = T)[[4]],
                     ics.cd8 = pheatmap(ist.scale$ics.cd8, 
                                        color = colorRampPalette(rev(RColorBrewer::brewer.pal(n=7, name="RdYlBu")), bias=3)(100), 
                                        #annotation_col = ist.md$ics, 
                                        cluster_rows = FALSE, 
                                        cluster_cols = FALSE, 
                                        show_colnames = F, 
                                        annotation_legend = F, 
                                        gaps_col = unname(cumsum(table(ist.md$ics))),
                                        silent = T)[[4]],
                     tphe.cd4 = pheatmap(ist.scale$tphe.cd4, 
                                         color = colorRampPalette(rev(RColorBrewer::brewer.pal(n=7, name="RdYlBu")), bias=3)(100), 
                                         #annotation_col = ist.md$tphe, 
                                         cluster_rows = FALSE, 
                                         cluster_cols = FALSE, 
                                         show_colnames = F, 
                                         annotation_legend = F, 
                                         gaps_col = unname(cumsum(table(ist.md$tphe))),
                                         silent = T)[[4]],
                     tphe.cd8 = pheatmap(ist.scale$tphe.cd8, 
                                         color = colorRampPalette(rev(RColorBrewer::brewer.pal(n=7, name="RdYlBu")), bias=3)(100), 
                                         #annotation_col = ist.md$tphe, 
                                         cluster_rows = FALSE, 
                                         cluster_cols = FALSE, 
                                         show_colnames = F, 
                                         annotation_legend = F, 
                                         gaps_col = unname(cumsum(table(ist.md$tphe))),
                                         silent = T)[[4]]
)
