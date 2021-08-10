source("flowcytometry/source_R/fSOM.functions_NGL.R")
library(grid)
library(gridExtra)
library(dplyr)
library(pheatmap)
library(viridis)

#### fsom (metacluster) heatmaps
## need to reorder based on original heatmaps

col_order <- list(ics.cd4 = c("CD45RA", "IL8", "IL2", "TNFA", "IL4", "IL17", "IFNG"),
                  ics.cd8 = c("CD45RA", "IL8", "IL2", "TNFA", "IL4", "IFNG", "CD107A"),
                  tphe.cd4 = c("CD127", "FOXP3", "CD31", "CD28", "CD45RO", "CD197", "CD185", "CD57" ),
                  tphe.cd8 = c("CD122", "CD127", "FOXP3", "CD31", "CD28", "CD45RO", "CD197", "CD185", "CD57", "GRZB", "PERFORIN")
)

row_order <- list(ics.cd4 = read.csv("flowcytometry/intermediates/meta_cluster.row_order.ics_CD4.csv")[, 1],
                  ics.cd8 = read.csv("flowcytometry/intermediates/meta_cluster.row_order.ics_CD8.csv")[, 1],
                  tphe.cd4 = read.csv("flowcytometry/intermediates/meta_cluster.row_order.tphe_CD4.csv")[, 1],
                  tphe.cd8 = read.csv("flowcytometry/intermediates/meta_cluster.row_order.tphe_CD8.csv")[, 1]
)
row_order$tphe.cd4 <- row_order$tphe.cd4[-which(row_order$tphe.cd4 %in% 7)]

## read in fsoms to generate median expression for heatmaps
fsom <- list(ics.cd4 = readRDS("flowcytometry/intermediates/03_CD4_fSOM.CD4Cytokines.rds"),
              ics.cd8 = readRDS("flowcytometry/intermediates/03B_CD8_fSOM.CD8Cytokines.RDS"),
              tphe.cd4 = readRDS("flowcytometry/intermediates/2018-11-14.fSOM_PASS3_CD4p_rerun.rds"),
              tphe.cd8 = readRDS("flowcytometry/intermediates/2019-01-18.fSOM_PASS3_CD8p.rds")
)
fsom$ics.cd8 <- nodes.to.meta_fSOM(fsom$ics.cd8, grep("^12$|^13$", fsom$ics.cd8$metaclustering))

fsom_expr <- list(descaled = lapply(fsom, function(i) data.frame(descale.fsom(i)[[1]], cell_clustering = as.numeric(i$metaclustering)) 
                                    %>% group_by(cell_clustering) %>% summarize_all(list(median))),
                  descaled01 = lapply(fsom, function(i) data.frame(descale.fsom(i)[[3]], cell_clustering = as.numeric(i$metaclustering)) 
                                      %>% group_by(cell_clustering) %>% summarize_all(list(median)))
)

saveRDS(fsom_expr, 'flowcytometry/intermediates/metacluster_mfi.rds')
reorder.expr <- function(expr, row.order, col.order){
  expr <- as.data.frame(expr[row.order, col.order])
  rownames(expr) <- row.order
  expr
}

fsom_expr_reordered <- list(descaled = lapply(names(fsom_expr$descaled), function(i) reorder.expr(fsom_expr$descaled[[i]], row_order[[i]], col_order[[i]])),
                            descaled01 = lapply(names(fsom_expr$descaled01), function(i) reorder.expr(fsom_expr$descaled01[[i]], row_order[[i]], col_order[[i]]))
)
names(fsom_expr_reordered$descaled) <- names(fsom_expr$descaled)
names(fsom_expr_reordered$descaled01) <- names(fsom_expr$descaled01)

## match columns; per panel -- TPHE and ICS

fsom_expr_reordered$descaled01$tphe.cd4[, colnames(fsom_expr_reordered$descaled01$tphe.cd8)[!colnames(fsom_expr_reordered$descaled01$tphe.cd8) %in% colnames(fsom_expr_reordered$descaled01$tphe.cd4)]] <- NA
fsom_expr_reordered$descaled01$tphe.cd8[, colnames(fsom_expr_reordered$descaled01$tphe.cd4)[!colnames(fsom_expr_reordered$descaled01$tphe.cd4) %in% colnames(fsom_expr_reordered$descaled01$tphe.cd8)]] <- NA


fsom_expr_reordered$descaled$ics.cd4[, colnames(fsom_expr_reordered$descaled$ics.cd8)[!colnames(fsom_expr_reordered$descaled$ics.cd8) %in% colnames(fsom_expr_reordered$descaled$ics.cd4)]] <- NA
fsom_expr_reordered$descaled$ics.cd8[, colnames(fsom_expr_reordered$descaled$ics.cd4)[!colnames(fsom_expr_reordered$descaled$ics.cd4) %in% colnames(fsom_expr_reordered$descaled$ics.cd8)]] <- NA

##
tphe.combined.col_order <- c("CD122", "CD127", "FOXP3", "CD31", "CD28", "CD45RO", "CD197", "CD185", "CD57", "GRZB", "PERFORIN")
tphe.figure.col_names <- c(expression(paste("IL2r", beta)), expression(paste("IL7r", alpha)), "FOXP3", "CD31", "CD28", "CD45RO", "CCR7", "CXCR5", "CD57", "Granzyme", "Perforin")
ics.combined.col_order <- c("CD45RA", "IL8", "IL2", "TNFA", "IL4", "IL17", "IFNG", "CD107A")
ics.figure.col_names <- c("CD45RA", "IL-8", "IL-2", expression(paste("TNF", alpha)), "IL-4", "IL-17", expression(paste("IFN", gamma)), "CD107a")
##

pheatmap(fsom_expr_reordered$descaled$ics.cd4[, ics.combined.col_order], color = viridis(100), cluster_cols = F, cluster_rows = F, border_color = NA, legend = F)
pheatmap(fsom_expr_reordered$descaled$ics.cd8[rev(seq(row_order$ics.cd8)), ics.combined.col_order], color = viridis(100), cluster_cols = F, cluster_rows = F, border_color = NA, legend = F)
pheatmap(fsom_expr_reordered$descaled01$tphe.cd4[rev(seq(row_order$tphe.cd4)), tphe.combined.col_order], color = viridis(100), cluster_cols = F, cluster_rows = F, border_color = NA, legend = F)
pheatmap(fsom_expr_reordered$descaled01$tphe.cd8[rev(seq(row_order$tphe.cd8)), tphe.combined.col_order], color = viridis(100), cluster_cols = F, cluster_rows = F, border_color = NA, legend = F)

## generate list of heatmpas for grid arrange

fsom.heatmaps <- list(ics.cd4 = pheatmap(fsom_expr_reordered$descaled$ics.cd4[, ics.combined.col_order],
                                         #color = viridis(100), 
                                         cluster_cols = F, 
                                         cluster_rows = F, 
                                         border_color = NA, 
                                         legend = F,
                                         show_rownames = F,
                                         show_colnames = F,
                                         silent = T),
                      ics.cd8 = pheatmap(fsom_expr_reordered$descaled$ics.cd8[rev(seq(row_order$ics.cd8)), ics.combined.col_order], 
                                         #color = viridis(100), 
                                         cluster_cols = F, 
                                         cluster_rows = F, 
                                         border_color = NA, 
                                         legend = F,
                                         show_rownames = F,
                                         show_colnames = F,
                                         silent = T),
                      ics.cols = pheatmap(fsom_expr_reordered$descaled$ics.cd4[, ics.combined.col_order],
                                          #color = viridis(100), 
                                          cluster_cols = F, 
                                          cluster_rows = F, 
                                          border_color = NA, 
                                          legend = F,
                                          show_rownames = F,
                                          show_colnames = T,
                                          labels_col = ics.figure.col_names,
                                          hjust = 0,
                                          fontsize_col = 14,
                                          silent = T)[[4]][5], #store only the column labels
                      tphe.cd4 = pheatmap(fsom_expr_reordered$descaled01$tphe.cd4[rev(seq(row_order$tphe.cd4)), tphe.combined.col_order], 
                                          #color = viridis(100), 
                                          cluster_cols = F, 
                                          cluster_rows = F, 
                                          border_color = NA, 
                                          legend = F,
                                          show_rownames = F,
                                          show_colnames = F,
                                          silent = T),
                      tphe.cd8 = pheatmap(fsom_expr_reordered$descaled01$tphe.cd8[rev(seq(row_order$tphe.cd8)), tphe.combined.col_order], 
                                          #color = viridis(100), 
                                          cluster_cols = F, 
                                          cluster_rows = F, 
                                          border_color = NA, 
                                          legend = F,
                                          show_rownames = F,
                                          show_colnames = F,
                                          silent = T),
                      tphe.cols = pheatmap(fsom_expr_reordered$descaled01$tphe.cd4[rev(seq(row_order$tphe.cd4)), tphe.combined.col_order], 
                                           #color = viridis(100), 
                                           cluster_cols = F, 
                                           cluster_rows = F, 
                                           border_color = NA, 
                                           legend = F,
                                           show_rownames = F,
                                           show_colnames = T,
                                           labels_col = tphe.figure.col_names,
                                           fontsize_col = 14,
                                           silent = T)[[4]][5] #store only the column labels
)

grid.arrange(arrangeGrob(grobs = list(fsom.heatmaps$ics.cols)),
             arrangeGrob(grobs = lapply(fsom.heatmaps[grep("ics.cd", names(fsom.heatmaps))], `[[`, "gtable")),
             heights = c(1,10))

grid.arrange(arrangeGrob(grobs = list(fsom.heatmaps$tphe.cols)),
             arrangeGrob(grobs = lapply(fsom.heatmaps[grep("tphe.cd", names(fsom.heatmaps))], `[[`, "gtable")), 
             heights = c(1,10))

pdf(file.path("./results_R/plots/TPHE_ICS", paste0(Sys.Date(), "_", "Heatmaps_ICS_CD4.CD8.combined.pdf")), width = 4, height = 6)
grid.arrange(arrangeGrob(grobs = list(fsom.heatmaps$ics.cols)),
             arrangeGrob(grobs = lapply(fsom.heatmaps[grep("ics.cd", names(fsom.heatmaps))], `[[`, "gtable")), 
             heights = c(10,11))
dev.off()

png(file.path("./results_R/plots/TPHE_ICS", paste0(Sys.Date(), "_", "Heatmaps_ICS_CD4.CD8.combined_RdYlBu.png")),    # create PNG for the heat map        
    width = 4*600,        
    height = 6*600,
    res = 600,            # pixels per inch
    pointsize = 8)
grid.arrange(arrangeGrob(grobs = list(fsom.heatmaps$ics.cols)),
             arrangeGrob(grobs = lapply(fsom.heatmaps[grep("ics.cd", names(fsom.heatmaps))], `[[`, "gtable")), 
             heights = c(1,11))
dev.off()

pdf(file.path("./results_R/plots/TPHE_ICS", paste0(Sys.Date(), "_", "Heatmaps_TPHE_CD4.CD8.combined.pdf")), width = 4, height = 6)
grid.arrange(arrangeGrob(grobs = list(fsom.heatmaps$tphe.cols)),
             arrangeGrob(grobs = lapply(fsom.heatmaps[grep("tphe.cd", names(fsom.heatmaps))], `[[`, "gtable")), 
             heights = c(1,11))
dev.off()

png(file.path("./results_R/plots/TPHE_ICS", paste0(Sys.Date(), "_", "Heatmaps_TPHE_CD4.CD8.combined_RdYlBu.png")),    # create PNG for the heat map        
    width = 4*600,        
    height = 6*600,
    res = 600,            # pixels per inch
    pointsize = 8)
grid.arrange(arrangeGrob(grobs = list(fsom.heatmaps$tphe.cols)),
             arrangeGrob(grobs = lapply(fsom.heatmaps[grep("tphe.cd", names(fsom.heatmaps))], `[[`, "gtable")), 
             heights = c(1,11))
dev.off()

#### IST composition heatmaps

ist.md <- list(tphe = read.csv("intermediates/dmn/tphe_md.txt", sep = "\t", row.names = 1),
               ics = read.csv("intermediates/dmn/ics_md.txt", sep = "\t", row.names = 1))

ist.md <- lapply(ist.md, function(i) i[order(i), , drop = F])

ist <- list(tphe.cd4 = read.csv("intermediates/dmn/tphe_cd4_comp.txt", sep = "\t", row.names = 1),
            tphe.cd8 = read.csv("intermediates/dmn/tphe_cd8_comp.txt", sep = "\t", row.names = 1),
            ics.cd4 = read.csv("intermediates/dmn/ics_cd4_comp.txt", sep = "\t", row.names = 1),
            ics.cd8 = read.csv("intermediates/dmn/ics_cd8_comp.txt", sep = "\t", row.names = 1))

ist[grep("tphe", names(ist))] <- lapply(ist[grep("tphe", names(ist))], function(i) i[, rownames(ist.md[[grep("tphe", names(ist.md))]])]) # order by IST
ist[grep("ics", names(ist))] <- lapply(ist[grep("ics", names(ist))], function(i) i[, rownames(ist.md[[grep("ics", names(ist.md))]])]) # order by IST

rownames(ist$tphe.cd4) <- sub("TPHE CD4 C", "", rownames(ist$tphe.cd4))
rownames(ist$tphe.cd8) <- sub("TPHE CD8 C", "", rownames(ist$tphe.cd8))
rownames(ist$ics.cd4) <- sub("ICS CD4 C", "", rownames(ist$ics.cd4))
rownames(ist$ics.cd8) <- sub("ICS CD8 C", "", rownames(ist$ics.cd8))

ist.scale <- lapply(ist, function(i) t(apply(i, 1L, scales::rescale)))

ist.heatmaps <- list(ics.cd4 = pheatmap(ist.scale$ics.cd4, 
                                        color = colorRampPalette(rev(RColorBrewer::brewer.pal(n=7, name="RdYlBu")), bias=3)(100), 
                                        #annotation_col = ist.md$ics,
                                        cluster_rows = FALSE, 
                                        cluster_cols = FALSE,
                                        fontsize_row = 14,
                                        show_colnames = F,
                                        legend = F,
                                        annotation_legend = F, 
                                        gaps_col = unname(cumsum(table(ist.md$ics))),
                                        silent = T),
                     ics.cd8 = pheatmap(ist.scale$ics.cd8[nrow(ist.scale$ics.cd8):1, ], 
                                        color = colorRampPalette(rev(RColorBrewer::brewer.pal(n=7, name="RdYlBu")), bias=3)(100), 
                                        #annotation_col = ist.md$ics, 
                                        cluster_rows = FALSE, 
                                        cluster_cols = FALSE,
                                        fontsize_row = 14,
                                        show_colnames = F,
                                        legend = F,
                                        annotation_legend = F, 
                                        gaps_col = unname(cumsum(table(ist.md$ics))),
                                        silent = T),
                     ics.cols = grid.text(levels(factor(1:8)), x = c(0.060, 0.175, 0.30, 0.412, 0.525, 0.65, 0.77, 0.89), y = rep(0.15, 8), gp=gpar(fontsize = 16), draw = F),
                     tphe.cd4 = pheatmap(ist.scale$tphe.cd4[nrow(ist.scale$tphe.cd4):1, ], 
                                         color = colorRampPalette(rev(RColorBrewer::brewer.pal(n=7, name="RdYlBu")), bias=3)(100), 
                                         #annotation_col = ist.md$tphe, 
                                         cluster_rows = FALSE, 
                                         cluster_cols = FALSE,
                                         fontsize_row = 14,
                                         show_colnames = F,
                                         legend = F,
                                         annotation_legend = F, 
                                         gaps_col = unname(cumsum(table(ist.md$tphe))),
                                         silent = T),
                     tphe.cd8 = pheatmap(ist.scale$tphe.cd8[nrow(ist.scale$tphe.cd8):1, ], 
                                         color = colorRampPalette(rev(RColorBrewer::brewer.pal(n=7, name="RdYlBu")), bias=3)(100), 
                                         #annotation_col = ist.md$tphe, 
                                         cluster_rows = FALSE, 
                                         cluster_cols = FALSE,
                                         fontsize_row = 14,
                                         show_colnames = F,
                                         legend = F,
                                         annotation_legend = F, 
                                         gaps_col = unname(cumsum(table(ist.md$tphe))),
                                         silent = T),
                     #tphe.cols = NA
                     tphe.cols = grid.text(levels(factor(1:7)), x = c(0.1, 0.27, 0.42, 0.55, 0.65, 0.74, 0.85), y = rep(0.15, 7), gp=gpar(fontsize = 16), draw = F)
)


#### combined plots
## ICS subjects by metacluster proportions (?)
grid.arrange(arrangeGrob(grobs = list(ist.heatmaps$ics.cols)),
             arrangeGrob(grobs = lapply(ist.heatmaps[grep("ics.cd", names(ist.heatmaps))], `[[`, "gtable")), 
             heights = c(1,11))
## ICS Marker composition
grid.arrange(arrangeGrob(grobs = list(fsom.heatmaps$ics.cols)),
             arrangeGrob(grobs = lapply(fsom.heatmaps[grep("ics.cd", names(fsom.heatmaps))], `[[`, "gtable")), 
             heights = c(1,11))

grid.arrange(arrangeGrob(grobs = list(ist.heatmaps$tphe.cols)),
             arrangeGrob(grobs = lapply(ist.heatmaps[grep("tphe.cd", names(ist.heatmaps))], `[[`, "gtable")), 
             heights = c(1,11))
grid.arrange(arrangeGrob(grobs = list(fsom.heatmaps$tphe.cols)),
             arrangeGrob(grobs = lapply(fsom.heatmaps[grep("tphe.cd", names(fsom.heatmaps))], `[[`, "gtable")), 
             heights = c(1,11))


grid.arrange(arrangeGrob(grobs = c(list(ist.heatmaps$ics.cols), list(fsom.heatmaps$ics.cols),
                                   lapply(c(ist.heatmaps, fsom.heatmaps)[grep("ics.cd4", names(c(ist.heatmaps, fsom.heatmaps)))], `[[`, "gtable"),
                                   lapply(c(ist.heatmaps, fsom.heatmaps)[grep("ics.cd8", names(c(ist.heatmaps, fsom.heatmaps)))], `[[`, "gtable")),
                         nrow = 3, ncol = 2, heights = c(1, 4, 4), widths = 2:1))

grid.arrange(arrangeGrob(grobs = c(list(ist.heatmaps$tphe.cols), list(fsom.heatmaps$tphe.cols),
                                   lapply(c(ist.heatmaps, fsom.heatmaps)[grep("tphe.cd4", names(c(ist.heatmaps, fsom.heatmaps)))], `[[`, "gtable"),
                                   lapply(c(ist.heatmaps, fsom.heatmaps)[grep("tphe.cd8", names(c(ist.heatmaps, fsom.heatmaps)))], `[[`, "gtable")),
                         nrow = 3, ncol = 2, heights = c(1, 4, 4), widths = 2:1))

pdf(file.path("./results_R/plots/TPHE_ICS", paste0(Sys.Date(), "_", "Heatmaps_ICS_CD4.CD8.combined.pdf")), width = 11, height = 8.5)
grid.arrange(arrangeGrob(grobs = c(list(ist.heatmaps$ics.cols), list(fsom.heatmaps$ics.cols),
                                   lapply(c(ist.heatmaps, fsom.heatmaps)[grep("ics.cd4", names(c(ist.heatmaps, fsom.heatmaps)))], `[[`, "gtable"),
                                   lapply(c(ist.heatmaps, fsom.heatmaps)[grep("ics.cd8", names(c(ist.heatmaps, fsom.heatmaps)))], `[[`, "gtable")),
                         nrow = 3, ncol = 2, heights = c(1, 4, 4), widths = 2:1))
dev.off()

pdf(file.path("./results_R/plots/TPHE_ICS", paste0(Sys.Date(), "_", "Heatmaps_TPHE_CD4.CD8.combined.pdf")), width = 11, height = 8.5)
grid.arrange(arrangeGrob(grobs = c(list(ist.heatmaps$tphe.cols), list(fsom.heatmaps$tphe.cols),
                                   lapply(c(ist.heatmaps, fsom.heatmaps)[grep("tphe.cd4", names(c(ist.heatmaps, fsom.heatmaps)))], `[[`, "gtable"),
                                   lapply(c(ist.heatmaps, fsom.heatmaps)[grep("tphe.cd8", names(c(ist.heatmaps, fsom.heatmaps)))], `[[`, "gtable")),
                         nrow = 3, ncol = 2, heights = c(1, 4, 4), widths = 2:1))
dev.off()





grid.arrange(arrangeGrob(grobs = c(lapply(c(ist.heatmaps, fsom.heatmaps)[grep("tphe.cd4", names(c(ist.heatmaps, fsom.heatmaps)))], `[[`, "gtable"),
                                   lapply(c(ist.heatmaps, fsom.heatmaps)[grep("tphe.cd8", names(c(ist.heatmaps, fsom.heatmaps)))], `[[`, "gtable"),
                                   lapply(c(ist.heatmaps, fsom.heatmaps)[grep("ics.cd4", names(c(ist.heatmaps, fsom.heatmaps)))], `[[`, "gtable"),
                                   lapply(c(ist.heatmaps, fsom.heatmaps)[grep("ics.cd8", names(c(ist.heatmaps, fsom.heatmaps)))], `[[`, "gtable")), 
                         nrow = 2, ncol= 4))
