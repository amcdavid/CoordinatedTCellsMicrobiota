## fSOM Build using subsampled files; CD4

## file paths
tphe.warped.files <- list.files("./data_modified/TPHE_warping/fSOM_CD4p_PASS1_scatter.warped_SUBSAMPLED/", full.names = TRUE)

input.set <- read.flowSet(tphe.warped.files, transformation = FALSE)

fSOM <- ReadInput(input.set, scale = TRUE)
fSOM$markers <- get.markers(input.set[[1]])

mset <- match(c("FSC.A", "FSC.W", "SSC.A", "SSC.W"), fSOM$markers)

fSOM <- BuildSOM(fSOM, colsToUse = mset, xdim = 5, ydim = 5) # PASS 1
fSOM <- BuildMST(fSOM)

nbclust <- 10
fSOM[["metaclustering"]] <- as.factor(metaClustering_consensus(fSOM$map$codes, k = nbclust))
PlotStars(fSOM, backgroundValues = fSOM$metaclustering)

fSOM$heatmap <- pheatmap(scale(fsom.clusters.heatmap(fSOM)), scale = "row")
fcs.explorer(fsom.dat(fSOM), "FSC.A", "SSC.A", 20000)

fSOM <- nodes.to.meta_fSOM(fSOM, grep("3|4|10", fSOM$metaclustering))  # stochastic

fSOM$heatmap <- pheatmap(scale(fsom.clusters.heatmap(fSOM)), scale = "row")
fcs.explorer(fsom.dat(fSOM), "FSC.A", "SSC.A", 20000)

fSOM <- nodes.to.meta_fSOM(fSOM, grep("5|8", fSOM$metaclustering))  # stochastic

fSOM$heatmap <- pheatmap(scale(fsom.clusters.heatmap(fSOM)), scale = "row")
fcs.explorer(fsom.dat(fSOM), "FSC.A", "SSC.A", 20000)

fSOM$metaclustering.anno <- vector(mode = "list", length = length(levels(fSOM$metaclustering)))
fSOM$metaclustering.anno[[7]] <- "singlets"

savefsom(fSOM, file.path("./results_R/fSOMs/TPHE", paste(Sys.Date(), "fSOM_SINGLETS.PASS2.rds", sep = ".")))

####
## fSOM Build using subsampled files; CD8

## file paths
tphe.warped.files <- list.files("./data_modified/TPHE_warping/fSOM_CD8p_PASS1_scatter.warped_SUBSAMPLED/", full.names = TRUE)

input.set <- read.flowSet(tphe.warped.files, transformation = FALSE)

fSOM <- ReadInput(input.set, scale = TRUE)
fSOM$markers <- get.markers(input.set[[1]])

mset <- match(c("FSC.A", "FSC.W", "SSC.A", "SSC.W"), fSOM$markers)

fSOM <- BuildSOM(fSOM, colsToUse = mset, xdim = 4, ydim = 4)
fSOM <- BuildMST(fSOM)

nbclust <- 8
fSOM[["metaclustering"]] <- as.factor(metaClustering_consensus(fSOM$map$codes, k = nbclust))
PlotStars(fSOM, backgroundValues = fSOM$metaclustering)

library(pheatmap)
fSOM$heatmap <- pheatmap(scale(fsom.clusters.heatmap(fSOM)), scale = "row")
library(shiny)
library(ggplot2)
library(viridis)
fcs.explorer(fsom.dat(fSOM), "FSC.A", "SSC.A", 20000)

fSOM <- nodes.to.meta_fSOM(fSOM, grep("1|2|3", fSOM$metaclustering))  # stochastic

fSOM$heatmap <- pheatmap(scale(fsom.clusters.heatmap(fSOM)), scale = "row")
fcs.explorer(fsom.dat(fSOM), "FSC.A", "SSC.A", 20000)

fSOM <- nodes.to.meta_fSOM(fSOM, grep("1|6", fSOM$metaclustering))  # stochastic

fSOM$heatmap <- pheatmap(scale(fsom.clusters.heatmap(fSOM)), scale = "row")
fcs.explorer(fsom.dat(fSOM), "FSC.A", "SSC.A", 20000)

fSOM$metaclustering.anno <- vector(mode = "list", length = length(levels(fSOM$metaclustering)))
fSOM$metaclustering.anno[[5]] <- "singlets"

savefsom(fSOM, file.path("./results_R/fSOMs/TPHE", paste(Sys.Date(), "fSOM_CD8p_SINGLETS.PASS2.rds", sep = ".")))
