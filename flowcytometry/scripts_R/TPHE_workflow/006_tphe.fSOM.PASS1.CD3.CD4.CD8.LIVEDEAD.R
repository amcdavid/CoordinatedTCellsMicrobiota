## fSOM Build using subsampled files
## !!!flowSOM clustering is stochastic, so can't be 100% reproduced (can it be 'seeded' ?)!!!
library(flowCore)
library(FlowSOM)

## file paths
tphe.warped.files <- list.files("./data_modified/TPHE_warping/CD3.CD4.CD8.LIVEDEAD_warped_SUBSAMPLED/", full.names = TRUE)

input.set <- fset.compensate(read.flowSet(tphe.warped.files, transformation = FALSE))

trans.obj <- readRDS("./results_R/FCS.transforms/TPHE/asinh.transforms.TPHE.rds")
trans.obj <- trans.obj[grep('Violet G|Violet E|Violet A|Green D', names(trans.obj))]

input.set.transformed <- transform.set(input.set)

fSOM <- ReadInput(input.set.transformed, scale = TRUE)
fSOM$markers <- get.markers(input.set.transformed[[1]])

mset <- match(c("CD3", "CD4", "CD8A", "LIVEDEAD"), fSOM$markers)

fSOM <- BuildSOM(fSOM, colsToUse = mset, xdim = 8, ydim = 8) # PASS 1
fSOM <- BuildMST(fSOM)

nbclust <- 12
fSOM[["metaclustering"]] <- as.factor(metaClustering_consensus(fSOM$map$codes, k = nbclust))
PlotStars(fSOM, backgroundValues = fSOM$metaclustering)

fSOM$heatmap <- pheatmap(scale(fsom.clusters.heatmap(fSOM)), scale = "row")
fcs.explorer(fsom.dat(fSOM), "FSC.A", "LIVEDEAD", 20000)

## merge nodes based on visual inspection
fSOM <- nodes.to.meta_fSOM(fSOM, 8)
fSOM <- nodes.to.meta_fSOM(fSOM, 29) 
fSOM <- nodes.to.meta_fSOM(fSOM, 54)

fSOM$heatmap <- pheatmap(scale(fsom.clusters.heatmap(fSOM)), scale = "row")
fcs.explorer(fsom.dat(fSOM), "CD4", "CD8A", 20000)

fSOM$metaclustering.anno <- vector(mode = "list", length = length(levels(fSOM$metaclustering)))
fSOM$metaclustering.anno[[3]] <- "CD4+"
fSOM$metaclustering.anno[[10]] <- "CD8+"

fSOM$transform.object <- trans.obj #!!! need to store this for downstream functions to work correctly

savefsom(fSOM, file.path("./results_R/fSOMs/TPHE", paste(Sys.Date(), "fSOM_CD3.CD4.CD8.LIVEDEAD.PASS1.rds", sep = ".")))
