## source flow-related libraries and functions
source("./source_R/flow.src.R", chdir = TRUE)

## file paths
tphe.warped.files <- list.files("./data_modified/TPHE_warping/targets_HD/scatter.VG.VE.GD.VA.warped/", recursive = TRUE, full.names = TRUE, pattern = ".fcs")

input.set <- fset.compensate(fset.trim(read.flowSet(tphe.warped.files)))

trans.obj <- readRDS("./results_R/FCS.transforms/TPHE/asinh.transforms.TPHE.rds")

input.set.trans <- input.set

input.set.trans <- transform(input.set.trans, trans.obj$`Violet G 550_40-A`$asinh.translist)
input.set.trans <- transform(input.set.trans, trans.obj$`Violet E 585_42-A`$asinh.translist)
input.set.trans <- transform(input.set.trans, trans.obj$`Green D 610_20-A`$asinh.translist)
input.set.trans <- transform(input.set.trans, trans.obj$`Violet A 780_60-A`$asinh.translist)

dat <- sample_n(as.data.frame(exprs(as(input.set.trans, "flowFrame"))), 100000)
markers <- get.markers(input.set.trans[[1]])
colnames(dat)[1:length(markers)] <- markers

fcs.explorer(dat, "CD4", "CD8A", 20000)

library(flowViz)
sub.set <- fset.downsample(input.set.trans, 10000)
channel = "Violet A 780_60-A"
densityplot(~ ., sub.set[, channel], filter = curv1Filter(channel, bwFac = 1))

fSOM <- ReadInput(input.set.trans, scale = TRUE)
fSOM$markers <- get.markers(input.set.trans[[1]])

mset <- match(c("LIVEDEAD", "CD3", "CD4", "CD8A"), fSOM$markers)

fSOM <- BuildSOM(fSOM, colsToUse = mset, xdim= 7, ydim = 7, rlen = 10) # PASS 1
fSOM <- BuildMST(fSOM)

nbclust <- 15
fSOM[["metaclustering"]] <- as.factor(metaClustering_consensus(fSOM$map$codes, k = nbclust))
PlotStars(fSOM, backgroundValues = fSOM$metaclustering)

fSOM$heatmap <- pheatmap(fsom.clusters.heatmap(fSOM), scale = "row")

fcs.explorer(fsom.dat(fSOM), "CD4", "CD8A", 20000)
