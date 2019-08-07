## fSOM Build using subsampled files

## file paths: warped samples
tphe.warped.files <- list.files("./data_modified/TPHE_warping/PERFORIN.FOXP3.CD185.CD31.GRZB_warped_CD4p_SUBSAMPLED/", full.names = TRUE)
bad.CD28 <- c("RPRC0530020_7_TPHE|RPRC0530091_19_TPHE|RPRC0540079_1_TPHE|RPRC0540075_7_TPHE") #(identified during QC)
tphe.warped.files <- tphe.warped.files[-grep(bad.CD28, tphe.warped.files)]
tphe.warped.files <- tphe.warped.files[-grep("HD0191", tphe.warped.files)]

input.set <- fset.compensate(read.flowSet(tphe.warped.files, transformation = FALSE))

trans.obj <- readRDS("./results_R/FCS.transforms/TPHE/asinh.transforms.TPHE.rds")

input.set.transformed <- transform.set(input.set)

fSOM <- ReadInput(input.set.transformed, scale = TRUE)
fSOM$markers <- get.markers(input.set.transformed[[1]])

mset <- match(c("CD57", "CD28", "FOXP3", "CD197", "CD185", "CD45RO", "CD127", "CD31"), fSOM$markers)

fSOM <- BuildSOM(fSOM, colsToUse = mset, xdim= 10, ydim = 10, rlen = 20) # PASS 1
fSOM <- BuildMST(fSOM)

nbclust <- 30
fSOM[["metaclustering"]] <- as.factor(metaClustering_consensus(fSOM$map$codes, k = nbclust))
PlotStars(fSOM, backgroundValues = fSOM$metaclustering)

fSOM$heatmap <- pheatmap(scale(fsom.clusters.heatmap(fSOM)), scale = "row")
fcs.explorer(fsom.dat(fSOM), "CD31", "CD197", 20000)

fSOM <- nodes.to.meta_fSOM(fSOM, grep("^4|30", fSOM$metaclustering))  # stochastic
fSOM$heatmap <- pheatmap(scale(fsom.clusters.heatmap(fSOM)), scale = "row")
fcs.explorer(fsom.dat(fSOM), "CD31", "CD197", 20000)

fSOM <- nodes.to.meta_fSOM(fSOM, c(28, 97))  # stochastic
fSOM$heatmap <- pheatmap(scale(fsom.clusters.heatmap(fSOM)), scale = "row")
fcs.explorer(fsom.dat(fSOM), "CD31", "CD197", 20000)

fSOM <- nodes.to.meta_fSOM(fSOM, c(2, 19))  # stochastic
fSOM$heatmap <- pheatmap(scale(fsom.clusters.heatmap(fSOM)), scale = "row")
fcs.explorer(fsom.dat(fSOM), "CD31", "CD197", 20000)

fSOM <- nodes.to.meta_fSOM(fSOM, grep("^10|^13", fSOM$metaclustering))  # stochastic
fSOM$heatmap <- pheatmap(scale(fsom.clusters.heatmap(fSOM)), scale = "row")
fcs.explorer(fsom.dat(fSOM), "CD31", "CD197", 20000)

fSOM <- nodes.to.meta_fSOM(fSOM, grep("^18|^25", fSOM$metaclustering))  # stochastic
fSOM$heatmap <- pheatmap(scale(fsom.clusters.heatmap(fSOM)), scale = "row")
fcs.explorer(fsom.dat(fSOM), "CD31", "CD197", 20000)

fSOM <- nodes.to.meta_fSOM(fSOM, grep("^26|^29", fSOM$metaclustering))  # stochastic
fSOM$heatmap <- pheatmap(scale(fsom.clusters.heatmap(fSOM)), scale = "row")
fcs.explorer(fsom.dat(fSOM), "CD31", "CD197", 20000)

fSOM <- nodes.to.meta_fSOM(fSOM, grep("^12|^16", fSOM$metaclustering))  # stochastic
fSOM$heatmap <- pheatmap(scale(fsom.clusters.heatmap(fSOM)), scale = "row")
fcs.explorer(fsom.dat(fSOM), "CD31", "CD197", 20000)


fSOM$metaclustering.anno <- vector(mode = "list", length = length(levels(fSOM$metaclustering)))

fSOM$transform.object <- trans.obj

## fSOM with 19 elements
savefsom(fSOM, file.path("./results_R/fSOMs/TPHE", paste(Sys.Date(), "fSOM_PASS3_CD4p_rerun.rds", sep = ".")))
saveRDS(fSOM, file.path("./results_R/fSOMs/TPHE", paste(Sys.Date(), "fSOM_PASS3_CD4p_rerun_ALLDATA.rds", sep = ".")))
