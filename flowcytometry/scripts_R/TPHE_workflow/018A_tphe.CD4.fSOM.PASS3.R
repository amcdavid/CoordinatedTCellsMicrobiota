## fSOM Build using subsampled files

## file paths: warped samples
tphe.warped.files <- list.files("./data_modified/TPHE_warping/PERFORIN.FOXP3.CD185.CD31.GRZB_warped_CD4p_SUBSAMPLED/", full.names = TRUE)

input.set <- fset.compensate(read.flowSet(tphe.warped.files, transformation = FALSE))

trans.obj <- readRDS("./results_R/FCS.transforms/TPHE/asinh.transforms.TPHE.rds")

input.set.transformed <- transform.set(input.set)

fSOM <- ReadInput(input.set.transformed, scale = TRUE)
fSOM$markers <- get.markers(input.set.transformed[[1]])

saveRDS(fSOM, "fSOM_TPHE_PASS3_temp.RDS")

mset <- match(c("CD57", "CD28", "FOXP3", "CD197", "CD185", "CD45RO", "CD127", "CD31"), fSOM$markers)

fSOM <- BuildSOM(fSOM, colsToUse = mset, xdim= 10, ydim = 10, rlen = 20) # PASS 3
fSOM <- BuildMST(fSOM)

nbclust <- 30
fSOM[["metaclustering"]] <- as.factor(metaClustering_consensus(fSOM$map$codes, k = nbclust))
PlotStars(fSOM, backgroundValues = fSOM$metaclustering)

fSOM$heatmap <- pheatmap(scale(fsom.clusters.heatmap(fSOM)), scale = "row")
fcs.explorer(fsom.dat(fSOM), "CD31", "CD197", 20000)

fSOM <- nodes.to.meta_fSOM(fSOM, "#")  # stochastic

fSOM$metaclustering.anno <- vector(mode = "list", length = length(levels(fSOM$metaclustering)))
fSOM$metaclustering.anno[[3]] <- "CD4+"
fSOM$metaclustering.anno[[10]] <- "CD8+"

fSOM$transform.object <- trans.obj

savefsom(fSOM, file.path("./results_R/fSOMs/TPHE", paste(Sys.Date(), "fSOM_PASS3_CD4p.rds", sep = ".")))
