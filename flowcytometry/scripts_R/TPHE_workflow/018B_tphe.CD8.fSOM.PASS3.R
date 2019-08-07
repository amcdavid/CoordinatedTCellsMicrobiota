## fSOM Build using warped files
library(flowCore); library(FlowSOM); library(pheatmap); library(shiny); library(ggplot2); library(viridis)

## file paths: warped samples
tphe.warped.files <- list.files("./data_modified/TPHE_warping/EXPERIMENTS/", recursive = TRUE, full.names = TRUE, pattern = "CD8p.SINGLETS_.fcs")
tphe.warped.files <- tphe.warped.files[grep("PERFORIN.CD57.CD28.FOXP3.CD197.CD185.CD31.CD122.GRZB", dirname(tphe.warped.files))]

samples.exclude.from.CD8fsom <- c("RPRC0530043_7_|RPRC0530043_1_|RPRC0530043_19_|RPRC0530042_7_|RPRC0530042_1_|RPRC0530042_19_|RPRC0530041_7_|RPRC0530041_1_|RPRC0530041_19_|RPRC0530059_7")
experiments.exclude.from.CD8fsom <- c("NGL047|NGL049|NGL051|NGL057|NGL059")

tphe.warped.files <- tphe.warped.files[-grep(samples.exclude.from.CD8fsom, tphe.warped.files)]
tphe.warped.files <- tphe.warped.files[-grep(experiments.exclude.from.CD8fsom, tphe.warped.files)]

input.set <- fset.compensate(read.flowSet(tphe.warped.files, transformation = FALSE))

trans.obj <- readRDS("./results_R/FCS.transforms/TPHE/asinh.transforms.TPHE.rds")

input.set.transformed <- transform.set(input.set)

fSOM <- ReadInput(input.set.transformed, scale = TRUE)
fSOM$markers <- get.markers(input.set.transformed[[1]])

saveRDS(fSOM, "./temp/fSOM_CD8.TPHE_PASS3_temp.RDS")

mset <- match(c("PERFORIN", "CD122", "CD57", "CD28", "FOXP3", "CD197", "CD185", "CD45RO", "CD127", "CD31", "GRZB"), fSOM$markers)

fSOM <- BuildSOM(fSOM, colsToUse = mset, xdim= 6, ydim = 6, rlen = 20) # PASS 3
fSOM <- BuildMST(fSOM)

nbclust <- 20
fSOM[["metaclustering"]] <- as.factor(metaClustering_consensus(fSOM$map$codes, k = nbclust))
PlotStars(fSOM, backgroundValues = fSOM$metaclustering)

fSOM$heatmap <- pheatmap(scale(fsom.clusters.heatmap(fSOM)), scale = "row")
fcs.explorer(fsom.dat(fSOM), "CD31", "CD197", 20000)

fSOM <- nodes.to.meta_fSOM(fSOM, grep("29|30", fSOM$metaclustering))
fSOM$heatmap <- pheatmap(scale(fsom.clusters.heatmap(fSOM)), scale = "row")
fcs.explorer(fsom.dat(fSOM), "CD57", "PERFORIN", 20000)

fSOM <- nodes.to.meta_fSOM(fSOM, grep("^1$|11", fSOM$metaclustering))
fSOM$heatmap <- pheatmap(scale(fsom.clusters.heatmap(fSOM)), scale = "row")
fcs.explorer(fsom.dat(fSOM), "CD57", "PERFORIN", 20000)

fSOM <- nodes.to.meta_fSOM(fSOM, c(15, 26))
fSOM$heatmap <- pheatmap(scale(fsom.clusters.heatmap(fSOM)), scale = "row")
fcs.explorer(fsom.dat(fSOM), "CD57", "PERFORIN", 20000)

fSOM <- nodes.to.meta_fSOM(fSOM, grep("^7$|26", fSOM$metaclustering))
fSOM$heatmap <- pheatmap(scale(fsom.clusters.heatmap(fSOM)), scale = "row")
fcs.explorer(fsom.dat(fSOM), "CD57", "PERFORIN", 20000)

fSOM <- nodes.to.meta_fSOM(fSOM, grep("19|28", fSOM$metaclustering))
fSOM$heatmap <- pheatmap(scale(fsom.clusters.heatmap(fSOM)), scale = "row")
fcs.explorer(fsom.dat(fSOM), "CD57", "PERFORIN", 20000)

fSOM <- nodes.to.meta_fSOM(fSOM, grep("19|28", fSOM$metaclustering))
fSOM$heatmap <- pheatmap(scale(fsom.clusters.heatmap(fSOM)), scale = "row")
fcs.explorer(fsom.dat(fSOM), "CD57", "PERFORIN", 20000)

fSOM <- nodes.to.meta_fSOM(fSOM, grep("^6$|20", fSOM$metaclustering))
fSOM$heatmap <- pheatmap(scale(fsom.clusters.heatmap(fSOM)), scale = "row")
fcs.explorer(fsom.dat(fSOM), "GRZB", "CD31", 20000)

fSOM <- nodes.to.meta_fSOM(fSOM, grep("10|26", fSOM$metaclustering))
fSOM$heatmap <- pheatmap(scale(fsom.clusters.heatmap(fSOM)), scale = "row")
fcs.explorer(fsom.dat(fSOM), "GRZB", "CD28", 20000)

fSOM <- nodes.to.meta_fSOM(fSOM, grep("24|25", fSOM$metaclustering))
fSOM$heatmap <- pheatmap(scale(fsom.clusters.heatmap(fSOM)), scale = "row")
fcs.explorer(fsom.dat(fSOM), "GRZB", "CD28", 20000)


fSOM$metaclustering.anno <- vector(mode = "list", length = length(levels(fSOM$metaclustering)))

fSOM$transform.object <- trans.obj

savefsom(fSOM, file.path("./results_R/fSOMs/TPHE", paste(Sys.Date(), "fSOM_PASS3_CD8p.rds", sep = ".")))
saveRDS(fSOM, "./temp/fSOM_PASS3_CD8p_ALLDATA.rds")
