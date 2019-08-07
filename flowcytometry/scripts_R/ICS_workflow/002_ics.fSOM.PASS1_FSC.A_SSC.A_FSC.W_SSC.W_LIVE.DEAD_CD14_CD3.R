## fSOM PASS1 build using subsampled files
## !!!flowSOM clustering is stochastic, so can't be 100% reproduced!!!
library(flowCore)
library(FlowSOM)
library(pheatmap)
library(shiny)
library(ggplot2)
library(viridis)

files.withoutHW <- paste(sub(".fcs", "", list.files("./data_source/ICS2_with.spillover/NGL046/", pattern = "*SEB.fcs$")), collapse = "|")# NGL046 has missing parameters
input.files <- list.files("./data_modified/ICS/SEB_subsampled/", full.names = TRUE)
input.files <- input.files[-grep(files.withoutHW, input.files)]

input.flowset <- read.flowSet(input.files, transformation = FALSE)
input.flowset <- fset.trim(input.flowset)
input.flowset <- fset.compensate(input.flowset)

elgcl <- elgcl.transform(input.flowset, 500000)

input.flowset <- fset.transform(input.flowset)# requires elgcl in the environment

##
fSOM <- ReadInput(input.flowset, scale = TRUE)# can take some time due to scaling and creation of $metaData
fSOM$markers <- get.markers(input.flowset[[1]])
fSOM$markers[grep("CD14", fSOM$markers)] <- "LIVE.DEAD_CD14"
fSOM$transform.object <- elgcl

fSOM <- BuildSOM(fSOM, 
                 colsToUse = match(c("FSC.A", "SSC.A", "FSC.W", "SSC.W", "LIVE.DEAD_CD14", "CD3", "CD69"), fSOM$markers), 
                 xdim = 10, 
                 ydim = 10
                 )# fSOM PASS 1: set columns to cluster for CD3+ LIVE singlets
fSOM <- BuildMST(fSOM)

nbclust <- 24
fSOM[["metaclustering"]] <- as.factor(metaClustering_consensus(fSOM$map$codes, k = nbclust))
PlotStars(fSOM, backgroundValues = fSOM$metaclustering)

fSOM$heatmap <- pheatmap(scale(fsom.clusters.heatmap(fSOM)), scale = "row")
fcs.explorer(fsom.dat(fSOM), "LIVE.DEAD_CD14", "FSC.A", 20000)

## merge nodes based on visual inspection for CD3+ CD69-/+ LIVE; these metacluster numbers may change due to seed
fSOM <- nodes.to.meta_fSOM(fSOM, grep("22|23", fSOM$metaclustering))

fSOM$heatmap <- pheatmap(scale(fsom.clusters.heatmap(fSOM)), scale = "row")
fcs.explorer(fsom.dat(fSOM), "LIVE.DEAD_CD14", "FSC.A", 20000)

##annotate
fSOM$metaclustering.anno <- vector(mode = "list", length = length(levels(fSOM$metaclustering)))
fSOM$metaclustering.anno[[23]] <- "CD3+ LIVE"

##save work
saveRDS(fSOM, file.path("./results_R/fSOMs/ICS", paste0("ICS_fSOM_PASS1_", paste(fSOM$markers[fSOM$map$colsUsed], collapse = "_"), "_ALLDATA",".rds")))
savefsom(fSOM, file.path("./results_R/fSOMs/ICS", paste0("ICS_fSOM_PASS1_", paste(fSOM$markers[fSOM$map$colsUsed], collapse = "_"), "_", paste(Sys.Date()),".rds")))
