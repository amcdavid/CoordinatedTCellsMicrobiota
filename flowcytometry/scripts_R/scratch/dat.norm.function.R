## source flow-related libraries and functions
source("./source_R/flow.src.R", chdir = TRUE)

## file paths
tphe.original.files <- list.files("./data_source/TPHE2_with.spillover", recursive = TRUE, full.names = TRUE, pattern = "_TPHE.fcs")

batches <- batch.list(tphe.original.files)

## warping targets
targets <- list.files("./data_modified/TPHE_warping/targets_HD/scatter.VG.VE.GD.VA.warped/", full.names = TRUE)

## transforms
trans.obj <- readRDS("./results_R/FCS.transforms/TPHE/asinh.transforms.TPHE.rds")

##
dat.norm <- function(batch.name){
  
  input.set <- fset.compensate(fset.trim(read.flowSet(c(batches[batch.name][[1]], targets[grep(batch.name, targets)]), 
                                                    transformation = FALSE)))
  
  input.set <- transform(input.set, trans.obj$`Violet G 550_40-A`$asinh.translist)
  input.set <- transform(input.set, trans.obj$`Violet E 585_42-A`$asinh.translist)
  input.set <- transform(input.set, trans.obj$`Green D 610_20-A`$asinh.translist)
  input.set <- transform(input.set, trans.obj$`Violet A 780_60-A`$asinh.translist)

  input.set <- fset.rangefix(input.set, "Violet G 550_40-A")
  input.set <- fset.rangefix(input.set, "Violet E 585_42-A")
  input.set <- fset.rangefix(input.set, "Green D 610_20-A")
  input.set <- fset.rangefix(input.set, "Violet A 780_60-A")
  
  input.set <- warpSet(input.set, c("Violet G 550_40-A"), monwrd = FALSE, peakNr = 2, clipRange = 0.01, bwFac = 1.5, 
                              target = sampleNames(input.set)[grep("HD", sampleNames(input.set))])
  input.set <- warpSet(input.set, c("Violet E 585_42-A"), monwrd = FALSE, peakNr = 2, clipRange = 0.01, bwFac = 1.5, 
                              target = sampleNames(input.set)[grep("HD", sampleNames(input.set))])
  input.set <- warpSet(input.set, c("Green D 610_20-A"), monwrd = FALSE, peakNr = 2, clipRange = 0.01, bwFac = 1.5, 
                              target = sampleNames(input.set)[grep("HD", sampleNames(input.set))])
  input.set <- warpSet(input.set, c("Violet A 780_60-A"), monwrd = FALSE, peakNr = 2, clipRange = 0.2, bwFac = 2.5, 
                              target = sampleNames(input.set)[grep("HD", sampleNames(input.set))])
  
  input.set
}

batch.name <- "NGL069"

input.set.warped <- dat.norm(batch.name)

## warp QC
channel <- "Violet G 550_40-A"
densityplot(~ ., fset.downsample(input.set.warped, 5000)[, channel], filter = curv1Filter(channel, bwFac = 1.5))

channel <- "Violet E 585_42-A"
densityplot(~ ., fset.downsample(input.set.warped, 5000)[, channel], filter = curv1Filter(channel, bwFac = 1.5))

channel <- "Green D 610_20-A"
densityplot(~ ., fset.downsample(input.set.warped, 5000)[, channel], filter = curv1Filter(channel, bwFac = 1.5))

channel <- "Violet A 780_60-A"
densityplot(~ ., fset.downsample(input.set.warped, 5000)[, channel], filter = curv1Filter(channel, bwFac = 2))

## inverse transform warped set
dat.deprocess <- function(transformed.flowSet){
  
  channel <- "Violet G 550_40-A"
  sinhTrans <- sinhTransform(transformationId = "inverse.arcsinh-transformation",
                             a = trans.obj[[channel]]$asinh.vals$a,
                             b = trans.obj[[channel]]$asinh.vals$b, 
                             c = trans.obj[[channel]]$asinh.vals$c)
  sinh.translist <- transformList(channel, sinhTrans)

    transformed.flowSet <- transform(transformed.flowSet, sinh.translist)

  channel <- "Violet E 585_42-A"
  sinhTrans <- sinhTransform(transformationId = "inverse.arcsinh-transformation", 
                             a = trans.obj[[channel]]$asinh.vals$a, 
                             b = trans.obj[[channel]]$asinh.vals$b, 
                             c = trans.obj[[channel]]$asinh.vals$c)
  sinh.translist <- transformList(channel, sinhTrans)

    transformed.flowSet <- transform(transformed.flowSet, sinh.translist)

  channel <- "Green D 610_20-A"
  sinhTrans <- sinhTransform(transformationId = "inverse.arcsinh-transformation", 
                             a = trans.obj[[channel]]$asinh.vals$a, 
                             b = trans.obj[[channel]]$asinh.vals$b, 
                             c = trans.obj[[channel]]$asinh.vals$c)
  sinh.translist <- transformList(channel, sinhTrans)

    transformed.flowSet <- transform(transformed.flowSet, sinh.translist)

  channel <- "Violet A 780_60-A"
  sinhTrans <- sinhTransform(transformationId = "inverse.arcsinh-transformation", 
                             a = trans.obj[[channel]]$asinh.vals$a, 
                             b = trans.obj[[channel]]$asinh.vals$b, 
                             c = trans.obj[[channel]]$asinh.vals$c)
  sinh.translist <- transformList(channel, sinhTrans)

    transformed.flowSet <- transform(transformed.flowSet, sinh.translist)
    
    
    transformed.flowSet <- fset.decompensate(transformed.flowSet)
    
    transformed.flowSet
}

input.set.warped.inverse.decomp <- dat.deprocess(input.set.warped)

## modify transformed files to include new identifiers
## change/add keywords
input.set.warped.inverse.decomp.mod <- fsApply(input.set.warped.inverse.decomp, function (i) {
  i@description$`EXPERIMENT NAME` <- unlist(strsplit(i@description$`EXPERIMENT NAME`, " "))[[1]]
  i@description$`WARPED CHANNELS` <- "VG_VE_VA_GD"
  i@description$`WARP TARGET` <- sampleNames(input.set.warped.inverse.decomp)[grep("HD", sampleNames(input.set.warped.inverse.decomp))]
  i@description$`WARPED` <- "fda_WARPED"
  i@description$`$FIL` <- paste(i@description$`$SRC`, 
                                i@description$`TUBE NAME`, 
                                i@description$`EXPERIMENT NAME`, 
                                i@description$WARPED, 
                                i@description$`WARPED CHANNELS`,
                                ".fcs",
                                sep = "_")
  i
})

## write modified files
outDir = file.path("./data_modified/TPHE_warping", batch.name, "VG_VE_VA_GD")

fsApply(input.set.warped.inverse.decomp.mod[-grep("HD", sampleNames(input.set.warped.inverse.decomp.mod))], function(i) write.FCS(x = i, filename = file.path(outDir, i@description$`$FIL`)))

