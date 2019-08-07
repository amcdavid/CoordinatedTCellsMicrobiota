## warp files using previously defined parameters

#### required objects/arguments
##
trans.obj <- readRDS("./results_R/FCS.transforms/TPHE/asinh.transforms.TPHE.rds") # arcsinh transfomation parameters
trans.obj <- trans.obj[grep('Violet G|Violet E|Violet A|Green D', names(trans.obj))]# channels to warp
warp.markers <- c("LIVEDEAD", "CD3", "CD4", "CD8A") # channels to be warped
warp.parameters <- readRDS("./results_R/FCS.transforms/TPHE/TPHE.warp.parameters.rds") # channels-of-interest
warp.parameters <- warp.parameters[names(warp.parameters) %in% warp.markers] # warping parameters per channels-of-interest

## file paths
tphe.files <- list.files("./data_source/TPHE2_with.spillover", recursive = TRUE, full.names = TRUE, pattern = ".fcs")
tphe.files <- tphe.files[grep("HD|_TPHE.fcs", tphe.files)]

batches <- batch.list.edit(tphe.files[-grep("NGL063", tphe.files)])

## warping target
target.file <- tphe.files[grep("HD0189_TPHE_NGL091", tphe.files)]
target <- "HD0189_TPHE_NGL091"

#############################################################################################################################################################################################
#############################################################################################################################################################################################
## test-bed for QC-ing/visualizing warping per batch

library(flowCore)
library(flowStats)
library(flowViz)

## read in batched set
batch.name <- "NGL091"

input.set <- fset.compensate(fset.trim(read.flowSet(c(batches[batch.name][[1]], target.file), transformation = FALSE))) # read in flowSet; compensate
input.set <- transform.set(input.set) # transform using 'trans.obj'
input.set.warped <- warp.set(input.set, target) # warp using 'warp.parameters'

## warp QC plots
warp.qc.plots(input.set, input.set.warped)

## inverse transform warped set
input.set.warped.inverse <- inverse.transform.set(input.set.warped)

## decomp
input.set.warped.inverse.decomp <- fset.decompensate(input.set.warped.inverse)

## modify transformed files to include new identifiers
## change/add keywords
input.set.warped.inverse.decomp.mod <- fsApply(input.set.warped.inverse.decomp, function (i) {
  i@description$`EXPERIMENT NAME` <- unlist(strsplit(i@description$`EXPERIMENT NAME`, " "))[[1]]
  i@description$`WARPED CHANNELS` <- paste(i@description$`WARPED CHANNELS`, paste(names(warp.parameters), collapse = "."), sep = ".")
  i@description$`WARP TARGET` <- sampleNames(input.set.warped.inverse.decomp)[grep(target, sampleNames(input.set.warped.inverse.decomp))]
  i@description$`WARPED` <- "fda_WARPED"
  i@description$`$FIL` <- paste(i@description$`$SRC`, 
                                i@description$`TUBE NAME`, 
                                i@description$`EXPERIMENT NAME`, 
                                i@description$WARPED,
                                ".fcs",
                                sep = "_")
  i
})

## write modified files
outDir = file.path("./data_modified/TPHE_warping", batch.name, "CD4.CD8A.CD3.LIVEDEAD_warped")

fsApply(input.set.warped.inverse.decomp.mod[-grep(target, sampleNames(input.set.warped.inverse.decomp.mod))], function(i) write.FCS(x = i, filename = file.path(outDir, i@description$`$FIL`)))

#############################################################################################################################################################################################
#############################################################################################################################################################################################
## looping function for warping per-batch/per-channel

for(i in seq(names(batches))){
  
  batch.name <- names(batches)[i]
  print(batch.name)
  
  input.set <- fset.compensate(fset.trim(read.flowSet(c(batches[batch.name][[1]], target.file), transformation = FALSE))) # read in flowSet; compensate
  input.set <- transform.set(input.set) # transform using 'trans.obj'
  input.set.warped <- warp.set(input.set, target) # warp using 'warp.parameters'
  
  warp.qc.plots(input.set, input.set.warped)
  
  ## inverse transform warped set
  input.set.warped <- inverse.transform.set(input.set.warped)
  
  ## decomp
  input.set.warped <- fset.decompensate(input.set.warped)
  
  ## modify transformed files to include new identifiers
  ## change/add keywords
  input.set.warped <- fsApply(input.set.warped, function (i) {
    i@description$`EXPERIMENT NAME` <- unlist(strsplit(i@description$`EXPERIMENT NAME`, " "))[[1]]
    i@description$`WARPED CHANNELS` <- paste(i@description$`WARPED CHANNELS`, paste(names(warp.parameters), collapse = "."), sep = ".")
    i@description$`WARP TARGET` <- sampleNames(input.set.warped)[grep(target, sampleNames(input.set.warped))]
    i@description$`WARPED` <- "fda_WARPED"
    i@description$`$FIL` <- paste(i@description$`$SRC`, 
                                  i@description$`TUBE NAME`, 
                                  i@description$`EXPERIMENT NAME`, 
                                  i@description$WARPED,
                                  ".fcs",
                                  sep = "_")
    i
  })
  
  ## write modified files
  outDir = file.path("./data_modified/TPHE_warping", batch.name, "CD4.CD8A.CD3.LIVEDEAD_warped")
  print("Writing warped FCS files")
  fsApply(input.set.warped[-grep(target, sampleNames(input.set.warped))], function(i) write.FCS(x = i, filename = file.path(outDir, i@description$`$FIL`)))
}
