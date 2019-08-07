## source flow-related libraries and functions
source("./source_R/flow.src.R", chdir = TRUE)

## required objects/arguments
trans.obj <- readRDS("./results_R/FCS.transforms/TPHE/asinh.transforms.TPHE.rds") # arcsinh transfomation parameters
trans.obj <- trans.obj[grep('Violet G|Violet E|Violet A|Green D', names(trans.obj))]
warp.markers <- c("LIVEDEAD", "CD3", "CD4", "CD8A") # channels to be warped
warp.parameters <- readRDS("./results_R/FCS.transforms/TPHE/TPHE.warp.parameters.rds") ; warp.parameters <- warp.parameters[names(warp.parameters) %in% warp.markers] # warping parameters per-channel

## file paths
tphe.files <- list.files("./data_source/TPHE2_with.spillover/", recursive = TRUE, full.names = TRUE, pattern = "_TPHE.fcs")
batches <- batch.list.edit(tphe.files)
batches$NGL063_CD4CD31samechannel <- NULL

## warping target
target <- list.files("./data_modified/TPHE_warping/targets_HD/scatter.VG.VE.GD.VA.warped/", full.names = TRUE) ; target <- target[grep("HD0189_TPHE_NGL091", target)]

#############################################################################################################################################################################################
#############################################################################################################################################################################################

## read in batched set
batch.name <- "NGL061"

input.set <- fset.compensate(fset.trim(read.flowSet(c(batches[batch.name][[1]], target), transformation = FALSE))) # read in flowSet; compensate
input.set <- transform.set(input.set) # transform using 'trans.obj'
input.set.warped <- warp.set(input.set) # warp using 'warp.parameters'

## warp QC plots
qc.plots <- warp.qc.plots(input.set, input.set.warped)

## inverse transform warped set
input.set.warped.inverse <- inverse.transform.set(input.set.warped)

## decomp
input.set.warped.inverse.decomp <- fset.decompensate(input.set.warped.inverse)

## modify transformed files to include new identifiers
## change/add keywords
input.set.warped.inverse.decomp.mod <- fsApply(input.set.warped.inverse.decomp, function (i) {
  i@description$`EXPERIMENT NAME` <- unlist(strsplit(i@description$`EXPERIMENT NAME`, " "))[[1]]
  i@description$`WARPED CHANNELS` <- paste(names(warp.parameters), collapse = ".")
  i@description$`WARP TARGET` <- sampleNames(input.set.warped.inverse.decomp)[grep("HD", sampleNames(input.set.warped.inverse.decomp))]
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

fsApply(input.set.warped.inverse.decomp.mod[-grep("HD", sampleNames(input.set.warped.inverse.decomp.mod))], function(i) write.FCS(x = i, filename = file.path(outDir, i@description$`$FIL`)))

#############################################################################################################################################################################################
#############################################################################################################################################################################################
## good luck

for(i in seq(names(batches))){
  
  batch.name <- names(batches)[i]
  print(batch.name)
  
  input.set <- fset.compensate(fset.trim(read.flowSet(c(batches[batch.name][[1]], target), transformation = FALSE))) # read in flowSet; compensate
  input.set <- transform.set(input.set) # transform using 'trans.obj'
  input.set.warped <- warp.set(input.set) # warp using 'warp.parameters'
  
  qc.plots <- warp.qc.plots(input.set, input.set.warped)
  
  ## inverse transform warped set
  input.set.warped <- inverse.transform.set(input.set.warped)
  
  ## decomp
  input.set.warped <- fset.decompensate(input.set.warped)
  
  ## modify transformed files to include new identifiers
  ## change/add keywords
  input.set.warped <- fsApply(input.set.warped, function (i) {
    i@description$`EXPERIMENT NAME` <- unlist(strsplit(i@description$`EXPERIMENT NAME`, " "))[[1]]
    i@description$`WARPED CHANNELS` <- paste(names(warp.parameters), collapse = ".")
    i@description$`WARP TARGET` <- sampleNames(input.set.warped)[grep("HD", sampleNames(input.set.warped))]
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
  fsApply(input.set.warped[-grep("HD", sampleNames(input.set.warped))], function(i) write.FCS(x = i, filename = file.path(outDir, i@description$`$FIL`)))
}
