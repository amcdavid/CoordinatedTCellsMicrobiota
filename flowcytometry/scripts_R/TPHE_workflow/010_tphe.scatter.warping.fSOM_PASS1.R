library(flowCore)
library(flowStats)
library(flowViz)
##

#### required objects/arguments
##
warp.parameters <- readRDS("./results_R/FCS.transforms/TPHE/TPHE.warp.scatter.parameters.rds") # channels-of-interest

## file paths; set pattern argument for either CD4 or CD8
fsom.files <- list.files("./data_modified/TPHE_warping/fSOM_CD3.CD4.CD8.LIVEDEAD_PASS1/", full.names = TRUE, pattern = "fSOM_CD8p_.fcs")
batch.names <- paste0("NGL0", seq(47, 91, 2)) ; batch.names <- batch.names[-grep("NGL063", batch.names)]

## warping target
target.file <- fsom.files[grep("HD0189_TPHE_NGL091", fsom.files)]
target <- "HD0189_TPHE_NGL091"; fsom.files <- fsom.files[-grep(target, fsom.files)]


#############################################################################################################################################################################################
#############################################################################################################################################################################################
## test-bed for QC-ing/visualizing warping per batch

## read in batched set
batch.name <- "NGL075"

input.set <- read.flowSet(c(fsom.files[grep(batch.name, fsom.files)], target.file), transformation = FALSE)
input.set.warped <- warp.set(input.set, target) # warp using 'warp.parameters'

## warp QC plots
warp.qc.plots(input.set, input.set.warped)

## modify warped files to include new identifiers
## change/add keywords
input.set.warped.mod <- fsApply(input.set.warped, function (i) {
  i@description$`EXPERIMENT NAME` <- unlist(strsplit(i@description$`EXPERIMENT NAME`, " "))[[1]]
  i@description$`WARPED CHANNELS` <- paste(i@description$`WARPED CHANNELS`, paste(names(warp.parameters), collapse = "_"), sep = ".")
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
outDir = file.path("./data_modified/TPHE_warping/fSOM_CD3.CD4.CD8.LIVEDEAD_PASS1/scatter_warped")

fsApply(input.set.warped.mod[-grep(target, sampleNames(input.set.warped.mod))], function(i) write.FCS(x = i, filename = file.path(outDir, i@description$`$FIL`)))


#############################################################################################################################################################################################
#############################################################################################################################################################################################
## looping function for warping per-batch/per-channel; set i@description$`WARPED` name to respective population; set 'outDir' path name

for(i in seq(batch.names)){
  
  batch.name <- batch.names[i]
  print(batch.name)
  
  input.set <- read.flowSet(c(fsom.files[grep(batch.name, fsom.files)], target.file), transformation = FALSE)
  input.set.warped <- warp.set(input.set, target) # warp using 'warp.parameters'
  
  warp.qc.plots(input.set, input.set.warped)

  ## modify warped files to include new identifiers
  ## change/add keywords
  input.set.warped <- fsApply(input.set.warped, function (i) {
    i@description$`EXPERIMENT NAME` <- unlist(strsplit(i@description$`EXPERIMENT NAME`, " "))[[1]]
    i@description$`WARPED CHANNELS` <- paste(i@description$`WARPED CHANNELS`, paste(names(warp.parameters), collapse = "_"), sep = ".")
    i@description$`WARP TARGET` <- sampleNames(input.set.warped)[grep(target, sampleNames(input.set.warped))]
    i@description$`WARPED` <- "fda_WARPED_fSOM_CD8p" # set name here
    i@description$`$FIL` <- paste(i@description$`$SRC`, 
                                  i@description$`TUBE NAME`, 
                                  i@description$`EXPERIMENT NAME`, 
                                  i@description$WARPED,
                                  ".fcs",
                                  sep = "_")
    i
  })
  
  ## write modified files
  outDir = file.path("./data_modified/TPHE_warping/fSOM_CD3.CD4.CD8.LIVEDEAD_PASS1/scatter_warped/CD8")
  print("Writing scatter-warped FCS files")
  fsApply(input.set.warped[-grep(target, sampleNames(input.set.warped))], function(i) write.FCS(x = i, filename = file.path(outDir, i@description$`$FIL`)))
}
