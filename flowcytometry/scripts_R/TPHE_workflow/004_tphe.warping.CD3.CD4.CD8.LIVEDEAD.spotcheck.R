## some files need fine-tuning of warping parameters (mainly bwfac value)
library(flowCore)
library(flowViz)
library(flowStats)

#### required objects/arguments
##
trans.obj <- readRDS("./results_R/FCS.transforms/TPHE/asinh.transforms.TPHE.rds") # arcsinh transfomation parameters
trans.obj <- trans.obj[grep('Violet G|Violet E|Violet A|Green D', names(trans.obj))]
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

## files to be 'fixed'; determined from reviewing 'warp.qc.plots' output

# CD8
fix.cd8 <- c("RPRC0540013_19_TPHE", "RPRC0540012_7_TPHE", "RPRC0530076_19_TPHE", "RPRC0530072_1_TPHE")

# CD3
fix.cd3 <- c("RPRC0540038_1_TPHE", "RPRC0540071_1_TPHE", "RPRC0540102_1_TPHE")

##
cd8.fix <- tphe.files[grep(paste(fix.cd8, collapse = "|"), tphe.files)]
cd3.fix <- tphe.files[grep(paste(fix.cd3, collapse = "|"), tphe.files)]

#############################################################################################################################################################################################
#############################################################################################################################################################################################
# modify the channel-specific bwfac value to appropriately normalize; 'warp.parameters$CD8A$bwFac <- new.val'

input.set <- fset.compensate(fset.trim(read.flowSet(c(cd8.fix, target.file), transformation = FALSE))) # read in flowSet; compensate
input.set <- transform.set(input.set) # transform using 'trans.obj'
input.set.warped <- warp.set(input.set, target) # warp using 'warp.parameters'

## warp QC plots
warp.qc.plots.edit <- function(input.set, input.set.warped){
  
  before.set <- input.set
  after.set <- input.set.warped
  
  qc.plots <- vector(mode = "list", length = length(warp.parameters)*2)
  
  pData(before.set)[[1]] <- fsApply(before.set, function(i) i@description$`$FIL` <- paste(i@description$`$SRC`, i@description$`EXPERIMENT NAME`, sep = "_"))
  pData(after.set)[[1]] <- fsApply(after.set, function(i) i@description$`$FIL` <- paste(i@description$`$SRC`, i@description$`EXPERIMENT NAME`, sep = "_"))
  
  before.set <- fset.downsample(before.set, 5000)
  after.set <- fset.downsample(after.set, 5000)
  
  for(i in seq(warp.parameters)){
    qc.plots[[i+(i-1)]] <- densityplot(~ ., 
                                       before.set[, warp.parameters[[i]]$channel], 
                                       filter = curv1Filter(warp.parameters[[i]]$channel, bwFac = warp.parameters[[i]]$bwFac),
                                       main = paste("Non-warped", names(warp.parameters)[i], sep = " "))
    qc.plots[[i+i]]     <- densityplot(~ ., 
                                       after.set[, warp.parameters[[i]]$channel], 
                                       filter = curv1Filter(warp.parameters[[i]]$channel, bwFac = warp.parameters[[i]]$bwFac),
                                       main = paste("Warped", names(warp.parameters)[i], sep = " "))
  }
  qc.plots
  
  doc.name <- paste("Warp.QC.densityplots", Sys.Date(), paste(names(warp.parameters), collapse = "."), "CD8.fix", "pdf", sep = ".")
  pdf(doc.name)
  print(qc.plots)
  dev.off()
}
warp.qc.plots.edit(input.set, input.set.warped)

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
outDir = file.path("./data_modified/TPHE_warping", "CD8.fix", "CD4.CD8A.CD3.LIVEDEAD_warped")
dir.create(dirname(outDir))
dir.create(outDir)
fsApply(input.set.warped.inverse.decomp.mod[-grep(target, sampleNames(input.set.warped.inverse.decomp.mod))], function(i) write.FCS(x = i, filename = file.path(outDir, i@description$`$FIL`)))


#############################################################################################################################################################################################
#############################################################################################################################################################################################
# modify the channel-specific bwfac value to appropriately normalize; 'warp.parameters$CD3$bwFac <- new.val'

input.set <- fset.compensate(fset.trim(read.flowSet(c(cd3.fix, target.file), transformation = FALSE))) # read in flowSet; compensate
input.set <- transform.set(input.set) # transform using 'trans.obj'
input.set.warped <- warp.set(input.set, target) # warp using 'warp.parameters'

## warp QC plots
warp.qc.plots.edit <- function(input.set, input.set.warped){
  
  before.set <- input.set
  after.set <- input.set.warped
  
  qc.plots <- vector(mode = "list", length = length(warp.parameters)*2)
  
  pData(before.set)[[1]] <- fsApply(before.set, function(i) i@description$`$FIL` <- paste(i@description$`$SRC`, i@description$`EXPERIMENT NAME`, sep = "_"))
  pData(after.set)[[1]] <- fsApply(after.set, function(i) i@description$`$FIL` <- paste(i@description$`$SRC`, i@description$`EXPERIMENT NAME`, sep = "_"))
  
  before.set <- fset.downsample(before.set, 5000)
  after.set <- fset.downsample(after.set, 5000)
  
  for(i in seq(warp.parameters)){
    qc.plots[[i+(i-1)]] <- densityplot(~ ., 
                                       before.set[, warp.parameters[[i]]$channel], 
                                       filter = curv1Filter(warp.parameters[[i]]$channel, bwFac = warp.parameters[[i]]$bwFac),
                                       main = paste("Non-warped", names(warp.parameters)[i], sep = " "))
    qc.plots[[i+i]]     <- densityplot(~ ., 
                                       after.set[, warp.parameters[[i]]$channel], 
                                       filter = curv1Filter(warp.parameters[[i]]$channel, bwFac = warp.parameters[[i]]$bwFac),
                                       main = paste("Warped", names(warp.parameters)[i], sep = " "))
  }
  qc.plots
  
  doc.name <- paste("Warp.QC.densityplots", Sys.Date(), paste(names(warp.parameters), collapse = "."), "CD3.fix", "pdf", sep = ".")
  pdf(doc.name)
  print(qc.plots)
  dev.off()
}
warp.qc.plots.edit(input.set, input.set.warped)

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
outDir = file.path("./data_modified/TPHE_warping", "CD3.fix", "CD4.CD8A.CD3.LIVEDEAD_warped")
dir.create(dirname(outDir))
dir.create(outDir)
fsApply(input.set.warped.inverse.decomp.mod[-grep(target, sampleNames(input.set.warped.inverse.decomp.mod))], function(i) write.FCS(x = i, filename = file.path(outDir, i@description$`$FIL`)))
