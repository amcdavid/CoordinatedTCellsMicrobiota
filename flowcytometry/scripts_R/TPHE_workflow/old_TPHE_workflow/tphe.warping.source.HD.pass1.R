## source flow-related libraries and functions
source("./source_R/flow.src.R", chdir = TRUE)

## required objects/arguments
trans.obj <- readRDS("./results_R/FCS.transforms/TPHE/asinh.transforms.TPHE.rds") # arcsinh transfomation parameters
trans.obj <- trans.obj[grep('Violet G|Violet E|Violet A|Green D', names(trans.obj))]
warp.markers <- c("LIVEDEAD", "CD3", "CD4", "CD8A") # channels to be warped
warp.parameters <- readRDS("./results_R/FCS.transforms/TPHE/TPHE.warp.parameters.rds") ; warp.parameters <- warp.parameters[names(warp.parameters) %in% warp.markers] # warping parameters per-channel

## file paths
tphe.HD.files <- list.files("./data_source/TPHE2_with.spillover/", recursive = TRUE, full.names = TRUE, pattern = "HD")
target <- tphe.HD.files[grep("HD0189_TPHE_NGL091", tphe.HD.files)]
tphe.HD.files <- tphe.HD.files[-grep("HD0189_TPHE_NGL091|NGL063", tphe.HD.files)]


#############################################################################################################################################################################################
#############################################################################################################################################################################################

## read in full set with target

input.set <- fset.compensate(fset.trim(read.flowSet(c(tphe.HD.files, target), transformation = FALSE))) # read in flowSet; compensate
input.set <- transform.set(input.set) # transform using 'trans.obj'
input.set.warped <- warp.set(input.set, "HD0189_TPHE_NGL091") # warp using 'warp.parameters'

## warp QC plots
warp.qc.HD.plots <- function(input.set, input.set.warped){
  
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
  
  doc.name <- paste("Warp.QC.densityplots", Sys.Date(), paste(names(warp.parameters), collapse = "."), "HD", "pdf", sep = ".")
  pdf(doc.name)
  print(qc.plots)
  dev.off()
} # needs 'warp.parameters'
warp.qc.HD.plots(input.set, input.set.warped)

## inverse transform warped set
input.set.warped.inverse <- inverse.transform.set(input.set.warped)

## decomp
input.set.warped.inverse.decomp <- fset.decompensate(input.set.warped.inverse)

## modify transformed files to include new identifiers
## change/add keywords
input.set.warped.inverse.decomp.mod <- fsApply(input.set.warped.inverse.decomp, function (i) {
  i@description$`EXPERIMENT NAME` <- unlist(strsplit(i@description$`EXPERIMENT NAME`, " "))[[1]]
  i@description$`WARPED CHANNELS` <- paste(names(warp.parameters), collapse = ".")
  i@description$`WARP TARGET` <- sampleNames(input.set.warped.inverse.decomp)[grep("HD0189_TPHE_NGL091", sampleNames(input.set.warped.inverse.decomp))]
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
outDir = file.path("./data_modified/TPHE_warping/targets_HD/", "CD4.CD8A.CD3.LIVEDEAD_warped")

fsApply(input.set.warped.inverse.decomp.mod[-grep("HD0189_TPHE_NGL091", sampleNames(input.set.warped.inverse.decomp.mod))], function(i) write.FCS(x = i, filename = file.path(outDir, i@description$`$FIL`)))
