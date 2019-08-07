## source flow-related libraries and functions
source("./source_R/flow.src.R", chdir = TRUE)

## script functions
##
batch.list.edit <- function(fcs.file.paths){
  experiments <- unique(lapply(seq(length(fcs.file.paths)), FUN = function(x) {
    (strsplit(fcs.file.paths, "/")[[x]][grep("NGL*", strsplit(fcs.file.paths, "/")[[x]])][1])}))
  batches <- lapply(seq(unique(experiments)), FUN = function(x) {fcs.file.paths[grep(experiments[x], fcs.file.paths)]})
  names(batches) <- unique(experiments)
  batches
}

transform.set <- function(input.set){
  
  for(j in seq(length(trans.obj))){
    print(paste0("transforming ", names(trans.obj)[j]))
    input.set <- transform(input.set, trans.obj[[j]]$asinh.translist)
    input.set <- fsApply(input.set, function(i){
      i@parameters@data$minRange[grep(names(trans.obj)[j], i@parameters@data$name)] <- min(exprs(i)[, names(trans.obj)[j]])
      i@parameters@data$maxRange[grep(names(trans.obj)[j], i@parameters@data$name)] <- max(exprs(i)[, names(trans.obj)[j]])
      i
    })
  }
  input.set
} # needs 'trans.obj'
  
warp.set <- function(input.set){
  
  for(i in seq(CD4.warp)){
    input.set <- warpSet(input.set, 
                       stains = CD4.warp[[i]][["channel"]], 
                       monwrd = CD4.warp[[i]][["monwrd"]], 
                       peakNr = CD4.warp[[i]][["peakNr"]], 
                       clipRange = CD4.warp[[i]][["clipRange"]], 
                       bwFac = CD4.warp[[i]][["bwFac"]], 
                       target = sampleNames(input.set)[grep("HD", sampleNames(input.set))])
  }
  input.set
} # needs 'CD4.warp' and HD target in flowSet

inverse.transform.set <- function(input.set){
  
  for(j in seq(length(trans.obj))){
    sinhTrans <- sinhTransform(transformationId = "inverse.arcsinh-transformation",
                               a = trans.obj[[j]]$asinh.vals$a,
                               b = trans.obj[[j]]$asinh.vals$b, 
                               c = trans.obj[[j]]$asinh.vals$c)
    sinh.translist <- transformList(names(trans.obj)[j], sinhTrans)
    
    input.set <- transform(input.set, sinh.translist)
  }
  input.set
} # needs 'trans.obj

warp.qc.plots <- function(input.set, input.set.warped){
  
  before.set <- input.set
  after.set <- input.set.warped
  
  qc.plots <- vector(mode = "list", length = length(CD4.warp)*2)
  
  pData(before.set)[[1]] <- fsApply(before.set, function(i) i@description$`$FIL` <- i@description$`$SRC`)
  pData(after.set)[[1]] <- fsApply(after.set, function(i) i@description$`$FIL` <- i@description$`$SRC`)
  
  before.set <- fset.downsample(before.set, 5000)
  after.set <- fset.downsample(after.set, 5000)
  
  for(i in seq(CD4.warp)){
    qc.plots[[i+(i-1)]] <- densityplot(~ ., 
                                       before.set[, CD4.warp[[i]]$channel], 
                                       filter = curv1Filter(CD4.warp[[i]]$channel, bwFac = CD4.warp[[i]]$bwFac),
                                       main = "Non-warped")
    qc.plots[[i+i]]     <- densityplot(~ ., 
                                       after.set[, CD4.warp[[i]]$channel], 
                                       filter = curv1Filter(CD4.warp[[i]]$channel, bwFac = CD4.warp[[i]]$bwFac),
                                       main = "Warped")
  }
  qc.plots
  
  doc.name <- paste("Warp.QC.densityplots", batch.name, "pdf", sep = ".")
  pdf(doc.name)
  print(qc.plots)
  dev.off()
} # needs 'CD.warp'
##
##

## required objects/arguments
trans.obj <- readRDS("./results_R/FCS.transforms/TPHE/asinh.transforms.TPHE.rds") # arcsinh transfomation parameters
warp.markers <- c("PERFORIN", "CD57", "CD28", "FOXP3", "CD197", "CD185", "CD45RO", "CD127", "CD31", "GRZB") # channels to be warped
CD4.warp <- readRDS("TPHE.CD4.warp.parameters.rds") ; CD4.warp <- CD4.warp[names(CD4.warp) %in% warp.markers] # warping parameters per-channel

## file paths
tphe.warped.files <- list.files("./data_modified/TPHE_warping", recursive = TRUE, full.names = TRUE, pattern = "fda_WARPED_FSC.AHW_SSC.AHW_VG_VE_VA_GD_fSOM_CD4pos_fSOM_CD4pos_scatter_.fcs")
batches <- batch.list.edit(tphe.warped.files)

## warping target
target <- list.files("./data_modified/TPHE_warping/targets_HD/scatter.VG.VE.GD.VA.warped/fSOM_CD3.CD4.CD8.LIVEDEAD/fSOM_CD4pos/fSOM_CD4pos_scatter/fSOM_CD4pos_scatter/", 
                      full.names = TRUE) ; target <- target[grep("HD0189_TPHE_NGL091", target)]

#############################################################################################################################################################################################
#############################################################################################################################################################################################

## read in batched set
batch.name <- "NGL049"

input.set <- fset.compensate(fset.trim(read.flowSet(c(batches[batch.name][[1]], target), transformation = FALSE))) # read in flowSet; compensate
input.set <- transform.set(input.set) # transform using 'trans.obj'
input.set.warped <- warp.set(input.set) # warp using 'CD4.warp'

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
  i@description$`WARPED CHANNELS` <- "FSC.AHW_SSC.AHW_BA_BB_GA_GC_GD_GE_RA_RB_RC_VA_VB_VC_VD_VE_VG_VH"
  i@description$`WARP TARGET` <- sampleNames(input.set.warped.inverse.decomp)[grep("HD", sampleNames(input.set.warped.inverse.decomp))]
  i@description$`WARPED` <- "fda_WARPED"
  i@description$`CLUSTERING` <- "fSOM"
  i@description$`METACLUSTER` <- "CD4pos"
  i@description$`$FIL` <- paste(i@description$`$SRC`, 
                                i@description$`TUBE NAME`, 
                                i@description$`EXPERIMENT NAME`, 
                                i@description$WARPED,
                                i@description$`CLUSTERING`,
                                i@description$`METACLUSTER`,
                                ".fcs",
                                sep = "_")
  i
})

## write modified files
outDir = file.path("./data_modified/TPHE_warping", batch.name, "fSOM_CD4pos_warped")

fsApply(input.set.warped.inverse.decomp.mod[-grep("HD", sampleNames(input.set.warped.inverse.decomp.mod))], function(i) write.FCS(x = i, filename = file.path(outDir, i@description$`$FIL`)))

#############################################################################################################################################################################################
#############################################################################################################################################################################################
## good luck

for(i in seq(names(batches))){
  
  batch.name <- names(batches)[i]
  print(batch.name)
  
  input.set <- fset.compensate(fset.trim(read.flowSet(c(batches[batch.name][[1]], target), transformation = FALSE))) # read in flowSet; compensate
  input.set <- transform.set(input.set) # transform using 'trans.obj'
  input.set.warped <- warp.set(input.set) # warp using 'CD4.warp'
  
  qc.plots <- warp.qc.plots(input.set, input.set.warped)
  
  ## inverse transform warped set
  input.set.warped.inverse <- inverse.transform.set(input.set.warped)
  
  ## decomp
  input.set.warped.inverse.decomp <- fset.decompensate(input.set.warped.inverse)
  
  ## modify transformed files to include new identifiers
  ## change/add keywords
  input.set.warped.inverse.decomp.mod <- fsApply(input.set.warped.inverse.decomp, function (i) {
    i@description$`EXPERIMENT NAME` <- unlist(strsplit(i@description$`EXPERIMENT NAME`, " "))[[1]]
    i@description$`WARPED CHANNELS` <- "FSC.AHW_SSC.AHW_BA_BB_GA_GC_GD_GE_RA_RB_RC_VA_VB_VC_VD_VE_VG_VH"
    i@description$`WARP TARGET` <- sampleNames(input.set.warped.inverse.decomp)[grep("HD", sampleNames(input.set.warped.inverse.decomp))]
    i@description$`WARPED` <- "fda_WARPED"
    i@description$`CLUSTERING` <- "fSOM"
    i@description$`METACLUSTER` <- "CD4pos"
    i@description$`$FIL` <- paste(i@description$`$SRC`, 
                                  i@description$`TUBE NAME`, 
                                  i@description$`EXPERIMENT NAME`, 
                                  i@description$WARPED,
                                  i@description$`CLUSTERING`,
                                  i@description$`METACLUSTER`,
                                  ".fcs",
                                  sep = "_")
    i
  })
  
  ## write modified files
  outDir = file.path("./data_modified/TPHE_warping", batch.name, "fSOM_CD4pos_warped")
  
  fsApply(input.set.warped.inverse.decomp.mod[-grep("HD", sampleNames(input.set.warped.inverse.decomp.mod))], function(i) write.FCS(x = i, filename = file.path(outDir, i@description$`$FIL`)))
  
}

