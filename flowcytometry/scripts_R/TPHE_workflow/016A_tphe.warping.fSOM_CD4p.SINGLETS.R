#### required objects/arguments
##
warp.parameters <- readRDS("./results_R/FCS.transforms/TPHE/TPHE.CD4.warp.parameters.rds") # channels-of-interest
warp.markers <- c("PERFORIN", "FOXP3", "CD185", "CD31", "GRZB")
warp.parameters <- warp.parameters[names(warp.parameters) %in% warp.markers] # warping parameters per channels-of-interest

trans.obj <- readRDS("./results_R/FCS.transforms/TPHE/asinh.transforms.TPHE.rds") # arcsinh transfomation parameters
trans.obj <- trans.obj[grep(paste(unlist(lapply(seq(warp.parameters), function(i) warp.parameters[[i]]$channel)), collapse = "|"), names(trans.obj))]
## file paths
fsom.files <- list.files("./data_modified/TPHE_warping/fSOM_CD4p.SINGLETS_PASS2/", 
                         recursive = FALSE, full.names = TRUE, pattern = ".fcs")
batch.names <- paste0("NGL0", seq(47, 91, 2)) ; batch.names <- batch.names[-grep("NGL063", batch.names)]

## warping target
target.file <- fsom.files[grep("HD0189_TPHE_NGL091", fsom.files)]
target <- "HD0189_TPHE_NGL091"; fsom.files <- fsom.files[-grep(target, fsom.files)]

#############################################################################################################################################################################################
#############################################################################################################################################################################################

## read in batched set
batch.name <- "NGL059"

input.set <- fset.compensate(read.flowSet(c(fsom.files[grep(batch.name, fsom.files)], target.file), transformation = FALSE)) # read in flowSet; compensate
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
outDir = file.path("./data_modified/TPHE_warping", batch.name, "PERFORIN.FOXP3.CD185.CD31.GRZB_warped_CD4p")

fsApply(input.set.warped.inverse.decomp.mod[-grep(target, sampleNames(input.set.warped.inverse.decomp.mod))], function(i) write.FCS(x = i, filename = file.path(outDir, i@description$`$FIL`)))

#############################################################################################################################################################################################
#############################################################################################################################################################################################
## good luck

for(i in seq(batch.names)){
  
  batch.name <- batch.names[i]
  print(batch.name)
  
  input.set <- fset.compensate(fset.trim(read.flowSet(c(fsom.files[grep(batch.name, fsom.files)], target.file), transformation = FALSE))) # read in flowSet; compensate
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
  outDir = file.path("./data_modified/TPHE_warping", batch.name, "PERFORIN.FOXP3.CD185.CD31.GRZB_warped_CD4p")
  print("Writing warped FCS files")
  fsApply(input.set.warped[-grep(target, sampleNames(input.set.warped))], function(i) write.FCS(x = i, filename = file.path(outDir, i@description$`$FIL`)))
}

#############################################################################################################################################################################################
#############################################################################################################################################################################################
## test bed
channels <- colnames(input.set)
markers <- get.markers(input.set[[1]])

## read in batched set
batch.name <- "NGL091"

input.set <- fset.compensate(read.flowSet(c(fsom.files[grep(batch.name, fsom.files)], target.file), transformation = FALSE)) # read in flowSet; compensate
input.set <- transform.set(input.set) # transform using 'trans.obj'

input.set <- input.set[grep("530053_1_|530048_19_|HD0189_TPHE_NGL091", sampleNames(input.set))] #NGL059
input.set <- input.set[grep("540034_7_|HD0189_TPHE_NGL091", sampleNames(input.set))] #NGL069
input.set <- input.set[grep("540071_19_|HD0189_TPHE_NGL091", sampleNames(input.set))] #NGL079
input.set <- input.set[grep("540096_19_|HD0189_TPHE_NGL091", sampleNames(input.set))] #NGL089
input.set <- input.set[grep("540102_1_|HD0189_TPHE_NGL091", sampleNames(input.set))] #NGL091

marker.name <- "GRZB"
bwFac.set <- 9
densityplot(~ ., fset.downsample(input.set, 5000)[, channels[grep(marker.name, markers)]], filter = curv1Filter(channels[grep(marker.name, markers)], bwFac = bwFac.set))
test.set <- warpSet(input.set, c(channels[grep(marker.name, markers)]), monwrd = FALSE, peakNr = 1, clipRange = 0, bwFac = bwFac.set, 
                    target = sampleNames(input.set)[grep(target, sampleNames(input.set))])
densityplot(~ ., fset.downsample(test.set, 5000)[, channels[grep(marker.name, markers)]], filter = curv1Filter(channels[grep(marker.name, markers)], bwFac = bwFac.set))
TPHE.CD4.warp.parameters[[marker.name]]["monwrd"] <- FALSE
TPHE.CD4.warp.parameters[[marker.name]]["peakNr"] <- 2
TPHE.CD4.warp.parameters[[marker.name]]["clipRange"] <- 0.01
TPHE.CD4.warp.parameters[[marker.name]]["bwFac"] <- bwFac.set