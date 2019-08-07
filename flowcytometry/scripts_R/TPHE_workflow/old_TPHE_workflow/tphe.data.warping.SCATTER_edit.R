## source flow-related libraries and functions
source("./source_R/flow.src.R", chdir = TRUE)

#### script function
## warp scatter channels
scatter.norm <- function(fcs.file.paths, experiment.name){
  
  input.set <- fset.trim(read.flowSet(c(fcs.file.paths[grep(experiment.name, fcs.file.paths)], targets[grep(experiment.name, targets)]),
                                      transformation = FALSE))
  
  input.set <- warpSet(input.set, c("FSC-A"), monwrd = FALSE, peakNr = 1, clipRange = 0.2, bwFac = bwfac.fscA, 
                       target = sampleNames(input.set)[grep("HD", sampleNames(input.set))])
  input.set <- warpSet(input.set, c("FSC-H"), monwrd = FALSE, peakNr = 1, clipRange = 0.1, bwFac = bwfac.fscH, 
                       target = sampleNames(input.set)[grep("HD", sampleNames(input.set))])
  input.set <- warpSet(input.set, c("FSC-W"), monwrd = FALSE, peakNr = 1, clipRange = 0.1, bwFac = bwfac.fscW, 
                       target = sampleNames(input.set)[grep("HD", sampleNames(input.set))])
  input.set <- warpSet(input.set, c("SSC-A"), monwrd = FALSE, peakNr = 1, clipRange = 0.1, bwFac = bwfac.sscA, 
                       target = sampleNames(input.set)[grep("HD", sampleNames(input.set))])
  input.set <- warpSet(input.set, c("SSC-H"), monwrd = FALSE, peakNr = 1, clipRange = 0.1, bwFac = bwfac.sscH, 
                       target = sampleNames(input.set)[grep("HD", sampleNames(input.set))])
  input.set <- warpSet(input.set, c("SSC-W"), monwrd = FALSE, peakNr = 1, clipRange = 0.1, bwFac = bwfac.sscW, 
                       target = sampleNames(input.set)[grep("HD", sampleNames(input.set))])
  input.set
}
####

## file paths
tphe.warped.files <- list.files("./data_modified/TPHE_warping", recursive = TRUE, full.names = TRUE, pattern = "fSOM_CD4pos_.fcs")

## warping targets
targets <- list.files("./data_modified/TPHE_warping/targets_HD/scatter.VG.VE.GD.VA.warped/fSOM_CD3.CD4.CD8.LIVEDEAD/fSOM_CD4pos", full.names = TRUE)

tphe.warped.files <- tphe.warped.files[-grep(paste(basename(targets), collapse = "|"), basename(tphe.warped.files))]

targets <- targets[-grep("HD0191", targets)]
## scatter channels bandwitch factors for fda landmark
bwfac.fscA <- 2
bwfac.fscH <- 2
bwfac.fscW <- 2
bwfac.sscA <- 2
bwfac.sscH <- 2
bwfac.sscW <- 2

## batch with respective target
experiment.name <- "NGL091" # set this argumanet

##
input.set.warped <- scatter.norm(tphe.warped.files, experiment.name)

## warp QC
densityplot(~ `FSC-A`, fset.downsample(input.set.warped, 1000), filter = curv1Filter("FSC-A", bwFac = bwfac.fscA))
densityplot(~ `FSC-H`, fset.downsample(input.set.warped, 1000), filter = curv1Filter("FSC-H", bwFac = bwfac.fscH))
densityplot(~ `FSC-W`, fset.downsample(input.set.warped, 1000), filter = curv1Filter("FSC-W", bwFac = bwfac.fscW), xlim = c(50000, 110000))
densityplot(~ `SSC-A`, fset.downsample(input.set.warped, 1000), filter = curv1Filter("SSC-A", bwFac = bwfac.sscA))
densityplot(~ `SSC-H`, fset.downsample(input.set.warped, 1000), filter = curv1Filter("SSC-H", bwFac = bwfac.sscH))
densityplot(~ `SSC-W`, fset.downsample(input.set.warped, 1000), filter = curv1Filter("SSC-W", bwFac = bwfac.sscW), xlim = c(50000, 110000))

## modify transformed files to include new identifiers
## change/add keywords
input.set.warped.mod <- fsApply(input.set.warped, function (i) {
  i@description$`EXPERIMENT NAME` <- unlist(strsplit(i@description$`EXPERIMENT NAME`, " "))[[1]]
  i@description$`WARPED CHANNELS` <- "FSC.AHW_SSC.AHW_VG_VE_VA_GD"
  i@description$`WARP TARGET` <- sampleNames(input.set.warped)[grep("HD", sampleNames(input.set.warped))]
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
outDir = file.path("./data_modified/TPHE_warping", experiment.name, "VG_VE_VA_GD/fSOM_CD3.CD4.CD8.LIVEDEAD/fSOM_CD4pos/FSC.AHW_SSC.AHW_VG_VE_VA_GD/")

fsApply(input.set.warped.mod[-grep("HD", sampleNames(input.set.warped.mod))], function(i) write.FCS(x = i, filename = file.path(outDir, i@description$`$FIL`)))

##QC

scatter.warped.files <- list.files("./data_modified/TPHE_warping/", recursive = TRUE, full.names = TRUE, pattern = "FSC.AHW_SSC.AHW_VG_VE_VA_GD")

med.qc <- median.QC.scatter(scatter.warped.files)

ggplot(med.qc, aes(x = random, y = FSC.A)) + geom_point()
ggplot(med.qc, aes(x = random, y = FSC.W)) + geom_point()
ggplot(med.qc, aes(x = random, y = FSC.H)) + geom_point()

ggplot(med.qc, aes(x = random, y = SSC.A)) + geom_point()
ggplot(med.qc, aes(x = random, y = SSC.W)) + geom_point()
ggplot(med.qc, aes(x = random, y = SSC.H)) + geom_point()
