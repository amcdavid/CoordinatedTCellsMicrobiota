## source flow-related libraries and functions
source("./source_R/flow.src.R", chdir = TRUE)

## file paths
tphe.original.files <- list.files("./data_source/TPHE2_with.spillover", recursive = TRUE, full.names = TRUE, pattern = "_TPHE.fcs")

batches <- batch.list(tphe.original.files)

## warping targets
targets <- list.files("./data_modified/TPHE_warping/targets_HD/scatter.warped", full.names = TRUE)

## batch with respective target
experiment.name <- "NGL069" # set this argumanet

# input.set <- fset.trim(read.flowSet(c(batches[experiment.name][[1]][-grep("530019_1_", batches[experiment.name][[1]])], targets[grep(experiment.name, targets)]), 
#                                    transformation = FALSE)) # problem sample removed; batch specific for NGL047

input.set <- fset.trim(read.flowSet(c(batches[experiment.name][[1]][-grep("540038_1_|540034_1", batches[experiment.name][[1]])], targets[grep(experiment.name, targets)]), 
                                    transformation = FALSE)) # problem samplea removed; batch specific for NGL069

input.set <- fset.trim(read.flowSet(c(batches[experiment.name][[1]], targets[grep(experiment.name, targets)]), 
                                    transformation = FALSE))

## downsample for quicker testing/plotting
sub.set <- fset.downsample(input.set, 10000)

## non-warped scatter channels
bwfac.fscA <- 1
bwfac.fscH <- 1
bwfac.fscW <- 2
bwfac.sscA <- 1.5
bwfac.sscH <- 1.5
bwfac.sscW <- 3.5


densityplot(~ `FSC-A`, sub.set, filter = curv1Filter("FSC-A", bwFac = bwfac.fscA))
densityplot(~ `FSC-H`, sub.set, filter = curv1Filter("FSC-H", bwFac = bwfac.fscH))
densityplot(~ `FSC-W`, sub.set, filter = curv1Filter("FSC-W", bwFac = bwfac.fscW), xlim = c(50000, 110000))
densityplot(~ `SSC-A`, sub.set, filter = curv1Filter("SSC-A", bwFac = bwfac.sscA))
densityplot(~ `SSC-H`, sub.set, filter = curv1Filter("SSC-H", bwFac = bwfac.sscH))
densityplot(~ `SSC-W`, sub.set, filter = curv1Filter("SSC-W", bwFac = bwfac.sscW), xlim = c(50000, 110000))

## warped scatter channels
sub.set.warped <- warpSet(sub.set, c("FSC-A"), monwrd = FALSE, peakNr = 2, clipRange = 0.01, bwFac = bwfac.fscA, target = sampleNames(sub.set)[grep("HD", sampleNames(sub.set))])
densityplot(~ `FSC-A`, sub.set.warped, filter = curv1Filter("FSC-A", bwFac = bwfac.fscA))

sub.set.warped <- warpSet(sub.set, c("FSC-H"), monwrd = FALSE, peakNr = 2, clipRange = 0.01, bwFac = bwfac.fscH, target = sampleNames(sub.set)[grep("HD", sampleNames(sub.set))])
densityplot(~ `FSC-H`, sub.set.warped, filter = curv1Filter("FSC-H", bwFac = bwfac.fscH))

sub.set.warped <- warpSet(sub.set, c("FSC-W"), monwrd = FALSE, peakNr = 1, clipRange = 0.1, bwFac = bwfac.fscW, target = sampleNames(sub.set)[grep("HD", sampleNames(sub.set))])
densityplot(~ `FSC-W`, sub.set.warped, filter = curv1Filter("FSC-W", bwFac = bwfac.fscW), xlim = c(50000, 110000))

sub.set.warped <- warpSet(sub.set, c("SSC-A"), monwrd = FALSE, peakNr = 2, clipRange = 0.1, bwFac = bwfac.sscA, target = sampleNames(sub.set)[grep("HD", sampleNames(sub.set))])
densityplot(~ `SSC-A`, sub.set.warped, filter = curv1Filter("SSC-A", bwFac = bwfac.sscA))

sub.set.warped <- warpSet(sub.set, c("SSC-H"), monwrd = FALSE, peakNr = 3, clipRange = 0.1, bwFac = bwfac.sscH, target = sampleNames(sub.set)[grep("HD", sampleNames(sub.set))])
densityplot(~ `SSC-H`, sub.set.warped, filter = curv1Filter("SSC-H", bwFac = bwfac.sscH))

sub.set.warped <- warpSet(sub.set, c("SSC-W"), monwrd = FALSE, peakNr = 1, clipRange = 0.1, bwFac = bwfac.sscW, target = sampleNames(sub.set)[grep("HD", sampleNames(sub.set))])
densityplot(~ `SSC-W`, sub.set.warped, filter = curv1Filter("SSC-W", bwFac = bwfac.sscW), xlim = c(50000, 110000))

## warp actual data
input.set.warped <- input.set

input.set.warped <- warpSet(input.set.warped, c("FSC-A"), monwrd = FALSE, peakNr = 2, clipRange = 0.1, bwFac = bwfac.fscA, target = sampleNames(input.set.warped)[grep("HD", sampleNames(input.set.warped))])
input.set.warped <- warpSet(input.set.warped, c("FSC-H"), monwrd = FALSE, peakNr = 2, clipRange = 0.1, bwFac = bwfac.fscH, target = sampleNames(input.set.warped)[grep("HD", sampleNames(input.set.warped))])
input.set.warped <- warpSet(input.set.warped, c("FSC-W"), monwrd = FALSE, peakNr = 1, clipRange = 0.1, bwFac = bwfac.fscW, target = sampleNames(input.set.warped)[grep("HD", sampleNames(input.set.warped))])
input.set.warped <- warpSet(input.set.warped, c("SSC-A"), monwrd = FALSE, peakNr = 2, clipRange = 0.1, bwFac = bwfac.sscA, target = sampleNames(input.set.warped)[grep("HD", sampleNames(input.set.warped))])
input.set.warped <- warpSet(input.set.warped, c("SSC-H"), monwrd = FALSE, peakNr = 1, clipRange = 0.2, bwFac = bwfac.sscH, target = sampleNames(input.set.warped)[grep("HD", sampleNames(input.set.warped))])
input.set.warped <- warpSet(input.set.warped, c("SSC-W"), monwrd = FALSE, peakNr = 1, clipRange = 0.2, bwFac = bwfac.sscW, target = sampleNames(input.set.warped)[grep("HD", sampleNames(input.set.warped))])

## warp QC
sub.set <- fset.downsample(input.set.warped, 5000)

densityplot(~ `FSC-A`, sub.set, filter = curv1Filter("FSC-A", bwFac = bwfac.fscA))
densityplot(~ `FSC-H`, sub.set, filter = curv1Filter("FSC-H", bwFac = bwfac.fscH))
densityplot(~ `FSC-W`, sub.set, filter = curv1Filter("FSC-W", bwFac = bwfac.fscW), xlim = c(50000, 110000))
densityplot(~ `SSC-A`, sub.set, filter = curv1Filter("SSC-A", bwFac = bwfac.sscA))
densityplot(~ `SSC-H`, sub.set, filter = curv1Filter("SSC-H", bwFac = bwfac.sscH))
densityplot(~ `SSC-W`, sub.set, filter = curv1Filter("SSC-W", bwFac = bwfac.sscW), xlim = c(50000, 110000))

## modify transformed files to include new identifiers
## change/add keywords
input.set.warped.mod <- fsApply(input.set.warped, function (i) {
  i@description$`EXPERIMENT NAME` <- unlist(strsplit(i@description$`EXPERIMENT NAME`, " "))[[1]]
  i@description$`WARPED CHANNELS` <- "FSC.AHW_SSC.AHW"
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
outDir = paste0("./data_modified/TPHE_warping/", experiment.name)
dir.create(outDir)

fsApply(input.set.warped.mod[-grep("HD", sampleNames(input.set.warped.mod))], function(i) write.FCS(x = i, filename = file.path(outDir, i@description$`$FIL`)))
