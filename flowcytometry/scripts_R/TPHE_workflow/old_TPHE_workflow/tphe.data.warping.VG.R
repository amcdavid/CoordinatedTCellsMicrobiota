## source flow-related libraries and functions
source("./source_R/flow.src.R", chdir = TRUE)

## file paths
tphe.warped.files <- list.files("./data_modified/TPHE_warping/targets_HD/scatter.VG.VE.GD.VA.warped/", recursive = TRUE, full.names = TRUE, pattern = ".fcs")

experiments <- folder.list(tphe.warped.files)

## warping targets
targets <- experiments$`./data_modified/TPHE_warping//targets_HD/scatter.VG.warped`
  
## batch with respective target
experiment.name <- "NGL047" # set this argumanet

input.set <- fset.compensate(fset.trim(read.flowSet(c(experiments[grep(experiment.name, names(experiments))][[1]], targets[grep(experiment.name, targets)]), 
                                    transformation = FALSE)))

trans.obj <- readRDS("./results_R/FCS.transforms/TPHE/asinh.transforms.TPHE.rds")

input.set.transformed <- transform(input.set, trans.obj$Violet.G$asinh.translist)
input.set.transformed.rangefix <- fset.rangefix(input.set.transformed, 'Violet G 550_40-A')

sub.set <- fset.downsample(input.set.transformed.rangefix, 10000)

library(flowViz)
bwfac.vg <- 1

densityplot(~ `Violet G 550_40-A`, sub.set, filter = curv1Filter("Violet G 550_40-A", bwFac = bwfac.vg))
sub.set.warped <- warpSet(sub.set, c('Violet G 550_40-A'), monwrd = FALSE, peakNr = 1, clipRange = 0.1, bwFac = bwfac.vg, target = sampleNames(sub.set)[grep("HD", sampleNames(sub.set))])
densityplot(~ `Violet G 550_40-A`, sub.set.warped, filter = curv1Filter("Violet G 550_40-A", bwFac = bwfac.vg))

input.set.transformed.rangefix.warped <- warpSet(input.set.transformed.rangefix, c('Violet G 550_40-A'), monwrd = FALSE, peakNr = 1, clipRange = 0.1, bwFac = bwfac.vg, 
                                                 target = sampleNames(input.set.transformed.rangefix)[grep("HD", sampleNames(input.set.transformed.rangefix))])
densityplot(~ `Violet G 550_40-A`, input.set.transformed.rangefix.warped, filter = curv1Filter("Violet G 550_40-A", bwFac = bwfac.vg))

#input.set.transformed.rangefix.warped <- input.set.transformed.rangefix #no warping needed

## inverse transform warped set
asinh.vals <- trans.obj$Violet.G$asinh.vals
sinhTrans <- sinhTransform(transformationId = "inverse.arcsinh-transformation", a = asinh.vals$a, b = asinh.vals$b, c = asinh.vals$c)
sinh.translist <- transformList('Violet G 550_40-A', sinhTrans)

input.set.transformed.rangefix.warped.inverse <- transform(input.set.transformed.rangefix.warped, sinh.translist)

## decompensate
input.set.transformed.rangefix.warped.inverse.decomp <- fset.decompensate(input.set.transformed.rangefix.warped.inverse)

## modify transformed files to include new identifiers
## change/add keywords
input.set.transformed.rangefix.warped.inverse.decomp.mod <- fsApply(input.set.transformed.rangefix.warped.inverse.decomp, function (i) {
  i@description$`EXPERIMENT NAME` <- unlist(strsplit(i@description$`EXPERIMENT NAME`, " "))[[1]]
  i@description$`WARPED CHANNELS` <- "FSC.AHW_SSC.AH_VG"
  i@description$`WARP TARGET` <- sampleNames(input.set)[grep("HD", sampleNames(input.set))]
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

fsApply(input.set.transformed.rangefix.warped.inverse.decomp.mod[-grep("HD", sampleNames(input.set.transformed.rangefix.warped.inverse.decomp.mod))], 
        function(i) write.FCS(x = i, filename = file.path(outDir, i@description$`$FIL`)))
