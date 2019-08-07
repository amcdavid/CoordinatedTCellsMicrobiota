## source flow-related libraries and functions
source("./source_R/flow.src.R", chdir = TRUE)

##warp Violet.G (LIVEDEAD); use scatter-warped HD samples w/ target
tphe.warped.hd <- list.files("./data_modified/TPHE_warping/", full.names = TRUE, pattern = "HD")
target <- "HD0189_TPHE_NGL091_fda_WARPED_FSC.AHW_SSC.AHW_.fcs"

input.set <- fset.compensate(fset.trim(read.flowSet(tphe.warped.hd)))

trans.obj <- readRDS("./results_R/FCS.transforms/TPHE/asinh.transforms.TPHE.rds")

input.set.transformed <- transform(input.set, trans.obj$Violet.G)
input.set.transformed.rangefix <- fset.rangefix(input.set.transformed, 'Violet G 550_40-A')

sub.set <- fset.downsample(input.set.transformed.rangefix, 20000)

library(flowViz)
densityplot(~ `Violet G 550_40-A`, sub.set, filter = curv1Filter("Violet G 550_40-A", bwFac = 1.5))
sub.set.warped <- warpSet(sub.set, c('Violet G 550_40-A'), monwrd = FALSE, peakNr = 1, clipRange = 0.01, bwFac = 1.5, target = target)
densityplot(~ `Violet G 550_40-A`, sub.set.warped, filter = curv1Filter("Violet G 550_40-A", bwFac = 1.5))

input.set.transformed.rangefix.warped <- warpSet(input.set.transformed.rangefix, c('Violet G 550_40-A'), monwrd = FALSE, peakNr = 1, clipRange = 0.01, bwFac = 1.5, target = target)
densityplot(~ `Violet G 550_40-A`, input.set.transformed.rangefix.warped, filter = curv1Filter("Violet G 550_40-A", bwFac = 1.5))

## inverse transform warped set
asinh.vals <- list(a = -3, b = 1/400, c = 2.5)
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
  i@description$`WARP TARGET` <- "HD0189_TPHE_NGL091"
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
outDir = "./data_modified/TPHE_warping"

fsApply(input.set.transformed.rangefix.warped.inverse.decomp.mod, function(i) write.FCS(x = i, filename = file.path(outDir, i@description$`$FIL`)))
