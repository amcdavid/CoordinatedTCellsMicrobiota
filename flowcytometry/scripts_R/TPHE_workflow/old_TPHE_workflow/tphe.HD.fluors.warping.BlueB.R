## source flow-related libraries and functions
source("./source_R/flow.src.R", chdir = TRUE)

##warp arguments; use fSOM subsets

channel.name <- "Blue B" # set
marker.name <- "CD122" # set

## subset data includes only CD4+ (low granzyme B) so use original data to optimize transforms
tphe.hd <- list.files("./data_modified/TPHE_warping/targets_HD/scatter.VG.VE.GD.VA.warped/", recursive = FALSE, full.names = TRUE, pattern = ".fcs")

target <- "HD0189_TPHE_NGL091"

input.set <- fset.compensate(fset.trim(read.flowSet(tphe.hd[grep(target, tphe.hd)])))

markers <- get.markers(input.set[[1]])

dat <- as.data.frame(exprs(input.set[[grep(target, sampleNames(input.set))]]))
colnames(dat) <- markers

tmp <- subset(dat, select = marker.name)
ggplot(tmp, aes_string(x = names(tmp))) + geom_density()

a = 0
b = 1/150
c = 2

asinh.vals <- list(a = a, b = b, c = c)
ggplot(asinh(asinh.vals$a + asinh.vals$b * tmp) + asinh.vals$c, aes_string(x = names(tmp))) + geom_histogram(bins = 400)

dat <- as.data.frame(exprs(input.set[[grep(target, sampleNames(input.set))]]))
colnames(dat) <- markers
dat[[marker.name]] <- asinh(asinh.vals$a + asinh.vals$b * dat[[marker.name]]) + asinh.vals$c

fcs.explorer(dat, marker.name, "FCS.A", 20000)

asinh.transforms.TPHE <- readRDS("./results_R/FCS.transforms/TPHE/asinh.transforms.TPHE.rds")

channel <- colnames(input.set[[1]])[grep(channel.name, colnames(input.set[[1]]))]

length(asinh.transforms.TPHE)

a = asinh.vals$a
b = asinh.vals$b
c = asinh.vals$c

asinh.transforms.TPHE[length(asinh.transforms.TPHE) + 1] <- list(channel = (list(asinh.translist = transformList(channel, 
                                                                                                                 arcsinhTransform(transformationId = "arcsinh-transformation", 
                                                                                                                                  a = a, 
                                                                                                                                  b = b, 
                                                                                                                                  c = c)),
                                                                                 asinh.vals = asinh.vals)))
names(asinh.transforms.TPHE)[length(asinh.transforms.TPHE)] <- channel
length(asinh.transforms.TPHE)
names(asinh.transforms.TPHE)

saveRDS(asinh.transforms.TPHE, "./results_R/FCS.transforms/TPHE/asinh.transforms.TPHE.rds")


##
tphe.warped.hd <- list.files("./data_modified/TPHE_warping/targets_HD/scatter.VG.VE.GD.VA.warped/fSOM_CD3.CD4.CD8.LIVEDEAD/fSOM_CD4pos/fSOM_CD4pos_scatter/fSOM_CD4pos_scatter/", 
                             full.names = TRUE, pattern = ".fcs")
target <- "HD0189_TPHE_NGL091"

input.set <- fset.compensate(fset.trim(read.flowSet(tphe.warped.hd)))

markers <- get.markers(input.set[[1]])

dat <- as.data.frame(exprs(input.set[[grep(target, sampleNames(input.set))]]))
colnames(dat) <- markers

tmp <- subset(dat, select = marker.name)
ggplot(tmp, aes_string(x = names(tmp))) + geom_density()

a = 0
b = 1/150
c = 2

asinh.vals <- list(a = a, b = b, c = c)
ggplot(asinh(asinh.vals$a + asinh.vals$b * tmp) + asinh.vals$c, aes_string(x = names(tmp))) + geom_histogram(bins = 400)

dat <- as.data.frame(exprs(input.set[[grep(target, sampleNames(input.set))]]))
colnames(dat) <- markers
dat[[marker.name]] <- asinh(asinh.vals$a + asinh.vals$b * dat[[marker.name]]) + asinh.vals$c

fcs.explorer(dat, marker.name, "FCS.A", 20000)

trans.obj <- readRDS("./results_R/FCS.transforms/TPHE/asinh.transforms.TPHE.rds")

input.set.transformed <- transform(input.set, trans.obj[[channel]]$asinh.translist)
input.set.transformed.rangefix <- fset.rangefix(input.set.transformed, channel)

library(flowViz)

bwFac.set <- 2

densityplot(~ ., fset.downsample(input.set.transformed.rangefix, 20000)[, channel], filter = curv1Filter(channel, bwFac = bwFac.set))
sub.set.warped <- warpSet(fset.downsample(input.set.transformed.rangefix, 20000)[, channel], c(channel), monwrd = FALSE, peakNr = 1, clipRange = 0.1, bwFac = bwFac.set, 
                          target = sampleNames(input.set.transformed.rangefix[grep(target, sampleNames(input.set.transformed.rangefix))]))
densityplot(~ ., sub.set.warped[, channel], filter = curv1Filter(channel, bwFac = bwFac.set))

input.set.transformed.rangefix.warped <- warpSet(input.set.transformed.rangefix, c(channel), monwrd = FALSE, peakNr = 1, clipRange = 0.01, bwFac = bwFac.set, 
                                                 target = sampleNames(input.set.transformed.rangefix[grep(target, sampleNames(input.set.transformed.rangefix))]))
densityplot(~ ., fset.downsample(input.set.transformed.rangefix.warped, 10000)[, channel], filter = curv1Filter(channel, bwFac = bwFac.set))

## inverse transform warped set
asinh.vals <- trans.obj[[channel]]$asinh.vals

sinhTrans <- sinhTransform(transformationId = "inverse.arcsinh-transformation", a = asinh.vals$a, b = asinh.vals$b, c = asinh.vals$c)
sinh.translist <- transformList(channel, sinhTrans)

input.set.transformed.rangefix.warped.inverse <- transform(input.set.transformed.rangefix.warped, sinh.translist)

## decompensate
input.set.transformed.rangefix.warped.inverse.decomp <- fset.decompensate(input.set.transformed.rangefix.warped.inverse)

## modify transformed files to include new identifiers
## change/add keywords
input.set.transformed.rangefix.warped.inverse.decomp.mod <- fsApply(input.set.transformed.rangefix.warped.inverse.decomp, function (i) {
  i@description$`WARPED CHANNELS` <- "FSC.AHW_SSC.AH_BA_BB_VH_VG_VE_VD_VC_VB_VA_RC_RB_RA_GE_GC_GD_GA"
  i
})

## write modified files
outDir = "./data_modified/TPHE_warping/targets_HD/scatter.VG.VE.GD.VA.warped/fSOM_CD3.CD4.CD8.LIVEDEAD/fSOM_CD4pos/fSOM_CD4pos_scatter/fSOM_CD4pos_scatter/"


fsApply(input.set.transformed.rangefix.warped.inverse.decomp.mod, function(i) write.FCS(x = i, filename = file.path(outDir, i@description$`$FIL`)))
