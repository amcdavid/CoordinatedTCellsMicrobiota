## source flow-related libraries and functions
source("./source_R/flow.src.R", chdir = TRUE)

##Violet.G (LIVEDEAD) warping parameters
tphe.hd <- list.files("./data_source/TPHE2_with.spillover/", full.names = TRUE, recursive = TRUE, pattern = "HD")
target <- "HD0189_TPHE_NGL091"

input.set <- fset.compensate(fset.trim(read.flowSet(tphe.hd)))

asinh.vals <- list(a = -3, b = 1/400, c = 2.5)

asinh.transforms.TPHE <- list(Violet.G = (list(asinh.translist = transformList('Violet G 550_40-A', 
                                                                               arcsinhTransform(transformationId = "arcsinh-transformation", 
                                                                                                a = -3, 
                                                                                                b = 1/400, 
                                                                                                c = 2.5)),
                                               asinh.vals = asinh.vals)))

saveRDS(asinh.transforms.TPHE, "./results_R/FCS.transforms/TPHE/asinh.transforms.TPHE.rds")

trans.obj <- readRDS("./results_R/FCS.transforms/TPHE/asinh.transforms.TPHE.rds")

input.set.transformed <- transform(input.set, trans.obj$`Violet G 550_40-A`$asinh.translist)
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
outDir = "./data_modified/TPHE_warping/targets_HD/scatter.VG.warped/"
dir.create(outDir)

fsApply(input.set.transformed.rangefix.warped.inverse.decomp.mod, function(i) write.FCS(x = i, filename = file.path(outDir, i@description$`$FIL`)))

#
#
#
#
#
#

##warp Violet.E (CD3); use scatter-VG warped HD samples w/ target
tphe.warped.hd <- list.files("./data_modified/TPHE_warping/targets_HD/scatter.VG.warped/", full.names = TRUE, pattern = ".fcs")
target <- "HD0189_TPHE_NGL091"

input.set <- fset.compensate(fset.trim(read.flowSet(tphe.warped.hd)))

markers <- get.markers(input.set[[1]])

dat <- as.data.frame(exprs(input.set[[grep(target, sampleNames(input.set))]]))
colnames(dat) <- markers

tmp <- subset(dat, select = "CD3")
ggplot(tmp, aes_string(x = names(tmp))) + geom_density()

asinh.vals <- list(a = -2, b = 1/400, c = 2.5)
ggplot(asinh(asinh.vals$a + asinh.vals$b * tmp) + asinh.vals$c, aes_string(x = names(tmp))) + geom_histogram(bins = 400)

asinh.transforms.TPHE <- readRDS("./results_R/FCS.transforms/TPHE/asinh.transforms.TPHE.rds")

colnames(input.set[[1]])[grep("Violet E", colnames(input.set[[1]]))]

asinh.transforms.TPHE[2] <- list(Violet.E = (list(asinh.translist = transformList('Violet E 585_42-A', 
                                                                               arcsinhTransform(transformationId = "arcsinh-transformation", 
                                                                                                a = -2, 
                                                                                                b = 1/400, 
                                                                                                c = 2.5)),
                                               asinh.vals = asinh.vals)))
names(asinh.transforms.TPHE)[2] <- "Violet.E"

saveRDS(asinh.transforms.TPHE, "./results_R/FCS.transforms/TPHE/asinh.transforms.TPHE.rds")

trans.obj <- readRDS("./results_R/FCS.transforms/TPHE/asinh.transforms.TPHE.rds")

input.set.transformed <- transform(input.set, trans.obj$Violet.E$asinh.translist)
input.set.transformed.rangefix <- fset.rangefix(input.set.transformed, 'Violet E 585_42-A')

sub.set <- fset.downsample(input.set.transformed.rangefix, 50000)

library(flowViz)

channel <- "Violet E 585_42-A"

densityplot(~ ., sub.set[, channel], filter = curv1Filter(channel, bwFac = 1))
sub.set.warped <- warpSet(sub.set[, channel], c(channel), monwrd = FALSE, peakNr = 2, clipRange = 0.01, bwFac = 1, target = sampleNames(sub.set[grep(target, sampleNames(sub.set))]))
densityplot(~ ., sub.set.warped[, channel], filter = curv1Filter(channel, bwFac = 1))

input.set.transformed.rangefix.warped <- warpSet(input.set.transformed.rangefix, c(channel), monwrd = FALSE, peakNr = 2, clipRange = 0.01, bwFac = 1, 
                                                 target = sampleNames(sub.set[grep(target, sampleNames(sub.set))]))
densityplot(~ ., input.set.transformed.rangefix.warped[, channel], filter = curv1Filter(channel, bwFac = 1))

## inverse transform warped set
asinh.vals <- trans.obj$Violet.E$asinh.vals

sinhTrans <- sinhTransform(transformationId = "inverse.arcsinh-transformation", a = asinh.vals$a, b = asinh.vals$b, c = asinh.vals$c)
sinh.translist <- transformList(channel, sinhTrans)

input.set.transformed.rangefix.warped.inverse <- transform(input.set.transformed.rangefix.warped, sinh.translist)

## decompensate
input.set.transformed.rangefix.warped.inverse.decomp <- fset.decompensate(input.set.transformed.rangefix.warped.inverse)

## modify transformed files to include new identifiers
## change/add keywords
input.set.transformed.rangefix.warped.inverse.decomp.mod <- fsApply(input.set.transformed.rangefix.warped.inverse.decomp, function (i) {
  i@description$`EXPERIMENT NAME` <- unlist(strsplit(i@description$`EXPERIMENT NAME`, " "))[[1]]
  i@description$`WARPED CHANNELS` <- "FSC.AHW_SSC.AH_VG_VE"
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
outDir = "./data_modified/TPHE_warping/targets_HD/scatter.VG.VE.warped/"
dir.create(outDir)

fsApply(input.set.transformed.rangefix.warped.inverse.decomp.mod, function(i) write.FCS(x = i, filename = file.path(outDir, i@description$`$FIL`)))

#
#
#
#
#
#

##warp Green.D (CD4); use scatter-VG.VE warped HD samples w/ target
tphe.warped.hd <- list.files("./data_modified/TPHE_warping/targets_HD/scatter.VG.VE.warped/", full.names = TRUE, pattern = ".fcs")
target <- "HD0189_TPHE_NGL091"

input.set <- fset.compensate(fset.trim(read.flowSet(tphe.warped.hd)))

markers <- get.markers(input.set[[1]])

dat <- as.data.frame(exprs(input.set[[grep(target, sampleNames(input.set))]]))
colnames(dat) <- markers

tmp <- subset(dat, select = "CD4")
ggplot(tmp, aes_string(x = names(tmp))) + geom_density()

asinh.vals <- list(a = -1, b = 1/450, c = 2)
ggplot(asinh(asinh.vals$a + asinh.vals$b * tmp) + asinh.vals$c, aes_string(x = names(tmp))) + geom_histogram(bins = 400)

asinh.transforms.TPHE <- readRDS("./results_R/FCS.transforms/TPHE/asinh.transforms.TPHE.rds")

channel <- colnames(input.set[[1]])[grep("Green D", colnames(input.set[[1]]))]

asinh.transforms.TPHE[3] <- list(Green.D = (list(asinh.translist = transformList(channel, 
                                                                                  arcsinhTransform(transformationId = "arcsinh-transformation", 
                                                                                                   a = -1, 
                                                                                                   b = 1/450, 
                                                                                                   c = 2)),
                                                  asinh.vals = asinh.vals)))
names(asinh.transforms.TPHE)[3] <- channel

saveRDS(asinh.transforms.TPHE, "./results_R/FCS.transforms/TPHE/asinh.transforms.TPHE.rds")

trans.obj <- readRDS("./results_R/FCS.transforms/TPHE/asinh.transforms.TPHE.rds")

input.set.transformed <- transform(input.set, trans.obj$`Green D 610_20-A`$asinh.translist)
input.set.transformed.rangefix <- fset.rangefix(input.set.transformed, channel)

sub.set <- fset.downsample(input.set.transformed.rangefix, 10000)

library(flowViz)

channel <- "Green D 610_20-A"

densityplot(~ ., sub.set[, channel], filter = curv1Filter(channel, bwFac = 1))
sub.set.warped <- warpSet(sub.set[, channel], c(channel), monwrd = FALSE, peakNr = 2, clipRange = 0.01, bwFac = 1, target = sampleNames(sub.set[grep(target, sampleNames(sub.set))]))
densityplot(~ ., sub.set.warped[, channel], filter = curv1Filter(channel, bwFac = 1))

input.set.transformed.rangefix.warped <- warpSet(input.set.transformed.rangefix, c(channel), monwrd = FALSE, peakNr = 2, clipRange = 0.01, bwFac = 1, 
                                                 target = sampleNames(sub.set[grep(target, sampleNames(sub.set))]))
densityplot(~ ., input.set.transformed.rangefix.warped[, channel], filter = curv1Filter(channel, bwFac = 1))

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
  i@description$`EXPERIMENT NAME` <- unlist(strsplit(i@description$`EXPERIMENT NAME`, " "))[[1]]
  i@description$`WARPED CHANNELS` <- "FSC.AHW_SSC.AH_VG_VE_GD"
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
outDir = "./data_modified/TPHE_warping/targets_HD/scatter.VG.VE.GD.warped/"
dir.create(outDir)

fsApply(input.set.transformed.rangefix.warped.inverse.decomp.mod, function(i) write.FCS(x = i, filename = file.path(outDir, i@description$`$FIL`)))

#
#
#
#
#
#

##warp Violet.A (CD8A); use scatter-VG.VE.GD warped HD samples w/ target
tphe.warped.hd <- list.files("./data_modified/TPHE_warping/targets_HD/scatter.VG.VE.GD.warped/", full.names = TRUE, pattern = ".fcs")
target <- "HD0189_TPHE_NGL091"

input.set <- fset.compensate(fset.trim(read.flowSet(tphe.warped.hd)))

markers <- get.markers(input.set[[1]])

dat <- as.data.frame(exprs(input.set[[grep(target, sampleNames(input.set))]]))
colnames(dat) <- markers

tmp <- subset(dat, select = "CD8A")
ggplot(tmp, aes_string(x = names(tmp))) + geom_density()

asinh.vals <- list(a = 0, b = 1/450, c = 1.5)
ggplot(asinh(asinh.vals$a + asinh.vals$b * tmp) + asinh.vals$c, aes_string(x = names(tmp))) + geom_histogram(bins = 400)

asinh.transforms.TPHE <- readRDS("./results_R/FCS.transforms/TPHE/asinh.transforms.TPHE.rds")

channel <- colnames(input.set[[1]])[grep("Violet A", colnames(input.set[[1]]))]

asinh.transforms.TPHE[4] <- list(Violet.A = (list(asinh.translist = transformList(channel, 
                                                                                 arcsinhTransform(transformationId = "arcsinh-transformation", 
                                                                                                  a = 0, 
                                                                                                  b = 1/450, 
                                                                                                  c = 1.5)),
                                                 asinh.vals = asinh.vals)))
names(asinh.transforms.TPHE)[4] <- channel

saveRDS(asinh.transforms.TPHE, "./results_R/FCS.transforms/TPHE/asinh.transforms.TPHE.rds")

trans.obj <- readRDS("./results_R/FCS.transforms/TPHE/asinh.transforms.TPHE.rds")

input.set.transformed <- transform(input.set, trans.obj[[channel]]$asinh.translist)
input.set.transformed.rangefix <- fset.rangefix(input.set.transformed, channel)

sub.set <- fset.downsample(input.set.transformed.rangefix, 50000)

library(flowViz)

densityplot(~ ., sub.set[, channel], filter = curv1Filter(channel, bwFac = 2))
sub.set.warped <- warpSet(sub.set[, channel], c(channel), monwrd = FALSE, peakNr = 2, clipRange = 0.01, bwFac = 2, target = sampleNames(sub.set[grep(target, sampleNames(sub.set))]))
densityplot(~ ., sub.set.warped[, channel], filter = curv1Filter(channel, bwFac = 1))

input.set.transformed.rangefix.warped <- warpSet(input.set.transformed.rangefix, c(channel), monwrd = FALSE, peakNr = 2, clipRange = 0.01, bwFac = 2, 
                                                 target = sampleNames(sub.set[grep(target, sampleNames(sub.set))]))
densityplot(~ ., input.set.transformed.rangefix.warped[, channel], filter = curv1Filter(channel, bwFac = 1))

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
  i@description$`EXPERIMENT NAME` <- unlist(strsplit(i@description$`EXPERIMENT NAME`, " "))[[1]]
  i@description$`WARPED CHANNELS` <- "FSC.AHW_SSC.AH_VG_VE_GD.VA"
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
outDir = "./data_modified/TPHE_warping/targets_HD/scatter.VG.VE.GD.VA.warped/"
dir.create(outDir)

fsApply(input.set.transformed.rangefix.warped.inverse.decomp.mod, function(i) write.FCS(x = i, filename = file.path(outDir, i@description$`$FIL`)))

#
#
#
#
#
#

##warp Violet.D (CD31); use fSOM subsets
tphe.warped.hd <- list.files("./data_modified/TPHE_warping/targets_HD/scatter.VG.VE.GD.VA.warped/fSOM_CD3.CD4.CD8.LIVEDEAD/fSOM_CD4pos/fSOM_CD4pos_scatter/fSOM_CD4pos_scatter/", 
                             full.names = TRUE, pattern = ".fcs")
target <- "HD0189_TPHE_NGL091"

input.set <- fset.compensate(fset.trim(read.flowSet(tphe.warped.hd)))

markers <- get.markers(input.set[[1]])

dat <- as.data.frame(exprs(input.set[[grep(target, sampleNames(input.set))]]))
colnames(dat) <- markers

tmp <- subset(dat, select = "CD31")
ggplot(tmp, aes_string(x = names(tmp))) + geom_density()

asinh.vals <- list(a = -2, b = 1/400, c = 2.5)
ggplot(asinh(asinh.vals$a + asinh.vals$b * tmp) + asinh.vals$c, aes_string(x = names(tmp))) + geom_histogram(bins = 400)

dat <- as.data.frame(exprs(input.set[[grep(target, sampleNames(input.set))]]))
colnames(dat) <- markers
dat$CD31 <- asinh(-2 + 1/400 * dat$CD31) + 2.5

fcs.explorer(dat, "CD31", "FCS.A", 20000)

asinh.transforms.TPHE <- readRDS("./results_R/FCS.transforms/TPHE/asinh.transforms.TPHE.rds")

channel <- colnames(input.set[[1]])[grep("Violet D", colnames(input.set[[1]]))]

asinh.transforms.TPHE[length(asinh.transforms.TPHE) + 1] <- list(channel = (list(asinh.translist = transformList(channel, 
                                                                                  arcsinhTransform(transformationId = "arcsinh-transformation", 
                                                                                                   a = -2, 
                                                                                                   b = 1/400, 
                                                                                                   c = 2.5)),
                                                  asinh.vals = asinh.vals)))
names(asinh.transforms.TPHE)[length(asinh.transforms.TPHE)] <- channel

saveRDS(asinh.transforms.TPHE, "./results_R/FCS.transforms/TPHE/asinh.transforms.TPHE.rds")

trans.obj <- readRDS("./results_R/FCS.transforms/TPHE/asinh.transforms.TPHE.rds")

input.set.transformed <- transform(input.set, trans.obj[[channel]]$asinh.translist)
input.set.transformed.rangefix <- fset.rangefix(input.set.transformed, channel)

library(flowViz)

densityplot(~ ., fset.downsample(input.set.transformed.rangefix, 10000)[, channel], filter = curv1Filter(channel, bwFac = 2))
sub.set.warped <- warpSet(fset.downsample(input.set.transformed.rangefix, 10000)[, channel], c(channel), monwrd = FALSE, peakNr = 2, clipRange = 0.01, bwFac = 2, 
                          target = sampleNames(input.set.transformed.rangefix[grep(target, sampleNames(input.set.transformed.rangefix))]))
densityplot(~ ., sub.set.warped[, channel], filter = curv1Filter(channel, bwFac = 2))

input.set.transformed.rangefix.warped <- warpSet(input.set.transformed.rangefix, c(channel), monwrd = FALSE, peakNr = 2, clipRange = 0.01, bwFac = 2, 
                                                 target = sampleNames(input.set.transformed.rangefix[grep(target, sampleNames(input.set.transformed.rangefix))]))
densityplot(~ ., input.set.transformed.rangefix.warped[, channel], filter = curv1Filter(channel, bwFac = 2))

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
  i@description$`$FIL` <- sub("FSC.AHW_SSC.AH_VG_VE_GD.VA", "FSC.AHW_SSC.AH_VG_VE_VD_VA_GD", i@description$`$FIL`)
  i
})

## write modified files
outDir = "./data_modified/TPHE_warping/targets_HD/scatter.VG.VE.GD.VA.warped/fSOM_CD3.CD4.CD8.LIVEDEAD/fSOM_CD4pos/fSOM_CD4pos_scatter/fSOM_CD4pos_scatter/"


fsApply(input.set.transformed.rangefix.warped.inverse.decomp.mod, function(i) write.FCS(x = i, filename = file.path(outDir, i@description$`$FIL`)))
