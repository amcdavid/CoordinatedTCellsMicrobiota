## define warp parameters on CD4+ subset

## file paths
fsom.files <- list.files("./data_modified/TPHE_warping/fSOM_CD4p.SINGLETS_PASS2/", 
                         recursive = FALSE, full.names = TRUE, pattern = ".fcs")
batch.names <- paste0("NGL0", seq(47, 91, 2)) ; batch.names <- batch.names[-grep("NGL063", batch.names)]

## warping target
target.file <- fsom.files[grep("HD0189_TPHE_NGL091", fsom.files)]
target <- "HD0189_TPHE_NGL091"; fsom.files <- fsom.files[-grep(target, fsom.files)]

## transforms
trans.obj <- readRDS("./results_R/FCS.transforms/TPHE/asinh.transforms.TPHE.rds")

##
batch.name <- "NGL047"

input.set <- fset.compensate(read.flowSet(c(fsom.files[grep(batch.name, fsom.files)], target.file), 
                                          transformation = FALSE))

input.set <- transform.set(input.set)

##
channels <- colnames(input.set)
markers <- get.markers(input.set[[1]])

##
TPHE.CD4.warp.parameters <- readRDS("./results_R/FCS.transforms/TPHE/TPHE.warp.parameters.edit.rds")

##
marker.name <- "PERFORIN"
bwFac.set <- 1.5
densityplot(~ ., fset.downsample(input.set, 5000)[, channels[grep(marker.name, markers)]], filter = curv1Filter(channels[grep(marker.name, markers)], bwFac = bwFac.set))
test.set <- warpSet(input.set, c(channels[grep(marker.name, markers)]), monwrd = FALSE, peakNr = 1, clipRange = 0.01, bwFac = bwFac.set, 
                    target = sampleNames(input.set)[grep(target, sampleNames(input.set))])
densityplot(~ ., fset.downsample(test.set, 5000)[, channels[grep(marker.name, markers)]], filter = curv1Filter(channels[grep(marker.name, markers)], bwFac = bwFac.set))
TPHE.CD4.warp.parameters[[marker.name]]["monwrd"] <- FALSE
TPHE.CD4.warp.parameters[[marker.name]]["peakNr"] <- 1
TPHE.CD4.warp.parameters[[marker.name]]["clipRange"] <- 0.01
TPHE.CD4.warp.parameters[[marker.name]]["bwFac"] <- bwFac.set


marker.name <- "CD122"
bwFac.set <- 1.5
densityplot(~ ., fset.downsample(input.set, 5000)[, channels[grep(marker.name, markers)]], filter = curv1Filter(channels[grep(marker.name, markers)], bwFac = bwFac.set))
test.set <- warpSet(input.set, c(channels[grep(marker.name, markers)]), monwrd = FALSE, peakNr = 1, clipRange = 0.01, bwFac = bwFac.set, 
                    target = sampleNames(input.set)[grep(target, sampleNames(input.set))])
densityplot(~ ., fset.downsample(test.set, 5000)[, channels[grep(marker.name, markers)]], filter = curv1Filter(channels[grep(marker.name, markers)], bwFac = bwFac.set))
TPHE.CD4.warp.parameters[[marker.name]]["monwrd"] <- FALSE
TPHE.CD4.warp.parameters[[marker.name]]["peakNr"] <- 1
TPHE.CD4.warp.parameters[[marker.name]]["clipRange"] <- 0.01
TPHE.CD4.warp.parameters[[marker.name]]["bwFac"] <- bwFac.set


marker.name <- "CD57"
bwFac.set <- 1
densityplot(~ ., fset.downsample(input.set, 5000)[, channels[grep(marker.name, markers)]], filter = curv1Filter(channels[grep(marker.name, markers)], bwFac = bwFac.set))
test.set <- warpSet(input.set, c(channels[grep(marker.name, markers)]), monwrd = FALSE, peakNr = 1, clipRange = 0.01, bwFac = bwFac.set, 
                    target = sampleNames(input.set)[grep(target, sampleNames(input.set))])
densityplot(~ ., fset.downsample(test.set, 5000)[, channels[grep(marker.name, markers)]], filter = curv1Filter(channels[grep(marker.name, markers)], bwFac = bwFac.set))
TPHE.CD4.warp.parameters[[marker.name]]["monwrd"] <- FALSE
TPHE.CD4.warp.parameters[[marker.name]]["peakNr"] <- 1
TPHE.CD4.warp.parameters[[marker.name]]["clipRange"] <- 0.01
TPHE.CD4.warp.parameters[[marker.name]]["bwFac"] <- bwFac.set


marker.name <- "CD28"
bwFac.set <- 1.5
densityplot(~ ., fset.downsample(input.set, 5000)[, channels[grep(marker.name, markers)]], filter = curv1Filter(channels[grep(marker.name, markers)], bwFac = bwFac.set))
test.set <- warpSet(input.set, c(channels[grep(marker.name, markers)]), monwrd = FALSE, peakNr = 1, clipRange = 0.01, bwFac = bwFac.set, 
                    target = sampleNames(input.set)[grep(target, sampleNames(input.set))])
densityplot(~ ., fset.downsample(test.set, 5000)[, channels[grep(marker.name, markers)]], filter = curv1Filter(channels[grep(marker.name, markers)], bwFac = bwFac.set))
TPHE.CD4.warp.parameters[[marker.name]]["monwrd"] <- FALSE
TPHE.CD4.warp.parameters[[marker.name]]["peakNr"] <- 1
TPHE.CD4.warp.parameters[[marker.name]]["clipRange"] <- 0.01
TPHE.CD4.warp.parameters[[marker.name]]["bwFac"] <- bwFac.set


marker.name <- "FOXP3"
bwFac.set <- 1.5
densityplot(~ ., fset.downsample(input.set, 5000)[, channels[grep(marker.name, markers)]], filter = curv1Filter(channels[grep(marker.name, markers)], bwFac = bwFac.set))
test.set <- warpSet(input.set, c(channels[grep(marker.name, markers)]), monwrd = FALSE, peakNr = 1, clipRange = 0.01, bwFac = bwFac.set, 
                    target = sampleNames(input.set)[grep(target, sampleNames(input.set))])
densityplot(~ ., fset.downsample(test.set, 5000)[, channels[grep(marker.name, markers)]], filter = curv1Filter(channels[grep(marker.name, markers)], bwFac = bwFac.set))
TPHE.CD4.warp.parameters[[marker.name]]["monwrd"] <- FALSE
TPHE.CD4.warp.parameters[[marker.name]]["peakNr"] <- 1
TPHE.CD4.warp.parameters[[marker.name]]["clipRange"] <- 0.01
TPHE.CD4.warp.parameters[[marker.name]]["bwFac"] <- bwFac.set


marker.name <- "CD197"
bwFac.set <- 1.5
densityplot(~ ., fset.downsample(input.set, 5000)[, channels[grep(marker.name, markers)]], filter = curv1Filter(channels[grep(marker.name, markers)], bwFac = bwFac.set))
test.set <- warpSet(input.set, c(channels[grep(marker.name, markers)]), monwrd = FALSE, peakNr = 2, clipRange = 0.01, bwFac = bwFac.set, 
                    target = sampleNames(input.set)[grep(target, sampleNames(input.set))])
densityplot(~ ., fset.downsample(test.set, 5000)[, channels[grep(marker.name, markers)]], filter = curv1Filter(channels[grep(marker.name, markers)], bwFac = bwFac.set))
TPHE.CD4.warp.parameters[[marker.name]]["monwrd"] <- FALSE
TPHE.CD4.warp.parameters[[marker.name]]["peakNr"] <- 2
TPHE.CD4.warp.parameters[[marker.name]]["clipRange"] <- 0.01
TPHE.CD4.warp.parameters[[marker.name]]["bwFac"] <- bwFac.set


marker.name <- "CD185"
bwFac.set <- 1.5
densityplot(~ ., fset.downsample(input.set, 5000)[, channels[grep(marker.name, markers)]], filter = curv1Filter(channels[grep(marker.name, markers)], bwFac = bwFac.set))
test.set <- warpSet(input.set, c(channels[grep(marker.name, markers)]), monwrd = FALSE, peakNr = 2, clipRange = 0.01, bwFac = bwFac.set, 
                    target = sampleNames(input.set)[grep(target, sampleNames(input.set))])
densityplot(~ ., fset.downsample(test.set, 5000)[, channels[grep(marker.name, markers)]], filter = curv1Filter(channels[grep(marker.name, markers)], bwFac = bwFac.set))
TPHE.CD4.warp.parameters[[marker.name]]["monwrd"] <- FALSE
TPHE.CD4.warp.parameters[[marker.name]]["peakNr"] <- 2
TPHE.CD4.warp.parameters[[marker.name]]["clipRange"] <- 0.01
TPHE.CD4.warp.parameters[[marker.name]]["bwFac"] <- bwFac.set


marker.name <- "KLRG1"
bwFac.set <- 1.5
densityplot(~ ., fset.downsample(input.set, 5000)[, channels[grep(marker.name, markers)]], filter = curv1Filter(channels[grep(marker.name, markers)], bwFac = bwFac.set))
test.set <- warpSet(input.set, c(channels[grep(marker.name, markers)]), monwrd = FALSE, peakNr = 1, clipRange = 0.01, bwFac = bwFac.set, 
                    target = sampleNames(input.set)[grep(target, sampleNames(input.set))])
densityplot(~ ., fset.downsample(test.set, 5000)[, channels[grep(marker.name, markers)]], filter = curv1Filter(channels[grep(marker.name, markers)], bwFac = bwFac.set))
TPHE.CD4.warp.parameters[[marker.name]]["monwrd"] <- FALSE
TPHE.CD4.warp.parameters[[marker.name]]["peakNr"] <- 1
TPHE.CD4.warp.parameters[[marker.name]]["clipRange"] <- 0.01
TPHE.CD4.warp.parameters[[marker.name]]["bwFac"] <- bwFac.set


marker.name <- "CD45RO"
bwFac.set <- 1.5
densityplot(~ ., fset.downsample(input.set, 5000)[, channels[grep(marker.name, markers)]], 
            filter = curv1Filter(channels[grep(marker.name, markers)], bwFac = bwFac.set))
test.set <- warpSet(input.set, c(channels[grep(marker.name, markers)]), monwrd = FALSE, peakNr = 2, clipRange = 0.3, bwFac = bwFac.set, 
                    target = sampleNames(input.set)[grep(target, sampleNames(input.set))])
densityplot(~ ., fset.downsample(test.set, 5000)[, channels[grep(marker.name, markers)]], filter = curv1Filter(channels[grep(marker.name, markers)], bwFac = bwFac.set))
TPHE.CD4.warp.parameters[[marker.name]]["monwrd"] <- FALSE
TPHE.CD4.warp.parameters[[marker.name]]["peakNr"] <- 2
TPHE.CD4.warp.parameters[[marker.name]]["clipRange"] <- 0.01
TPHE.CD4.warp.parameters[[marker.name]]["bwFac"] <- bwFac.set


marker.name <- "CD127"
bwFac.set <- 1.5
densityplot(~ ., fset.downsample(input.set, 5000)[, channels[grep(marker.name, markers)]], filter = curv1Filter(channels[grep(marker.name, markers)], bwFac = bwFac.set))
test.set <- warpSet(input.set, c(channels[grep(marker.name, markers)]), monwrd = FALSE, peakNr = 1, clipRange = 0.01, bwFac = bwFac.set, 
                    target = sampleNames(input.set)[grep(target, sampleNames(input.set))])
densityplot(~ ., fset.downsample(test.set, 5000)[, channels[grep(marker.name, markers)]], filter = curv1Filter(channels[grep(marker.name, markers)], bwFac = bwFac.set))
TPHE.CD4.warp.parameters[[marker.name]]["monwrd"] <- FALSE
TPHE.CD4.warp.parameters[[marker.name]]["peakNr"] <- 1
TPHE.CD4.warp.parameters[[marker.name]]["clipRange"] <- 0.01
TPHE.CD4.warp.parameters[[marker.name]]["bwFac"] <- bwFac.set


marker.name <- "CD31"
bwFac.set <- 1
densityplot(~ ., fset.downsample(input.set, 5000)[, channels[grep(marker.name, markers)]], 
            filter = curv1Filter(channels[grep(marker.name, markers)], bwFac = bwFac.set),
            main = paste("Non-warped", marker.name, sep = " "))
test.set <- warpSet(input.set, c(channels[grep(marker.name, markers)]), monwrd = FALSE, peakNr = 2, clipRange = 0.01, bwFac = bwFac.set, 
                    target = sampleNames(input.set)[grep(target, sampleNames(input.set))])
densityplot(~ ., fset.downsample(test.set, 5000)[, channels[grep(marker.name, markers)]], 
            filter = curv1Filter(channels[grep(marker.name, markers)], bwFac = bwFac.set),
            main = paste("Warped", marker.name, sep = " "))
TPHE.CD4.warp.parameters[[marker.name]]["monwrd"] <- FALSE
TPHE.CD4.warp.parameters[[marker.name]]["peakNr"] <- 2
TPHE.CD4.warp.parameters[[marker.name]]["clipRange"] <- 0.01
TPHE.CD4.warp.parameters[[marker.name]]["bwFac"] <- bwFac.set


marker.name <- "GRZB"
bwFac.set <- 1.5
densityplot(~ ., fset.downsample(input.set, 5000)[, channels[grep(marker.name, markers)]], 
            filter = curv1Filter(channels[grep(marker.name, markers)], bwFac = bwFac.set))
test.set <- warpSet(input.set, c(channels[grep(marker.name, markers)]), monwrd = FALSE, peakNr = 1, clipRange = 0.01, bwFac = bwFac.set, 
                    target = sampleNames(input.set)[grep(target, sampleNames(input.set))])
densityplot(~ ., fset.downsample(test.set, 5000)[, channels[grep(marker.name, markers)]], filter = curv1Filter(channels[grep(marker.name, markers)], bwFac = bwFac.set))
TPHE.CD4.warp.parameters[[marker.name]]["monwrd"] <- FALSE
TPHE.CD4.warp.parameters[[marker.name]]["peakNr"] <- 1
TPHE.CD4.warp.parameters[[marker.name]]["clipRange"] <- 0.01
TPHE.CD4.warp.parameters[[marker.name]]["bwFac"] <- bwFac.set

saveRDS(TPHE.CD4.warp.parameters, "TPHE.CD4.warp.parameters.rds")