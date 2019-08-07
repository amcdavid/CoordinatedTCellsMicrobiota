library(flowCore)
library(flowStats)
library(flowViz)

## define warping parameter for scatter channels using fSOM subset; set pattern argument for either CD4 or CD8

## file paths
fsom.files <- list.files("./data_modified/TPHE_warping/fSOM_CD3.CD4.CD8.LIVEDEAD_PASS1/", full.names = TRUE, pattern = "CD8p_.fcs")
batch.names <- paste0("NGL0", seq(47, 91, 2)) ; batch.names <- batch.names[-grep("NGL063", batch.names)]

## warping target
target.file <- fsom.files[grep("HD0189_TPHE_NGL091", fsom.files)]
target <- "HD0189_TPHE_NGL091"; fsom.files <- fsom.files[-grep(target, fsom.files)]

##
batch.name <- "NGL087"

input.set <- read.flowSet(c(fsom.files[grep(batch.name, fsom.files)], target.file), transformation = FALSE)

##
channels <- colnames(input.set)
markers <- get.markers(input.set[[1]])

## initialize warp.parameters object
TPHE.warp.scatter.parameters <- vector(mode = "list", length = length(markers[grep("FSC|SSC", markers)]))
names(TPHE.warp.scatter.parameters) <- markers[grep("FSC|SSC", markers)]
channel.names <- channels[grep("FSC|SSC", channels)]
for(i in seq(names(TPHE.warp.scatter.parameters))){
  TPHE.warp.scatter.parameters[[i]] <- list(channel = channel.names[i])
}

## define warp parameters per scatter channel; iteratively test by modifying warpSet parameters then store in list


##
marker.name <- "FSC.A"
bwFac.set <- 2
densityplot(~ ., fset.downsample(input.set, 5000)[, channels[grep(marker.name, markers)]], 
            filter = curv1Filter(channels[grep(marker.name, markers)], bwFac = bwFac.set))
test.set <- warpSet(input.set, c(channels[grep(marker.name, markers)]), monwrd = FALSE , peakNr = 1, clipRange = 0.01, bwFac = bwFac.set, 
                    target = sampleNames(input.set)[grep(target, sampleNames(input.set))])
densityplot(~ ., fset.downsample(test.set, 5000)[, channels[grep(marker.name, markers)]], 
            filter = curv1Filter(channels[grep(marker.name, markers)], bwFac = bwFac.set))
TPHE.warp.scatter.parameters[[marker.name]]["monwrd"] <- FALSE
TPHE.warp.scatter.parameters[[marker.name]]["peakNr"] <- 1
TPHE.warp.scatter.parameters[[marker.name]]["clipRange"] <- 0.01
TPHE.warp.scatter.parameters[[marker.name]]["bwFac"] <- bwFac.set

##
marker.name <- "FSC.H"
bwFac.set <- 2
densityplot(~ ., fset.downsample(input.set, 5000)[, channels[grep(marker.name, markers)]], filter = curv1Filter(channels[grep(marker.name, markers)], bwFac = bwFac.set))
test.set <- warpSet(input.set, c(channels[grep(marker.name, markers)]), monwrd = FALSE, peakNr = 1, clipRange = 0.2, bwFac = bwFac.set, 
                    target = sampleNames(input.set)[grep(target, sampleNames(input.set))])
densityplot(~ ., fset.downsample(test.set, 5000)[, channels[grep(marker.name, markers)]], filter = curv1Filter(channels[grep(marker.name, markers)], bwFac = bwFac.set))
TPHE.warp.scatter.parameters[[marker.name]]["monwrd"] <- FALSE
TPHE.warp.scatter.parameters[[marker.name]]["peakNr"] <- 1
TPHE.warp.scatter.parameters[[marker.name]]["clipRange"] <- 0.1
TPHE.warp.scatter.parameters[[marker.name]]["bwFac"] <- bwFac.set

##
marker.name <- "FSC.W"
bwFac.set <- 2
densityplot(~ ., fset.downsample(input.set, 5000)[, channels[grep(marker.name, markers)]], 
            filter = curv1Filter(channels[grep(marker.name, markers)], bwFac = bwFac.set),
            xlim = c(50000, 110000))
test.set <- warpSet(input.set, c(channels[grep(marker.name, markers)]), monwrd = FALSE, peakNr = 1, clipRange = 0.1, bwFac = bwFac.set, 
                    target = sampleNames(input.set)[grep(target, sampleNames(input.set))])
densityplot(~ ., fset.downsample(test.set, 5000)[, channels[grep(marker.name, markers)]], 
            filter = curv1Filter(channels[grep(marker.name, markers)], bwFac = bwFac.set),
            xlim = c(50000, 110000))
TPHE.warp.scatter.parameters[[marker.name]]["monwrd"] <- FALSE
TPHE.warp.scatter.parameters[[marker.name]]["peakNr"] <- 1
TPHE.warp.scatter.parameters[[marker.name]]["clipRange"] <- 0.1
TPHE.warp.scatter.parameters[[marker.name]]["bwFac"] <- bwFac.set

##
marker.name <- "SSC.A"
bwFac.set <- 2
densityplot(~ ., fset.downsample(input.set, 5000)[, channels[grep(marker.name, markers)]], 
            filter = curv1Filter(channels[grep(marker.name, markers)], bwFac = bwFac.set))
test.set <- warpSet(input.set, c(channels[grep(marker.name, markers)]), monwrd = FALSE, peakNr = 1, clipRange = 0.1, bwFac = bwFac.set, 
                    target = sampleNames(input.set)[grep(target, sampleNames(input.set))])
densityplot(~ ., fset.downsample(test.set, 5000)[, channels[grep(marker.name, markers)]], 
            filter = curv1Filter(channels[grep(marker.name, markers)], bwFac = bwFac.set))
TPHE.warp.scatter.parameters[[marker.name]]["monwrd"] <- FALSE
TPHE.warp.scatter.parameters[[marker.name]]["peakNr"] <- 1
TPHE.warp.scatter.parameters[[marker.name]]["clipRange"] <- 0.1
TPHE.warp.scatter.parameters[[marker.name]]["bwFac"] <- bwFac.set

##
marker.name <- "SSC.H"
bwFac.set <- 2
densityplot(~ ., fset.downsample(input.set, 5000)[, channels[grep(marker.name, markers)]], 
            filter = curv1Filter(channels[grep(marker.name, markers)], bwFac = bwFac.set))
test.set <- warpSet(input.set, c(channels[grep(marker.name, markers)]), monwrd = FALSE, peakNr = 1, clipRange = 0.1, bwFac = bwFac.set, 
                    target = sampleNames(input.set)[grep(target, sampleNames(input.set))])
densityplot(~ ., fset.downsample(test.set, 5000)[, channels[grep(marker.name, markers)]], 
            filter = curv1Filter(channels[grep(marker.name, markers)], bwFac = bwFac.set))
TPHE.warp.scatter.parameters[[marker.name]]["monwrd"] <- FALSE
TPHE.warp.scatter.parameters[[marker.name]]["peakNr"] <- 1
TPHE.warp.scatter.parameters[[marker.name]]["clipRange"] <- 0.1
TPHE.warp.scatter.parameters[[marker.name]]["bwFac"] <- bwFac.set

##
marker.name <- "SSC.W"
bwFac.set <- 2
densityplot(~ ., fset.downsample(input.set, 5000)[, channels[grep(marker.name, markers)]], 
            filter = curv1Filter(channels[grep(marker.name, markers)], bwFac = bwFac.set),
            xlim = c(50000, 110000))
test.set <- warpSet(input.set, c(channels[grep(marker.name, markers)]), monwrd = FALSE, peakNr = 1, clipRange = 0.1, bwFac = bwFac.set, 
                    target = sampleNames(input.set)[grep(target, sampleNames(input.set))])
densityplot(~ ., fset.downsample(test.set, 5000)[, channels[grep(marker.name, markers)]], 
            filter = curv1Filter(channels[grep(marker.name, markers)], bwFac = bwFac.set),
            xlim = c(50000, 110000))
TPHE.warp.scatter.parameters[[marker.name]]["monwrd"] <- FALSE
TPHE.warp.scatter.parameters[[marker.name]]["peakNr"] <- 1
TPHE.warp.scatter.parameters[[marker.name]]["clipRange"] <- 0.1
TPHE.warp.scatter.parameters[[marker.name]]["bwFac"] <- bwFac.set


saveRDS(TPHE.warp.scatter.parameters, "./results_R/FCS.transforms/TPHE/TPHE.warp.scatter.parameters.rds")
