## TPHE data needs to be normalized to reduce HTS-related batch effects in earlier experiments; 'warping' is defined and implemented per channel

## file paths
tphe.files <- list.files("./data_source/TPHE2_with.spillover", recursive = TRUE, full.names = TRUE, pattern = ".fcs")
tphe.files <- tphe.files[grep("HD|_TPHE.fcs", tphe.files)] # only experimental (subjects) or healthy donor (controls)

batches <- batch.list.edit(tphe.files[-grep("NGL063", tphe.files)]) # batched by experiment number

## warping target
target.file <- tphe.files[grep("HD0189_TPHE_NGL091", tphe.files)] # 'gold-standard' file -- well distributed
target <- "HD0189_TPHE_NGL091"

## transforms
trans.obj <- readRDS("./results_R/FCS.transforms/TPHE/asinh.transforms.TPHE.rds")
trans.obj <- trans.obj[grep('Violet G|Violet E|Violet A|Green D', names(trans.obj))]

## read in FCS files by batch to test out warping parameters
batch.name <- "NGL047"

input.set <- fset.compensate(fset.trim(read.flowSet(c(batches[batch.name][[1]], target.file), 
                                                    transformation = FALSE)))

input.set <- transform.set(input.set)

##
channels <- colnames(input.set)
markers <- get.markers(input.set[[1]])

## initialize warp.parameters object
TPHE.warp.parameters <- vector(mode = "list", length = length(markers[-grep("FSC|SSC|Time", markers)]))
names(TPHE.warp.parameters) <- markers[-grep("FSC|SSC|Time", markers)]
channel.names <- channels[-grep("FSC|SSC|Time", channels)]
for(i in seq(names(TPHE.warp.parameters))){
  TPHE.warp.parameters[[i]] <- list(channel = channel.names[i])
}

## define warp parameters per channel
library(flowViz)
library(flowStats)

##
marker.name <- "LIVEDEAD"
bwFac.set <- 1.5
densityplot(~ ., fset.downsample(input.set, 5000)[, channels[grep(marker.name, markers)]], filter = curv1Filter(channels[grep(marker.name, markers)], bwFac = bwFac.set))
test.set <- warpSet(input.set, c(channels[grep(marker.name, markers)]), monwrd = FALSE, peakNr = 1, clipRange = 0.01, bwFac = bwFac.set, 
                    target = sampleNames(input.set)[grep(target, sampleNames(input.set))])
densityplot(~ ., fset.downsample(test.set, 5000)[, channels[grep(marker.name, markers)]], filter = curv1Filter(channels[grep(marker.name, markers)], bwFac = bwFac.set))
TPHE.warp.parameters[[marker.name]]["monwrd"] <- FALSE
TPHE.warp.parameters[[marker.name]]["peakNr"] <- 1
TPHE.warp.parameters[[marker.name]]["clipRange"] <- 0.01
TPHE.warp.parameters[[marker.name]]["bwFac"] <- bwFac.set


marker.name <- "CD3"
bwFac.set <- 1
densityplot(~ ., fset.downsample(input.set, 5000)[, channels[grep(paste0(marker.name, "$"), markers)]], 
            filter = curv1Filter(channels[grep(paste0(marker.name, "$"), markers)], bwFac = bwFac.set))
test.set <- warpSet(input.set, c(channels[grep(paste0(marker.name, "$"), markers)]), monwrd = FALSE, peakNr = 2, clipRange = 0.01, bwFac = bwFac.set, 
                    target = sampleNames(input.set)[grep(target, sampleNames(input.set))])
densityplot(~ ., fset.downsample(test.set, 5000)[, channels[grep(paste0(marker.name, "$"), markers)]], 
            filter = curv1Filter(channels[grep(paste0(marker.name, "$"), markers)], bwFac = bwFac.set))
TPHE.warp.parameters[[marker.name]]["monwrd"] <- FALSE
TPHE.warp.parameters[[marker.name]]["peakNr"] <- 2
TPHE.warp.parameters[[marker.name]]["clipRange"] <- 0.01
TPHE.warp.parameters[[marker.name]]["bwFac"] <- bwFac.set


marker.name <- "CD4"
bwFac.set <- 1.25
densityplot(~ ., fset.downsample(input.set, 5000)[, channels[grep(paste0(marker.name, "$"), markers)]], 
            filter = curv1Filter(channels[grep(paste0(marker.name, "$"), markers)], bwFac = bwFac.set))
test.set <- warpSet(input.set, c(channels[grep(paste0(marker.name, "$"), markers)]), monwrd = FALSE, peakNr = 2, clipRange = 0.1, bwFac = bwFac.set, 
                    target = sampleNames(input.set)[grep(target, sampleNames(input.set))])
densityplot(~ ., fset.downsample(test.set, 5000)[, channels[grep(paste0(marker.name, "$"), markers)]], 
            filter = curv1Filter(channels[grep(paste0(marker.name, "$"), markers)], bwFac = bwFac.set))
TPHE.warp.parameters[[marker.name]]["monwrd"] <- FALSE
TPHE.warp.parameters[[marker.name]]["peakNr"] <- 2
TPHE.warp.parameters[[marker.name]]["clipRange"] <- 0.1
TPHE.warp.parameters[[marker.name]]["bwFac"] <- bwFac.set


marker.name <- "CD8A"
bwFac.set <- 2.5
densityplot(~ ., fset.downsample(input.set, 5000)[, channels[grep(paste0(marker.name, "$"), markers)]], 
            filter = curv1Filter(channels[grep(paste0(marker.name, "$"), markers)], bwFac = bwFac.set))
test.set <- warpSet(input.set, c(channels[grep(paste0(marker.name, "$"), markers)]), monwrd = FALSE, peakNr = 2, clipRange = 0.01, bwFac = bwFac.set, 
                    target = sampleNames(input.set)[grep(target, sampleNames(input.set))])
densityplot(~ ., fset.downsample(test.set, 5000)[, channels[grep(paste0(marker.name, "$"), markers)]], 
            filter = curv1Filter(channels[grep(paste0(marker.name, "$"), markers)], bwFac = bwFac.set))
TPHE.warp.parameters[[marker.name]]["monwrd"] <- FALSE
TPHE.warp.parameters[[marker.name]]["peakNr"] <- 2
TPHE.warp.parameters[[marker.name]]["clipRange"] <- 0.01
TPHE.warp.parameters[[marker.name]]["bwFac"] <- bwFac.set


saveRDS(TPHE.warp.parameters, "./results_R/FCS.transforms/TPHE/TPHE.warp.parameters.rds")