## source flow-related libraries and functions
source("./source_R/flow.src.R", chdir = TRUE)

## file paths
tphe.warped.files <- list.files("./data_modified/TPHE_warping/fSOM_CD3.CD4.CD8.LIVEDEAD_PASS1/", 
                                recursive = FALSE, 
                                full.names = TRUE, 
                                pattern = ".fcs")

batches <- paste0("NGL0", seq(47, 91, 2)) ; batches[-grep("NGL063", batches)]

## warping targets
targets <- list.files("./data_modified/TPHE_warping/fSOM_CD3.CD4.CD8.LIVEDEAD_PASS1/", full.names = TRUE)
target <- targets[grep("HD0189_TPHE_NGL091", targets)]

## transforms
trans.obj <- readRDS("./results_R/FCS.transforms/TPHE/asinh.transforms.TPHE.rds")

##
batch.name <- "NGL047"

input.set <- fset.compensate(fset.trim(read.flowSet(c(batches[batch.name][[1]], target), 
                                                    transformation = FALSE)))

input.set <- transform.set(input.set)

##
channels <- colnames(input.set)
markers <- get.markers(input.set[[1]])

TPHE.CD4.warp.parameters <- vector(mode = "list", length = length(markers[-grep("FSC|SSC|Time", markers)]))
names(TPHE.CD4.warp.parameters) <- markers[-grep("FSC|SSC|Time", markers)]
channel.names <- channels[-grep("FSC|SSC|Time", channels)]

for(i in seq(names(TPHE.CD4.warp.parameters))){
  TPHE.CD4.warp.parameters[[i]] <- list(channel = channel.names[i])
}

##
marker.name <- "PERFORIN"
bwFac.set <- 1
densityplot(~ ., fset.downsample(input.set, 5000)[, channels[grep(marker.name, markers)]], filter = curv1Filter(channels[grep(marker.name, markers)], bwFac = bwFac.set))
test.set <- warpSet(input.set, c(channels[grep(marker.name, markers)]), monwrd = FALSE, peakNr = 2, clipRange = 0.01, bwFac = bwFac.set, 
                    target = sampleNames(input.set)[grep("HD", sampleNames(input.set))])
densityplot(~ ., fset.downsample(test.set, 5000)[, channels[grep(marker.name, markers)]], filter = curv1Filter(channels[grep(marker.name, markers)], bwFac = bwFac.set))
TPHE.CD4.warp.parameters[[marker.name]]["monwrd"] <- FALSE
TPHE.CD4.warp.parameters[[marker.name]]["peakNr"] <- 2
TPHE.CD4.warp.parameters[[marker.name]]["clipRange"] <- 0.01
TPHE.CD4.warp.parameters[[marker.name]]["bwFac"] <- bwFac.set

saveRDS(TPHE.CD4.warp.parameters, "TPHE.CD4.warp.parameters.rds")
