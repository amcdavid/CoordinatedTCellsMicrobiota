## source flow-related libraries and functions
source("./source_R/flow.src.R", chdir = TRUE)

## file paths
tphe.warped.files <- list.files("./data_modified/TPHE_warping/", recursive = TRUE, full.names = TRUE, pattern = "fda_WARPED_FSC.AHW_SSC.AHW_VG_VE_VA_GD_fSOM_CD4pos_subsampled.fcs")

input.set <- fset.trim(read.flowSet(tphe.warped.files, transformation = FALSE))

fSOM <- ReadInput(input.set, scale = TRUE)
fSOM$markers <- get.markers(input.set[[1]])

mset <- match(c("FSC.A", "FSC.W", "SSC.A", "SSC.W"), fSOM$markers)

fSOM <- BuildSOM(fSOM, colsToUse = mset, xdim= 5, ydim = 5, rlen = 10) # PASS 1
fSOM <- BuildMST(fSOM)

nbclust <- 6
fSOM[["metaclustering"]] <- as.factor(metaClustering_consensus(fSOM$map$codes, k = nbclust))
PlotStars(fSOM, backgroundValues = fSOM$metaclustering)

fSOM$heatmap <- pheatmap(scale(fsom.clusters.heatmap(fSOM)), scale = "row")
fcs.explorer(fsom.dat(fSOM), "FSC.A", "SSC.A", 20000)

fSOM$metaclustering.anno <- vector(mode = "list", length = length(levels(fSOM$metaclustering)))
fSOM$metaclustering.anno[[4]] <- "lymphs - singlets"

savefsom(fSOM, "./results_R/fSOMs/TPHE/fSOM_PASS2_singlets.rds")

## mapping

## source flow-related libraries and functions
source("./source_R/flow.src.R", chdir = TRUE)

## file paths
tphe.warped.files <- list.files("./data_modified/TPHE_warping/", recursive = TRUE, full.names = TRUE, pattern = "fda_WARPED_FSC.AHW_SSC.AHW_VG_VE_VA_GD_fSOM_CD4pos_.fcs")

##load fSOM
fSOM <- readRDS("./results_R/fSOMs/TPHE/fSOM_PASS2_singlets.rds") # needs stored $transform.object

##
newfcs.fromsom.scatter(fSOM.object = fSOM, 
                    meta.cluster = grep("lymphs", fSOM$metaclustering.anno), 
                    FCS.paths = tphe.warped.files,
                    basefolder.name = "fSOM_CD4pos_scatter", 
                    write.count.files = TRUE, 
                    write.fcs.files = TRUE, 
                    fcsfolder.name = "fSOM_CD4pos_scatter")
