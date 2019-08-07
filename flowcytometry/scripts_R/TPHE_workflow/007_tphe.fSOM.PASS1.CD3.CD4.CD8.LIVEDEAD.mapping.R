## fSOM PASS1 mapping
library(flowCore)
library(FlowSOM)

## file paths
tphe.warped.files <- list.files("./data_modified/TPHE_warping", recursive = TRUE, full.names = TRUE, pattern = "WARPED_.fcs")
tphe.warped.files <- tphe.warped.files[grep("CD4.CD8A.CD3.LIVEDEAD_warped", dirname(tphe.warped.files))]

## load fSOM
fSOM <- readRDS("./results_R/fSOMs/TPHE/2018-11-11.fSOM_CD3.CD4.CD8.LIVEDEAD.PASS1.rds") # needs stored $transform.object

##
newfcs.fromfsom.singlefolder(fSOM.object = fSOM,
                             meta.cluster = grep("CD4", fSOM$metaclustering.anno),
                             FCS.paths = tphe.warped.files,
                             basefolder.name = "./data_modified/TPHE_warping/fSOM_CD3.CD4.CD8.LIVEDEAD_PASS1",
                             fcsfolder.name =  "fSOM_CD4p",
                             write.count.files = TRUE,
                             write.fcs.files = TRUE)

##
newfcs.fromfsom.singlefolder(fSOM.object = fSOM,
                             meta.cluster = grep("CD8", fSOM$metaclustering.anno),
                             FCS.paths = tphe.warped.files,
                             basefolder.name = "./data_modified/TPHE_warping/fSOM_CD3.CD4.CD8.LIVEDEAD_PASS1",
                             fcsfolder.name =  "fSOM_CD8p",
                             write.count.files = FALSE,
                             write.fcs.files = TRUE)



## 'target file'
target.file <- list.files("./data_source/TPHE2_with.spillover/NGL091/", full.names = TRUE, pattern = "HD0189")

newfcs.fromfsom.singlefolder(fSOM.object = fSOM,
                             meta.cluster = grep("CD4", fSOM$metaclustering.anno),
                             FCS.paths = target.file,
                             basefolder.name = "./data_modified/TPHE_warping/fSOM_CD3.CD4.CD8.LIVEDEAD_PASS1",
                             fcsfolder.name =  "fSOM_CD4p",
                             write.count.files = TRUE,
                             write.fcs.files = TRUE)

newfcs.fromfsom.singlefolder(fSOM.object = fSOM,
                             meta.cluster = grep("CD8", fSOM$metaclustering.anno),
                             FCS.paths = target.file,
                             basefolder.name = "./data_modified/TPHE_warping/fSOM_CD3.CD4.CD8.LIVEDEAD_PASS1",
                             fcsfolder.name =  "fSOM_CD8p",
                             write.count.files = TRUE,
                             write.fcs.files = TRUE)
