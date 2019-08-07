## fSOM PASS3 mapping; map fSOM_CD4p warped files

## file paths: warped samples
tphe.warped.files <- list.files("./data_modified/TPHE_warping/", recursive = TRUE, full.names = TRUE, 
                           pattern = "._fda_WARPED_.fcs")
tphe.warped.files <- tphe.warped.files[grep("PERFORIN.FOXP3.CD185.CD31.GRZB", dirname(tphe.warped.files))]

## load fSOM
fSOM <- readRDS("./results_R/fSOMs/TPHE/2018-11-13.fSOM_PASS3_CD4p.rds") # needs stored $transform.object

## generate counts only
newfcs.fromfsom.singlefolder(fSOM.object = fSOM,
                             meta.cluster = NULL,
                             FCS.paths = tphe.warped.files,
                             basefolder.name = "./data_modified/TPHE_warping/fSOM_PASS3_CD4p.clusters",
                             fcsfolder.name =  NULL,
                             write.count.files = TRUE,
                             write.fcs.files = FALSE)

## 'target file'
target.file <- list.files("./data_source/TPHE2_with.spillover/NGL091/", full.names = TRUE, pattern = "HD0189")

newfcs.fromfsom.singlefolder(fSOM.object = fSOM,
                             meta.cluster = grep("CD4", fSOM$metaclustering.anno),
                             FCS.paths = target.file,
                             basefolder.name = "./data_modified/TPHE_warping/fSOM_CD3.CD4.CD8.LIVEDEAD_PASS1",
                             fcsfolder.name =  "fSOM_CD4p",
                             write.count.files = TRUE,
                             write.fcs.files = TRUE)