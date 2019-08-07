## fSOM PASS2 mapping

## file paths
tphe.warped.files <- list.files("./data_modified/TPHE_warping/fSOM_CD3.CD4.CD8.LIVEDEAD_PASS1/scatter_warped/CD4/", 
                                recursive = FALSE, full.names = TRUE, pattern = ".fcs")

## load fSOM
fSOM <- readRDS("./results_R/fSOMs/TPHE/2018-11-11.fSOM_SINGLETS.PASS2.rds") # needs stored $transform.object

##
newfcs.fromfsom.scatter.singlefolder (fSOM.object = fSOM,
                                      meta.cluster = grep("singlets", fSOM$metaclustering.anno),
                                      FCS.paths = tphe.warped.files,
                                      basefolder.name = "./data_modified/TPHE_warping/fSOM_CD4p.SINGLETS_PASS2",
                                      fcsfolder.name =  "fSOM_CD4p.SINGLETS",
                                      write.count.files = TRUE,
                                      write.fcs.files = TRUE)

## 'target file'
target.file <- list.files("./data_modified/TPHE_warping/fSOM_CD3.CD4.CD8.LIVEDEAD_PASS1/", full.names = TRUE, pattern = "HD0189_TPHE_NGL091")

newfcs.fromfsom.scatter.singlefolder (fSOM.object = fSOM,
                                      meta.cluster = grep("singlets", fSOM$metaclustering.anno),
                                      FCS.paths = target.file,
                                      basefolder.name = "./data_modified/TPHE_warping/fSOM_CD4p.SINGLETS_PASS2",
                                      fcsfolder.name =  "fSOM_CD4p.SINGLETS",
                                      write.count.files = TRUE,
                                      write.fcs.files = TRUE)

####

## file paths
tphe.warped.files <- list.files("./data_modified/TPHE_warping/fSOM_CD3.CD4.CD8.LIVEDEAD_PASS1/scatter_warped/CD8/", 
                                recursive = FALSE, full.names = TRUE, pattern = ".fcs")

## load fSOM
fSOM <- readRDS("./results_R/fSOMs/TPHE/2019-01-16.fSOM_CD8p_SINGLETS.PASS2.rds") # needs stored $transform.object

##
newfcs.fromfsom.scatter.singlefolder (fSOM.object = fSOM,
                                      meta.cluster = grep("singlets", fSOM$metaclustering.anno),
                                      FCS.paths = tphe.warped.files,
                                      basefolder.name = "./data_modified/TPHE_warping/fSOM_CD8p.SINGLETS_PASS2",
                                      fcsfolder.name =  "fSOM_CD8p.SINGLETS",
                                      write.count.files = TRUE,
                                      write.fcs.files = TRUE)

## 'target file'
target.file <- list.files("./data_modified/TPHE_warping/fSOM_CD3.CD4.CD8.LIVEDEAD_PASS1/", full.names = TRUE, pattern = "HD0189_TPHE_NGL091 TPHE.fcs_fSOM_CD8p_.fcs")

newfcs.fromfsom.scatter.singlefolder (fSOM.object = fSOM,
                                      meta.cluster = grep("singlets", fSOM$metaclustering.anno),
                                      FCS.paths = target.file,
                                      basefolder.name = "./data_modified/TPHE_warping/fSOM_CD8p.SINGLETS_PASS2",
                                      fcsfolder.name =  "fSOM_CD8p.SINGLETS",
                                      write.count.files = TRUE,
                                      write.fcs.files = TRUE)
