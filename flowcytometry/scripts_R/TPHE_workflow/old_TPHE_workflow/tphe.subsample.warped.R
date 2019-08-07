## source flow-related libraries and functions
source("./source_R/flow.src.R", chdir = TRUE)

## file paths: warped samples
warped.files <- list.files("./data_modified/TPHE_warping/", recursive = TRUE, full.names = TRUE, pattern = "WARPED_.fcs")
fsom.files <- list.files("./data_modified/TPHE_warping/fSOM_CD3.CD4.CD8.LIVEDEAD_PASS1/", full.names = TRUE, pattern = ".fcs")

FCS.subsampled.edit <- function(FCS.paths, samp.size, newfolder.path.name){
  ##  Randomly subsamples FCS file(s)
  for (i in seq(FCS.paths)){
    ##read in and sample raw FCS (linear/non-transformed data); use emptyValue to avoild delimiter keyword error (?)
    FCS.raw <- read.FCS(FCS.paths[i], transformation = FALSE)

    print(paste("Subsampling ", FCS.raw@description$`$FIL`, i, "of ", length(FCS.paths)))
    
    samp <- sample(1:nrow(exprs(FCS.raw)), samp.size, replace = TRUE)
    exprs(FCS.raw) <- exprs(FCS.raw)[samp, ]
    FCS.raw@description$`$FIL` <- sub("_.fcs", "_subsampled_.fcs", FCS.raw@description$`$FIL`)

    ##write FCS to directory
    if (dir.exists(newfolder.path.name)) {
      outFile <- file.path(newfolder.path.name, FCS.raw@description$`$FIL`)
    } else {
      dir.create(newfolder.path.name)
      outFile <- file.path(newfolder.path.name, FCS.raw@description$`$FIL`)
    }
    write.FCS(FCS.raw, outFile)
  }
}

FCS.subsampled.edit(fsom.files, 10000, "./data_modified/TPHE_warping/fSOM_CD3.CD4.CD8.LIVEDEAD_PASS1_SUBSAMPLED")

##
subsampled.files <- list.files("./data_modified/TPHE_warping/fSOM_CD3.CD4.CD8.LIVEDEAD_PASS1_SUBSAMPLED/", full.names = TRUE)
trans.obj <- readRDS("./results_R/FCS.transforms/TPHE/asinh.transforms.TPHE.rds")
trans.obj <- trans.obj[grep('Violet G|Violet E|Violet A|Green D', names(trans.obj))]

med.qc <- median.QC.edit(subsampled.files, trans.obj)

ggplot(med.qc, aes(x = random, y = CD3)) + geom_point()
ggplot(med.qc, aes(x = random, y = CD4)) + geom_point()
ggplot(med.qc, aes(x = random, y = CD8A)) + geom_point()





## file paths: warped samples
warped.files <- list.files("./data_modified/TPHE_warping", recursive = TRUE, full.names = TRUE, pattern = "fda_WARPED_FSC.AHW_SSC.AHW_VG_VE_VA_GD_fSOM_CD4pos_.fcs")

FCS.subsampled.edit <- function(FCS.paths, samp.size, newfolder.name, post.name){
  ##  Randomly subsamples FCS file(s)
  for (i in seq(FCS.paths)){
    ##read in and sample raw FCS (linear/non-transformed data); use emptyValue to avoild delimiter keyword error (?)
    FCS.raw <- read.FCS(FCS.paths[i])
    #print(paste0("subsampling... ", FCS.raw@description$`$FIL`))
    print(paste("Subsampling ", FCS.raw@description$`$FIL`, i, "of ", length(FCS.paths)))
    if (nrow(exprs(FCS.raw)) > samp.size) {
      samp <- sample(1:nrow(exprs(FCS.raw)), samp.size, replace = TRUE)
      exprs(FCS.raw) <- exprs(FCS.raw)[samp, ]
      FCS.raw@description$`$FIL` <- paste(sub("_.fcs", "", FCS.raw@description$`$FIL`), post.name, sep = "_")
    } else {
      FCS.raw@description$`$FIL` <- paste(sub("_.fcs", "", FCS.raw@description$`$FIL`), post.name, sep = "_")
    }
    ##write FCS to directory
    if (dir.exists(file.path(dirname(FCS.paths[i]), newfolder.name))) {
      outFile <- file.path(dirname(FCS.paths[i]), newfolder.name, FCS.raw@description$`$FIL`)
    } else {
      dir.create(file.path(dirname(FCS.paths[i]), newfolder.name))
      outFile <- file.path(dirname(FCS.paths[i]), newfolder.name, FCS.raw@description$`$FIL`)
    }
    write.FCS(FCS.raw, outFile)
  }
}

FCS.subsampled.edit(warped.files[1], 10000, "subsampled", "subsampled.fcs")