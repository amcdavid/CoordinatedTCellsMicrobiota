## subsample PERFORIN, FOXP3, CD185, CD31, GRZB warped files for QC

## file paths: warped samples
warped.files <- list.files("./data_modified/TPHE_warping/", recursive = TRUE, full.names = TRUE, 
                           pattern = "._fda_WARPED_.fcs")
warped.files <- warped.files[grep("PERFORIN.FOXP3.CD185.CD31.GRZB", dirname(warped.files))]

FCS.subsampled.singlefolder <- function(FCS.paths, samp.size, newfolder.path.name){
  ##  Randomly subsamples events from FCS file(s)
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

FCS.subsampled.singlefolder(warped.files, 10000, "./data_modified/TPHE_warping/PERFORIN.FOXP3.CD185.CD31.GRZB_warped_CD4p_SUBSAMPLED")
