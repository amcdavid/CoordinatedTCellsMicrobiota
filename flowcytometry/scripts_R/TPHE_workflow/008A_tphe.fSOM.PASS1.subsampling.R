## subsample fSOM PASS 1 CD4+ for QC
library(flowCore)

## file paths: warped samples
fsom.files <- list.files("./data_modified/TPHE_warping/fSOM_CD3.CD4.CD8.LIVEDEAD_PASS1", 
                         recursive = FALSE, full.names = TRUE, pattern = ".fcs")

FCS.subsampled.singlefolder <- function(FCS.paths, samp.size, newfolder.path.name){
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

FCS.subsampled.singlefolder(fsom.files, 10000, "./data_modified/TPHE_warping/fSOM_CD4p_CD8p_PASS1_SUBSAMPLED")
