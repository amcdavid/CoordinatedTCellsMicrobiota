## subsample warped files for use in building fSOM; can directly read and cluster on full files but will take much longer

## file paths: warped samples
warped.files <- list.files("./data_modified/TPHE_warping/", recursive = TRUE, full.names = TRUE, pattern = "WARPED_.fcs")
warped.files <- warped.files[grep("CD4.CD8A.CD3.LIVEDEAD_warped", dirname(warped.files))]

FCS.subsampled.singlefolder <- function(FCS.paths, samp.size, newfolder.path.name){
  ##  subsamples FCS file(s); randomly selects rows per file
  for (i in seq(FCS.paths)){
    ##read in and sample raw FCS (linear/non-transformed data)
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

library(flowCore)

FCS.subsampled.singlefolder(warped.files, 10000, "./data_modified/TPHE_warping/CD3.CD4.CD8.LIVEDEAD_warped_SUBSAMPLED")
