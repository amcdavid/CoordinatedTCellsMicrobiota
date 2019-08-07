## subsample files for use in building fSOM PASS2; can directly read and cluster on full files but will take much longer
library(flowCore)

## file paths
input.files <- list.files("./data_modified/ICS/fSOM_PASS1/", recursive = TRUE, full.names = TRUE, pattern = "fSOM_CD3p.fcs$")

FCS.subsampled.singlefolder(input.files, 20000, "./data_modified/ICS/fSOM_PASS1_CD3p_subsampled")
