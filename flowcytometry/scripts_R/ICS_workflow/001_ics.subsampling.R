## subsample files for use in building fSOM PASS1; can directly read and cluster on full files but will take much longer
library(flowCore)

## file paths
input.files <- list.files("./data_source/ICS2_with.spillover/", recursive = TRUE, full.names = TRUE, pattern = "*SEB.fcs$")

FCS.subsampled.singlefolder(input.files, 20000, "./data_modified/ICS/SEB_subsampled")