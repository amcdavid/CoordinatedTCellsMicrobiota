tphe.warped.files <- list.files("./data_modified/TPHE_warping", recursive = TRUE, full.names = TRUE, pattern = "FSC.AHW_SSC.AHW_VG_VE_VA_GD")

name.fix <- function(fcs.file.paths) {
  for(i in seq(fcs.file.paths)){
    fcs.tmp <- read.FCS(fcs.file.paths[i], transformation = FALSE)
    fcs.tmp@description$`$FIL` <- paste(sub("_.fcs", "", fcs.tmp@description$`$FIL`), "fSOM_CD4pos", ".fcs", sep = "_")
    
    write.FCS(x = fcs.tmp, filename = file.path(dirname(fcs.file.paths[i]), fcs.tmp@description$`$FIL`))
  }
}