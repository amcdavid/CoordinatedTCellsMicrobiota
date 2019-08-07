folder.names <- paste0("NGL0", seq(47, 91, 2)) ; folder.names[-grep("NGL063", folder.names)]

for(i in seq(folder.names)) {
  if (!dir.exists(file.path("./data_modified/TPHE_warping", folder.names[i]))) { 
    dir.create(file.path(file.path("./data_modified/TPHE_warping", folder.names[i])))
  }
}

for(i in seq(folder.names)) {
  if (!dir.exists(file.path("./data_modified/TPHE_warping", folder.names[i], "VG_VE_VA_GD"))) { 
    dir.create(file.path(file.path("./data_modified/TPHE_warping", folder.names[i], "VG_VE_VA_GD")))
  }
}

for(i in seq(folder.names)) {
  if (!dir.exists(file.path("./data_modified/TPHE_warping", folder.names[i], "VG_VE_VA_GD/fSOM_CD3.CD4.CD8.LIVEDEAD/fSOM_CD4pos/FSC.AHW_SSC.AHW_VG_VE_VA_GD"))) { 
    dir.create(file.path(file.path("./data_modified/TPHE_warping", folder.names[i], "VG_VE_VA_GD/fSOM_CD3.CD4.CD8.LIVEDEAD/fSOM_CD4pos/FSC.AHW_SSC.AHW_VG_VE_VA_GD")))
  }
}

for(i in seq(folder.names)) {
  if (!dir.exists(file.path("./data_modified/TPHE_warping", folder.names[i], "fSOM_CD4pos_warped"))) { 
    dir.create(file.path(file.path("./data_modified/TPHE_warping", folder.names[i], "fSOM_CD4pos_warped")))
  }
}

for(i in seq(folder.names)) {
  if (!dir.exists(file.path("./data_modified/TPHE_warping", folder.names[i], "CD4.CD8A.CD3.LIVEDEAD_warped"))) { 
    dir.create(file.path(file.path("./data_modified/TPHE_warping", folder.names[i], "CD4.CD8A.CD3.LIVEDEAD_warped")))
  }
}

for(i in seq(folder.names)) {
  if (!dir.exists(file.path("./data_modified/TPHE_warping", folder.names[i], "PERFORIN.FOXP3.CD185.CD31.GRZB_warped_CD4p"))) { 
    dir.create(file.path(file.path("./data_modified/TPHE_warping", folder.names[i], "PERFORIN.FOXP3.CD185.CD31.GRZB_warped_CD4p")))
  }
}
