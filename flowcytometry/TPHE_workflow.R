# STEP 1; # define and store asinh transformation parameters for transformation of .fcs data
file.edit("./scripts_R/TPHE_workflow/001_tphe.asinh.transforms.R")

# STEP 2; # define and store data warping parameters - initial fluors for FlowSOM Pass 1
file.edit("./scripts_R/TPHE_workflow/002_tphe.warp.parameters.VG.VE.VA.GD.R")

# STEP 3; # warp source files
file.edit("./scripts_R/TPHE_workflow/003_tphe.warping.CD3.CD4.CD8.LIVEDEAD.R")

# STEP 4; # spot-check/fix individual samples/channels
file.edit("./scripts_R/TPHE_workflow/004_tphe.warping.CD3.CD4.CD8.LIVEDEAD.spotcheck.R")

# STEP 5A; # subsample warped files
file.edit("./scripts_R/TPHE_workflow/005A_tphe.subsampling.R")
# STEP 5B; # QC
file.edit("./scripts_R/TPHE_workflow/005B_tphe.warping.CD3.CD4.CD8.LIVEDEAD.QC.R")

# STEP 6; # FlowSOM using subsampled files
file.edit("./scripts_R/TPHE_workflow/006_tphe.fSOM.PASS1.CD3.CD4.CD8.LIVEDEAD.R")

# STEP 7; # FlowSOM mapping
file.edit("./scripts_R/TPHE_workflow/007_tphe.fSOM.PASS1.CD3.CD4.CD8.LIVEDEAD.mapping.R")

# STEP 8A; # subsample FlowSOM mapped files
file.edit("./scripts_R/TPHE_workflow/008A_tphe.fSOM.PASS1.subsampling.R")
# STEP 8B; # QC
file.edit("./scripts_R/TPHE_workflow/008B_tphe.fSOM.PASS1.QC.R")

# STEP 9; # define and store data warping parameters - scatter
file.edit("./scripts_R/TPHE_workflow/009_tphe.warp.parameters.fSOM_PASS1.scatter.R")

# STEP 10; # scatter warp fsom files
file.edit("./scripts_R/TPHE_workflow/010_tphe.scatter.warping.fSOM_PASS1.R")

# STEP 11A; # subsample scatter warped files
file.edit("./scripts_R/TPHE_workflow/011A_tphe.subsample.scatter.warped.R")
# STEP 11B; # QC
file.edit("./scripts_R/TPHE_workflow/011B_tphe.scatter.warping.fSOM_PASS1.QC.R")

# STEP 12; # FlowSOM using subsampled files - scatter parameters
file.edit("./scripts_R/TPHE_workflow/012_tphe.fSOM.PASS2.scatter.R")

# STEP 13; # FlowSOM mapping
file.edit("./scripts_R/TPHE_workflow/013_tphe.fSOM.PASS2.SINGLETS.mapping.R")

# STEP 14A; # subsample FlowSOM mapped files
file.edit("./scripts_R/TPHE_workflow/014A_tphe.fSOM.PASS2.SINGLETS.subsampling.R")
# STEP 14B; # QC
file.edit("./scripts_R/TPHE_workflow/014B_tphe.fSOM.PASS2.SINGLETS.QC.R")

# STEP 15A; # define and store data warping parameters - fluors
file.edit("./scripts_R/TPHE_workflow/015A_tphe.warp.parameters.fSOM_CD4p.SINGLETS.R")
# STEP 15B; # define and store data warping parameters - fluors
file.edit("./scripts_R/TPHE_workflow/015B_tphe.warp.parameters.fSOM_CD8p.SINGLETS.R")

# STEP 16A; # warp CD4 fsom files; PERFORIN, FOXP3, CD185, CD31, GRZB
file.edit("./scripts_R/TPHE_workflow/016A_tphe.warping.fSOM_CD4p.SINGLETS.R")
# STEP 16B; # warp CD8 fsom files; PERFORIN.CD57.CD28.FOXP3.CD197.CD185.CD31
file.edit("./scripts_R/TPHE_workflow/016B_tphe.warping.fSOM_CD8p.SINGLETS.R")
# STEP 16C; # warp CD8 fsom files; CD122.GRZB
file.edit("./scripts_R/TPHE_workflow/016C_tphe.warping.fSOM_CD8p.SINGLETS.PASS2.R")

# STEP 17A; # subsample PERFORIN, FOXP3, CD185, CD31, GRZB CD4 warped files
file.edit("./scripts_R/TPHE_workflow/017A_tphe.subsample.PERFORIN.FOXP3.CD185.CD31.GRZB_CD4.warped.R")
# STEP 17B; # QC
file.edit("./scripts_R/TPHE_workflow/017B_tphe.PERFORIN.FOXP3.CD185.CD31.GRZB_CD4.warping.QC.R")

# STEP 18A; # FlowSOM using CD4 subsampled files
file.edit("./scripts_R/TPHE_workflow/018A_tphe.CD4.fSOM.PASS3.R")
# STEP 18B; # FlowSOM using CD8 warped files
file.edit("./scripts_R/TPHE_workflow/018B_tphe.CD8.fSOM.PASS3.R")

# STEP 19A; # FlowSOM mapping CD4
file.edit("./scripts_R/TPHE_workflow/019A_tphe.PASS3_CD4p.mapping.R")
# STEP 19B; # FlowSOM mapping CD8
file.edit("./scripts_R/TPHE_workflow/019B_tphe.PASS3_CD8p.mapping.R")

# STEP 20A; # FlowSOM cluster counts/QC
file.edit("./scripts_R/TPHE_workflow/020A_tphe.PASS3_CD4p.counts_QC.R")
# STEP 20B; # FlowSOM cluster counts/QC
file.edit("./scripts_R/TPHE_workflow/020B_tphe.PASS3_CD8p.counts_QC.R")

# STEP 21; # FlowSOM using subsampled files
file.edit("./scripts_R/TPHE_workflow/021_tphe.CD4.fSOM.PASS3_rerun.R")

# STEP 22; # FlowSOM mapping
file.edit("./scripts_R/TPHE_workflow/022_tphe.PASS3_CD4p.mapping_rerun.R")

# STEP 23; # FlowSOM cluster counts/QC
file.edit("./scripts_R/TPHE_workflow/023_tphe.PASS3_CD4p.counts_QC_rerun.R")

# STEP 24; # FlowSOM + tSNE
file.edit("./scripts_R/TPHE_workflow/024_tphe.PASS3_CD4p.fSOM_tSNE.R")
