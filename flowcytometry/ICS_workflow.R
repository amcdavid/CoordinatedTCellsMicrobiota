# STEP 1; # subsample source files
file.edit("./scripts_R/ICS_workflow/001_ics.subsampling.R")

# STEP 2; # FlowSOM PASS1 using subsampled files
file.edit("./scripts_R/ICS_workflow/002_ics.fSOM.PASS1_FSC.A_SSC.A_FSC.W_SSC.W_LIVE.DEAD_CD14_CD3.R")

# STEP 3; # FlowSOM mapping PASS1
file.edit("./scripts_R/ICS_workflow/003_ics.fSOM.PASS1.mapping.R")

# STEP 4; # subsample PASS1: CD3+ files
file.edit("./scripts_R/ICS_workflow/004_ics.CD3p_subsampling.R")

# STEP 5; # FlowSOM PASS2 using subsampled files
file.edit("./scripts_R/ICS_workflow/005_ics.fSOM.PASS2.R")
