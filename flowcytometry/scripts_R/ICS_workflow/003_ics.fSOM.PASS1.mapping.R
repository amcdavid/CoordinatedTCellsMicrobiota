## fSOM PASS1 mapping
library(flowCore)
library(FlowSOM)

## file paths
files.withoutHW <- paste(sub(".fcs", "", list.files("./data_source/ICS2_with.spillover/NGL046/", pattern = "*SEB.fcs$")), collapse = "|")# NGL046 has missing parameters
input.files <- list.files("./data_source/ICS2_with.spillover/", recursive = TRUE, full.names = TRUE, pattern = "*SEB.fcs$")
input.files <- input.files[-grep(files.withoutHW, input.files)]

## load fSOM
fSOM <- readRDS("./results_R/fSOMs/ICS/ICS_fSOM_PASS1_FSC.A_SSC.A_FSC.W_SSC.W_LIVE.DEAD_CD14_CD3_CD69_2019-01-23.rds") # needs stored $transform.object

##
fsom.mapping.elgcl(fSOM.object = fSOM,
                   meta.cluster = grep("CD3", fSOM$metaclustering.anno),
                   FCS.paths =  input.files,
                   trim = TRUE,
                   basefolder.name = "./data_modified/ICS/fSOM_PASS1",
                   fcsfolder.name = "fSOM_CD3p",
                   write.count.files = TRUE,
                   write.fcs.files = TRUE)
