library(flowCore)
library(FlowSOM)
library(shiny)
library(ggplot2)
library(viridis)
library(pheatmap)

fSOM <- readRDS("./analyses_R/FSOM_analysis/03_CD4_fSOM.CD4Cytokines.rds")

tmp.fcs <- read.FCS("./data_source/ICS2_with.spillover/NGL092/fSOM.CD3LIVE/fSOM.CD4pCD69p/RPRC0540103_1_SEB_fSOM.CD4pCD69p_.fcs", transformation = FALSE)
tmp.fcs <- compensate(tmp.fcs, (keyword(tmp.fcs)[grep("SPILL", names(keyword(tmp.fcs)))][[1]]))
tmp.fcs <- transform(tmp.fcs, elgcl)

fSOM.sub <- NewData(fSOM, tmp.fcs)

dat <- fsom.dat(fSOM.sub)

fcs.explorer(dat, "CD69", "IL8", 30000)


fSOM <- readRDS("./analyses_R/FSOM_analysis/03B_CD8_fSOM.CD8Cytokines.RDS")

tmp.fcs <- read.FCS("./data_source/ICS2_with.spillover/NGL092/fSOM.CD3LIVE/fSOM.CD8pCD69p/RPRC0530106_1_SEB_fSOM.CD8pCD69p_.fcs", transformation = FALSE)
tmp.fcs <- compensate(tmp.fcs, (keyword(tmp.fcs)[grep("SPILL", names(keyword(tmp.fcs)))][[1]]))
tmp.fcs <- transform(tmp.fcs, elgcl)

fSOM.sub <- NewData(fSOM, tmp.fcs)

fcs.explorer(fsom.dat(fSOM.sub), "CD69", "IL8", 30000)


fSOM <- readRDS("./results_R/fSOMs/TPHE/2018-11-14.fSOM_PASS3_CD4p_rerun_ALLDATA.rds")

fcs.explorer(fsom.dat(fSOM), "CD31", "CD197", 30000)


fSOM <- readRDS("./results_R/fSOMs/TPHE/2019-01-18.fSOM_PASS3_CD8p.rds")

input.files <- list.files("./data_modified/TPHE_warping/EXPERIMENTS/", full.names = TRUE, recursive = TRUE, pattern = "CD8p.SINGLETS")
input.files <- input.files[grep("CD122.GRZB", dirname(input.files))]
input.files <- input.files[-grep("HD", input.files)]
tmp.fcs <- read.flowSet(input.files[grep("106_19_", input.files)], transformation = FALSE)
tmp.fcs <- fset.compensate(tmp.fcs)
trans.obj <- fSOM$transform.object
tmp.fcs <- transform.set(tmp.fcs)

fSOM.sub <- NewData(fSOM, tmp.fcs[[1]])

fcs.explorer(fsom.dat(fSOM.sub), "CD31", "CD197", 30000)

CD31p.clusters <- sort(unique(subset(fsom.dat(fSOM.sub), fsom.dat(fSOM.sub)$CD31 > 1)$MCluster))

