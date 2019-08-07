library(flowCore); library(FlowSOM)
## fSOM PASS3 mapping; map fSOM_CD8p warped files

## file paths: warped samples
tphe.warped.files <- list.files("./data_modified/TPHE_warping/EXPERIMENTS/", recursive = TRUE, full.names = TRUE, pattern = "CD8p.SINGLETS_.fcs")
tphe.warped.files <- tphe.warped.files[grep("PERFORIN.CD57.CD28.FOXP3.CD197.CD185.CD31.CD122.GRZB", dirname(tphe.warped.files))]

## load fSOM
fSOM <- readRDS("./results_R/fSOMs/TPHE/2019-01-18.fSOM_PASS3_CD8p.rds") # needs stored $transform.object

## generate counts only
newfcs.fromfsom.singlefolder(fSOM.object = fSOM,
                             meta.cluster = NULL,
                             FCS.paths = tphe.warped.files,
                             basefolder.name = "./data_modified/TPHE_warping/fSOM_PASS3_CD8p.clusters",
                             fcsfolder.name =  NULL,
                             write.count.files = TRUE,
                             write.fcs.files = FALSE)

## 'target file'
target.file <- list.files("./data_modified/TPHE_warping/fSOM_CD8p.SINGLETS_PASS2/", pattern = "HD0189_TPHE_NGL091", full.names = TRUE)

newfcs.fromfsom.singlefolder(fSOM.object = fSOM,
                             meta.cluster = grep("CD4", fSOM$metaclustering.anno),
                             FCS.paths = target.file,
                             basefolder.name = "./data_modified/TPHE_warping/fSOM_PASS3_CD8p.clusters",
                             fcsfolder.name =  NULL,
                             write.count.files = TRUE,
                             write.fcs.files = FALSE)

# ## write new parameters -- nodes and meta clusters -- to fcs files; not working...can't read them back in (but flowjo can...?)
# add.fsom.nodes.meta_FCS <- function(fsom.obj, fcs.paths, fcsfolder.name) {
#   
#   fsom <- fsom.obj
#   trans.obj <- fsom$transform.object
#   
#   for(i in seq(fcs.paths)) {
#     
#     fcs.raw <- read.flowSet(fcs.paths[i], transformation = FALSE)
#     
#     fcs.tmp <- compensate(fcs.raw, fcs.raw@description$`$SPILLOVER`)
#     
#     #fcs.tmp <- transform.set(fcs.tmp)
#     
#     for(j in seq(length(trans.obj))){
#       print(paste0("transforming ", names(trans.obj)[j]))
#       fcs.tmp <- transform(fcs.tmp, trans.obj[[j]]$asinh.translist)
#       # fcs.tmp <- fsApply(fcs.tmp, function(i){
#       #   i@parameters@data$minRange[grep(names(trans.obj)[j], i@parameters@data$name)] <- min(exprs(i)[, names(trans.obj)[j]])
#       #   i@parameters@data$maxRange[grep(names(trans.obj)[j], i@parameters@data$name)] <- max(exprs(i)[, names(trans.obj)[j]])
#       #   i
#       # })
#     }
# 
#     message(paste("Mapping", fcs.tmp[[1]]@description$`$FIL`, "to", "fSOM.object...",  i, "of ", length(fcs.paths), sep = " "))
#     
#     fsom.tmp <- NewData(fsom, fcs.tmp)
#     
#     fcs.raw <- fsApply(fcs.raw, function(i){
#       i <- fr_append_cols(i, matrix(fsom.tmp$map$mapping[, 1], dimnames = list(NULL, "node_PASS3")))
#       i <- fr_append_cols(i, matrix(as.integer(fsom.tmp$metaclustering[fsom.tmp$map$mapping[, 1]]), dimnames = list(NULL, "meta.cluster_PASS3")))
#     })
#     
#     outFile <- file.path(fcsfolder.name, fcs.raw[[1]]@description$`$FIL`)
#     write.FCS(fcs.raw[[1]], outFile)
#   }
# }
# 
# add.fsom.nodes.meta_FCS(fsom.obj = readRDS("./results_R/fSOMs/TPHE/2019-01-18.fSOM_PASS3_CD8p.rds"),
#                         fcs.paths = tphe.warped.files,
#                         fcsfolder.name = "./data_modified/TPHE_warping/fSOM_PASS3_CD8p.clusters/")

add.fsom.nodes.meta_FCS <- function(fsom.obj, fcs.paths) {
  
  frames <- list()
  
  fsom <- fsom.obj
  trans.obj <- fsom$transform.object
  
  for(i in seq(fcs.paths)) {
    
    fcs.tmp <- read.FCS(fcs.paths[i], transformation = FALSE)
    
    fcs.tmp <- compensate(fcs.tmp, fcs.tmp@description$`$SPILLOVER`)
    
    for(j in seq(length(trans.obj))){
      #print(paste0("transforming ", names(trans.obj)[j]))
      fcs.tmp <- transform(fcs.tmp, trans.obj[[j]]$asinh.translist)
      # fcs.tmp <- fsApply(fcs.tmp, function(i){
      #   i@parameters@data$minRange[grep(names(trans.obj)[j], i@parameters@data$name)] <- min(exprs(i)[, names(trans.obj)[j]])
      #   i@parameters@data$maxRange[grep(names(trans.obj)[j], i@parameters@data$name)] <- max(exprs(i)[, names(trans.obj)[j]])
      #   i
      # })
    }
    
    message(paste("Mapping", fcs.tmp@description$`$FIL`, "to", "fSOM.object...",  i, "of ", length(fcs.paths), sep = " "))
    
    fsom.tmp <- NewData(fsom, fcs.tmp)
    
    new_cols <- matrix(c(fsom.tmp$map$mapping[, 1], fsom.tmp$metaclustering[fsom.tmp$map$mapping[, 1]]),
                       ncol = 2,
                       dimnames = list(NULL, c("node", "meta.cluster")))
    
    fcs.tmp <- fr_append_cols(fcs.tmp, new_cols)
    
    frames[[i]] <- as.data.frame(exprs(fcs.tmp))
    names(frames)[i] <- fcs.tmp@description$`$SRC`
  }
  frames
}

tphe.cd8.fsom_frames <- add.fsom.nodes.meta_FCS(fsom.obj = readRDS("./results_R/fSOMs/TPHE/2019-01-18.fSOM_PASS3_CD8p.rds"),
                                                fcs.paths = tphe.warped.files[-grep("HD", tphe.warped.files)])

tphe.cd8.fsom_conc <- as.data.frame(data.table::rbindlist(tphe.cd8.fsom_frames))
markers <- readRDS("./results_R/fSOMs/TPHE/2019-01-18.fSOM_PASS3_CD8p.rds")$markers
colnames(tphe.cd8.fsom_conc)[1:length(markers)] <- readRDS("./results_R/fSOMs/TPHE/2019-01-18.fSOM_PASS3_CD8p.rds")$markers

library(shiny); library(ggplot2); library(viridis)
fcs.explorer(tphe.cd8.fsom_conc, "CD185", "CD57", 30000)
