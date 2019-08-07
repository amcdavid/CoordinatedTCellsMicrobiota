library(flowCore); library(FlowSOM)

## file paths
tphe.warped.files <- list.files("./data_modified/TPHE_warping", recursive = TRUE, full.names = TRUE, pattern = "WARPED_.fcs")
tphe.warped.files <- tphe.warped.files[grep("CD4.CD8A.CD3.LIVEDEAD_warped", dirname(tphe.warped.files))]

## load fSOM
fSOM <- readRDS("./results_R/fSOMs/TPHE/2018-11-11.fSOM_CD3.CD4.CD8.LIVEDEAD.PASS1.rds") # needs stored $transform.object

# ## write new parameters -- nodes and meta clusters -- to fcs files; can't get this when trying to read back in the appened fcs file
# add.fsom.nodes.meta_FCS <- function(fsom.obj, fcs.paths, fcsfolder.name) {
# 
#   fsom <- fsom.obj
#   trans.obj <- fsom$transform.object
# 
#   for(i in seq(fcs.paths)) {
# 
#     fcs.raw <- read.FCS(fcs.paths[i], transformation = FALSE)
# 
#     fcs.tmp <- compensate(fcs.raw, fcs.raw@description$`$SPILLOVER`)
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
#     message(paste("Mapping", fcs.tmp@description$`$FIL`, "to", "fSOM.object...",  i, "of ", length(fcs.paths), sep = " "))
# 
#     fsom.tmp <- NewData(fsom, fcs.tmp)
#     
#     new_cols <- matrix(c(fsom.tmp$map$mapping[, 1], fsom.tmp$metaclustering[fsom.tmp$map$mapping[, 1]]), 
#                        ncol = 2, 
#                        dimnames = list(NULL, c("nodes", "clusters")))
#    
#     fcs.raw <- fr_append_cols(fcs.raw, new_cols)
#  
#     outFile <- file.path(fcsfolder.name, fcs.raw@description$`$FIL`)
#     write.FCS(fcs.raw, outFile)
#   }
# }
# 
# add.fsom.nodes.meta_FCS(fsom.obj = readRDS("./results_R/fSOMs/TPHE/2018-11-11.fSOM_CD3.CD4.CD8.LIVEDEAD.PASS1.rds"),
#                         fcs.paths = tphe.warped.files[1],
#                         fcsfolder.name = "./data_modified/TPHE_warping/fSOM_PASS3_CD8p.clusters/")

## write new parameters -- nodes and meta clusters -- to fcs files; can't get this when trying to read back in the appened fcs file
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
                       dimnames = list(NULL, c("node", "cluster")))

    fcs.tmp <- fr_append_cols(fcs.tmp, new_cols)
    
    frames[[i]] <- as.data.frame(exprs(fcs.tmp))
  }
  frames
}

frames <- add.fsom.nodes.meta_FCS(readRDS("./results_R/fSOMs/TPHE/2018-11-11.fSOM_CD3.CD4.CD8.LIVEDEAD.PASS1.rds"),
                                  tphe.warped.files[1:3])
