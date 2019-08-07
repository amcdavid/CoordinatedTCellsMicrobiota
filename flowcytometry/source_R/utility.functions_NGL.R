##create batch list (experiment #); dataset specific
batch.list <- function(fcs.file.paths){
  experiments <- unique(lapply(seq(length(fcs.file.paths)), FUN = function(x) {
    (strsplit(fcs.file.paths, "/")[[x]][grep("NGL*", strsplit(fcs.file.paths, "/")[[x]])])}))
  batches <- lapply(seq(unique(experiments)), FUN = function(x) {fcs.file.paths[grep(experiments[x], fcs.file.paths)]})
  names(batches) <- unique(experiments)
  batches
}

batch.list.edit <- function(fcs.file.paths){
  experiments <- unique(lapply(seq(length(fcs.file.paths)), FUN = function(x) {
    (strsplit(fcs.file.paths, "/")[[x]][grep("NGL*", strsplit(fcs.file.paths, "/")[[x]])][1])}))
  batches <- lapply(seq(unique(experiments)), FUN = function(x) {fcs.file.paths[grep(experiments[x], fcs.file.paths)]})
  names(batches) <- unique(experiments)
  batches
}

folder.list <- function(fcs.file.paths){
  unique.paths <- unique(dirname(fcs.file.paths))
  batches <- lapply(seq(unique.paths), FUN = function(x) {list.files(path = unique.paths[x], full.names = TRUE, pattern = ".fcs")})
  names(batches) <- unique.paths
  batches
}

##flowViz function for assessing 'quality' of median values obtained using median.QC function; needs 'median.dat'
density.QC.fromfiles <- function(fcs.file.paths, channel, max.files, decreasing = TRUE) {
  
  unique.sub <- data.frame("unique" = sub("_fSOM.*", "", basename(fcs.file.paths)))
  
  median.dat.sub <- merge(median.dat, unique.sub, by = "unique")
  
  median.files <- fcs.file.paths[grep(paste(median.dat.sub$unique[order(median.dat.sub[, channel], decreasing = decreasing)[1:max.files]], collapse = "|"), fcs.file.paths)]
  
  input.flowset <- fset.transform(fset.compensate(read.flowSet(median.files)))
  channel.name <- colnames(input.flowset)[grep(channel, input.flowset[[1]]@parameters@data$desc)]
  print(channel.name)
  
  densityplot(~ ., input.flowset,  main = "original",   filter= curv1Filter(channel.name), channel = channel.name)
  
}

density.QC.fromflowsets <- function(input.flowset, channel) {
  
  channel.name <- colnames(input.flowset)[grep(channel, input.flowset[[1]]@parameters@data$desc)]
  
  densityplot(~ ., input.flowset,  main = "original",   filter= curv1Filter(channel.name), channel = channel.name)
  
}




