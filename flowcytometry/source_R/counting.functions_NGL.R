generate.counts <- function(count.paths){
  ##read .csv count files and create named list
  counts <- lapply(seq(count.paths), function(i) fread(count.paths[i]))
  names(counts) <- lapply(seq(count.paths), function(i) gsub("_counts.csv","",basename(count.paths[i])))
  ##
  counts <- as.data.frame(t(mapply(`[[`,counts, 2)))
  colnames(counts) <- sub("V","Meta.Cluster_",colnames(counts))
  counts
}