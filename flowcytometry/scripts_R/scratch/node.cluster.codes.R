cols <- matrix(c(as.integer(fSOM.tmp$map$mapping[, 1]), as.integer(fSOM.tmp$metaclustering[fSOM.tmp$map$mapping[ ,1]])),
               ncol = 2,
               dimnames = list(NULL, c("nodes", "clusters")))

cols[, 1]

node.code <- sapply(cols[, 1], function(x) {
  if(nchar(x) <= 1) {
    paste0(0, 0, x)
  } else if(nchar(x) <= 2) {
    paste0(0, x)
  } else {
    paste0(x)
  }
})
nodes.1E2 <- as.integer(sapply(node.code, function(x) strsplit(x, "")[[1]][1]))
nodes.1E1 <- as.integer(sapply(node.code, function(x) strsplit(x, "")[[1]][2]))
nodes.1E0 <- as.integer(sapply(node.code, function(x) strsplit(x, "")[[1]][3]))


cluster.code <- sapply(cols[, 2], function(x) {
  if(nchar(x) <= 1) {
    paste0(0, x)
  } else if(nchar(x) <= 2) {
    paste0(x)
  }
})
clusters.1E1 <- as.integer(sapply(cluster.code, function(x) strsplit(x, "")[[1]][1]))
clusters.1E0 <- as.integer(sapply(cluster.code, function(x) strsplit(x, "")[[1]][2]))


cols <- matrix(c(as.integer(fSOM.tmp$map$mapping[, 1]), 
                 as.integer(fSOM.tmp$metaclustering[fSOM.tmp$map$mapping[ ,1]]),
                 nodes.1E2,
                 nodes.1E1,
                 nodes.1E0,
                 clusters.1E1,
                 clusters.1E0),
               ncol = 7,
               dimnames = list(NULL, 
                               c("nodes", 
                                 "clusters",
                                 "nodes.1E2",
                                 "nodes.1E1",
                                 "nodes.1E0",
                                 "clusters.1E1",
                                 "clusters.1E0")))

fcs.tmp.newcols <- fr_append_cols(fcs.raw, cols)

write.FCS(fcs.tmp.newcols, "test_tmp.fcs")
