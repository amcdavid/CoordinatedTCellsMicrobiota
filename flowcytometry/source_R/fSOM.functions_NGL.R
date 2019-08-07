savefsom <- function(fSOM.object,fSOM.name){
  tmp <- fSOM.object
  tmp$data <- head(tmp$data)
  tmp$map$mapping <- head(tmp$map$mapping)
  saveRDS(tmp, file = fSOM.name)
}

fsom.mapping.elgcl <- function(fSOM.object, meta.cluster = NULL, FCS.paths, trim = TRUE, basefolder.name, fcsfolder.name, write.count.files = FALSE, write.fcs.files = FALSE) {
  
  if(!dir.exists(file.path(basefolder.name))) {
    dir.create(file.path(basefolder.name))
  }
  
  if(!dir.exists(file.path(basefolder.name, "counts"))) {
    dir.create(file.path(basefolder.name, "counts"))
  }
  
  for(i in seq(FCS.paths)) {
    
    fcs.raw <- read.FCS(FCS.paths[i], transformation = FALSE)
    
    if(trim) {
      exprs(fcs.raw) <- subset(exprs(fcs.raw),
                               exprs(fcs.raw)[, 'FSC-A'] > 0 
                               & exprs(fcs.raw)[, 'FSC-A'] < 250000
                               & exprs(fcs.raw)[, 'SSC-A'] > 0 
                               & exprs(fcs.raw)[, 'SSC-A'] < 250000 
                               & exprs(fcs.raw)[, 'FSC-W'] > 50000
                               & exprs(fcs.raw)[, 'FSC-W'] < 150000
                               & exprs(fcs.raw)[, 'SSC-W'] > 50000
                               & exprs(fcs.raw)[, 'SSC-W'] < 150000)
    }
    
    fcs.tmp <- compensate(fcs.raw, (keyword(fcs.raw)[grep("SPILL", names(keyword(fcs.raw)))][[1]]))
    
    fcs.tmp <- transform(fcs.tmp, fSOM.object$transform.object)
    
    message(paste("Mapping", fcs.tmp@description$`$FIL`, "to", "fSOM.object...",  i, "of ", length(FCS.paths), sep = " "))
    
    fSOM.tmp <- NewData(fSOM.object, fcs.tmp)
    
    if(write.count.files) {
      message("Generating tabled counts for each meta-cluster\n")
      counts <- as.data.frame(table(fSOM.tmp$metaclustering[fSOM.tmp$map$mapping[,1]]))
      names(counts) <- c("Meta.Cluster","Freq")
      
      counts.fname <- paste(sub(".fcs", "", fcs.raw@description$`$FIL`), 
                            basename(basefolder.name), 
                            "counts.csv", 
                            sep = "_")
      
      outCounts <- file.path(basefolder.name, "counts", counts.fname)
      
      write.csv(counts, outCounts, row.names = FALSE)
    }
    
    if(write.fcs.files) {
      message(paste0(paste("Subsetting fSOM-mapped", fcs.raw@description$`$FIL`, "by chosen meta-cluster(s)", sep = " ")), "\n")
      index <- which(!is.na(match(fSOM.tmp$metaclustering[fSOM.tmp$map$mapping[,1]], meta.cluster)))
      
      exprs(fcs.raw) <- exprs(fcs.raw)[index, ]
      
      fcs.raw@description$`$FIL` <- paste0(sub(".fcs", "", fcs.raw@description$`$FIL`),
                                           "_",
                                           fcsfolder.name,
                                           ".fcs"
      )
      
      outFile <- file.path(basefolder.name, fcs.raw@description$`$FIL`)
      write.FCS(fcs.raw, outFile)
    }
  }
}

newfcs.fromsom <- function(fSOM.object, meta.cluster = NULL, FCS.paths, basefolder.name, write.count.files = FALSE, write.fcs.files = FALSE, fcsfolder.name = NULL) {
  for(i in seq(FCS.paths)) {
    if (!dir.exists(file.path(dirname(FCS.paths[i]), basefolder.name))) { 
         dir.create(file.path(dirname(FCS.paths[i]), basefolder.name))
    }
    if (write.count.files) {
      if (!dir.exists(file.path(dirname(FCS.paths[i]), basefolder.name, paste(basefolder.name, "counts", sep = "_")))) {
           dir.create(file.path(dirname(FCS.paths[i]), basefolder.name, paste(basefolder.name, "counts", sep = "_")))
      }
    }
    if (write.fcs.files) {
      if (!dir.exists(file.path(dirname(FCS.paths[i]), basefolder.name, fcsfolder.name))) {
           dir.create(file.path(dirname(FCS.paths[i]), basefolder.name, fcsfolder.name))
      }
    }
    FCS.raw <- read.FCS(FCS.paths[i])
    exprs(FCS.raw) <- subset(exprs(FCS.raw), 
                           exprs(FCS.raw)[, 'FSC-A'] > 0 & exprs(FCS.raw)[, 'FSC-A'] < 250000 & 
                           exprs(FCS.raw)[, 'SSC-A'] > 0 & exprs(FCS.raw)[, 'SSC-A'] < 250000)
    FCS.COMP <- compensate(FCS.raw, (keyword(FCS.raw)[grep("SPILL", names(keyword(FCS.raw)))][[1]]))
    FCS.COMPTRANS <- transform(FCS.COMP, elgcl)
    message(paste("Mapping", FCS.raw@description$`$FIL`, "to", "fSOM.object...",  i, "of ", length(FCS.paths), sep = " "))
    fSOM.tmp <- NewData(fSOM.object, FCS.COMPTRANS)
    if (write.count.files) {
      message("Generating tabled counts for each meta-cluster\n")
      counts <- as.data.frame(table(fSOM.tmp$metaclustering[fSOM.tmp$map$mapping[,1]]))
      names(counts) <- c("Meta.Cluster","Freq")
      counts.fname <- paste(FCS.raw@description$`$SRC`, FCS.raw@description$`TUBE NAME`, basefolder.name, "counts.csv", sep="_")
      outCounts <- file.path(dirname(FCS.paths[i]), basefolder.name, paste(basefolder.name, "counts", sep = "_"), counts.fname)
      write.csv(counts, outCounts, row.names = FALSE)
    }      
    if (write.fcs.files) {
      message(paste0(paste("Subsetting fSOM-mapped", FCS.raw@description$`$FIL`, "by chosen meta-cluster(s)", sep = " ")), "\n")
      index <- which(!is.na(match(fSOM.tmp$metaclustering[fSOM.tmp$map$mapping[,1]], meta.cluster)))
      exprs(FCS.raw) <- exprs(FCS.raw)[index, ]
      FCS.raw@description$`$FIL` <- paste(FCS.raw@description$`$SRC`, FCS.raw@description$`TUBE NAME`, fcsfolder.name, ".fcs", sep="_")
      outFile <- file.path(dirname(FCS.paths[i]), basefolder.name, fcsfolder.name, FCS.raw@description$`$FIL`)
      write.FCS(FCS.raw, outFile)
    }
  }
}

get.elgcl <- function(fSOM.object.rds){
  fSOM.tmp <- readRDS(fSOM.object.rds)
  elgcl <- fSOM.tmp$elgcl
}

nodes.to.meta_fSOM <- function(fSOM.object, nodes) {
  levels(fSOM.object$metaclustering) <- c(levels(fSOM.object$metaclustering), length(levels(fSOM.object$metaclustering))+1)
  fSOM.object$metaclustering[nodes] <- length(levels(fSOM.object$metaclustering))
  fSOM.object$metaclustering <- factor(fSOM.object$metaclustering)
  levels(fSOM.object$metaclustering) <- c(1:length(levels(fSOM.object$metaclustering)))
  return(fSOM.object)
}

newfcs.fromsom.edit <- function(fSOM.object, 
                                meta.cluster = NULL, 
                                FCS.paths, 
                                basefolder.name, 
                                write.count.files = FALSE, 
                                write.fcs.files = FALSE, 
                                fcsfolder.name = NULL) {
  for(i in seq(FCS.paths)) {
    if (!dir.exists(file.path(dirname(FCS.paths[i]), basefolder.name))) { 
      dir.create(file.path(dirname(FCS.paths[i]), basefolder.name))
    }
    if (write.count.files) {
      if (!dir.exists(file.path(dirname(FCS.paths[i]), basefolder.name, paste(basefolder.name, "counts", sep = "_")))) {
        dir.create(file.path(dirname(FCS.paths[i]), basefolder.name, paste(basefolder.name, "counts", sep = "_")))
      }
    }
    if (write.fcs.files) {
      if (!dir.exists(file.path(dirname(FCS.paths[i]), basefolder.name, fcsfolder.name))) {
        dir.create(file.path(dirname(FCS.paths[i]), basefolder.name, fcsfolder.name))
      }
    }
    
    fcs.raw <- read.FCS(FCS.paths[i])
    
    exprs(fcs.raw) <- subset(exprs(fcs.raw),
                             exprs(fcs.raw)[, 'FSC-A'] > 0 & exprs(fcs.raw)[, 'FSC-A'] < 250000 
                             & exprs(fcs.raw)[, 'SSC-A'] > 0 & exprs(fcs.raw)[, 'SSC-A'] < 250000
                             & exprs(fcs.raw)[, 'FSC-W'] > 50000
                             & exprs(fcs.raw)[, 'FSC-W'] < 150000
                             & exprs(fcs.raw)[, 'SSC-W'] > 50000
                             & exprs(fcs.raw)[, 'SSC-W'] < 150000)
    
    fcs.tmp <- compensate(fcs.raw, (keyword(fcs.raw)[grep("SPILL", names(keyword(fcs.raw)))][[1]]))
    
    
    for(j in seq(length(fSOM.object$transform.object))){
      #print(paste0("transforming ", names(fSOM.object$transform.object)[i]))
      fcs.tmp <- transform(fcs.tmp, fSOM.object$transform.object[[j]]$asinh.translist)
      fcs.tmp@parameters@data$minRange[grep(names(fSOM.object$transform.object)[j], fcs.tmp@parameters@data$name)] <- min(exprs(fcs.tmp)[, names(fSOM.object$transform.object)[j]])
      fcs.tmp@parameters@data$maxRange[grep(names(fSOM.object$transform.object)[j], fcs.tmp@parameters@data$name)] <- max(exprs(fcs.tmp)[, names(fSOM.object$transform.object)[j]])
      fcs.tmp
    }
    
    message(paste("Mapping", fcs.tmp@description$`$FIL`, "to", "fSOM.object...",  i, "of ", length(FCS.paths), sep = " "))
    
    fSOM.tmp <- NewData(fSOM.object, fcs.tmp)
    
    if (write.count.files) {
      message("Generating tabled counts for each meta-cluster\n")
      counts <- as.data.frame(table(fSOM.tmp$metaclustering[fSOM.tmp$map$mapping[,1]]))
      names(counts) <- c("Meta.Cluster","Freq")
      if(grepl("HD", fcs.raw@description$`$FIL`)) {
        counts.fname <- paste(sub("_.fcs", "", fcs.raw@description$`$FIL`), 
                              basefolder.name, 
                              "counts.csv", 
                              sep="_")
      } else {
        counts.fname <- paste(sub("_.fcs", "", fcs.raw@description$`$FIL`), 
                              basefolder.name, 
                              "counts.csv", 
                              sep="_")
      }
      outCounts <- file.path(dirname(FCS.paths[i]), 
                             basefolder.name, 
                             paste(basefolder.name, 
                                   "counts", 
                                   sep = "_"), 
                             counts.fname)
      write.csv(counts, outCounts, row.names = FALSE)
    }      
    if (write.fcs.files) {
      message(paste0(paste("Subsetting fSOM-mapped", fcs.raw@description$`$FIL`, "by chosen meta-cluster(s)", sep = " ")), "\n")
      index <- which(!is.na(match(fSOM.tmp$metaclustering[fSOM.tmp$map$mapping[,1]], meta.cluster)))
      exprs(fcs.raw) <- exprs(fcs.raw)[index, ]
      if(grepl("HD", fcs.raw@description$`$FIL`)) {
        fcs.raw@description$`$FIL` <- paste(sub("_.fcs", "", fcs.raw@description$`$FIL`),
                                            fcsfolder.name,
                                            ".fcs",
                                            sep = "_")
      } else {
        fcs.raw@description$`$FIL` <- paste(sub("_.fcs", "", fcs.raw@description$`$FIL`),
                                            fcsfolder.name,
                                            ".fcs",
                                            sep = "_")
      }
      outFile <- file.path(dirname(FCS.paths[i]), basefolder.name, fcsfolder.name, fcs.raw@description$`$FIL`)
      write.FCS(fcs.raw, outFile)
    }
  }
}

newfcs.fromsom.scatter <- function(fSOM.object, 
                                meta.cluster = NULL, 
                                FCS.paths, 
                                basefolder.name, 
                                write.count.files = FALSE, 
                                write.fcs.files = FALSE, 
                                fcsfolder.name = NULL) {
  for(i in seq(FCS.paths)) {
    if (!dir.exists(file.path(dirname(FCS.paths[i]), basefolder.name))) { 
      dir.create(file.path(dirname(FCS.paths[i]), basefolder.name))
    }
    if (write.count.files) {
      if (!dir.exists(file.path(dirname(FCS.paths[i]), basefolder.name, paste(basefolder.name, "counts", sep = "_")))) {
        dir.create(file.path(dirname(FCS.paths[i]), basefolder.name, paste(basefolder.name, "counts", sep = "_")))
      }
    }
    if (write.fcs.files) {
      if (!dir.exists(file.path(dirname(FCS.paths[i]), basefolder.name, fcsfolder.name))) {
        dir.create(file.path(dirname(FCS.paths[i]), basefolder.name, fcsfolder.name))
      }
    }
    
    fcs.raw <- read.FCS(FCS.paths[i], transformation = FALSE)
    
    exprs(fcs.raw) <- subset(exprs(fcs.raw),
                             exprs(fcs.raw)[, 'FSC-A'] > 0 & exprs(fcs.raw)[, 'FSC-A'] < 250000 
                             & exprs(fcs.raw)[, 'SSC-A'] > 0 & exprs(fcs.raw)[, 'SSC-A'] < 250000
                             & exprs(fcs.raw)[, 'FSC-W'] > 50000
                             & exprs(fcs.raw)[, 'FSC-W'] < 150000
                             & exprs(fcs.raw)[, 'SSC-W'] > 50000
                             & exprs(fcs.raw)[, 'SSC-W'] < 150000)
    
    fcs.tmp <- fcs.raw
    
    message(paste("Mapping", fcs.tmp@description$`$FIL`, "to", "fSOM.object...",  i, "of ", length(FCS.paths), sep = " "))
    
    fSOM.tmp <- NewData(fSOM.object, fcs.tmp)
    
    if (write.count.files) {
      message("Generating tabled counts for each meta-cluster\n")
      counts <- as.data.frame(table(fSOM.tmp$metaclustering[fSOM.tmp$map$mapping[,1]]))
      names(counts) <- c("Meta.Cluster","Freq")
      if(grepl("HD", fcs.raw@description$`$FIL`)) {
        counts.fname <- paste(sub("_.fcs", "", fcs.raw@description$`$FIL`), 
                              basefolder.name, 
                              "counts.csv", 
                              sep="_")
      } else {
        counts.fname <- paste(sub("_.fcs", "", fcs.raw@description$`$FIL`), 
                              basefolder.name, 
                              "counts.csv", 
                              sep="_")
      }
      outCounts <- file.path(dirname(FCS.paths[i]), 
                             basefolder.name, 
                             paste(basefolder.name, 
                                   "counts", 
                                   sep = "_"), 
                             counts.fname)
      write.csv(counts, outCounts, row.names = FALSE)
    }      
    if (write.fcs.files) {
      message(paste0(paste("Subsetting fSOM-mapped", fcs.raw@description$`$FIL`, "by chosen meta-cluster(s)", sep = " ")), "\n")
      index <- which(!is.na(match(fSOM.tmp$metaclustering[fSOM.tmp$map$mapping[,1]], meta.cluster)))
      exprs(fcs.raw) <- exprs(fcs.raw)[index, ]
      if(grepl("HD", fcs.raw@description$`$FIL`)) {
        fcs.raw@description$`$FIL` <- paste(sub("_.fcs", "", fcs.raw@description$`$FIL`),
                                            fcsfolder.name,
                                            ".fcs",
                                            sep = "_")
      } else {
        fcs.raw@description$`$FIL` <- paste(sub("_.fcs", "", fcs.raw@description$`$FIL`),
                                            fcsfolder.name,
                                            ".fcs",
                                            sep = "_")
      }
      outFile <- file.path(dirname(FCS.paths[i]), basefolder.name, fcsfolder.name, fcs.raw@description$`$FIL`)
      write.FCS(fcs.raw, outFile)
    }
  }
}

newfcs.fromfsom.singlefolder <- function(fSOM.object, meta.cluster = NULL, FCS.paths, basefolder.name, fcsfolder.name, write.count.files = FALSE, write.fcs.files = FALSE) {
  
  if(!dir.exists(file.path(basefolder.name))) {
    dir.create(file.path(basefolder.name))
  }
  
  if(!dir.exists(file.path(basefolder.name, "counts"))) {
    dir.create(file.path(basefolder.name, "counts"))
  }
  
  for(i in seq(FCS.paths)) {
    
    fcs.raw <- read.FCS(FCS.paths[i], transformation = FALSE)
    
    fcs.tmp <- compensate(fcs.raw, (keyword(fcs.raw)[grep("SPILL", names(keyword(fcs.raw)))][[1]]))
    
    for(j in seq(fSOM.object$transform.object)){
      #print(paste0("transforming ", names(fSOM.object$transform.object)[i]))
      fcs.tmp <- transform(fcs.tmp, fSOM.object$transform.object[[j]]$asinh.translist)
      fcs.tmp@parameters@data$minRange[grep(names(fSOM.object$transform.object)[j], fcs.tmp@parameters@data$name)] <- min(exprs(fcs.tmp)[, names(fSOM.object$transform.object)[j]])
      fcs.tmp@parameters@data$maxRange[grep(names(fSOM.object$transform.object)[j], fcs.tmp@parameters@data$name)] <- max(exprs(fcs.tmp)[, names(fSOM.object$transform.object)[j]])
      fcs.tmp
    }
    
    message(paste("Mapping", fcs.tmp@description$`$FIL`, "to", "fSOM.object...",  i, "of ", length(FCS.paths), sep = " "))
    
    fSOM.tmp <- NewData(fSOM.object, fcs.tmp)
    
    if (write.count.files) {
      message("Generating tabled counts for each meta-cluster\n")
      counts <- as.data.frame(table(fSOM.tmp$metaclustering[fSOM.tmp$map$mapping[,1]]))
      names(counts) <- c("Meta.Cluster","Freq")
      
      counts.fname <- paste(sub("_.fcs", "", fcs.raw@description$`$FIL`), 
                            basename(basefolder.name), 
                            "counts.csv", 
                            sep="_")
      
      outCounts <- file.path(basefolder.name, "counts", counts.fname)
      
      write.csv(counts, outCounts, row.names = FALSE)
    }      
    if (write.fcs.files) {
      message(paste0(paste("Subsetting fSOM-mapped", fcs.raw@description$`$FIL`, "by chosen meta-cluster(s)", sep = " ")), "\n")
      index <- which(!is.na(match(fSOM.tmp$metaclustering[fSOM.tmp$map$mapping[,1]], meta.cluster)))
      exprs(fcs.raw) <- exprs(fcs.raw)[index, ]
      
      fcs.raw@description$`$FIL` <- paste(sub("_.fcs", "", fcs.raw@description$`$FIL`),
                                          fcsfolder.name,
                                          ".fcs",
                                          sep = "_")
      
      outFile <- file.path(basefolder.name, fcs.raw@description$`$FIL`)
      write.FCS(fcs.raw, outFile)
    }
  }
}

newfcs.fromfsom.scatter.singlefolder <- function(fSOM.object, meta.cluster = NULL, FCS.paths, comp.trans = TRUE, basefolder.name, fcsfolder.name, write.count.files = FALSE, write.fcs.files = FALSE) {
  
  if(!dir.exists(file.path(basefolder.name))) {
    dir.create(file.path(basefolder.name))
  }
  
  if(!dir.exists(file.path(basefolder.name, "counts"))) {
    dir.create(file.path(basefolder.name, "counts"))
  }
  
  for(i in seq(FCS.paths)) {
    
    fcs.raw <- read.FCS(FCS.paths[i], transformation = FALSE)
    
    fcs.tmp <- fcs.raw
    
    message(paste("Mapping", fcs.tmp@description$`$FIL`, "to", "fSOM.object...",  i, "of ", length(FCS.paths), sep = " "))
    
    fSOM.tmp <- NewData(fSOM.object, fcs.tmp)
    
    if (write.count.files) {
      message("Generating tabled counts for each meta-cluster\n")
      counts <- as.data.frame(table(fSOM.tmp$metaclustering[fSOM.tmp$map$mapping[,1]]))
      names(counts) <- c("Meta.Cluster","Freq")
      
      counts.fname <- paste(sub("_.fcs", "", fcs.raw@description$`$FIL`), 
                            basename(basefolder.name), 
                            "counts.csv", 
                            sep="_")
      
      outCounts <- file.path(basefolder.name, "counts", counts.fname)
      
      write.csv(counts, outCounts, row.names = FALSE)
    }      
    if (write.fcs.files) {
      message(paste0(paste("Subsetting fSOM-mapped", fcs.raw@description$`$FIL`, "by chosen meta-cluster(s)", sep = " ")), "\n")
      index <- which(!is.na(match(fSOM.tmp$metaclustering[fSOM.tmp$map$mapping[,1]], meta.cluster)))
      exprs(fcs.raw) <- exprs(fcs.raw)[index, ]
      
      fcs.raw@description$`$FIL` <- paste(sub("_.fcs", "", fcs.raw@description$`$FIL`),
                                          fcsfolder.name,
                                          ".fcs",
                                          sep = "_")
      
      outFile <- file.path(basefolder.name, fcs.raw@description$`$FIL`)
      write.FCS(fcs.raw, outFile)
    }
  }
}
