fset.trim <- function(flowSet){
  ## Trims extreme light-scatter area/width values in FCS file(s).
  ## Args:
  ##  flowSet: a 'class' flowSet created using flowCore::read.flowSet()
  ##
  ## Returns:
  ##  FCS file(s) with a modified expression matrix.
  fsApply(flowSet, FUN = function(frame) {
    
    print(paste0("trimming scatter: area/width ", frame@description$`$FIL`))
    
    exprs(frame) <- subset(exprs(frame),
                           exprs(frame)[, 'FSC-A'] > 0 
                           & exprs(frame)[, 'FSC-A'] < 250000
                           & exprs(frame)[, 'SSC-A'] > 0 
                           & exprs(frame)[, 'SSC-A'] < 250000 
                           & exprs(frame)[, 'FSC-W'] > 50000
                           & exprs(frame)[, 'FSC-W'] < 150000
                           & exprs(frame)[, 'SSC-W'] > 50000
                           & exprs(frame)[, 'SSC-W'] < 150000)
    frame
    })
}

fset.trim.channel <- function(flowset.samples, channel, cutoff) {
  fsApply(flowset.samples, FUN = function(frame) {
    print(paste0("trimming ", channel, "in ", frame@description$`$FIL`))
    exprs(frame) <- subset(exprs(frame), exprs(frame)[, channel] > cutoff)
    frame
  })
}

fset.rangefix <- function(flowSet, channel) {
  fsApply(flowSet, function(frame) {
    frame@parameters@data$minRange[grep(channel, frame@parameters@data$name)] <- min(exprs(frame)[, channel])
    frame@parameters@data$maxRange[grep(channel, frame@parameters@data$name)] <- max(exprs(frame)[, channel])
    frame
  })
}

fset.downsample <- function(flowSet, size) {
  fsApply(flowSet, function(frame) {
    samp <- sample(1:nrow(exprs(frame)), size, replace = TRUE)
    exprs(frame) <- exprs(frame)[samp, ]
    frame
    })
}

FCS.subsampled <- function(FCS.paths, samp.size, newfolder.name, post.name){
  ##  Randomly subsamples FCS file(s)
  for (i in seq(FCS.paths)){
    ##read in and sample raw FCS (linear/non-transformed data); use emptyValue to avoild delimiter keyword error (?)
    FCS.raw <- read.FCS(FCS.paths[i])
    #print(paste0("subsampling... ", FCS.raw@description$`$FIL`))
    print(paste("Subsampling ", FCS.raw@description$`$FIL`, i, "of ", length(FCS.paths)))
    if (nrow(exprs(FCS.raw)) > samp.size) {
      samp <- sample(1:nrow(exprs(FCS.raw)), samp.size, replace = TRUE)
      exprs(FCS.raw) <- exprs(FCS.raw)[samp, ]
      FCS.raw@description$`$FIL` <- paste(FCS.raw@description$`$SRC`, FCS.raw@description$`TUBE NAME`, post.name, sep ="_")
    } else {
      FCS.raw@description$`$FIL` <- paste(FCS.raw@description$`$SRC`, FCS.raw@description$`TUBE NAME`, post.name, sep ="_")
    }
    ##write FCS to directory
    if (dir.exists(file.path(dirname(FCS.paths[i]), newfolder.name))) {
      outFile <- file.path(dirname(FCS.paths[i]), newfolder.name, FCS.raw@description$`$FIL`)
    } else {
      dir.create(file.path(dirname(FCS.paths[i]), newfolder.name))
      outFile <- file.path(dirname(FCS.paths[i]), newfolder.name, FCS.raw@description$`$FIL`)
    }
    write.FCS(FCS.raw, outFile)
  }
}

FCS.subsampled.singlefolder <- function(FCS.paths, samp.size, newfolder.path.name){
  ##  subsamples FCS file(s); randomly selects rows per file
  for (i in seq(FCS.paths)){
    ##read in and sample raw FCS (linear/non-transformed data)
    FCS.raw <- read.FCS(FCS.paths[i], transformation = FALSE)
    
    print(paste("Subsampling ", FCS.raw@description$`$FIL`, i, "of ", length(FCS.paths)))
    
    samp <- sample(1:nrow(exprs(FCS.raw)), samp.size, replace = TRUE)
    exprs(FCS.raw) <- exprs(FCS.raw)[samp, ]
    FCS.raw@description$`$FIL` <- sub(".fcs", "_subsampled.fcs", FCS.raw@description$`$FIL`)
    
    ##write FCS to directory
    if (dir.exists(newfolder.path.name)) {
      outFile <- file.path(newfolder.path.name, FCS.raw@description$`$FIL`)
    } else {
      dir.create(newfolder.path.name)
      outFile <- file.path(newfolder.path.name, FCS.raw@description$`$FIL`)
    }
    write.FCS(FCS.raw, outFile)
  }
}

fset.compensate <- function(flowSet.samples){fsApply(flowSet.samples, function(frame) {
  ## compensate each 'frame'(sample) using stored compensation matrix (FCS header; keyword = $SPILL or SPILLOVER)
  print(paste0("compensating ", frame@description$`$FIL`))
  comp <- keyword(frame)[grep("SPILL", names(keyword(frame)))][[1]]
  frame_comped <- compensate(frame, comp)
  frame_comped
  })
}


fset.transform <- function(flowSet.samples) {
  ## transform fcs data using ^ estimated logicile
  fsApply(flowSet.samples, function(frame) {
  f.name <- frame@description$`$FIL`
  print(paste0("transforming ", frame@description$`$FIL`))
  frame_trans <- transform(frame, elgcl)
  frame_trans
  })
}


get.markers <- function(frame) {gsub("-",".",unlist(lapply(1:length(frame@parameters@data$desc), function(x) {
  if(is.na(frame@parameters@data$desc[[x]])){
    (strsplit(frame@parameters@data$name," ")[[x]][1])
  } else {
    (strsplit(frame@parameters@data$desc," ")[[x]][1])
  }
  })))
}

elgcl.transform <- function(compensated.set, concatenate.size) {
  concatenate <- as(compensated.set, "flowFrame")
  exprs(concatenate) <- exprs(concatenate)[sample(1:nrow(exprs(concatenate)), concatenate.size), ]
  channels.fluor <- as.vector(concatenate@parameters$name[-grep(c("FSC|SSC|Time|Original"), (concatenate@parameters$name))])
  elgcl <- estimateLogicle(concatenate, channel = channels.fluor)
  elgcl
}

median.QC <- function(fcs.file.paths) {
  median.QC.list <- vector(mode = "list", length = length(fcs.file.paths))
  for (i in seq(fcs.file.paths)) {
    fcs.tmp <- fset.transform(fset.compensate(fset.trim(read.flowSet(fcs.file.paths[i]))))
    median.QC.list[[i]] <- apply(exprs(fcs.tmp[[1]]), 2, median)
    names(median.QC.list)[i] <- fcs.tmp[[1]]@description$`$FIL`
  }
  dat <- as.data.frame(t(data.frame(median.QC.list)))
  colnames(dat) <- get.markers(read.FCS(fcs.file.paths[1]))
  dat$random <- runif(length(row.names(dat)), 1, length(row.names(dat)))
  dat
}

median.QC.edit <- function(fcs.file.paths, transform.object) {
  median.QC.list <- vector(mode = "list", length = length(fcs.file.paths))
  
  for (i in seq(fcs.file.paths)) {
    
    print(paste("Calculating transformed medians...", i, "of", length(fcs.file.paths), sep = " "))
    
    fcs.tmp <- read.FCS(fcs.file.paths[i], transformation = FALSE)
    
    # exprs(fcs.tmp) <- subset(exprs(fcs.tmp),
    #                          exprs(fcs.tmp)[, 'FSC-A'] > 0 & exprs(fcs.tmp)[, 'FSC-A'] < 250000 
    #                          & exprs(fcs.tmp)[, 'SSC-A'] > 0 & exprs(fcs.tmp)[, 'SSC-A'] < 250000
    #                          & exprs(fcs.tmp)[, 'FSC-W'] > 50000
    #                          & exprs(fcs.tmp)[, 'FSC-W'] < 150000
    #                          & exprs(fcs.tmp)[, 'SSC-W'] > 50000
    #                          & exprs(fcs.tmp)[, 'SSC-W'] < 150000)
    
    fcs.tmp <- compensate(fcs.tmp, (keyword(fcs.tmp)[grep("SPILL", names(keyword(fcs.tmp)))][[1]]))
    fcs.tmp <- fcs.tmp[, colnames(fcs.tmp)[grep(paste(names(transform.object), collapse = "|"), colnames(fcs.tmp))]]
    
    for(j in seq(length(transform.object))){
      #print(paste0("transforming ", names(fSOM.object$transform.object)[i]))
      fcs.tmp <- transform(fcs.tmp, transform.object[[j]]$asinh.translist)
      fcs.tmp@parameters@data$minRange[grep(names(transform.object)[j], fcs.tmp@parameters@data$name)] <- min(exprs(fcs.tmp)[, names(transform.object)[j]])
      fcs.tmp@parameters@data$maxRange[grep(names(transform.object)[j], fcs.tmp@parameters@data$name)] <- max(exprs(fcs.tmp)[, names(transform.object)[j]])
      fcs.tmp
    }
    
    median.QC.list[[i]] <- apply(exprs(fcs.tmp), 2, median)
    names(median.QC.list)[i] <- fcs.tmp@description$`$FIL`
  }
  dat <- as.data.frame(t(data.frame(median.QC.list)))
  colnames(dat) <- get.markers(read.FCS(fcs.file.paths[1]))
  dat$random <- runif(length(row.names(dat)), 1, length(row.names(dat)))
  dat
}

median.QC.scatter <- function(fcs.file.paths) {
  median.QC.list <- vector(mode = "list", length = length(fcs.file.paths))
  
  for (i in seq(fcs.file.paths)) {
    
    print(paste(i, "of", length(fcs.file.paths), sep = " "))
    
    fcs.tmp <- read.FCS(fcs.file.paths[i])
    
    exprs(fcs.tmp) <- exprs(fcs.tmp)[, 1:6]
    
    exprs(fcs.tmp) <- subset(exprs(fcs.tmp),
                             exprs(fcs.tmp)[, 'FSC-A'] > 0 & exprs(fcs.tmp)[, 'FSC-A'] < 250000 
                             & exprs(fcs.tmp)[, 'SSC-A'] > 0 & exprs(fcs.tmp)[, 'SSC-A'] < 250000
                             & exprs(fcs.tmp)[, 'FSC-W'] > 50000
                             & exprs(fcs.tmp)[, 'FSC-W'] < 150000
                             & exprs(fcs.tmp)[, 'SSC-W'] > 50000
                             & exprs(fcs.tmp)[, 'SSC-W'] < 150000)
    
    median.QC.list[[i]] <- apply(exprs(fcs.tmp), 2, median)
    names(median.QC.list)[i] <- fcs.tmp@description$`$FIL`
  }
  dat <- as.data.frame(t(data.frame(median.QC.list)))
  colnames(dat) <- get.markers(read.FCS(fcs.file.paths[1]))[1:length(colnames(dat))]
  dat$random <- runif(length(row.names(dat)), 1, length(row.names(dat)))
  dat
}

##inverse of flowCore::arcsinhTransform
sinhTransform <- function (transformationId = "defaultInverseArcsinhTransform", a = 1, b = 1, c = 0){
  t <- new("transform", .Data = function(x) (sinh(x - c) - a) * 1/(b))
  t@transformationId <- transformationId
  t
}

fset.decompensate <- function(flowSet.samples){fsApply(flowSet.samples, function(frame) {
  ## decompensate each 'frame'(sample) using stored compensation matrix (FCS header; keyword = $SPILL or SPILLOVER)
  print(paste0("de-compensating ", frame@description$`$FIL`))
  cols <- colnames(keyword(frame)[grep("SPILL", names(keyword(frame)))][[1]])
  sel <- cols %in% colnames(frame)
  if(!all(sel)) {
    stop(keyword(frame)[["FILENAME"]], 
         "\\nThe following parameters in the spillover matrix are not present in the flowFrame:\\n",
         paste(cols[!sel], collapse=", "), call.=FALSE)
  }
  e <- exprs(frame)
  e[, cols] <- e[, cols] %*% keyword(frame)[grep("SPILL", names(keyword(frame)))][[1]]
  exprs(frame) = e
  frame

  # comp <- keyword(frame)[grep("SPILL", names(keyword(frame)))][[1]]
  # frame_decomped <- decompensate(frame, as.matrix(comp))
  # frame_decomped
})
}

fset.equality <- function(fset.1, fset.2, channels) {
  sapply(seq(length(fset.1)), function(i) all.equal(fset.1[[i]]@exprs[, channels], fset.2[[i]]@exprs[, channels]), simplify = TRUE)
}

transform.set <- function(input.set){
  
  for(j in seq(length(trans.obj))){
    print(paste0("transforming ", names(trans.obj)[j]))
    input.set <- transform(input.set, trans.obj[[j]]$asinh.translist)
    input.set <- fsApply(input.set, function(i){
      i@parameters@data$minRange[grep(names(trans.obj)[j], i@parameters@data$name)] <- min(exprs(i)[, names(trans.obj)[j]])
      i@parameters@data$maxRange[grep(names(trans.obj)[j], i@parameters@data$name)] <- max(exprs(i)[, names(trans.obj)[j]])
      i
    })
  }
  input.set
} # needs 'trans.obj'

inverse.transform.set <- function(input.set){
  
  for(j in seq(length(trans.obj))){
    sinhTrans <- sinhTransform(transformationId = "inverse.arcsinh-transformation",
                               a = trans.obj[[j]]$asinh.vals$a,
                               b = trans.obj[[j]]$asinh.vals$b, 
                               c = trans.obj[[j]]$asinh.vals$c)
    sinh.translist <- transformList(names(trans.obj)[j], sinhTrans)
    
    input.set <- transform(input.set, sinh.translist)
  }
  input.set
} # needs 'trans.obj

warp.set <- function(input.set, target.name){
  
  for(i in seq(warp.parameters)){
    input.set <- warpSet(input.set, 
                         stains    = warp.parameters[[i]][["channel"]], 
                         monwrd    = warp.parameters[[i]][["monwrd"]], 
                         peakNr    = warp.parameters[[i]][["peakNr"]], 
                         clipRange = warp.parameters[[i]][["clipRange"]], 
                         bwFac     = warp.parameters[[i]][["bwFac"]], 
                         target    = sampleNames(input.set)[grep(target.name, sampleNames(input.set))])
  }
  input.set
} # needs 'warp.parameters' and HD target in flowSet

warp.qc.plots <- function(input.set, input.set.warped){
  
  before.set <- input.set
  after.set <- input.set.warped
  
  qc.plots <- vector(mode = "list", length = length(warp.parameters)*2)
  
  pData(before.set)[[1]] <- fsApply(before.set, function(i) i@description$`$FIL` <- paste(i@description$`$SRC`, i@description$`EXPERIMENT NAME`, sep = "_"))
  pData(after.set)[[1]] <- fsApply(after.set, function(i) i@description$`$FIL` <- paste(i@description$`$SRC`, i@description$`EXPERIMENT NAME`, sep = "_"))
  
  before.set <- fset.downsample(before.set, 5000)
  after.set <- fset.downsample(after.set, 5000)
  
  for(i in seq(warp.parameters)){
    qc.plots[[i+(i-1)]] <- densityplot(~ ., 
                                       before.set[, warp.parameters[[i]]$channel], 
                                       filter = curv1Filter(warp.parameters[[i]]$channel, bwFac = warp.parameters[[i]]$bwFac),
                                       main = paste("Non-warped", names(warp.parameters)[i], sep = " "))
    qc.plots[[i+i]]     <- densityplot(~ ., 
                                       after.set[, warp.parameters[[i]]$channel], 
                                       filter = curv1Filter(warp.parameters[[i]]$channel, bwFac = warp.parameters[[i]]$bwFac),
                                       main = paste("Warped", names(warp.parameters)[i], sep = " "))
  }
  qc.plots
  
  doc.name <- paste("Warp.QC.densityplots", Sys.Date(), paste(names(warp.parameters), collapse = "."), batch.name, "pdf", sep = ".")
  pdf(doc.name)
  print(qc.plots)
  dev.off()
} # needs 'warp.parameters'
