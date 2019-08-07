source("/Volumes/nlaniewski/R_source/workflow_20190403.R")

## FCS files: path to, read in, and trim FSC/SSC values; write to new folder

fcs.paths <- list.files("./data_source/ICS2_with.spillover/NGL092", full.names = T, pattern = ".fcs")
fcs.paths <- fcs.paths[intersect(grep("SEB|(HD&SEB)", basename(fcs.paths)), grep("swift|merge|split|cct", basename(fcs.paths), ignore.case = T, invert = T))]

input.flowset <- read.flowSet(fcs.paths, transformation = F, truncate_max_range = F)
before.trim <- fsApply(input.flowset, function(frame) nrow(frame)) # event counts before trim

input.flowset <- scatter.trim(input.flowset)
after.trim <- fsApply(input.flowset, function(frame) nrow(frame)) # event counts after trim

before.trim - after.trim # difference

if(!dir.exists("./data_modified/ICS/SEB_trimmed")){
  dir.create("./data_modified/ICS/SEB_trimmed", recursive = T)
}
fsApply(input.flowset, function(frame) write.FCS(frame, filename = file.path("./data_modified/ICS/SEB_trimmed", frame@description$`$FIL`))) # write trimmed files

## FCS files: path to, read in, compensate, and transform trimmed files

fcs.paths <- list.files("./data_modified/ICS/SEB_trimmed/", full.names = T, pattern = ".fcs")

input.flowset <- read.flowSet(fcs.paths, transformation = F, truncate_max_range = F)
# input.flowset <- input.flowset[, -grep('Time', colnames(input.flowset))] # drop time column
input.flowset <- fset.compensate(input.flowset)

elgcl <- elgcl.pars(as(input.flowset, "flowFrame")) # estimate logicle transform parameters using concatenate of all files

input.flowset <- fset.transform(input.flowset, elgcl) # transform files using estimated parameters

channel.density(input.flowset[[1]], 10000) # spot-check transforms

## FlowSOM: cluster using transformed FCS files
fsom <- ReadInput(input.flowset, scale = T, transform = F) # generate object using 'ReadInput'
fsom$elgcl <- elgcl
fsom$markers <- sub("CD14", "CD14_LIVE.DEAD", sub("-", ".", sapply(seq(fsom$prettyColnames), function(i) strsplit(fsom$prettyColnames[i], " ")[[1]][1]))) # marker names
fsom$dims.used <- c("FSC.A", "FSC.W", "SSC.A", "SSC.W", "CD3", "CD14_LIVE.DEAD")
fsom$seed <- 4242 # define/store seed
set.seed(fsom$seed)
fsom <- BuildSOM(fsom, colsToUse = sort(match(fsom$dims.used, fsom$markers)), xdim = 10, ydim = 10, rlen = 20) # generate SOMs
fsom <- BuildMST(fsom) # build a minimal-spanning tree
PlotStars(fsom) # spot-check tree

fsom$nbclust <- 30 # number of meta-clusters to define
fsom$clustering <- ConsensusClusterPlus(t(fsom$map$codes), 
                                        maxK = fsom$nbclust, 
                                        reps = 100, 
                                        pItem = 0.9, 
                                        pFeature = 1, 
                                        title = "Clustering Consensus", 
                                        plot = "png",
                                        clusterAlg = "hc", 
                                        innerLinkage = "average", 
                                        finalLinkage = "average",
                                        distance = "euclidean", 
                                        seed = 4242)
fsom[["metaclustering"]] <- as.factor(fsom$clustering[[fsom$nbclust]]$consensusClass)
PlotStars(fsom, backgroundValues = fsom$metaclustering)
fsom.heatmap(fsom, fsom$dims.used, display_numbers = F)

fcs.explorer(fsom.dat(fsom), "CD14_LIVE.DEAD", "FSC.A", "CD3", "CD14_LIVE.DEAD", 10000)

dead.nodes <- which(fsom$map$medianValues[, grep("LIVE.DEAD", fsom$markers)] > 0) # vector of nodes(SOMs) that are CD14+/dead
fsom <- nodes.to.meta_fSOM(fsom, dead.nodes)
fsom.heatmap(fsom, fsom$dims.used, display_numbers = F)
fcs.explorer(fsom.dat(fsom), "CD14_LIVE.DEAD", "FSC.A", "CD3", "CD14_LIVE.DEAD", 10000)

fsom <- nodes.to.meta_fSOM(fsom, grep(paste0("^", c(4,1,3,9,5,7), "$", collapse = "|"), fsom$metaclustering))
fsom.heatmap(fsom, fsom$dims.used, display_numbers = F)
fcs.explorer(fsom.dat(fsom), "CD14_LIVE.DEAD", "FSC.A", "CD3", "CD14_LIVE.DEAD", 10000)

fsom <- nodes.to.meta_fSOM(fsom, node.median.bycluster(fsom, 1, "CD3", "greater", 0))
fsom.heatmap(fsom, fsom$dims.used, display_numbers = F)
fcs.explorer(fsom.dat(fsom), "CD14_LIVE.DEAD", "FSC.A", "CD3", "CD14_LIVE.DEAD", 10000)

fsom$metaclustering.anno <- vector(mode = "list", length = length(levels(fsom$metaclustering)))
fsom$metaclustering.anno[[9]] <- "CD3+"

for(i in seq(fcs.paths)){
  fcs.add_fsom(fsom, fcs.paths[i], "./data_modified/ICS/20190404_fsom_pass1", "fsom_pass1")
}

## Subset mapped FCS files based on cluster-of-interest
fcs.paths <- list.files("./data_modified/ICS/20190404_fsom_pass1/", full.names = T, pattern = ".fcs")

for(i in seq(fcs.paths)){
  subset.by.fsomcluster(fcs.paths[i], 9, "./data_modified/ICS/20190404_fsom_pass1/CD3pos", "fsom_pass1_CD3pos")
}

## FCS files: path to, read in, compensate, and transform CD3+ files

fcs.paths <- list.files("./data_modified/ICS/20190404_fsom_pass1/CD3pos/", full.names = T, pattern = ".fcs")

input.flowset <- read.flowSet(fcs.paths, transformation = F, truncate_max_range = F)
input.flowset <- fset.compensate(input.flowset)

elgcl <- elgcl.pars(as(input.flowset, "flowFrame")) # estimate logicle transform parameters using concatenate of all files

input.flowset <- fset.transform(input.flowset, elgcl) # transform files using estimated parameters

channel.density(input.flowset[[1]], 10000) # spot-check transforms

## FlowSOM: cluster CD3+ using transformed FCS files
fsom <- ReadInput(input.flowset, scale = F, transform = F) # generate object using 'ReadInput'
fsom$elgcl <- elgcl
fsom$markers <- sub("CD14", "CD14_LIVE.DEAD", sub("-", ".", sapply(seq(fsom$prettyColnames), function(i) strsplit(fsom$prettyColnames[i], " ")[[1]][1]))) # marker names
fsom$dims.used <- c("CD4", "CD8A", "CD69")
fsom$seed <- 4242 # define/store seed
set.seed(fsom$seed)
fsom <- BuildSOM(fsom, colsToUse = sort(match(fsom$dims.used, fsom$markers)), xdim = 10, ydim = 10, rlen = 20) # generate SOMs
fsom <- BuildMST(fsom) # build a minimal-spanning tree
PlotStars(fsom) # spot-check tree

fsom$nbclust <- 15 # number of meta-clusters to define
fsom$clustering <- ConsensusClusterPlus(t(fsom$map$codes), 
                                        maxK = fsom$nbclust, 
                                        reps = 100, 
                                        pItem = 0.9, 
                                        pFeature = 1, 
                                        title = "Clustering Consensus", 
                                        plot = "png",
                                        clusterAlg = "hc", 
                                        innerLinkage = "average", 
                                        finalLinkage = "average",
                                        distance = "euclidean", 
                                        seed = 4242)
fsom[["metaclustering"]] <- as.factor(fsom$clustering[[fsom$nbclust]]$consensusClass)
PlotStars(fsom, backgroundValues = fsom$metaclustering)
fsom.heatmap(fsom, fsom$dims.used, display_numbers = F)

fcs.explorer(fsom.dat(fsom), "CD4", "CD8A", "CD3", "CD69", 10000)

fsom$metaclustering.anno <- vector(mode = "list", length = length(levels(fsom$metaclustering)))
fsom$metaclustering.anno[[3]] <- "CD4+ CD69+"
fsom$metaclustering.anno[[15]] <- "CD8+ CD69+"

for(i in seq(fcs.paths)){
  fcs.add_fsom(fsom, fcs.paths[i], "./data_modified/ICS/20190404_fsom_pass2/", "fsom_pass2")
}

## Subset mapped FCS files based on cluster-of-interest
fcs.paths <- list.files("./data_modified/ICS/20190404_fsom_pass2/", full.names = T, pattern = ".fcs")

for(i in seq(fcs.paths)){
  subset.by.fsomcluster(fcs.paths[i], 3, "./data_modified/ICS/20190404_fsom_pass2/CD4posCD69pos", "fsom_pass2_CD4posCD69pos")
}

for(i in seq(fcs.paths)){
  subset.by.fsomcluster(fcs.paths[i], 15, "./data_modified/ICS/20190404_fsom_pass2/CD8posCD69pos", "fsom_pass2_CD8posCD69pos")
}

## FCS files: path to, read in, compensate, and transform CD4+CD69+ files
fcs.paths <- list.files("./data_modified/ICS/20190404_fsom_pass2/CD4posCD69pos/", full.names = T, pattern = ".fcs")

input.flowset <- read.flowSet(fcs.paths, transformation = F, truncate_max_range = F)
input.flowset <- input.flowset[, -which(colnames(input.flowset) %in% c("node", "cluster"))]
input.flowset <- fset.compensate(input.flowset)

elgcl <- elgcl.pars(as(input.flowset, "flowFrame")) # estimate logicle transform parameters using concatenate of all files

input.flowset <- fset.transform(input.flowset, elgcl) # transform files using estimated parameters

channel.density(input.flowset[[1]], 10000) # spot-check transforms

flowset.cols <- sub("CD14", "CD14_LIVE.DEAD", get.markers(input.flowset[[1]], " "))
dim.vec <- c(7,8,10,11,14,15,16,17,18,21)
dims <- flowset.cols[dim.vec]

meta.data <- data.frame(file = sampleNames(input.flowset),
                        sample.ID = sub("_SEB.*", "", meta.data$file),
                        visit = NA,
                        term = NA)
meta.data$visit[grep("_19", meta.data$sample.ID)] <- "12-month"
meta.data$visit[grep("_7", meta.data$sample.ID)] <- "Discharge"
meta.data$visit[grep("_1$", meta.data$sample.ID)] <- "Birth"
meta.data$visit[which(is.na(meta.data$visit))] <- "Adult"
meta.data$term[grep("RPRC053", meta.data$sample.ID)] <- "Full-term"
meta.data$term[grep("RPRC054", meta.data$sample.ID)] <- "Pre-term"
meta.data$term[which(is.na(meta.data$term))] <- "Adult"
meta.data$term.visit <- paste(meta.data$term, meta.data$visit, sep = "_")

mds.markers.median(input.flowset, dims, factor(meta.data$visit))
mds.markers.median(input.flowset[3:19], dims, factor(meta.data$visit[3:19]))
mds.markers.median(input.flowset[3:19], dims, factor(meta.data$term[3:19]))
mds.markers.median(input.flowset[3:19], dims, factor(meta.data$term.visit[3:19]))

## Define a function that calculates the NRS per sample
NRS <- function(x, ncomp = 3){
  pr <- prcomp(x, center = TRUE, scale. = FALSE)
  score <- rowSums(outer(rep(1, ncol(x)),
                         pr$sdev[1:ncomp]^2) * abs(pr$rotation[,1:ncomp]))
  return(score)
}

## Calculate the score
nrs_sample <- fsApply(input.flowset[, dim.vec], NRS, use.exprs = TRUE)
rownames(nrs_sample) <- sub("_SEB.*", "", sampleNames(input.flowset))
colnames(nrs_sample) <- flowset.cols[dim.vec]
nrs <- colMeans(nrs_sample, na.rm = TRUE)

## Plot the NRS for ordered markers
markers_ord <- names(sort(nrs, decreasing = TRUE))
nrs_sample <- data.frame(nrs_sample)
nrs_sample$sample_id <- rownames(nrs_sample)

ggdf <- melt(nrs_sample, id.var = "sample_id",
             value.name = "nrs", variable.name = "marker")

ggdf$marker <- factor(ggdf$marker, levels = markers_ord)
ggdf$sample.type <- meta.data$term.visit[match(ggdf$sample_id, meta.data$sample.ID)]

ggplot(ggdf, aes(x = marker, y = nrs)) +
  geom_point(aes(color = sample.type), alpha = 0.9,
             position = position_jitter(width = 0.3, height = 0)) +
  geom_boxplot(outlier.color = NA, fill = NA) +
  stat_summary(fun.y = "mean", geom = "point", shape = 21, fill = "white") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

## FlowSOM: cluster CD4+CD69+ using transformed FCS files
fsom <- ReadInput(input.flowset, scale = F, transform = F) # generate object using 'ReadInput'
fsom$elgcl <- elgcl
fsom$markers <- sub("CD14", "CD14_LIVE.DEAD", sub("-", ".", sapply(seq(fsom$prettyColnames), function(i) strsplit(fsom$prettyColnames[i], " ")[[1]][1]))) # marker names
fsom$dims.used <- c("IL8", "IFNG", "IL4", "TNFA", "CD45RA", "IL2", "IL17")
fsom$seed <- 4242 # define/store seed
set.seed(fsom$seed)
fsom <- BuildSOM(fsom, colsToUse = sort(match(fsom$dims.used, fsom$markers)), xdim = 11, ydim = 11, rlen = 20) # generate SOMs
fsom <- BuildMST(fsom) # build a minimal-spanning tree
PlotStars(fsom) # spot-check tree

fsom$nbclust <- 25 # number of meta-clusters to define
fsom$clustering <- ConsensusClusterPlus(t(fsom$map$codes), 
                                        maxK = fsom$nbclust, 
                                        reps = 100, 
                                        pItem = 0.9, 
                                        pFeature = 1, 
                                        title = "Clustering Consensus", 
                                        plot = "png",
                                        clusterAlg = "hc", 
                                        innerLinkage = "average", 
                                        finalLinkage = "average",
                                        distance = "euclidean", 
                                        seed = 4242)
fsom[["metaclustering"]] <- as.factor(fsom$clustering[[fsom$nbclust]]$consensusClass)
PlotStars(fsom, backgroundValues = fsom$metaclustering)
fsom.heatmap(fsom, fsom$dims.used, display_numbers = F)

fcs.explorer(fsom.dat(fsom), "IL8", "IL2", "TNFA", "IL2", 15000)

fsom <- nodes.to.meta_fSOM(fsom, grep("^1$|^3$", fsom$metaclustering))
fsom <- nodes.to.meta_fSOM(fsom, grep("^2$|^7$", fsom$metaclustering))
fsom <- nodes.to.meta_fSOM(fsom, grep("^3$|^23$", fsom$metaclustering))
fsom <- nodes.to.meta_fSOM(fsom, grep("^13$|^19$", fsom$metaclustering))

fsom.heatmap(fsom, fsom$dims.used, display_numbers = F, use.scaled.expr = F)
fcs.explorer(fsom.dat(fsom), "IL8", "IL2", "TNFA", "IL2", 15000)

fsom$metaclustering.anno <- vector(mode = "list", length = length(levels(fsom$metaclustering)))

counts <- as.data.frame.matrix(table(cell_clustering1, rep(meta.data$sample.ID, fsApply(input.flowset, nrow))))
props <- as.data.frame.matrix(t(t(counts_table) / colSums(counts_table)) * 100)

ggdf <- melt(data.frame(cluster = rownames(props), props), 
             id.vars = "cluster", value.name = "proportion", variable.name = "sample_id")
ggdf$condition <- meta.data$visit[ggdf$sample_id]

ggplot(ggdf, aes(x = sample_id, y = proportion, fill = cluster)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ condition, scales = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = color_clusters) 

for(i in seq(fcs.paths)){
  fcs.add_fsom(fsom, fcs.paths[i], "./data_modified/ICS/20190404_fsom_pass2/", "fsom_pass2")
}

## Subset mapped FCS files based on cluster-of-interest
fcs.paths <- list.files("./data_modified/ICS/20190404_fsom_pass2/", full.names = T, pattern = ".fcs")

for(i in seq(fcs.paths)){
  subset.by.fsomcluster(fcs.paths[i], 3, "./data_modified/ICS/20190404_fsom_pass2/CD4posCD69pos", "fsom_pass2_CD4posCD69pos")
}

for(i in seq(fcs.paths)){
  subset.by.fsomcluster(fcs.paths[i], 15, "./data_modified/ICS/20190404_fsom_pass2/CD8posCD69pos", "fsom_pass2_CD8posCD69pos")
}

##