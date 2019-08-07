library(flowCore); library(FlowSOM); library(pheatmap); library(shiny); library(ggplot2); library(viridis); library(Rtsne); library(gganimate)
###
input.files <- list.files("./data_source/ICS2_with.spillover/", recursive = TRUE, full.names = TRUE, pattern = "*SEB.fcs$")
input.files[grep("RPRC0530106", input.files)]

tmp <- read.flowSet(input.files[grep("RPRC0530106", input.files)], transformation = F)

tmp <- fsApply(tmp, function(i) compensate(i, i@description$`$SPILLOVER`))

fsom <- readRDS("./analyses_R/FSOM_analysis/01_fSOM.CD3LIVE.rds")
fsom$markers <- get.markers(tmp[[1]])
fsom$markers[grep("CD14", fsom$markers)] <- "CD14_LIVEDEAD"

tmp <- transform(tmp, fsom$elgcl)

# tmp.inverse <- transform(tmp.trans, inverseLogicleTransform(estimateLogicle(tmp.comp[[1]], colnames(tmp.comp[[1]])[-grep("FSC|SSC|Time|Original", colnames(tmp.comp[[1]]))])))
# 
# tmp.comp[[1]]@exprs[, 7][1] == tmp.inverse[[1]]@exprs[, 7][1]                
# 
# all(round(tmp.comp[[1]]@exprs, digits = 5) == round(tmp.inverse[[1]]@exprs, digits = 5))

fsom.tmp <- NewData(fsom, tmp[[1]])
dat <- fsom.dat(fsom.tmp)
index <- sample(1:nrow(dat), 10000)

tsne.out <- Rtsne(dat[index, ][, fsom.tmp$map$colsUsed], perplexity = 50, verbose = TRUE)

dat.tsne <- cbind(dat[index, ], data.frame(tsne.x = tsne.out$Y[, 1], tsne.y = tsne.out$Y[, 2], frame = 1))
vec <- sapply(unique(dat.tsne$MCluster), function(i) nrow(subset(dat.tsne, MCluster == i)))

tmp.dat <- as.data.frame(aggregate(dat.tsne, list(meta.clusters = dat.tsne$MCluster), median)[, c("MCluster", "tsne.x", "tsne.y")])
frame2 <- list()
frame2$tsne.x <- unlist(sapply(seq(vec), function(i) rep(tmp.dat$tsne.x[i], vec[i])))
frame2$tsne.y <- unlist(sapply(seq(vec), function(i) rep(tmp.dat$tsne.y[i], vec[i])))
frame2$MCluster <- unlist(sapply(seq(vec), function(i) rep(tmp.dat$MCluster[i], vec[i])))

tmp.frame1 <- dat.tsne[, -grep("tsne", colnames(dat.tsne))]
tmp.frame2 <- data.frame(frame2)

dat.tsne <- rbind(dat.tsne, cbind(dat.tsne[, -grep("tsne|MCluster|frame", colnames(dat.tsne))], frame2, frame = 2))


x.lims <- c(min(dat.tsne$tsne.x), max(dat.tsne$tsne.x))
y.lims <- c(min(dat.tsne$tsne.y), max(dat.tsne$tsne.y))

ggplot(dat.tsne, aes(tsne.x, tsne.y)) +
  geom_point() +
  transition_states(frame)


tmp.fsomed <- fsApply(tmp, function(i){
  fsom.tmp <- NewData(fsom, i)
  i <- fr_append_cols(i, matrix(fsom.tmp$map$mapping[, 1], dimnames = list(NULL, "node_PASS1")))
  i <- fr_append_cols(i, matrix(as.integer(fsom.tmp$metaclustering[fsom.tmp$map$mapping[, 1]]), dimnames = list(NULL, "meta.cluster_PASS1")))
})

fsom.heat <- function(fsom.object){
  heat <- aggregate(fsom.object$map$codes, list(meta.cluster = fsom.object$metaclustering), median)[, -1]
  rownames(heat) <- paste0(rownames(heat), "_meta.cluster")
  colnames(heat) <- fsom.object$markers[fsom.object$map$colsUsed]
  heat
}

fsom$heatmap <- pheatmap(scale(fsom.heat(fsom)), scale = "row")
fcs.explorer(fsom.dat(fsom.tmp), "CD14_LIVEDEAD", "FSC.A", 20000)

fsom$metaclustering.anno <- vector(mode = "list", length = length(levels(fsom$metaclustering)))
fsom$metaclustering.anno[[17]] <- "CD3+ LIVE"

tmp[[1]] <- fsApply(tmp, function(i){
  i <- fr_append_cols(i, matrix(fsom$map$mapping[, 1], dimnames = list(NULL, "node_PASS1")))
  i <- fr_append_cols(i, matrix(as.integer(fsom$metaclustering[fsom$map$mapping[, 1]]), dimnames = list(NULL, "meta.cluster_PASS1")))
})

write.FCS(tmp[[1]], "./temp/test.fcs")

tmp.trans <- fsApply(tmp.trans, function(i){
  i <- fr_append_cols(i, matrix(fsom$map$mapping[, 1], dimnames = list(NULL, "node_PASS1")))
  i <- fr_append_cols(i, matrix(as.integer(fsom$metaclustering[fsom$map$mapping[, 1]]), dimnames = list(NULL, "meta.cluster_PASS1")))
})

tmp.trans.cd3_live <-fsApply(tmp.trans, function(i){
  exprs(i) <- exprs(i)[which(exprs(i)[, "meta.cluster_PASS1"] == 16), ]
  i
})

fsom <- ReadInput(tmp.trans.cd3_live, scale = T)
fsom$markers <- get.markers(tmp.trans.cd3_live[[1]])
fsom$markers[grep("CD14", fsom$markers)] <- "LIVE.DEAD_CD14"
fsom$transform.object <- fsom_PASS1$transform.object
fsom <- BuildSOM(fsom, 
                 colsToUse = match(c("CD3", "CD4", "CD8A", "CD69"), fsom$markers), 
                 xdim = 10, 
                 ydim = 10
)
fsom <- BuildMST(fsom)
nbclust <- 15
fsom[["metaclustering"]] <- as.factor(metaClustering_consensus(fsom$map$codes, k = nbclust))
PlotStars(fsom, backgroundValues = fsom$metaclustering)

fsom$heatmap <- pheatmap(scale(fsom.heat(fsom)), scale = "row")
fcs.explorer(fsom.dat(fsom), "CD4", "CD8A", 20000)

fsom <- nodes.to.meta_fSOM(fsom, grep("^7$|^8$", fsom$metaclustering))
fsom$heatmap <- pheatmap(scale(fsom.heat(fsom)), scale = "row")
fcs.explorer(fsom.dat(fsom), "CD4", "CD8A", 20000)

fsom <- nodes.to.meta_fSOM(fsom, grep("^1$|^2$", fsom$metaclustering))
fsom$heatmap <- pheatmap(scale(fsom.heat(fsom)), scale = "row")
fcs.explorer(fsom.dat(fsom), "CD4", "CD8A", 20000)

fsom$metaclustering.anno <- vector(mode = "list", length = length(levels(fsom$metaclustering)))
fsom$metaclustering.anno[[13]] <- "CD4+ CD69+"
fsom$metaclustering.anno[[12]] <- "CD8+ CD69+"

fsom_PASS2 <- fsom

total <- nrow(tmp[[1]])
cd3.total <- length(which(exprs(tmp[[1]])[, "meta.cluster_PASS1"] == 16))

x.nodes <- vector("integer", total)
x.nodes[cd3.index] <- fsom$map$mapping[, 1]

x.clusters <- vector("integer", total)
x.clusters[cd3.index] <- as.integer(fsom$metaclustering[fsom$map$mapping[, 1]])

tmp <- fsApply(tmp, function(i){
  i <- fr_append_cols(i, matrix(x.nodes, dimnames = list(NULL, "node_PASS2")))
  i <- fr_append_cols(i, matrix(x.clusters, dimnames = list(NULL, "meta.cluster_PASS2")))
})

write.FCS(tmp[[1]], "./temp/test2.fcs")
