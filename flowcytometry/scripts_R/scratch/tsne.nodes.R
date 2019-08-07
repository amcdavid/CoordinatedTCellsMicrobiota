## fSOM PASS3 mapping; map fSOM_CD4p warped files to rerun

## file paths: warped samples
tphe.warped.files <- list.files("./data_modified/TPHE_warping/", recursive = TRUE, full.names = TRUE, 
                                pattern = "._fda_WARPED_.fcs")
tphe.warped.files <- tphe.warped.files[grep("PERFORIN.FOXP3.CD185.CD31.GRZB", dirname(tphe.warped.files))]

hd <- tphe.warped.files[grep("HD", tphe.warped.files)]

trans.obj <- readRDS("./results_R/FCS.transforms/TPHE/asinh.transforms.TPHE.rds")
fSOM <- readRDS("./results_R/fSOMs/TPHE/2018-11-14.fSOM_PASS3_CD4p_rerun.rds")

fsom.mapping.concatenate <- function(fcs.paths, fsom.obj) {
  input.set <- fset.compensate(read.flowSet(fcs.paths, transformation = FALSE))
  input.set <- transform.set(input.set)
  input.set <- as(input.set, "flowFrame")
  exprs(input.set) <- exprs(input.set)[, -grep("Original", colnames(exprs(input.set)))]
  fSOM <- NewData(fsom.obj, input.set)
  fSOM
}

fSOM.fullV1 <- fsom.mapping.concatenate(tphe.warped.files[grep("RPRC053.*_1_", tphe.warped.files)], fSOM)
fSOM.fullV7 <- fsom.mapping.concatenate(tphe.warped.files[grep("RPRC053.*_7_", tphe.warped.files)], fSOM)
fSOM.fullV19 <- fsom.mapping.concatenate(tphe.warped.files[grep("RPRC053.*_19_", tphe.warped.files)], fSOM)

fSOM.preV1 <- fsom.mapping.concatenate(tphe.warped.files[grep("RPRC054.*_1_|PROP.*_1_", tphe.warped.files)], fSOM)
fSOM.preV7 <- fsom.mapping.concatenate(tphe.warped.files[grep("RPRC054.*_7_|PROP.*_7_", tphe.warped.files)], fSOM)
fSOM.preV19 <- fsom.mapping.concatenate(tphe.warped.files[grep("RPRC054.*_19_|PROP.*_19_", tphe.warped.files)], fSOM)

##readRDSs and generate median.by node
fsom.median.by.node <- function(path.to.fsom) {
  dat <- readRDS(path.to.fsom)
  dat <- fsom.dat(dat)
  dat <- cbind(aggregate(. ~ nodes, data = dat, FUN = median), node.counts = as.data.frame(table(dat$nodes))$Freq)
  dat
}

fsoms.to.read <- list.files("./temp/", pattern = "full|pre", full.names = TRUE)
median.by.node <- vector(mode = "list", length = length(fsoms.to.read))

for(i in seq(fsoms.to.read)) {
  dat <- fsom.median.by.node(fsoms.to.read[i])
  median.by.node[[i]] <- dat
  names(median.by.node)[i] <- sub(".rds", "", basename(fsoms.to.read[i]))
}

for(i in seq(median.by.node)) {
  median.by.node[[i]]$cohort <- unlist(strsplit(sub("fSOM.", "", names(median.by.node[i])), "V"))[1]
  median.by.node[[i]]$visit <- unlist(strsplit(sub("fSOM.", "", names(median.by.node[i])), "V"))[2]
}

saveRDS(median.by.node, "./temp/median.by.node.RDS")

median.by.node <- do.call(rbind.data.frame, readRDS("./temp/median.by.node.RDS"))

library(Rtsne)

dims.used <- c("CD57", "CD28", "FOXP3", "CD197", "CD185", "CD45RO", "CD127", "CD31")
set.seed(42)
nodes.tsne <- Rtsne(median.by.node[, dims.used], perplexity = 20, max_iter = 2000, verbose = TRUE, theta = 0.0)
# nodes.tsne <- Rtsne(aggregate(. ~ nodes, median.by.node[, -grep("cohort|visit", colnames(median.by.node))], median)[, dims.used], 
#                     perplexity = 30, max_iter = 3000, verbose = TRUE, theta = 0.5)

dat.nodes.tsne <- cbind(median.by.node, tsne.x = nodes.tsne$Y[, 1], tsne.y = nodes.tsne$Y[, 2])
mcluster.agg <- aggregate(. ~ MCluster, dat.nodes.tsne[, -grep("cohort|visit", colnames(dat.nodes.tsne))], median)
mc.df <- merge(dat.nodes.tsne[1:100, ], aggregate(cbind(median.x = tsne.x, median.y = tsne.y) ~ MCluster, dat.nodes.tsne, median), by = "MCluster")

ggplot(dat.nodes.tsne[1:100, ], aes(x = tsne.x, y = tsne.y, size = node.counts, color = factor(MCluster))) + 
  geom_point(show.legend = FALSE) +
  geom_segment(data = mc.df, aes(x = tsne.x, y = tsne.y, xend = median.x, yend = median.y), size = 0.4, show.legend = FALSE) +
  geom_text(data = mcluster.agg, inherit.aes = FALSE, aes(x = tsne.x, y = tsne.y, label = rownames(mcluster.agg)), nudge_x = 1.2, nudge_y = 1.2, size = 4.5) +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  labs(title = "T-cell Phenotyping Panel: tSNE embedding of 100 CD4+ fSOM nodes; 
       colored by meta-cluster assignment (nodes connected to meta-cluster centroid", align = "justify")

mcluster.df <- merge(aggregate(cbind(tsne.x, tsne.y) ~ MCluster, dat.nodes.tsne, median), aggregate(node.counts ~ MCluster, dat.nodes.tsne, sum), by = "MCluster")

ggplot(mcluster.df, aes(x = tsne.x, y = tsne.y, size = node.counts, color = factor(MCluster))) +
  geom_point() +
  scale_size(range = c(0.5, 20))
