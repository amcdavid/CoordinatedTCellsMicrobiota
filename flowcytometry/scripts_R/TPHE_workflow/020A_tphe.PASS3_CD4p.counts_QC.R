count.paths <- list.files("./data_modified/TPHE_warping/fSOM_PASS3_CD4p.clusters/counts/", full.names=TRUE, pattern=".csv")

raw.counts <- generate.counts(count.paths)

counts <- as.data.frame(t(apply(raw.counts, 1, function(i) i/sum(i)))) * 100 # sum100 (percent)

counts$unique <- sub(" ", "_", sub("_fda.*", "", row.names(counts)))
counts$unique[-grep("HD", counts$unique)] <- sub("_NGL.*", "", counts$unique[-grep("HD", counts$unique)])

TPHE_meta <- readRDS("./data_modified/TPHE_meta.RDS")

counts.HD <- counts[grep("HD", counts$unique), ]
counts.HD$ExperimentALIAS <- unlist(lapply(seq(counts.HD$unique), function(i) strsplit(counts.HD$unique, "_")[[i]][3]))

meta.explorer(counts.HD)


counts.meta <- merge(counts, TPHE_meta, by = "unique")
counts.meta$random <- runif(length(counts.meta$unique), 1, length(counts.meta$unique))

meta.explorer(counts.meta, "ExperimentVALUE", "Meta.Cluster_29")

bad.CD28 <- counts.meta$unique[order(counts.meta$Meta.Cluster_30, decreasing = TRUE)[1:4]]


ggord(prcomp(counts.meta[,grep("Meta.Cluster",colnames(counts.meta))]), counts.meta$VisitALIAS, ellipse = FALSE, arrow = NULL, txt = NULL)
ggord(prcomp(counts.meta[,grep("Meta.Cluster",colnames(counts.meta))]), counts.meta$VisitTermALIAS, ellipse = TRUE, ellipse_pro = .5, arrow = NULL, txt = NULL)
ggord(prcomp(counts.meta[,grep("Meta.Cluster",colnames(counts.meta))]), counts.meta$ExperimentALIAS, ellipse = FALSE, arrow = NULL, txt = NULL)
biplot(prcomp(counts.meta[,grep("Meta.Cluster",colnames(counts.meta))]))
