library(data.table); library(shiny); library(ggplot2); library(viridis); library(ggord)

count.paths <- list.files("./data_modified/TPHE_warping/fSOM_PASS3_CD4p_rerun.clusters/counts/", full.names=TRUE, pattern=".csv")

raw.counts <- generate.counts(count.paths)

counts <- as.data.frame(t(apply(raw.counts, 1, function(i) i/sum(i)))) * 100 # sum100 (percent)

counts$unique <- sub(" ", "_", sub("_fda.*", "", row.names(counts)))
counts$unique[-grep("HD", counts$unique)] <- sub("_NGL.*", "", counts$unique[-grep("HD", counts$unique)])

TPHE_meta <- readRDS("./data_modified/TPHE_meta.RDS")

counts.HD <- counts[grep("HD", counts$unique), ]
counts.HD$donor <- sub("_TPHE.*", "", counts.HD$unique)
counts.HD$ExperimentALIAS <- unlist(lapply(seq(counts.HD$unique), function(i) strsplit(counts.HD$unique, "_")[[i]][3]))
counts.HD$ExperimentVALUE <- sub("NGL0", "", counts.HD$ExperimentALIAS)

meta.explorer(counts.HD, "ExperimentVALUE", "Meta.Cluster_1")

ggord(prcomp(counts.HD[,grep("Meta.Cluster", colnames(counts.HD))]), counts.HD$donor, ellipse = FALSE, arrow = NULL, txt = NULL)
ggord(prcomp(counts.HD[,grep("Meta.Cluster", colnames(counts.HD))]), counts.HD$ExperimentVALUE, ellipse = FALSE, arrow = NULL, txt = NULL)

##

counts.meta <- merge(counts, TPHE_meta, by = "unique")
counts.meta$random <- runif(length(counts.meta$unique), 1, length(counts.meta$unique))

meta.explorer(counts.meta, "ExperimentVALUE", "Meta.Cluster_7")

ggord(prcomp(counts.meta[,grep("Meta.Cluster",colnames(counts.meta))]), counts.meta$VisitALIAS, ellipse = FALSE, arrow = NULL, txt = NULL)
ggord(prcomp(counts.meta[,grep("Meta.Cluster",colnames(counts.meta))]), counts.meta$VisitTermALIAS, ellipse = TRUE, ellipse_pro = .85, arrow = NULL, txt = NULL)
ggord(prcomp(counts.meta[,grep("Meta.Cluster",colnames(counts.meta))]), counts.meta$ExperimentALIAS, ellipse = FALSE, arrow = NULL, txt = NULL)
biplot(prcomp(counts.meta[,grep("Meta.Cluster",colnames(counts.meta))]))

