library(data.table)
library(shiny)
library(ggplot2)
library(viridis)

count.paths <- list.files("./data_modified/TPHE_warping/fSOM_PASS3_CD8p.clusters/counts/", full.names=TRUE, pattern=".csv")

raw.counts <- generate.counts(count.paths)
row.names(raw.counts) <- sub(" ", "_", sub("_fda.*", "", row.names(raw.counts)))
raw.counts <- raw.counts[-grep("HD", row.names(raw.counts)), ]
row.names(raw.counts) <- sub("_NGL.*", "", row.names(raw.counts))

counts <- as.data.frame(t(apply(raw.counts, 1, function(i) i/sum(i)))) * 100 # sum100 (percent)

counts$unique <- sub(" ", "_", sub("_fda.*", "", row.names(counts)))
counts$unique[grep("HD0189_TPHE_NGL091", counts$unique)] <- "HD0189_TPHE_NGL091"
counts$unique[-grep("HD", counts$unique)] <- sub("_NGL.*", "", counts$unique[-grep("HD", counts$unique)])

counts.HD <- counts[grep("HD", counts$unique), ]
counts.HD$ExperimentALIAS <- unlist(lapply(seq(counts.HD$unique), function(i) strsplit(counts.HD$unique, "_")[[i]][3]))


TPHE_meta <- readRDS("./data_modified/TPHE_meta.RDS")

counts.meta <- merge(counts, TPHE_meta, by = "unique")
counts.meta$random <- runif(length(counts.meta$unique), 1, length(counts.meta$unique))

meta.explorer(counts.meta, "ExperimentVALUE", "Meta.Cluster_24")


library(ggord)
ggord(prcomp(scale(counts.meta[,grep("Meta.Cluster",colnames(counts.meta))])), counts.meta$VisitALIAS, ellipse = FALSE, arrow = NULL, txt = NULL)
ggord(prcomp(counts.meta[,grep("Meta.Cluster",colnames(counts.meta))]), counts.meta$VisitTermALIAS, ellipse = TRUE, ellipse_pro = .5, arrow = NULL, txt = NULL)
ggord(prcomp(counts.meta[,grep("Meta.Cluster",colnames(counts.meta))]), counts.meta$ExperimentALIAS, ellipse = FALSE, arrow = NULL, txt = NULL)
biplot(prcomp(counts.meta[,grep("Meta.Cluster",colnames(counts.meta))]))
