library(umap); library(ggplot2); library(ggrepel)

counts.merged <- read.csv("./results_R/fSOMs/TPHE_ICS/Cluster.counts.ALLmerged.csv")
counts.merged$term.visit <- sub("_", " ", sub("V19", "12-months", sub("V7", "Discharge", sub("V1$", "Birth", counts.merged$term.visit))))
counts.merged$term.visit <- factor(counts.merged$term.visit)
counts.merged$term.visit <- factor(counts.merged$term.visit, levels(counts.merged$term.visit)[c(2,3,1,5,6,4)])
counts.merged$SubID <- as.character(counts.merged$SubID)
colnames(counts.merged)[grep("subject", colnames(counts.merged))] <- "SampleID"

ics <- readr::read_tsv("./results_R/MI2/ics_basic.txt")
ics$Subj <- forcats::fct_reorder(ics$Subject, ics$GAB)

tphe <- readr::read_tsv("./results_R/MI2/tphe_basic.txt")
tphe$Subj <- forcats::fct_reorder(tphe$Subject, tphe$GAB)

counts.merged.ist.ics <- merge(counts.merged, ics, by = "SampleID")
counts.merged.ist.tphe <- merge(counts.merged, tphe, by = "SampleID")

ist.tphe.umap <- umap(scale(counts.merged.ist.tphe[, grep("cluster", colnames(counts.merged.ist.tphe), ignore.case = T)]))
ist.ics.umap <- umap(scale(counts.merged.ist.ics[, grep("cluster", colnames(counts.merged.ist.ics), ignore.case = T)]))

counts.merged.ist.tphe.umap <- cbind(counts.merged.ist.tphe, 
                                data.frame(umap1 = ist.tphe.umap$layout[, 1], 
                                           umap2 = ist.tphe.umap$layout[, 2]))

counts.merged.ist.ics.umap <- cbind(counts.merged.ist.ics, 
                                data.frame(umap1 = ist.ics.umap$layout[, 1], 
                                           umap2 = ist.ics.umap$layout[, 2]))

####
tmp <- counts.merged.ist.tphe.umap
term.visit.centroids <- merge(tmp[grep("term.visit|umap", colnames(tmp))],
                              aggregate(cbind(umap.centroid.x = umap1, umap.centroid.y = umap2) ~ term.visit, tmp, median),
                              by = "term.visit")

ist.centroids <- merge(tmp[grep("Renamed_IST|umap", colnames(tmp))], 
                       aggregate(cbind(umap.centroid.x = umap1, umap.centroid.y = umap2) ~ Renamed_IST, tmp, median), 
                       by = "Renamed_IST")

ggplot(counts.merged.ist.tphe.umap, aes(umap1, umap2, color = term.visit)) + 
  geom_point() +
  geom_segment(data = term.visit.centroids, aes(xend = umap.centroid.x, yend = umap.centroid.y)) + 
  theme_void()

ggplot(counts.merged.ist.tphe.umap, aes(umap1, umap2, color = term.visit)) + 
  geom_point(show.legend = F) +
  geom_segment(data = term.visit.centroids, aes(xend = umap.centroid.x, yend = umap.centroid.y), alpha = 1.0, show.legend = F) +
  geom_label_repel(data = aggregate(cbind(umap.centroid.x = umap1, umap.centroid.y = umap2) ~ term.visit, tmp, median), 
                   aes(x = umap.centroid.x, y = umap.centroid.y, label = term.visit), inherit.aes = F) +
  theme_void()

ggplot(counts.merged.ist.ics.umap, aes(umap1, umap2, color = Renamed_IST)) + 
  geom_point(show.legend = F) +
  geom_segment(data = ist.centroids, aes(xend = umap.centroid.x, yend = umap.centroid.y), alpha = 1.0, show.legend = F) +
  geom_label_repel(data = aggregate(cbind(umap.centroid.x = umap1, umap.centroid.y = umap2) ~ Renamed_IST, tmp, median), 
                   aes(x = umap.centroid.x, y = umap.centroid.y, label = Renamed_IST), inherit.aes = F) +
  theme_void()

ggplot(counts.merged.ist.tphe.umap, aes(umap1, umap2, color = term.visit)) + 
  geom_point() +
  #geom_segment(data = term.visit.centroids, aes(xend = umap.centroid.x, yend = umap.centroid.y)) + 
  theme_void() +
  scale_color_discrete()
