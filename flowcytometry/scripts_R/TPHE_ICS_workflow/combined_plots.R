library(ggplot2); library(ggrepel); library(viridis)

counts.merged <- read.csv("./results_R/fSOMs/TPHE_ICS/Cluster.counts.ALLmerged.csv")
counts.merged$term.visit <- sub("_", " ", sub("V19", "12-months", sub("V7", "Discharge", sub("V1$", "Birth", counts.merged$term.visit))))
counts.merged$term.visit <- factor(counts.merged$term.visit)
counts.merged$term.visit <- factor(counts.merged$term.visit, levels(counts.merged$term.visit)[c(2,3,1,5,6,4)])
counts.merged$SubID <- as.character(counts.merged$SubID)
colnames(counts.merged)[grep("subject", colnames(counts.merged))] <- "SampleID"

ggplot(counts.merged, aes(tsne.y, (tsne.x * -1), color = term.visit)) + 
  geom_point(show.legend = T, size = 2) +
  geom_segment(aes(xend = centroid.y, yend = (centroid.x * -1)), size = 0.7, alpha = 0.3, show.legend = FALSE) +
  # geom_text(data = aggregate(cbind(centroid.x = tsne.x, centroid.y = tsne.y) ~ term.visit, counts.merged, median), 
  #           aes(x = centroid.x, y = centroid.y, label = term.visit), nudge_x = 2, nudge_y = 2, size = 6) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~ExperimentALIAS)
# labs(title = "tSNE embedding of 385 CD4+ fSOM nodes; 
#      colored by meta-cluster assignment (nodes connected to group centroid)", align = "justify")


ggplot(subset(counts.merged, counts.merged$ExperimentVALUE == 59), aes(tsne.x, tsne.y, color = SubID)) + 
  geom_point(show.legend = FALSE, size = 2) +
  geom_text_repel(aes(label = subject), size = 4, show.legend = FALSE) +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5))
# labs(title = "tSNE embedding of 385 CD4+ fSOM nodes; 
#      colored by meta-cluster assignment (nodes connected to group centroid)", align = "justify")

##color-by these clusters:
ics.cd4.clusters <- sapply(c(1, 2, 4, 5, 13, 17, 18), function(i) paste0("ICS.CD4_Meta.Cluster_", i))
ics.cd8.clusters <- sapply(c(3, 4, 8, 11), function(i) paste0("ICS.CD8_Meta.Cluster_", i))
tphe.cd4.clusters <- sapply(c(2, 6, 8, 16, 23), function(i) paste0("TPHE.CD4_Meta.Cluster_", i))
tphe.cd8.clusters <- sapply(c(1, 2, 6, 15, 22, 23), function(i) paste0("TPHE.CD8_Meta.Cluster_", i))
clusters.of.interest <- c(ics.cd4.clusters, ics.cd8.clusters, tphe.cd4.clusters, tphe.cd8.clusters)

ggplot(counts.merged, aes(ICS.CD4_Meta.Cluster_18)) + geom_density(bw = "SJ")

i <- clusters.of.interest[1]
ggplot(counts.merged, aes_string("tsne.y", "(tsne.x * -1)", color = i, size = i)) +
  scale_size_area() +
  geom_point(show.legend = T) +
  scale_color_viridis(option = "B") +
  geom_segment(aes(xend = centroid.y, yend = (centroid.x * -1)), size = 0.7, alpha = 0.3, show.legend = FALSE) +
  # geom_text(data = aggregate(cbind(centroid.x = tsne.x, centroid.y = tsne.y) ~ term.visit, counts.merged, median), 
  #           aes(x = centroid.x, y = centroid.y, label = term.visit), nudge_x = 2, nudge_y = 2, size = 6) +
  theme_dark() +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5))

###
plots <- list()
for(i in seq(clusters.of.interest)){
  plots[[i]] <- ggplot(counts.merged, aes_string("tsne.y", "(tsne.x * -1)", color = clusters.of.interest[i], size = clusters.of.interest[i])) +
    scale_size_area() +
    geom_point(show.legend = T) +
    scale_color_viridis(option = "B") +
    geom_segment(aes(xend = centroid.y, yend = (centroid.x * -1)), size = 0.7, alpha = 0.3, show.legend = FALSE) +
    # geom_text(data = aggregate(cbind(centroid.x = tsne.x, centroid.y = tsne.y) ~ term.visit, counts.merged, median), 
    #           aes(x = centroid.x, y = centroid.y, label = term.visit), nudge_x = 2, nudge_y = 2, size = 6) +
    theme_dark() +
    theme(panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size = 12),
          plot.title = element_text(hjust = 0.5),
          legend.position = c(0.96, 0.2),
          legend.background = element_rect(fill = "dark grey")) +
    labs(title = paste0(clusters.of.interest[i]))
}

pdf("clusters.of.interest_overlay.pdf", width = 11, height = 8.5)
plots
dev.off()
####

ics <- readr::read_tsv("./results_R/MI2/ics_basic.txt")
ics$Subj <- forcats::fct_reorder(ics$Subject, ics$GAB)

tphe <- readr::read_tsv("./results_R/MI2/tphe_basic.txt")
tphe$Subj <- forcats::fct_reorder(tphe$Subject, tphe$GAB)

tmp.ics <- merge(counts.merged, ics, by = "SampleID")
tmp.tphe <- merge(counts.merged, tphe, by = "SampleID")

ggplot(tmp.tphe, aes(CGA, Subj, color = Renamed_IST)) + 
  geom_point() + 
  theme(axis.text.y=element_blank()) +
  facet_wrap(c("TermALIAS", "Visit"))

ggplot(tmp.tphe, aes(CGA, GAatBirth, color = Renamed_IST)) + 
  scale_color_discrete(name = "Immune State Type") +
  geom_point(size = 2) + 
  facet_wrap(c("TermALIAS", "Visit")) +
  theme(strip.text = element_text(size = 14)) +
  xlab("Corrected Gestational Age (weeks)") +
  ylab("Gestational Age at Birth (weeks)")


ist.tphe.facets <- ggplot(tmp.tphe, aes(tsne.y, (tsne.x * -1), color = term.visit)) +
  geom_point(show.legend = T, size = 2) +
  scale_color_discrete() +
  geom_segment(aes(xend = centroid.y, yend = (centroid.x * -1)), size = 0.7, alpha = 0.3, show.legend = FALSE) +
  # geom_text(data = aggregate(cbind(centroid.x = tsne.x, centroid.y = tsne.y) ~ term.visit, counts.merged, median), 
  #           aes(x = centroid.x, y = centroid.y, label = term.visit), nudge_x = 2, nudge_y = 2, size = 6) +
  theme_dark() +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~Renamed_IST)

ist.ics.facets <- ggplot(tmp.ics, aes(tsne.y, (tsne.x * -1), color = term.visit)) +
  geom_point(show.legend = T, size = 2) +
  scale_color_discrete() +
  geom_segment(aes(xend = centroid.y, yend = (centroid.x * -1)), size = 0.7, alpha = 0.3, show.legend = FALSE) +
  # geom_text(data = aggregate(cbind(centroid.x = tsne.x, centroid.y = tsne.y) ~ term.visit, counts.merged, median), 
  #           aes(x = centroid.x, y = centroid.y, label = term.visit), nudge_x = 2, nudge_y = 2, size = 6) +
  theme_dark() +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~Renamed_IST)

bysample.map <- ggplot(tmp.ics, aes(tsne.y, (tsne.x * -1), color = term.visit)) +
  geom_point(show.legend = T, size = 2) +
  scale_color_discrete() +
  geom_segment(aes(xend = centroid.y, yend = (centroid.x * -1)), size = 0.7, alpha = 0.3, show.legend = FALSE) +
  # geom_text(data = aggregate(cbind(centroid.x = tsne.x, centroid.y = tsne.y) ~ term.visit, counts.merged, median), 
  #           aes(x = centroid.x, y = centroid.y, label = term.visit), nudge_x = 2, nudge_y = 2, size = 6) +
  theme_dark() +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5))

pdf("bysample_ISToverlay.pdf", width = 11, height = 8.5)
bysample.map
ist.tphe.facets
bysample.map
ist.ics.facets
dev.off()

ggplot(tmp.tphe, aes(tsne.y, (tsne.x * -1), color = term.visit)) +
  geom_point(show.legend = T, size = 2) +
  scale_color_discrete() +
  geom_segment(aes(xend = centroid.y, yend = (centroid.x * -1)), size = 0.7, alpha = 0.3, show.legend = FALSE) +
  # geom_text(data = aggregate(cbind(centroid.x = tsne.x, centroid.y = tsne.y) ~ term.visit, counts.merged, median), 
  #           aes(x = centroid.x, y = centroid.y, label = term.visit), nudge_x = 2, nudge_y = 2, size = 6) +
  theme_dark() +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~MOD)

ics.exp.counts <- as.data.frame(sapply(sort(unique(tmp.ics$Renamed_IST)), function(i) table(tmp.ics$ExperimentALIAS[tmp.ics$Renamed_IST == paste(i)])))
tphe.exp.counts <- as.data.frame(sapply(sort(unique(tmp.tphe$Renamed_IST)), function(i) table(tmp.tphe$ExperimentALIAS[tmp.tphe$Renamed_IST == paste(i)])))

plots.ics <- list()
for(i in seq(colnames(ics.exp.counts))) {
  plots.ics[[i]] <- ggplot(ics.exp.counts, aes_string(x = factor(rownames(ics.exp.counts)), y = paste0("ICS_", i))) + geom_point(size = 3)
}

plots.tphe <- list()
for(i in seq(colnames(tphe.exp.counts))) {
  plots.tphe[[i]] <- ggplot(tphe.exp.counts, aes_string(x = factor(rownames(tphe.exp.counts)), y = paste0("TPHE_", i))) + geom_point(size = 3)
}

##CMV
cmv.status <- read.csv("./data_source/TPHE2_with.spillover/CMVstatus.edit.csv")
cmv.status <- cmv.status[, c("Patient.ID", "Positive")]
colnames(cmv.status) <- c("ParticipantID", "CMV.status")       

cmv.status$ParticipantID %in% tmp.tphe$ParticipantID

tmp.tphe.cmv <- merge(tmp.tphe, cmv.status, by = "ParticipantID")

ggplot(subset(tmp.tphe.cmv, tmp.tphe.cmv$Visit == "OneYear"), aes(Renamed_IST, CMV.status)) + geom_point()

ggplot(subset(tmp.tphe.cmv, tmp.tphe.cmv$CMV.status == 1), aes(TermALIAS, TPHE.CD8_Meta.Cluster_22, color = Renamed_IST)) + geom_point()

##
ist.centroids <- merge(tmp.tphe[grep("Renamed_IST|tsne", colnames(tmp.tphe))], 
                   aggregate(cbind(centroid.x = tsne.x, centroid.y = tsne.y) ~ Renamed_IST, tmp.tphe, median), 
                   by = "Renamed_IST")

ggplot(tmp.tphe, aes(tsne.y, (tsne.x * -1), color = Renamed_IST)) + 
  geom_point() +
  geom_segment(data = ist.centroids, aes(xend = centroid.y, yend = (centroid.x * -1)), size = 0.7, alpha = 0.3, show.legend = FALSE) +
  stat_ellipse(level = .9)

ist.centroids <- merge(tmp.ics[grep("Renamed_IST|tsne", colnames(tmp.ics))], 
                       aggregate(cbind(centroid.x = tsne.x, centroid.y = tsne.y) ~ Renamed_IST, tmp.ics, median), 
                       by = "Renamed_IST")

ggplot(tmp.ics, aes(tsne.y, (tsne.x * -1), color = Renamed_IST)) + 
  geom_point() +
  geom_segment(data = ist.centroids, aes(xend = centroid.y, yend = (centroid.x * -1)), size = 0.7, alpha = 0.3, show.legend = FALSE) +
  stat_ellipse()

##

term.visit.centroids <- merge(tmp.tphe[grep("term.visit|umap", colnames(tmp.tphe))],
                              aggregate(cbind(umap.centroid.x = umap.x, umap.centroid.y = umap.y) ~ term.visit, tmp.tphe, median),
                              by = "term.visit")

ist.centroids <- merge(tmp.tphe[grep("Renamed_IST|umap", colnames(tmp.tphe))], 
                       aggregate(cbind(umap.centroid.x = umap.x, umap.centroid.y = umap.y) ~ Renamed_IST, tmp.tphe, median), 
                       by = "Renamed_IST")

ggplot(tmp.tphe, aes(umap.x, umap.y, color = Renamed_IST)) + 
  geom_point() +
  geom_segment(data = ist.centroids, aes(xend = umap.centroid.x, yend = umap.centroid.y)) +
  facet_wrap(~Renamed_IST)

ggplot(tmp.tphe, aes(umap.x, umap.y, color = term.visit)) + 
  geom_point() +
  geom_segment(data = term.visit.centroids, aes(xend = umap.centroid.x, yend = umap.centroid.y)) +
  stat_ellipse(level = .8)

