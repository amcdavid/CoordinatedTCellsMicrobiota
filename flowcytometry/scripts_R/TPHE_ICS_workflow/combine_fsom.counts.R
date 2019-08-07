## list counts datasets
library(Rtsne)
library(ggplot2)

counts <- list(CD4.TPHE = read.csv("./results_R/fSOMs/TPHE/FINAL/TPHE_fSOM.CD4_MetaClusters_Counts.csv", stringsAsFactors = F),
               CD4.ICS  = read.csv("./results_R/fSOMs/ICS/FINAL/ICS_fSOM.CD4cytokines_MetaClusters_Counts.csv", stringsAsFactors = F),
               CD8.TPHE = read.csv("./results_R/fSOMs/TPHE/FINAL/TPHE_fSOM.CD8_MetaClusters_Counts.csv", stringsAsFactors = F),
               CD8.ICS  = read.csv("./results_R/fSOMs/ICS/FINAL/ICS_fSOM.CD8cytokines_MetaClusters_Counts.csv", stringsAsFactors = F)
               )

## rename clusters to be unique; create subject column
counts <- sapply(names(counts), function(i) {
  colnames(counts[[i]])[grep("Cluster", colnames(counts[[i]]))] <- paste0(i, "_", colnames(counts[[i]])[grep("Cluster", colnames(counts[[i]]))])
  counts[[i]]$subject <- sapply(strsplit(counts[[i]]$unique, "_"), function(x) paste((x)[1:2], collapse = "_"))
  return(counts[[i]])
})

## raw counts to frequency
counts <- sapply(names(counts), function(i) {
  counts[[i]][, grep("Cluster", colnames(counts[[i]]))] <- as.data.frame(t(apply(counts[[i]][grep("Cluster", colnames(counts[[i]]))], 1, function(x) x/sum(x) * 100)))
  return(counts[[i]])
})

## merge
counts.merged <- Reduce(function(x,y) merge(x = x, y = y, by = "subject"), 
                        list(counts$CD4.TPHE, counts$CD4.ICS, counts$CD8.TPHE, counts$CD8.ICS))
counts.merged <- counts.merged[-grep("unique", colnames(counts.merged))]
rowSums(counts.merged[grep("Cluster", colnames(counts.merged))]) # equals 4; sum 1 for each of the 4 datasets


## meta values
{
counts.merged$term[grep("RPRC054|PROP", counts.merged$subject)] <- "Pre-term"
counts.merged$term[grep("RPRC053", counts.merged$subject)] <- "Full-term"
counts.merged$term <- as.factor(counts.merged$term)

counts.merged$visit[grep("_1$", counts.merged$subject)] <- "V1"
counts.merged$visit[grep("_19", counts.merged$subject)] <- "V19"
counts.merged$visit[grep("_7", counts.merged$subject)] <- "V7"
counts.merged$visit <- as.factor(counts.merged$visit)

counts.merged$term.visit <- factor(paste(counts.merged$term, counts.merged$visit, sep = "_"),
                                      levels = c("Full-term_V1", 
                                                 "Full-term_V7", 
                                                 "Full-term_V19",
                                                 "Pre-term_V1",
                                                 "Pre-term_V7",
                                                 "Pre-term_V19"))
}

## z-score
counts.merged[grep("Cluster", colnames(counts.merged))] <- sapply(grep("Cluster", colnames(counts.merged), value = TRUE), 
                                                                  function(i) (counts.merged[, i] - mean(counts.merged[, i])) / sd(counts.merged[, i]))
## dimensional reduction
tmp <- Rtsne(scale(counts.merged[grep("Cluster", colnames(counts.merged))]), verbose = TRUE, pca = FALSE, theta = 0.5, perplexity = 40, max_iter = 1000)

counts.merged$tsne.x  <- tmp$Y[, 1]
counts.merged$tsne.y  <- tmp$Y[, 2]

ggplot(counts.merged, aes(tsne.x, tsne.y, color = term.visit)) + 
  geom_point(size = 3) + 
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())

centroids <- merge(counts.merged[grep("term.visit|tsne", colnames(counts.merged))], 
                   aggregate(cbind(centroid.x = tsne.x, centroid.y = tsne.y) ~ term.visit, counts.merged, median), 
                   by = "term.visit")

ggplot(counts.merged, aes(tsne.x, tsne.y, color = term.visit)) + 
  geom_point(show.legend = TRUE, size = 2) +
  geom_segment(data = centroids, aes(xend = centroid.x, yend = centroid.y), size = 0.5, alpha = 0.3, show.legend = FALSE) +
  # geom_text(data = aggregate(cbind(centroid.x = tsne.x, centroid.y = tsne.y) ~ term.visit, counts.merged, median), 
  #           aes(x = centroid.x, y = centroid.y, label = term.visit), nudge_x = 2, nudge_y = 2, size = 6) +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5))
# labs(title = "tSNE embedding of 385 CD4+ fSOM nodes; 
#      colored by meta-cluster assignment (nodes connected to group centroid)", align = "justify")





meta <- readRDS("./data_modified/TPHE_meta.RDS")
colnames(meta)[grep("dolSRC", colnames(meta))] <- "subject"

counts.merged.meta <- merge(counts.merged, meta, by = "subject")

write.csv(counts.merged.meta, "./results_R/fSOMs/TPHE_ICS/counts.merged.meta.csv")
##
ggplot(counts.merged.meta, aes(tsne.x, tsne.y, color = term.visit)) + 
  geom_point(size = 3) + 
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())

ggplot(counts.merged.meta, aes(tsne.x, tsne.y, color = term.visit)) + 
  geom_point(show.legend = TRUE, size = 2) +
  geom_segment(aes(xend = centroid.x, yend = centroid.y), size = 0.5, alpha = 0.3, show.legend = FALSE) +
  # geom_text(data = aggregate(cbind(centroid.x = tsne.x, centroid.y = tsne.y) ~ term.visit, counts.merged, median), 
  #           aes(x = centroid.x, y = centroid.y, label = term.visit), nudge_x = 2, nudge_y = 2, size = 6) +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5))
  # labs(title = "tSNE embedding of 385 CD4+ fSOM nodes; 
  #      colored by meta-cluster assignment (nodes connected to group centroid)", align = "justify")
