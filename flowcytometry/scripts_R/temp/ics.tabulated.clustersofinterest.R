ics.cluster.counts <- list(cd4 = read.csv("./results_R/fSOMs/ICS/FINAL/ICS_fSOM.CD4cytokines_MetaClusters_Counts.csv"),
                           cd8 = read.csv("./results_R/fSOMs/ICS/FINAL/ICS_fSOM.CD8cytokines_MetaClusters_Counts.csv"))

ics.cluster.counts <- sapply(names(ics.cluster.counts), function(i) {
  ics.cluster.counts[[i]]$total <- rowSums(ics.cluster.counts[[i]][, grep("cluster", colnames(ics.cluster.counts[[i]]), ignore.case = T)])
  ics.cluster.counts[[i]]$subject.id <- as.factor(sub("_SEB", "", ics.cluster.counts[[i]]$unique))
  ics.cluster.counts[[i]]$visit.timepoint <- as.factor(sub("_SEB", "", sub("^[^_]*_", "", ics.cluster.counts[[i]]$unique)))
  ics.cluster.counts[[i]]$stimulation <- as.factor("SEB")
  ics.cluster.counts[[i]]
})

clusters.of.interest <- list(cd4 = list(IFNg = c(5),
                                        IL4 = c(4,7,8),
                                        IL8 = c(9,10,11,12,13),
                                        IL17 = c(1)),
                             cd8 = list(CD107a = c(2,3,5),
                                        IFNg = c(2,3,4,5),
                                        TNFa = c(2,3,6,7)))

clusters.of.interest$cd4 <- sapply(names(clusters.of.interest[["cd4"]]), function(i) paste0("Meta.Cluster_", clusters.of.interest[["cd4"]][[i]]))
clusters.of.interest$cd8 <- sapply(names(clusters.of.interest[["cd8"]]), function(i) paste0("Meta.Cluster_", clusters.of.interest[["cd8"]][[i]]))


ics.cluster.counts$cd4$IFNg.cluster.proportion <- (ics.cluster.counts$cd4[, clusters.of.interest$cd4$IFNg]/ics.cluster.counts$cd4$total) * 100
ics.cluster.counts$cd4$IL4.cluster.proportion <- (rowSums(ics.cluster.counts$cd4[, clusters.of.interest$cd4$IL4])/ics.cluster.counts$cd4$total) * 100
ics.cluster.counts$cd4$IL8.cluster.proportion <- (rowSums(ics.cluster.counts$cd4[, clusters.of.interest$cd4$IL8])/ics.cluster.counts$cd4$total) * 100
ics.cluster.counts$cd4$IL17.cluster.proportion <- (ics.cluster.counts$cd4[, clusters.of.interest$cd4$IL17]/ics.cluster.counts$cd4$total) * 100

ics.cluster.counts$cd8$CD107a.cluster.proportion <- (rowSums(ics.cluster.counts$cd8[, clusters.of.interest$cd8$CD107a])/ics.cluster.counts$cd8$total) * 100
ics.cluster.counts$cd8$IFNg.cluster.proportion <- (rowSums(ics.cluster.counts$cd8[, clusters.of.interest$cd8$IFNg])/ics.cluster.counts$cd8$total) * 100
ics.cluster.counts$cd8$TNFa.cluster.proportion <- (rowSums(ics.cluster.counts$cd8[, clusters.of.interest$cd8$TNFa])/ics.cluster.counts$cd8$total) * 100

write.csv(ics.cluster.counts$cd4, "./results_R/fSOMs/ICS/FINAL/ICS_fSOM.CD4cytokines_MetaClusters_Counts_EDIT.csv", row.names = FALSE)
write.csv(ics.cluster.counts$cd8, "./results_R/fSOMs/ICS/FINAL/ICS_fSOM.CD8cytokines_MetaClusters_Counts_EDIT.csv", row.names = FALSE)
