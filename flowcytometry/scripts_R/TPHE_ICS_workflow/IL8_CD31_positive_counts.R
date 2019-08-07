library(ggplot2)
#ICS IL8
#identified 9, 10, 11, 12, 13 as IL8+ in CD4+CD69+

counts <- read.csv("./results_R/fSOMs/ICS/FINAL/ICS_fSOM.CD4cytokines_MetaClusters_Counts.csv")

CD4p.CD69p.IL8p <- data.frame(unique = as.character(counts$unique),
                              CD4p.CD69p.IL8p_Meta.Clusters = (counts$Meta.Cluster_9 + counts$Meta.Cluster_10 + counts$Meta.Cluster_11 + counts$Meta.Cluster_12 + counts$Meta.Cluster_13),
                              CD4p.CD69p.TOTAL = rowSums(counts[, grep("Cluster", colnames(counts))])
                              )

CD4p.CD69p.IL8p$CD4p.CD69p.IL8p_frequency <- CD4p.CD69p.IL8p$CD4p.CD69p.IL8p_Meta.Clusters / CD4p.CD69p.IL8p$CD4p.CD69p.TOTAL


#identified 8, 10, as IL8+ in CD8+CD69+

counts <- read.csv("./results_R/fSOMs/ICS/FINAL/ICS_fSOM.CD8cytokines_MetaClusters_Counts.csv")

CD8p.CD69p.IL8p <- data.frame(unique = as.character(counts$unique),
                              CD8p.CD69p.IL8p_Meta.Clusters = (counts$Meta.Cluster_8 + counts$Meta.Cluster_10),
                              CD8p.CD69p.TOTAL = rowSums(counts[, grep("Cluster", colnames(counts))])
)

CD8p.CD69p.IL8p$CD8p.CD69p.IL8p_frequency <- CD8p.CD69p.IL8p$CD8p.CD69p.IL8p_Meta.Clusters / CD8p.CD69p.IL8p$CD8p.CD69p.TOTAL

CD69p_CD4p.CD8p_IL8p <- merge(CD4p.CD69p.IL8p, CD8p.CD69p.IL8p, by = "unique")

ICS.meta <- readRDS("./data_source/ICS_SEB.meta.rds")

CD69p_CD4p.CD8p_IL8p_meta <- merge(CD69p_CD4p.CD8p_IL8p, ICS.meta, by = "unique")

ggplot(CD69p_CD4p.CD8p_IL8p_meta, aes(TermALIAS, CD4p.CD69p.IL8p_frequency)) +
  geom_jitter(width = 0.2) +
  facet_wrap("VisitALIAS")

##

#TPHE CD31
#identified 2, 13, 22, 25 as CD31+ in CD4+ TPHE

counts <- read.csv("./results_R/fSOMs/TPHE/FINAL/TPHE_fSOM.CD4_MetaClusters_Counts.csv")

CD4p.CD31p <- data.frame(unique = as.character(counts$unique),
                         CD4p.CD31p_Meta.Clusters = (counts$Meta.Cluster_2 + counts$Meta.Cluster_13 + counts$Meta.Cluster_22 + counts$Meta.Cluster_25),
                         CD4p.TOTAL = rowSums(counts[, grep("Cluster", colnames(counts))])
)

CD4p.CD31p$CD4p.CD31p_frequency <- CD4p.CD31p$CD4p.CD31p_Meta.Clusters / CD4p.CD31p$CD4p.TOTAL

#identified "2,3,6,7,9,11,14,15,16,18,19,20,21,23,24" as CD31+ in CD8+ TPHE
CD31p.clusters <- c(2,3,6,7,9,11,14,15,16,18,19,20,21,23,24)

counts <- read.csv("./results_R/fSOMs/TPHE/FINAL/TPHE_fSOM.CD8_MetaClusters_Counts.csv")

CD8p.CD31p <- data.frame(unique = as.character(counts$unique),
                         CD8p.CD31p_Meta.Clusters = rowSums(counts[, grep("Cluster", colnames(counts))][, CD31p.clusters]),
                         CD8p.TOTAL = rowSums(counts[, grep("Cluster", colnames(counts))])
)

CD8p.CD31p$CD8p.CD31p_frequency <- CD8p.CD31p$CD8p.CD31p_Meta.Clusters / CD8p.CD31p$CD8p.TOTAL

CD4p.CD8p_CD31p <- merge(CD4p.CD31p, CD8p.CD31p, by = "unique")

CD4p.CD8p_CD31p$dolSRC <- sub("_TPHE", "", CD4p.CD8p_CD31p$unique)

CD31.IL8.meta <- merge(CD4p.CD8p_CD31p, CD69p_CD4p.CD8p_IL8p_meta, by = "dolSRC")

ggplot(CD31.IL8.meta, aes(TermALIAS, CD4p.CD69p.IL8p_frequency)) +
  geom_jitter(width = 0.2) +
  facet_wrap("VisitALIAS")

ggplot(CD31.IL8.meta, aes(TermALIAS, CD8p.CD31p_frequency)) +
  geom_jitter(width = 0.2) +
  facet_wrap("VisitALIAS")

write.csv(CD31.IL8.meta, "./results_R/fSOMs/TPHE_ICS/CD31.IL8.csv", row.names = FALSE)
