TPHE.CD4p.clusters <- readRDS("./results_R/fSOMs/TPHE/TPHE.CD4pos.unified.RDS")

counts <- TPHE.CD4p.clusters$counts.meta

ggplot(subset(counts, counts$TermALIAS == "Full-Term"), aes(x = VisitALIAS, y = Meta.Cluster_25, group = ParticipantID)) + geom_line() + geom_point()
ggplot(subset(counts, counts$TermALIAS == "Pre-Term"), aes(x = VisitALIAS, y = Meta.Cluster_25, group = ParticipantID)) + geom_line() + geom_point()

ggplot(counts, aes(x = VisitALIAS, y = Meta.Cluster_25, group = TermALIAS, color = TermALIAS)) + geom_line() + geom_point()
ggplot(counts, aes(x = VisitALIAS, y = Meta.Cluster_7, group = SubID)) + geom_line() + geom_point() + facet_wrap("TermALIAS")


plots <- vector(mode = "list", length = 27)
for(i in seq(1:27)) {
  plots[[i]] <- ggplot(counts, aes_string(x = "VisitALIAS", y = paste0("Meta.Cluster_", i), group = "SubID")) + geom_line() + geom_point() + facet_wrap("TermALIAS")
}
