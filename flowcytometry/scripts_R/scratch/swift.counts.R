####start build
##paths to merged median cluster counts
count.paths <- dir(path = "data_modified/SWIFT_results/NGL092_V19.SEB/SWIFT_OUTPUT/", 
                                       full.names = TRUE, 
                                       recursive = TRUE, 
                                       pattern = "Merge_Cluster_Sizes.txt")

##build counts list
counts <- vector(mode = "list")
counts$input <- lapply(seq(count.paths), function(i) fread(count.paths[i]))
names(counts$input) <- gsub(".fcs.*", "", basename(count.paths))

##add specific cluster id - 'cid' - indicies to list
counts$index$CD4pCD69p <- as.numeric(unlist(
  data.frame(fread("data_modified/SWIFT_results/NGL092_V19.SEB/Files_generated/CD4+_CD69+_CIDs.csv"))))

saveRDS(counts, "data_modified/SWIFT_results/NGL092_V19.SEB/NGL092_V19.SEB.merged.cluster.counts.RDS")  
####end build

####start with built object
counts <- readRDS("data_modified/SWIFT_results/NGL092_V19.SEB/NGL092_V19.SEB.merged.cluster.counts.RDS")

## index to parent population; sum 1 or 100
counts.frame <- as.data.frame(t(mapply(`[[`,counts$input, 2)))[, counts$index$CD4pCD69p]
colnames(counts.frame) <- gsub("V", "C.", colnames(counts.frame))
counts.frame <- counts.frame/rowSums(counts.frame) * 100
counts.frame$unique <- rownames(counts.frame)

counts.meta <- merge(counts.frame, readRDS("data_source/ICS_SEB.meta.rds"), by = "unique")

library(ggrepel)

ggplot(counts.meta, aes(TermALIAS, C.100)) + 
  geom_jitter(width = 0.1) + 
  facet_wrap(~VisitALIAS) +
  theme(axis.text.x = element_blank(),
        axis.title.y = element_text(size = 14, vjust = 2.5),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 14)) + 
  geom_text_repel(data = subset(counts.meta, C.100 > 0.1e-03), inherit.aes = TRUE,
                  aes(label = ExperimentALIAS, color = ExperimentALIAS), force = 10, show.legend = FALSE) +
  labs(x = "Pre-term VS Full-Term",
       y = paste("% of Cluster.100 in CD4+ CD69+"))

ggplot(counts.meta, aes(TermALIAS, C.87)) + 
  geom_jitter(width = 0.1) + 
  facet_wrap(~VisitALIAS) +
  theme(axis.text.x = element_blank(),
        axis.title.y = element_text(size = 14, vjust = 2.5),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 14)) + 
  geom_text_repel(data = subset(counts.meta, C.87 > 0.1e-03), inherit.aes = TRUE,
                  aes(label = ExperimentALIAS, color = ExperimentALIAS), force = 10, show.legend = FALSE) +
  labs(x = "Pre-term VS Full-Term",
       y = paste("% of Cluster.87 in CD4+ CD69+"))

ggplot(counts.meta, aes(TermALIAS, C.86)) + 
  geom_jitter(width = 0.1) + 
  facet_wrap(~VisitALIAS) +
  theme(axis.text.x = element_blank(),
        axis.title.y = element_text(size = 14, vjust = 2.5),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 14)) + 
  geom_text_repel(data = subset(counts.meta, C.86 > 0.1e-03), inherit.aes = TRUE,
                  aes(label = ExperimentALIAS, color = ExperimentALIAS), force = 10, show.legend = FALSE) +
  labs(x = "Pre-term VS Full-Term",
       y = paste("% of Cluster.86 in CD4+ CD69+"))

ggplot(counts.meta, aes(TermALIAS, C.166, color = ExperimentALIAS)) + 
  geom_jitter(width = 0.1) + 
  facet_wrap(~VisitALIAS) +
  theme(axis.text.x = element_blank(),
        axis.title.y = element_text(size = 14, vjust = 2.5),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 14))

meta.explorer(counts.meta)
saveRDS(counts.meta, "data_modified/SWIFT_results/NGL092_V19.SEB/NGL092_V19.SEB.merged.cluster.counts.meta.RDS")
