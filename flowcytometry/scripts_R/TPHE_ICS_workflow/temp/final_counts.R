counts <- list(tphe.cd4 = read.csv("./results_R/fSOMs/TPHE/FINAL/TPHE_fSOM.CD4_MetaClusters_Counts.csv", stringsAsFactors = F),
               ics.cd4  = read.csv("./results_R/fSOMs/ICS/FINAL/ICS_fSOM.CD4cytokines_MetaClusters_Counts.csv", stringsAsFactors = F),
               tphe.cd8 = read.csv("./results_R/fSOMs/TPHE/FINAL/TPHE_fSOM.CD8_MetaClusters_Counts.csv", stringsAsFactors = F),
               ics.cd8  = read.csv("./results_R/fSOMs/ICS/FINAL/ICS_fSOM.CD8cytokines_MetaClusters_Counts.csv", stringsAsFactors = F)
)

## rename clusters to be unique; create subject column
counts <- sapply(names(counts), function(i) {
  colnames(counts[[i]])[grep("Cluster", colnames(counts[[i]]))] <- paste0(i, "_", colnames(counts[[i]])[grep("Cluster", colnames(counts[[i]]))])
  counts[[i]]$subject <- sapply(strsplit(counts[[i]]$unique, "_"), function(x) x[1])
  counts[[i]]$visit <- sapply(strsplit(counts[[i]]$unique, "_"), function(x) x[2])
  counts[[i]]$id <- sapply(strsplit(counts[[i]]$unique, "_"), function(x) paste((x)[1:2], collapse = "_"))
  return(counts[[i]])
})

## raw counts to frequency
counts <- sapply(names(counts), function(i) {
  counts[[i]][, grep("Cluster", colnames(counts[[i]]))] <- as.data.frame(t(apply(counts[[i]][grep("Cluster", colnames(counts[[i]]))], 1, function(x) x/sum(x) * 100)))
  return(counts[[i]])
})
sapply(names(counts), function(i) rowSums(counts[[i]][grep("cluster", colnames(counts[[i]]), ignore.case = TRUE)]))


ist.md <- list(tphe = read.csv("./ISTs/tphe_md.txt", sep = "\t", row.names = 1),
               ics = read.csv("./ISTs/ics_md.txt", sep = "\t", row.names = 1))

counts$tphe.cd4$IST <- ist.md$tphe$TPHE.IST[match(counts$tphe.cd4$id, rownames(ist.md$tphe))]
counts$tphe.cd8$IST <- ist.md$tphe$TPHE.IST[match(counts$tphe.cd8$id, rownames(ist.md$tphe))]
counts$ics.cd4$IST <- ist.md$ics$ICS.IST[match(counts$ics.cd4$id, rownames(ist.md$ics))]
counts$ics.cd8$IST <- ist.md$ics$ICS.IST[match(counts$ics.cd8$id, rownames(ist.md$ics))]

tphe.meta.data <- read.csv("./data_source/TPHE2_with.spillover/PRISM.TPHE2_MetaData.csv")
colnames(tphe.ics.clusters)[grep("subject", colnames(tphe.ics.clusters), ignore.case = T)] <- "id"

counts$tphe.cd4 <- merge(counts$tphe.cd4, tphe.meta.data, by = "unique")
counts$tphe.cd8 <- merge(counts$tphe.cd8, tphe.meta.data, by = "unique")

ggplot(counts$tphe.cd4, aes(x = ExperimentALIAS, y = IST))
ggplot(counts$tphe.cd4, aes(x = IST, y = ExperimentALIAS, color = ExperimentVALUE, size = unique)) + geom_point()
