ics.cd4.counts <- read.csv("./results_R/fSOMs/ICS/FINAL/ICS_fSOM.CD4cytokines_MetaClusters_Counts.csv")
ics.cd8.counts <- read.csv("./results_R/fSOMs/ICS/FINAL/ICS_fSOM.CD8cytokines_MetaClusters_Counts.csv")

tphe.cd4.counts <- read.csv("./results_R/fSOMs/TPHE/TPHE_fSOM.CD4_MetaClusters_Counts.csv")
tphe.cd8.counts <- read.csv("./results_R/fSOMs/TPHE/TPHE_fSOM.CD8_MetaClusters_Counts.csv")

combinded.counts <- read.csv("./results_R/fSOMs/TPHE_ICS/Cluster.counts.ALLmerged.csv")

##
ics.subjects <- sapply(as.character(ics.cd4.counts$unique), function(i) strsplit(as.character(i), "_")[[1]][1])
ics.visits <- sapply(as.character(ics.cd4.counts$unique), function(i) strsplit(as.character(i), "_")[[1]][2])

ics <- data.frame(subject = ics.subjects,
                  visit = ics.visits,
                  term = NA)


ics$term[grep("53", ics$subject)] <- "full-term"
ics$term[-grep("53", ics$subject)] <- "pre-term"

length(unique(ics$subject))
##

##
tphe.subjects <- sapply(as.character(tphe.cd4.counts$unique), function(i) strsplit(as.character(i), "_")[[1]][1])
tphe.visits <- sapply(as.character(tphe.cd4.counts$unique), function(i) strsplit(as.character(i), "_")[[1]][2])

tphe <- data.frame(subject = tphe.subjects,
                  visit = tphe.visits,
                  term = NA)


tphe$term[grep("53", tphe$subject)] <- "full-term"
tphe$term[-grep("53", tphe$subject)] <- "pre-term"

length(unique(tphe$subject))
##

##
length(unique(ics$subject[which(ics$term == "full-term")]))
length(unique(ics$subject[which(ics$term == "pre-term")]))

length(unique(tphe$subject[which(tphe$term == "full-term")]))
length(unique(tphe$subject[which(tphe$term == "pre-term")]))
##

##
combined.subjects <- sapply(as.character(combinded.counts$unique), function(i) strsplit(as.character(i), "_")[[1]][1])
combined.visits <- sapply(as.character(combinded.counts$unique), function(i) strsplit(as.character(i), "_")[[1]][2])

combined <- data.frame(subject = combined.subjects,
                   visit = combined.visits,
                   term = NA)


combined$term[grep("53", combined$subject)] <- "full-term"
combined$term[-grep("53", combined$subject)] <- "pre-term"

length(unique(combined$subject))
length(unique(combined$subject[which(combined$term == "full-term")]))
length(unique(combined$subject[which(combined$term == "pre-term")]))
##

##
subject.numbers <- data.frame(matrix(nrow = 2, 
                                     ncol = 3, 
                                     dimnames = list(c("full-term", "pre-term"), 
                                                     c("ics", "tphe", "combined"))))

subject.numbers[grep("full-term", row.names(subject.numbers)), ] <- c(length(unique(ics$subject[which(ics$term == "full-term")])),
                                                                      length(unique(tphe$subject[which(tphe$term == "full-term")])),
                                                                      length(unique(combined$subject[which(combined$term == "full-term")])))

subject.numbers[grep("pre-term", row.names(subject.numbers)), ] <- c(length(unique(ics$subject[which(ics$term == "pre-term")])),
                                                                      length(unique(tphe$subject[which(tphe$term == "pre-term")])),
                                                                      length(unique(combined$subject[which(combined$term == "pre-term")])))

