library(data.table)
library(ggplot2)
cmv.dat <- fread("./data_source/TPHE2_with.spillover/CD8posCD28negCD57pos_V19.csv")
cmv.status <- fread("./data_source/TPHE2_with.spillover/CMVstatus.edit.csv")

colnames(cmv.status)[grep("Positive", colnames(cmv.status))] <- "CMV.PCR_status.Saliva"
cmv.status$CMV.PCR_status.Saliva[grep("1", cmv.status$CMV.PCR_status.Saliva)] <- "Positive"
cmv.status$CMV.PCR_status.Saliva[grep("0", cmv.status$CMV.PCR_status.Saliva)] <- "Negative"
cmv.status$CMV.PCR_status.Saliva[is.na(cmv.status$CMV.PCR_status.Saliva)] <- "Plasma sample but no Saliva?"

cmv.skit <- fread("./data_source/TPHE2_with.spillover/CMVserologykit.csv")

cmv <- merge(cmv.dat, cmv.status, by = "Patient.ID")

cmv$serology.kit <- 0

skit.positive <- cmv$Patient.ID[cmv.skit$Serology.Kit == 1]
cmv$Patient.ID  %in% skit.positive

ggplot(cmv, aes(x = CMV.PCR_status.Saliva, y = CD8posCD28negCD57pos)) + 
  labs(title = "Visit 19 Samples") +
  geom_jitter(width = 0.15) +
  ylab("% CD28- CD57+ of CD8+ T-cells") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16, vjust = +2))

ggplot(skit, aes(x = Serology.Kit, y = CD8posCD28negCD57pos)) + 
  labs(title = "Visit 19 Samples") +
  geom_jitter(width = 0.15) +
  ylab("% CD28- CD57+ of CD8+ T-cells") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16, vjust = +2))

####
library(data.table); library(shiny); library(ggplot2); library(viridis); library(ggord)

counts.merged <- read.csv("./results_R/fSOMs/TPHE_ICS/Cluster.counts.ALLmerged.csv")
counts.merged$term.visit <- sub("_", " ", sub("V19", "12-months", sub("V7", "Discharge", sub("V1$", "Birth", counts.merged$term.visit))))
counts.merged$term.visit <- factor(counts.merged$term.visit)
counts.merged$term.visit <- factor(counts.merged$term.visit, levels(counts.merged$term.visit)[c(2,3,1,5,6,4)])
counts.merged$SubID <- as.character(counts.merged$SubID)
counts.merged$VisitALIAS <- factor(counts.merged$VisitALIAS, levels(counts.merged$VisitALIAS)[c(2, 3, 1)])
colnames(counts.merged)[grep("subject", colnames(counts.merged))] <- "SampleID"

tphe <- readr::read_tsv("./results_R/MI2/tphe_basic.txt")
tphe$Subj <- forcats::fct_reorder(tphe$Subject, tphe$GAB)

counts.merged.tpheISTs <- merge(counts.merged, tphe, by = "SampleID")

ggplot(counts.merged.tpheISTs, aes(TermALIAS, TPHE.CD8_Meta.Cluster_22, color = Renamed_IST)) + geom_point() + facet_wrap(~VisitALIAS)

##
cmv.status <- read.csv("./data_source/TPHE2_with.spillover/CMVstatus.edit.csv")
cmv.status <- cmv.status[, c("Patient.ID", "CMV.status", "CMV.status.ALIAS")]
colnames(cmv.status) <- c("ParticipantID", "CMV.status", "CMV.status.ALIAS")       

as.character(subset(cmv.status, cmv.status$CMV.status == 1)$ParticipantID) %in% as.character(counts.merged.tpheISTs$ParticipantID)

counts.merged.tpheISTs.cmv <- merge(counts.merged.tpheISTs, cmv.status, by = "ParticipantID")
ggplot(counts.merged.tpheISTs.cmv, aes(CMV.status.ALIAS, TPHE.CD8_Meta.Cluster_22, color = Renamed_IST)) + geom_jitter(width = .2) + facet_wrap(~VisitALIAS+Renamed_IST)

cmv.12month.dat <- subset(counts.merged.tpheISTs.cmv, counts.merged.tpheISTs.cmv$VisitALIAS == "12-month")
cmv.12month.dat <- subset(cmv.12month.dat, cmv.12month.dat$CMV.status != 2)

cmv.plot <- ggplot(subset(cmv.12month.dat, 
                          cmv.12month.dat$Renamed_IST == "TPHE_6" | cmv.12month.dat$Renamed_IST == "TPHE_7"), 
                   aes(Renamed_IST, CMV.status.ALIAS, color = TPHE.CD8_Meta.Cluster_22)) + 
  geom_jitter(width = .08, size = 5) +
  scale_color_viridis("TPHE CD8\nMetacluster 22\n(% of CD8+)", option = "inferno") +
  labs(title = "12-month samples tested for CMV (by PCR)") +
  xlab("Immune State Type") +
  ylab("CMV Infection (PCR-confirmed)") +
  theme_dark() +
  geom_hline(aes(yintercept = 1.5), color = "white", size = 1) +
  theme(legend.key.size = unit(1, "cm"),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 14),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 18))

pdf("TPHE.IST_12month_CMV.pdf", width = 8, height = 6)
cmv.plot
dev.off()

dat <- tphe.cd8.fsom_frames$RPRC0540076_19
markers <- readRDS("./results_R/fSOMs/TPHE/2019-01-18.fSOM_PASS3_CD8p.rds")$markers
colnames(dat)[1:length(markers)] <- markers

dat <- subset(dat, dat$CD57 > 0)
dat <- subset(dat, dat$GRZB > 0)
dat <- subset(dat, dat$CD28 > 0)

pdf("TPHE.CD8.metacluster22_overlays.pdf", width = 6, height = 6)

ggplot(dat, aes(GRZB, PERFORIN)) + 
  geom_point(shape = ".") +
  geom_point(data = subset(dat, dat$meta.cluster == 22), size = .4, color = "red") +
  theme_minimal() +
  xlab("Granzyme B") +
  ylab("Perforin") +
  labs(title = "TPHE CD8 Meta Cluster 22 Overlay") +
  theme(axis.title = element_text(size = 18),
        plot.title = element_text(size = 18))

ggplot(dat, aes(CD57, CD28)) + 
  geom_point(shape = ".") +
  geom_point(data = subset(dat, dat$meta.cluster == 22), size = .4, color = "red") +
  theme_minimal() +
  xlab("CD57") +
  ylab("CD28") +
  labs(title = "TPHE CD8 Meta Cluster 22 Overlay") +
  theme(axis.title = element_text(size = 18),
        plot.title = element_text(size = 18))

dev.off()
