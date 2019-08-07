##
source("source_R/flow.src.R")

##

swift.template.fcs.path <- "/Volumes/nlaniewski/SWIFT_results/NGL092_V19.SEB/Consensus_Template/export_RPRC0530105_19_SEB_CD4+CD69+_MM.fcs"

swift.template.fcs <- read.FCS(swift.template.fcs.path, transformation = FALSE)

dat <- data.frame(exprs(swift.template.fcs))

cluster.ids <- c("Merge_1e0", "Merge_1e1", "Merge_1e2")

dat$cid.1   <- (dat[, "Merge_1e0"] - 25000)/25000 * 1
dat$cid.10  <- (dat[, "Merge_1e1"] - 25000)/25000 * 10
dat$cid.100 <- (dat[, "Merge_1e2"] - 25000)/25000 * 100

unique(dat$cid.1)
unique(dat$cid.10)
unique(dat$cid.100)

dat$cid <- dat$cid.1 + dat$cid.10 + dat$cid.100

cid.index <- data.frame(cid.index = sort(unique(dat$cid)))

dat.sub <- sample_n(dat, 200000)
dat.sub.index <- as.numeric(rownames(dat.sub))

dat.sub.tsne <- dat[, ]

dims.original <- c("FSC.A", "FSC.H", "FSC.W", "SSC.A", "SSC.H", "SSC.W",
                   "Blue.B.515_20.A",
                   "Green.A.780_40.A", "Green.B.710_50.A", "Green.D.610_20.A", "Green.E.575_25.A",
                   "Red.A.780_60.A", "Red.B.710_50.A", "Red.C.660_20.A",
                   "Violet.A.780_60.A", "Violet.B.705_70.A", "Violet.C.660_40.A", "Violet.D.605_40.A", "Violet.E.585_42.A", "Violet.G.550_40.A", "Violet.H.450_50.A")

swift.template.fcs.original <- fset.compensate(read.flowSet(swift.template.fcs.path, transformation = FALSE, column.pattern = paste(dims.original, collapse = "|")))
elgcl <- elgcl.transform(swift.template.fcs.original, 370000)
swift.template.fcs.original <- fset.transform(swift.template.fcs.original)

dat.original <- data.frame(exprs(swift.template.fcs.original[[1]]))
markers <- get.markers(swift.template.fcs.original[[1]])
markers[grep("CD14", markers)] <- "LIVE.DEAD"
colnames(dat.original) <- markers

dat.original.sub <- dat.original[dat.sub.index, ]

dims.used <- c("IL8", "IFNG", "IL4", "TNFA", "CD45RA", "IL2", "IL17")

dat.original.sub.tsne <- dat.original.sub[, dims.used]

write.table(dat.original.sub.tsne, "tsne.swift.csv", row.names = FALSE, col.names = FALSE, sep = ",")

saveRDS(dat.original.sub.tsne, "dat.original.sub.tsne.RDS")
saveRDS(dat.original.sub, "dat.original.sub.RDS")
saveRDS(dat.sub.index, "dat.sub.index.RDS")
saveRDS(dat.sub, "dat.sub.RDS")

dat.original.sub$cid <- dat.sub$cid
fcs.explorer(dat.original.sub, "IL8", "IL2", 20000)

final.embedding <- fread("final.embedding.csv", col.names = c("tsne.X", "tsne.Y"))

dat.original.sub <- cbind(dat.original.sub, final.embedding)
fcs.explorer(dat.original.sub, "tsne.X", "tsne.Y", 20000)

dat.original.sub.aggregate <- aggregate(dat.original.sub, list(dat.original.sub$cid), median)
heatmap <- data.frame(row.names = dat.original.sub.aggregate$cid, dat.original.sub.aggregate[, dims.used])

fcs.tsne.swift <- list(input.data = dat.original.sub,
                       agg.data   = dat.original.sub.aggregate,
                       hm.data    = data.frame(row.names = dat.original.sub.aggregate$cid, dat.original.sub.aggregate[, dims.used]))

fcs.explorer(fcs.tsne.swift$input.data, "tsne.X", "tsne.Y", 20000)
pheatmap(fcs.tsne.swift.frame$hm.data)

saveRDS(fcs.tsne.swift, "fcs.tsne.swift.RDS")
