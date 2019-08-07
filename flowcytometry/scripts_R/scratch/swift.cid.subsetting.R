##
source("./source_R/flow.src.R")

##paths to merged median cluster counts
fcs.paths <- dir(path = "data_modified/SWIFT_results/NGL092_V19.SEB/Consensus_Template/", 
                   full.names = TRUE, 
                   recursive = TRUE, 
                   pattern = "_SWIFT.fcs")

fcs.template <- read.FCS(fcs.paths[1])

dat <- data.frame(exprs(fcs.template))

cid <- "100"

nrow(subset(dat, dat$Merge_1e2 == (as.numeric(strsplit(as.character(cid), "")[[1]][1]) + 1) * 25000
            &dat$Merge_1e1 == (as.numeric(strsplit(as.character(cid), "")[[1]][2]) + 1) * 25000
            &dat$Merge_1e0 == (as.numeric(strsplit(as.character(cid), "")[[1]][3]) + 1) * 25000))

  dat$cid_1e2 <- dat$Merge_1e2/25000 - 1
  dat$cid_1e1 <- dat$Merge_1e1/25000 - 1
  dat$cid_1e0 <- dat$Merge_1e0/25000 - 1

dat$cid <- paste0(dat$cid_1e2, dat$cid_1e1, dat$cid_1e0)  

nrow(dat[dat$cid == "203", ])


for (i in seq(cids)) {
  if(nchar(cids[i]) < 3) {
    cids[i] <- paste0("0", cids[i])
  }
}

