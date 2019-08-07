##stacked densities per batch for visual QC

subsampled.fsom.files <- list.files("./data_modified/TPHE_warping/PERFORIN.FOXP3.CD185.CD31.GRZB_warped_CD4p_SUBSAMPLED/", full.names = TRUE)
batch.names <- paste0("NGL0", seq(47, 91, 2)) ; batch.names <- batch.names[-grep("NGL063", batch.names)]

sub.set <- read.flowSet(subsampled.fsom.files, transformation = FALSE)
sub.set <- fset.compensate(sub.set)

##
warp.parameters <- readRDS("./results_R/FCS.transforms/TPHE/TPHE.CD4.warp.parameters.rds") # channels-of-interest
warp.markers <- c("PERFORIN", "FOXP3", "CD185", "CD31", "GRZB")
warp.parameters <- warp.parameters[names(warp.parameters) %in% warp.markers] # warping parameters per channels-of-interest

trans.obj <- readRDS("./results_R/FCS.transforms/TPHE/asinh.transforms.TPHE.rds") # arcsinh transfomation parameters
trans.obj <- trans.obj[grep(paste(unlist(lapply(seq(warp.parameters), function(i) warp.parameters[[i]]$channel)), collapse = "|"), names(trans.obj))]

sub.set <- transform.set(sub.set)

library(gridExtra)

##
for(i in seq(batch.names)){
  dat <- data.frame(exprs(as(sub.set[sampleNames(sub.set)[grep(batch.names[i], sampleNames(sub.set))]], "flowFrame")))
  qc.plots <- vector(mode = "list", length = length(trans.obj))
  for(j in seq(trans.obj)){
    qc.plots[[j]] <- ggplot(dat, aes_string(x = gsub("-", ".", gsub(" ", ".", names(trans.obj)))[j], group = "Original", fill = "Original")) + 
      geom_density(adjust = 1.5, alpha = 0.1, show.legend = FALSE)
  }
  doc.name <- paste("PERFORIN.FOXP3.CD185.CD31.GRZB_warped_CD4p.QC", Sys.Date(), batch.names[i], "pdf", sep = ".")
  pdf(doc.name)
  grid.arrange(qc.plots[[1]], 
               qc.plots[[2]], 
               qc.plots[[3]],
               qc.plots[[4]],
               qc.plots[[5]],
               nrow = 2)
  dev.off()
}
