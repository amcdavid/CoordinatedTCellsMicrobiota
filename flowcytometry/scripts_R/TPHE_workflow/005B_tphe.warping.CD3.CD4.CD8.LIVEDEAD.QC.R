##stacked densities per batch for visual QC
library(flowCore)
library(ggplot2)

subsampled.warped.files <- list.files("./data_modified/TPHE_warping/CD3.CD4.CD8.LIVEDEAD_warped_SUBSAMPLED/", full.names = TRUE)
batch.names <- paste0("NGL0", seq(47, 91, 2)) ; batch.names <- batch.names[-grep("NGL063", batch.names)]

trans.obj <- readRDS("./results_R/FCS.transforms/TPHE/asinh.transforms.TPHE.rds")
trans.obj <- trans.obj[grep('Violet G|Violet E|Violet A|Green D', names(trans.obj))]

sub.set <- fset.compensate(read.flowSet(subsampled.warped.files, transformation = FALSE))
sub.set <- transform.set(sub.set)

library(gridExtra)

for(i in seq(batch.names)){
  dat <- data.frame(exprs(as(sub.set[sampleNames(sub.set)[grep(batch.names[i], sampleNames(sub.set))]], "flowFrame")))
  qc.plots <- vector(mode = "list", length = length(trans.obj))
  for(j in seq(trans.obj)){
    qc.plots[[j]] <- ggplot(dat, aes_string(x = gsub("-", ".", gsub(" ", ".", names(trans.obj)))[j], group = "Original", fill = "Original")) + 
      geom_density(adjust = 1.5, alpha = 0.1, show.legend = FALSE)
  }
  grid.arrange(qc.plots[[1]], 
               qc.plots[[2]], 
               qc.plots[[3]],
               qc.plots[[4]],
               nrow = 2)
  doc.name <- paste("Warped.QC.stacked.densityplots", Sys.Date(), batch.names[i], "pdf", sep = ".")
  pdf(doc.name)
  grid.arrange(qc.plots[[1]], 
               qc.plots[[2]], 
               qc.plots[[3]],
               qc.plots[[4]],
               nrow = 2)
  dev.off()
}