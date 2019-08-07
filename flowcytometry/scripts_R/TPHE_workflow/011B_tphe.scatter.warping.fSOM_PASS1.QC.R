library(flowCore)
library(ggplot2)
library(gridExtra)

##stacked densities per batch for visual QC; set folder to population

subsampled.fsom.files <- list.files("./data_modified/TPHE_warping/fSOM_CD8p_PASS1_scatter.warped_SUBSAMPLED/", full.names = TRUE)
batch.names <- paste0("NGL0", seq(47, 91, 2)) ; batch.names <- batch.names[-grep("NGL063", batch.names)]

sub.set <- read.flowSet(subsampled.fsom.files, transformation = FALSE)

## scatter
scatter <- c("FSC.A", "FSC.H", "FSC.W", "SSC.A", "SSC.H", "SSC.W")
for(i in seq(batch.names)){
  dat <- data.frame(exprs(as(sub.set[sampleNames(sub.set)[grep(batch.names[i], sampleNames(sub.set))]], "flowFrame")))
  qc.plots <- vector(mode = "list", length = length(scatter))
  for(j in seq(scatter)){
    qc.plots[[j]] <- ggplot(dat, aes_string(x = scatter[j], group = "Original", fill = "Original")) + 
      geom_density(adjust = 1.5, alpha = 0.1, show.legend = FALSE)
  }
  doc.name <- paste("fSOM_CD8p.QC.stacked.scatter.densityplots", Sys.Date(), batch.names[i], "pdf", sep = ".") # set name
  pdf(doc.name)
  grid.arrange(qc.plots[[1]], 
               qc.plots[[2]], 
               qc.plots[[3]],
               qc.plots[[4]],
               qc.plots[[5]],
               qc.plots[[6]],
               nrow = 2)
  dev.off()
}
