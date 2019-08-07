library(flowCore)
library(ggplot2)
library(gridExtra)

##stacked densities per batch for visual QC

subsampled.fsom.files <- list.files("./data_modified/TPHE_warping/fSOM_CD4p_CD8p_PASS1_SUBSAMPLED/", full.names = TRUE)
batch.names <- paste0("NGL0", seq(47, 91, 2)) ; batch.names <- batch.names[-grep("NGL063", batch.names)]

trans.obj <- readRDS("./results_R/FCS.transforms/TPHE/asinh.transforms.TPHE.rds")
trans.obj <- trans.obj[grep('Violet G|Violet E|Violet A|Green D', names(trans.obj))]

## set grep argument to equal subset of either CD4 or CD8
sub.set <- fset.compensate(read.flowSet(subsampled.fsom.files[grep("CD8p_", basename(subsampled.fsom.files))], transformation = FALSE))
sub.set <- transform.set(sub.set) # needs 'trans.obj'

## fluors; set doc.name to equal population subset
for(i in seq(batch.names)){
  dat <- data.frame(exprs(as(sub.set[sampleNames(sub.set)[grep(batch.names[i], sampleNames(sub.set))]], "flowFrame")))
  qc.plots <- vector(mode = "list", length = length(trans.obj))
  for(j in seq(trans.obj)){
    qc.plots[[j]] <- ggplot(dat, aes_string(x = gsub("-", ".", gsub(" ", ".", names(trans.obj)))[j], group = "Original", fill = "Original")) + 
      geom_density(adjust = 1.5, alpha = 0.1, show.legend = FALSE)
  }
  doc.name <- paste("fSOM_CD8p.QC.stacked.fluors.densityplots", Sys.Date(), batch.names[i], "pdf", sep = ".")
  pdf(doc.name)
  grid.arrange(qc.plots[[1]], 
               qc.plots[[2]], 
               qc.plots[[3]],
               qc.plots[[4]],
               nrow = 2)
  dev.off()
}

## scatter
scatter <- c("FSC.A", "FSC.H", "FSC.W", "SSC.A", "SSC.H", "SSC.W")
for(i in seq(batch.names)){
  dat <- data.frame(exprs(as(sub.set[sampleNames(sub.set)[grep(batch.names[i], sampleNames(sub.set))]], "flowFrame")))
  qc.plots <- vector(mode = "list", length = length(scatter))
  for(j in seq(scatter)){
    qc.plots[[j]] <- ggplot(dat, aes_string(x = scatter[j], group = "Original", fill = "Original")) + 
      geom_density(adjust = 1.5, alpha = 0.1, show.legend = FALSE)
  }
  doc.name <- paste("fSOM_CD8p.QC.stacked.scatter.densityplots", Sys.Date(), batch.names[i], "pdf", sep = ".")
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
