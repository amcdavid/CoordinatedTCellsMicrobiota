## fSOM.object with full dataset used to generate clusters
fSOM <- readRDS("./results_R/fSOMs/TPHE/2018-11-14.fSOM_PASS3_CD4p_rerun_ALLDATA.rds")

## extract data.frame
dat <- fsom.dat(fSOM)

## subset to smaller frame for passing to fastTSNE
dat <- sample_n(dat, 200000)

## write-to-csv dimensions used in flowSOM
write.table(dat[, fSOM$map$colsUsed], file = "./results_R/fSOMs/TPHE/input.forTSNE.TPHE_CD4p_PASS3.csv", row.names = FALSE, col.names = FALSE, sep = ",")

## run fastTSNE to generate final embeddings; read-in
file.edit("./scripts_R/TPHE/fastTSNE_PYTHON.R")
final.embeddings <- fread("./results_R/fSOMs/TPHE/embeddings.TSNE.TPHE_CD4p_PASS3.csv", col.names = c("tsne.X", "tsne.Y"))
ggplot(final.embeddings, aes(tsne.X, tsne.Y)) + geom_point(size = 0.5)

## combine input dat with tSNE XY
dat.fsom.tsne <- cbind(dat, final.embeddings)

fcs.explorer(dat.fsom.tsne, "tsne.X", "tsne.Y", 50000)

## unified frame for shiny

# TPHE CD4+ clusters
strip.fsom <- function(fSOM){
  tmp <- fSOM
  tmp$data <- head(tmp$data)
  tmp$map$mapping <- head(tmp$map$mapping)
  tmp
}

TPHE.CD4pos <- list(input.data  = dat.fsom.tsne,
                    counts.meta = counts.meta,
                    dims.used   = fSOM$markers[fSOM$map$colsUsed],
                    fsom        = strip.fsom(fSOM),
                    fsom.mclust.agg = aggregate(dat.fsom.tsne, list(dat.fsom.tsne$MCluster), median)
)
saveRDS(TPHE.CD4pos, "./results_R/fSOMs/TPHE/TPHE.CD4pos.unified.RDS")

file.edit("./scripts_R/TPHE/TPHE_Explorer.app.R")
