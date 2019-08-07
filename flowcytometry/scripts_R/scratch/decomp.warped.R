deCompensateFlowFrame <- function(x, spillover) 
{
  cols <- colnames(spillover)
  sel <- cols %in% colnames(x)
  if(!all(sel)) {
    stop(keyword(x)[["FILENAME"]], 
         "\\nThe following parameters in the spillover matrix are not present in the flowFrame:\\n",
         paste(cols[!sel], collapse=", "), call.=FALSE)
  }
  e <- exprs(x)
  e[, cols] <- e[, cols] %*% spillover
  exprs(x) = e
  x
}

cols <- colnames(set.transformed.rangefix.warped.inverse[[1]]@description$`$SPILLOVER`)
cols <- cols[-grep("Violet G 550_40-A", cols)]
sel <- cols %in% colnames(dat)
dat <- exprs(set.transformed.rangefix.warped.inverse[[1]])
dat[, cols] <- dat[, cols] %*% set.transformed.rangefix.warped.inverse[[1]]@description$`$SPILLOVER`[-4, -4]
dat[, -21] - set.original[[1]]@exprs[, -grep("Violet G 550_40-A", colnames(set.original[[1]]@exprs))]

cols2 <- colnames(set.transformed.rangefix.warped.inverse[[2]]@description$`$SPILLOVER`)
sel2 <- cols2 %in% colnames(set.transformed.rangefix.warped.inverse[[2]])
dat2 <- exprs(set.transformed.rangefix.warped.inverse[[2]])
dat2[, cols2] <- dat2[, cols2] %*% set.transformed.rangefix.warped.inverse[[2]]@description$`$SPILLOVER`
dat2 - set.original[[2]]@exprs

dat2.comp <- exprs(set.comp[[2]])


cols <- colnames(set.transformed.rangefix.warped.inverse[[1]]@description$`$SPILLOVER`)
cols <- cols[1]
sel <- cols %in% colnames(dat)
dat <- exprs(set.transformed.rangefix.warped.inverse[[1]])
dat[, cols] <- dat[, cols] %*% set.transformed.rangefix.warped.inverse[[1]]@description$`$SPILLOVER`[1, ]
dat[, -21] - set.original[[1]]@exprs[, -grep("Violet G 550_40-A", colnames(set.original[[1]]@exprs))]

warped.decomp <- set.transformed.rangefix.warped.inverse.decomp[[1]]@exprs[, 21]
set.original.mod <- set.original[[1]]
set.original.mod@exprs[, 21] <- warped.decomp

set.original.mod@exprs - set.original[[1]]@exprs

range(set.original[[1]])

write.FCS(set.transformed.rangefix.warped.inverse.decomp[[1]], "test3.fcs")
