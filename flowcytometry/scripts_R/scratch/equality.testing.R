## source flow-related libraries and functions
source("./source_R/flow.src.R", chdir = TRUE)

## original -> compensation -> transformation -> warping -> inverse transformation -> decompensation -> original
## test with two files; warping
tphe.warped.files = list.files(path = "data_modified/TPHE_warping/", recursive = TRUE, full.names = TRUE, pattern = "WARPED")

set.original <- read.flowSet(tphe.warped.files[c(2, 21)], column.pattern = "Time", invert.pattern = TRUE)

## compensate
set.comp <- fset.compensate(set.original)

## transform
asinhTrans <- arcsinhTransform(transformationId = "arcsinh-transformation", a = -1, b = 1/400, c = 2)
asinh.translist <- transformList('Violet G 550_40-A', asinhTrans)

set.transformed <- transform(set.comp, asinh.translist)

## inverse transformation
sinhTrans <- sinhTransform(transformationId = "inverse.arcsinh-transformation", a = -1, b = 1/400, c = 2)
sinh.translist <- transformList('Violet G 550_40-A', sinhTrans)

set.inverse <- transform(set.transformed, sinh.translist)

## decompensate
set.decomp <- fset.decompensate(set.inverse)

## test equalites
fset.equality(set.original, set.decomp, c(1:22)) # TRUE
fset.equality(set.comp, set.inverse, c(1:22)) # TRUE


## data warping; try range.fix
library(flowStats)

fsApply(set.transformed, range)[[1]][21]

set.transformed.rangefix <- fset.rangefix(set.transformed, 'Violet G 550_40-A')

fsApply(set.transformed.rangefix, range)[[1]][21]

set.transformed.rangefix.warped <- warpSet(set.transformed.rangefix
                                           , c("Violet G 550_40-A")
                                           , monwrd = FALSE
                                           , peakNr = 1
                                           , clipRange = 0.01
                                           , bwFac = 1.3
                                           , target = sampleNames(set.transformed)[grep("HD0189_TPHE_NGL091", sampleNames(set.transformed))])

## plot/QC results
densityplot(~ `Violet G 550_40-A`, set.transformed.rangefix, filter = curv1Filter("Violet G 550_40-A", bwFac = 1.3))
densityplot(~ `Violet G 550_40-A`, set.transformed.rangefix.warped, filter = curv1Filter("Violet G 550_40-A", bwFac = 1.3))

## inverse transformation of warped+range.fixed data
set.transformed.rangefix.warped.inverse <- transform(set.transformed.rangefix.warped, sinh.translist)

fset.equality(set.comp, set.transformed.rangefix.warped.inverse, c(1:20, 22)) # TRUE for non-warped columns
fset.equality(set.comp, set.transformed.rangefix.warped.inverse, c(21)) # TRUE for target; difference in warped sample

  set.transformed.rangefix.warped.inverse[[1]]@exprs - set.comp[[1]]@exprs
  set.transformed.rangefix.warped.inverse[[2]]@exprs - set.comp[[2]]@exprs

## decompensate

set.transformed.rangefix.warped.inverse.decomp <- fset.decompensate(set.transformed.rangefix.warped.inverse)

fset.equality(set.decomp, set.original, c(1:length(colnames(set.original)))) # no issues full-circle without warping
fset.equality(set.transformed.rangefix.warped.inverse.decomp, set.original, c(1:22)) # target unafffected; check off-target for column differences
fset.equality(set.transformed.rangefix.warped.inverse.decomp, set.original, c(1:6)) # scatter unaffected
fset.equality(set.transformed.rangefix.warped.inverse.decomp, set.original, c(7))

set.transformed.rangefix.warped.inverse.decomp[[1]]@exprs - set.original[[1]]@exprs
set.transformed.rangefix.warped.inverse.decomp[[2]]@exprs - set.original[[2]]@exprs

set.transformed.warped.inverse[[1]]@exprs - set.inverse[[1]]@exprs # fine here
set.transformed.warped.inverse[[2]]@exprs - set.inverse[[2]]@exprs # fine here

set.transformed.warped[[2]]@exprs - set.transformed[[2]]@exprs # fine here

set.inverse[[2]]@exprs - set.comp[[2]]@exprs # fine here

set.inverse.decomp[[2]]@exprs - set.original[[2]]@exprs # fine here
set.inverse.decomp[[1]]@exprs - set.original[[1]]@exprs # fine here

# final decomp for off-targets is the issue....

colnames(set.transformed.rangefix.warped.inverse.decomp[[1]])
colnames(set.transformed.rangefix.warped.inverse.decomp[[1]]@exprs)
colnames(set.original[[1]]@exprs)
colnames(set.transformed.rangefix.warped.inverse.decomp[[1]]@exprs) == colnames(set.original[[1]]@exprs)

colnames(set.transformed.rangefix.warped.inverse.decomp[[2]])
colnames(set.transformed.rangefix.warped.inverse.decomp[[2]]@exprs) # problem here
colnames(set.original[[2]]@exprs)
colnames(set.transformed.rangefix.warped.inverse.decomp[[2]]@exprs) == colnames(set.original[[2]]@exprs)

colnames(set.inverse[[2]]@exprs)
colnames(set.transformed[[2]]@exprs)
colnames(set.transformed.rangefix.warped[[2]]@exprs) # FUUUUUUU....
colnames(set.transformed.rangefix.warped[[1]]@exprs) # okay here; need to revisit warping function
