NGL091 <- fset.trim(read.flowSet("./data/TPHE2_compensated/NGL091/HD0189_TPHE_NGL091 TPHE.fcs"))
NGL091.comp <- fset.compensate(NGL091)
markers <- get.markers(NGL091.comp[[1]])
dat <- data.frame(exprs(as(NGL091.comp, "flowFrame"))) %>% sample_n(50000)
colnames(dat)[1:length(markers)] <- markers
dat.sub <- subset(dat, select = c(FSC.A, SSC.A))
dat.sub.trans <- asinh(0 + 1 * dat.sub) + 0
ggplot(dat.sub.trans, aes(x = FSC.A, y = SSC.A)) + geom_point(size = 0.1)

NGL091.trans <- fset.transform(NGL091.comp)

NGL091.FSC <- data.frame(exprs(NGL091[[1]]))[1]
NGL091.comp.FSC <- data.frame(exprs(NGL091.comp[[1]]))[1]
NGL091.trans.FSC <- data.frame(exprs(NGL091.trans[[1]]))[1]


asinhTrans <- arcsinhTransform(transformationId="ln-transformation", a=5, b=1, c=1)
translist <- transformList(c('FSC-A','SSC-A'), asinhTrans)
NGL091.asinhtrans <- transform(NGL091, translist)

dat <- data.frame(exprs(as(NGL091.comp.asinhtrans, "flowFrame"))) %>% sample_n(50000)
colnames(dat)[1:length(markers)] <- markers
fcs.explorer(dat, "FSC.A", "SSC.A", 25000)

####
sample.raw <- read.FCS(tphe.original.files.good.HD.67.91[2])
sample.comp <- compensate(sample.raw, sample.raw@description$`$SPILLOVER`)
markers <- get.markers(sample.comp)

channels <- colnames(sample.comp)

exprs(sample.comp) <- exprs(sample.comp)[exprs(sample.comp)[,"FSC-A"] <= 250000, ]
exprs(sample.comp) <- exprs(sample.comp)[exprs(sample.comp)[,"SSC-A"] <= 250000, ]

dat <- as.data.frame(exprs(sample.comp))
names(dat)[1:length(markers)] <- markers

dat.sub <- sample_n(dat,50000)

cofactors <- list()

asinh.histo <- function(input,marker.name,a,b,c){
  vals <- as.data.frame(input[,marker.name])
  names(vals) <- marker.name
  vals.asinh <- asinh(a + b * vals) + c
  p <- ggplot(data = vals.asinh, aes_string(x=marker.name)) + geom_density()
  print(p)
  print(marker.name)
  cofactors[[marker.name]] <<- list("a" = a, "b" = b, "c" = c)
  # if(is.null(cofactors[[marker.name]])){
  # cofactors$marker.name <<- list("a" = a, "b" = b, "c" = c)
  # names(cofactors)[grep("marker.name",names(cofactors))] <<- marker.name
  # }
  # else{
  #   cofactors[[marker.name]] <<- list("a" = a, "b" = b, "c" = c)
  # }
}

asinh.histo(dat.sub,"SSC.A",1,1/2,1)
