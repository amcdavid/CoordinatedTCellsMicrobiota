## workflow
source("flowsom.libfun.R")
##

### 01: define FCS file paths

## set 'file.path' to a named directory

  fp <- file.path("./data_source/ICS2_with.spillover/")
  
## string-match for data-set specific FCS files

FCS_FILES <- list(
  SEB.original              = dir(fp, recursive = TRUE, full.names = TRUE, pattern = "*SEB.fcs$"),
  SEB.original.subsampled   = dir(fp, recursive = TRUE, full.names = TRUE, pattern = "*SEB_subsampled.fcs$"),
  fSOM.CD3LIVE              = dir(fp, recursive = TRUE, full.names = TRUE, pattern = "*fSOM.CD3LIVE_.fcs$"),
  fSOM.CD3LIVE_subsampled   = dir(fp, recursive = TRUE, full.names = TRUE, pattern = "*fSOM.CD3LIVE_subsampled.fcs$"),
  fSOM.CD4pCD69p            = dir(fp, recursive = TRUE, full.names = TRUE, pattern = "*fSOM.CD4pCD69p_.fcs$"),
  fSOM.CD4pCD69p_subsampled = dir(fp, recursive = TRUE, full.names = TRUE, pattern = "*fSOM.CD4pCD69p_subsampled.fcs$"),
  fSOM.CD8pCD69p            = dir(fp, recursive = TRUE, full.names = TRUE, pattern = "*fSOM.CD8pCD69p_.fcs$")
                  )
str(FCS_FILES)

fSOM.objects <- list(
  fSOM.CD3LIVE_Pass1      = dir(pattern = "01_fSOM.CD3LIVE.rds"),
  fSOM.CD4CD8CD69_Pass2   = dir(pattern = "02_fSOM.CD4CD8CD69.rds"),
  fSOM.CD4cytokines_Pass3 = dir(pattern = "03_CD4_fSOM.CD4Cytokines.rds")
)
str(fSOM.objects)  

##create batch list (experiment #); dataset specific
batch.list <- function(fcs.file.paths){
  experiments <- unique(lapply(seq(length(fcs.file.paths)), FUN = function(x){
    (strsplit(fcs.file.paths,"/")[[x]][grep("NGL*",strsplit(fcs.file.paths,"/")[[x]])])}))
  batches <- lapply(seq(unique(experiments)), FUN = function(x){fcs.file.paths[grep(experiments[x],fcs.file.paths)]})
  names(batches) <- unique(experiments)
  batches
}
  
  batches <- batch.list(FCS_FILES$SEB.original)
  batches.sub <- batch.list(FCS_FILES$SEB.original.subsampled)
  
  ##'good' batches NGL056-NGL092; 'bad'/'high' batches NGL048-52;54
  bad.blue      <- unlist(batches[names(batches)[grep("NGL048|NGL050|NGL052", names(batches))]], use.names = FALSE)
  bad.blue.sub  <- unlist(batches.sub[names(batches.sub)[grep("NGL048|NGL050|NGL052", names(batches.sub))]], use.names = FALSE)
  
  high.blue     <- unlist(batches[names(batches)[grep("NGL054", names(batches))]], use.names = FALSE)
  high.blue.sub <- unlist(batches.sub[names(batches.sub)[grep("NGL054", names(batches.sub))]], use.names = FALSE)
  
  good.batches  <- unlist(batches[names(batches)[-grep("NGL046|NGL048|NGL050|NGL052|NGL054", names(batches))]], use.names = FALSE)
  
######

### 02A: read in .fcs files as a flowSet; no compensation or transformation applied; access FCS_FILES list using '$'

##dataset specific
  fcs.file.paths <- FCS_FILES$SEB.original.subsampled

  raw.set <- fset.trim(read.flowSet(fcs.file.paths))
  
  event.counts <- lapply(seq(raw.set), function(i) nrow(raw.set[[i]]))
  sum(unlist(event.counts))
  summary(unlist(event.counts))

### 02B: compensate  
  compensated.set <- fset.compensate(raw.set)
  
  rm(raw.set) #remove; some sets stress computer memory


### 02C: transform

## estimate and store the logicle transform parameters (per fluorescense channel) from the full range (concatenate; as(...,"flowFrame")) of data
  channels.fluor <- as.vector(compensated.set[[1]]@parameters@data$name[!is.na(compensated.set[[1]]@parameters@data$desc)])
  concatenate <- as(compensated.set, "flowFrame")
  elgcl <- estimateLogicle(concatenate, channel = channels.fluor)
  rm(concatenate)# remove; concatenate will be quite large...
  
  
  ## Or extract from fSOM (if stored)
  get.elgcl <- function(fSOM.object){
    fSOM.tmp <- readRDS(fSOM.object)
    elgcl <- fSOM.tmp$elgcl
  }
  elgcl <- elgcl()
  rm(concatenate)# remove; concatenate will be quite large...
  
  transformed.set <- fset.transform(compensated.set)

### plot/visualize
  ## create 'plot-friendly' marker names; get marker names from a single flowFrame
  ## dataset-specific renaming
  markers <- sub("CD14", "LIVE.DEAD.CD14", get.markers(transformed.set[[1]]))
  (markers)

  ## place-holder data.frame for plotting
  dat <- as.data.frame(exprs(as(transformed.set, "flowFrame")))
  colnames(dat)[1:length(markers)] <- markers
  
  fcs.explorer(dat)
  
### clean up enviroment
  rm(compensated.set)
  input.flowset <- transformed.set
  rm(transformed.set)

###    
## link it all together (02A,B,C); must assign file paths and 'elgcl' to environment first in order to read and transform
  input.flowset <- fset.transform(fset.compensate(fset.trim(read.flowSet(good.batches))))  

######
  
  
####should now have a 'fully-processed' flowSet that can be passed to flowSOM
####end FCS; begin flowSOM
  
##01: create initial fSOM object from flowset; uses data scaling
fSOM <- ReadInput(input.flowset, scale = TRUE)

##create 'plot-friendly' marker names from fSOM$prettyColnames; dataset specific
fSOM$markers <- sub("CD14","LIVE.DEAD_CD14",(sub("-",".",lapply(seq(fSOM$prettyColnames), function(i) (strsplit(fSOM$prettyColnames," "))[[i]][1]))))
fSOM$markers   
##02-A: create matched column index based on user-defined input; must assign 'markers'
  ##fSOM PASS 1: set columns to cluster for CD3+ LIVE singlets
  mset <- match(c("FSC.A","SSC.A","FSC.W","SSC.W","LIVE.DEAD_CD14","CD3"), fSOM$markers)
  head(fSOM$data[,mset], n = 20)
  
  ##fSOM PASS 2: set columns to cluster for CD4/8+ CD69+
  mset <- match(c("LIVE.DEAD_CD14", "CD3", "CD4", "CD8A", "CD69"), fSOM$markers)
  head(fSOM$data[,mset], n = 20)
  
  ##fSOM PASS 3: set columns to cluster for functional/cytokine+ T-cells
  mset <- match(c("IL8","IFNG","IL4","TNFA","CD45RA","IL2","IL17"), fSOM$markers)
  head(fSOM$data[,mset], n = 20)
  
##02-B: create self-organizing maps (SOM) from user-defined grid-size (xdim;ydim) and pre-determined markers
  ##PASS 1
  fSOM <- BuildSOM(fSOM, colsToUse = mset, xdim=7, ydim = 7, rlen = 10)
  ##PASS 2
  fSOM <- BuildSOM(fSOM, colsToUse = mset, xdim=5, ydim = 5, rlen = 10)
  ##PASS 3
  fSOM <- BuildSOM(fSOM, colsToUse = mset, xdim=10, ydim = 10, rlen = 10)
  
##03: build a minimal-spanning tree from SOMs
  fSOM <- BuildMST(fSOM)  

###plot results
  PlotStars(fSOM, view="grid")
  PlotNumbers(fSOM, view = "grid")
  PlotNode(fSOM, 24)
  PlotStars(fSOM)
  PlotNumbers(UpdateNodeSize(fSOM,reset=TRUE))
  PlotMarker(fSOM, channels.fluor[14])
  PlotMarker(fSOM, channels.fluor[7])

##04: meta-cluster SOM nodes based on user-defined input; iterate through meta-clustering 'nbclust' to identify/isolate 'pure' clusters
nbclust <- 4
{fSOM[["metaclustering"]] <- as.factor(metaClustering_consensus(fSOM$map$codes, k = nbclust))
  PlotStars(fSOM, backgroundValues = as.factor(fSOM$metaclustering))
}
PlotStars(fSOM, backgroundValues = as.factor(fSOM$metaclustering))

##05: fSOM$data into data.frame appended with node/meta-cluster assignment for use with shiny.FCSplots
dat <- as.data.frame(cbind(fSOM$data,
                           nodes =fSOM$map$mapping[,1], 
                           MCluster = fSOM$metaclustering[fSOM$map$mapping[,1]]))
names(dat)[1:length(fSOM$markers)] <- fSOM$markers

####pass dat to shiny.FCSplots.R

##heatmaps

####based on meta-clustering results: merge/split manually using factor levels
length(levels(fSOM$metaclustering))
##add new factor level
levels(fSOM$metaclustering) <- c(levels(fSOM$metaclustering),length(levels(fSOM$metaclustering))+1)
levels(fSOM$metaclustering)
##assign nodes-to-be-split to new factor level
fSOM$metaclustering[c(78,80)] <- 21
fSOM$metaclustering[c(86)] <- 6
fSOM$metaclustering

##for saving fSOM.object; reduce file size by trimming out extraneous rows/data; store estimated transform 'elgcl'
fSOM$elgcl <- elgcl
fSOM$metaclusters.anno <- vector("list", nbclust)
fSOM$metaclusters.anno[[1]] <- "CD8+CD69+"
savefsom <- function(fSOM.object,fSOM.name){
  tmp <- fSOM.object
  tmp$data <- head(tmp$data)
  tmp$map$mapping <- head(tmp$map$mapping)
  saveRDS(tmp, file = fSOM.name)
}

####should now have a fSOM object with meta-clusters
####end flowSOM; begin fcs subsetting

##01: subset original (raw) FCS file by identified meta-cluster number; index by row
##based on fSOM meta-clustering result; reads in original `raw` FCS file (no comp, no transform) and subsets based on chosen fSOM meta-cluster
##also generates meta-cluster counts
elgcl <- fSOM$elgcl
newfcs.fromsom <- function(fSOM.object, meta.cluster, FCS.paths, newfolder.name, post.name){
  for(i in seq(FCS.paths)){
    ##read in and scatter-trim raw FCS (linear/non-transformed data)
    FCS.raw <- read.FCS(FCS.paths[i])
    exprs(FCS.raw) <- exprs(FCS.raw)[exprs(FCS.raw)[,"FSC-A"] > 0,]
    exprs(FCS.raw) <- exprs(FCS.raw)[exprs(FCS.raw)[,"FSC-A"] < 250000,]
    exprs(FCS.raw) <- exprs(FCS.raw)[exprs(FCS.raw)[,"SSC-A"] > 0,]
    exprs(FCS.raw) <- exprs(FCS.raw)[exprs(FCS.raw)[,"SSC-A"] < 250000,]
    
    ##compensate using stored SPILL/SPILLOVER matrix  
    FCS.COMP <- compensate(FCS.raw,(keyword(FCS.raw)[grep("SPILL", names(keyword(FCS.raw)))][[1]]))
    ##transform using estimate logicle (parameter/channel dependent)
    FCS.COMPTRANS <- transform(FCS.COMP, elgcl)
    #FCS.COMPTRANS <- transform(FCS.COMP, estimateLogicle(FCS.COMP, channels = colnames(keyword(FCS.COMP)[grep("SPILL", names(keyword(FCS.COMP)))][[1]])))
    ##map compensated and transformed data to existing fSOM object (NewData function scales existing data)
    print(paste0("Mapping new data to fSOM.object...",FCS.COMPTRANS@description$`$FIL`))
    fSOM.tmp <- NewData(fSOM.object,FCS.COMPTRANS)
    ##generate tabled meta-cluster counts
    counts <- as.data.frame(table(fSOM.tmp$metaclustering[fSOM.tmp$map$mapping[,1]]))
    names(counts) <- c("Meta.Cluster","Freq")
    counts.fname <- paste(FCS.raw@description$`$SRC`,FCS.raw@description$`TUBE NAME`,post.name,"counts.csv",sep="_")
    ##generates a row index where rows match assigned meta.cluster; used to subset raw FCS
    #index <- grep(TRUE,fSOM.tmp$metaclustering[fSOM.tmp$map$mapping[,1]]==meta.cluster)
    index <- which(!is.na(match(fSOM.tmp$metaclustering[fSOM.tmp$map$mapping[,1]], meta.cluster)))
    ##subset raw FCS based on row index
    exprs(FCS.raw) <- exprs(FCS.raw)[index,]
    ##modify/update $FIL name
    FCS.raw@description$`$FIL` <- paste(FCS.raw@description$`$SRC`,FCS.raw@description$`TUBE NAME`,post.name,".fcs",sep="_")
    ##write FCS to directory
    if (dir.exists(file.path(dirname(FCS.paths[i]),newfolder.name))){
      outFile <- file.path(dirname(FCS.paths[i]),newfolder.name,FCS.raw@description$`$FIL`)
      #outCounts <- file.path(dirname(FCS.paths[i]),newfolder.name,counts.fname)
    }
    else{
      dir.create(file.path(dirname(FCS.paths[i]),newfolder.name))
      outFile <- file.path(dirname(FCS.paths[i]),newfolder.name,FCS.raw@description$`$FIL`)
      #outCounts <- file.path(dirname(FCS.paths[i]),newfolder.name,counts.fname)
    }
    if (dir.exists(file.path(dirname(FCS.paths[i]),newfolder.name,"counts"))){
      outCounts <- file.path(dirname(FCS.paths[i]),newfolder.name,"counts",counts.fname)
    }
    else{
      dir.create(file.path(dirname(FCS.paths[i]),newfolder.name,"counts"))
      outCounts <- file.path(dirname(FCS.paths[i]),newfolder.name,"counts",counts.fname)
    }
    write.FCS(FCS.raw,outFile)
    write.csv(counts,outCounts,row.names = FALSE)
  }}

##
newfcs.fromsom(fSOM.CD3LIVE,17,good.batches,"fSOM.CD3LIVE","fSOM.CD3LIVE")
newfcs.fromsom(fSOM,c(3,4,5),FCS_FILES$fSOM.CD3LIVE_cluster17,"fSOM.CD4pCD69p","fSOM.CD4pCD69p")
##
##
newfcs.fromsom(fSOM,17,bad.blue,"fSOM.CD3LIVE","fSOM.CD3LIVE")
newfcs.fromsom(fSOM.CD4CD8CD69,c(3,4,5),FCS_FILES$fSOM.CD3LIVE[1:42],"fSOM.CD4pCD69p","fSOM.CD4pCD69p")
newfcs.fromsom(fSOM.CD4CD8CD69,c(1,13),FCS_FILES$fSOM.CD3LIVE[1:42],"fSOM.CD8pCD69p","fSOM.CD8pCD69p")
##

COUNTS_PATHS <- dir(fp, recursive=TRUE, full.names=TRUE, pattern="*fSOM.CD3LIVE_counts.csv$")
COUNTS_PATHS <- dir(fp, recursive=TRUE, full.names=TRUE, pattern="*fSOM.CD4pCD69p_counts.csv$")
COUNTS_PATHS <- dir(fp, recursive=TRUE, full.names=TRUE, pattern="*CD4.Cytokines")
(COUNTS_PATHS)

library(data.table)
counts <- lapply(seq(COUNTS_PATHS), function(i) fread(COUNTS_PATHS[i]))
counts.names <- lapply(seq(COUNTS_PATHS), function(i) gsub("_counts.csv","",basename(COUNTS_PATHS[i])))
counts.names
names(counts) <- counts.names
counts <- mapply(`[[`, counts, 2)
counts <- t(counts)
counts <- as.data.frame(counts)
##
##as frequency (sum1)
counts <- as.data.frame(t(apply(counts, 1, function(i) i/sum(i))))
##as percent (sum100)
counts <- as.data.frame(t(apply(counts, 1, function(i) i/sum(i)))) * 100
##
##
rowSums(counts)
class(counts)
colnames(counts) <- sub("V","Meta.Cluster_",colnames(counts))
colnames(counts)
row.names(counts)
counts$dolSRC_Stim <- sub("_CD4.*", "", row.names(counts))
counts$dolSRC_Stim
meta.data <- fread("PRISM.ICS_MetaData.csv")
meta.data <- meta.data[meta.data$StimulationALIAS == "SEB"]
meta.data <- meta.data[meta.data$TermALIAS != "Healthy Adult"]
meta.data$TermALIAS <- factor(meta.data$TermALIAS, c("Pre-term", "Full-term"))
meta.data$VisitALIAS <- factor(meta.data$VisitALIAS,c("Birth (cord)","Discharge","12-month"))
meta.data$Cohort = factor(meta.data$Cohort)
meta.data$Cohort = factor(meta.data$Cohort,levels(meta.data$Cohort)[c(2:7,1)])
meta.data$Cohort
counts.meta <- merge(counts,meta.data,by = "dolSRC_Stim")
class(counts.meta)

ggplot(data = counts.meta, mapping = aes(x = TermALIAS)) + geom_jitter(mapping = aes(y = Meta.Cluster_17), width = 0.1)
ggplot(data = counts.meta, mapping = aes(x = TermALIAS)) + geom_jitter(mapping = aes(y = Meta.Cluster_17), width = 0.1) + facet_wrap("VisitALIAS")
ggplot(data = counts.meta, mapping = aes(x = TermALIAS)) + geom_jitter(mapping = aes(y = Meta.Cluster_17, color = GenderALIAS), width = 0.1) + facet_wrap("VisitALIAS")


##CD8 Cytokines
elgcl <- get.elgcl("02_fSOM.CD4CD8CD69.rds")

input.flowsetCD8CD69 <- fset.transform(fset.compensate(read.flowSet(FCS_FILES$fSOM.CD8pCD69p)))
markers <- sub("CD14", "LIVE.DEAD.CD14", get.markers(input.flowsetCD8CD69[[1]]))

fSOM <- ReadInput(input.flowsetCD8CD69, scale = TRUE)
fSOM$markers <- markers

##fSOM PASS 3: set columns to cluster for CD8 functional/cytokine+ T-cells
mset <- match(c("IL8", "IFNG", "CD107A", "IL4", "TNFA", "CD45RA", "IL2", "IL17"), fSOM$markers)

fSOM <- BuildSOM(fSOM, colsToUse = mset, xdim = 10, ydim = 10)
fSOM <- BuildMST(fSOM)
PlotStars(fSOM)

nbclust <- 20
{fSOM[["metaclustering"]] <- as.factor(metaClustering_consensus(fSOM$map$codes, k = nbclust))
  PlotStars(fSOM, backgroundValues = as.factor(fSOM$metaclustering))
}

pheatmap(fsom.clusters.heatmap(fSOM), scale = "row")

fSOM <- nodes.to.meta_fSOM(fSOM, grep("^1$|^11$", fSOM$metaclustering))
fSOM <- nodes.to.meta_fSOM(fSOM, 51)
fSOM <- nodes.to.meta_fSOM(fSOM, grep("^14$|^15$", fSOM$metaclustering))
fSOM <- nodes.to.meta_fSOM(fSOM, grep("^2$|^9$", fSOM$metaclustering))
fSOM <- nodes.to.meta_fSOM(fSOM, grep("^11$|^18$", fSOM$metaclustering))
fSOM <- nodes.to.meta_fSOM(fSOM, 39)
fSOM <- nodes.to.meta_fSOM(fSOM, grep("^1$|^12$", fSOM$metaclustering))
fSOM <- nodes.to.meta_fSOM(fSOM, grep("^12$|^17$", fSOM$metaclustering))
fSOM <- nodes.to.meta_fSOM(fSOM, grep("^15$|^16$", fSOM$metaclustering))
fSOM <- nodes.to.meta_fSOM(fSOM, grep("^3$|^11$", fSOM$metaclustering))
fSOM <- nodes.to.meta_fSOM(fSOM, grep("^2$|^6$", fSOM$metaclustering))
fSOM <- nodes.to.meta_fSOM(fSOM, 70)
fSOM <- nodes.to.meta_fSOM(fSOM, grep("^10$|^11$", fSOM$metaclustering))
fSOM <- nodes.to.meta_fSOM(fSOM, grep("^1$|^11$", fSOM$metaclustering))


pheatmap(fsom.clusters.heatmap(fSOM), scale = "row")
fcs.explorer(fsom.dat(fSOM), "IL8", "IL2", 5000)
PlotStars(fSOM, backgroundValues = as.factor(fSOM$metaclustering))

hm <- fsom.clusters.heatmap(fSOM)
hm <- hm[, -grep("IL17", colnames(hm))]
fSOM$heatmap <- pheatmap(hm, scale = "row")
fSOM$elgcl <- elgcl
savefsom(fSOM, "03B_CD8_fSOM.CD8Cytokines.RDS")

fSOM.CD8Cytokines <- readRDS("03B_CD8_fSOM.CD8Cytokines.RDS")

newfcs.fromsom(fSOM.CD8Cytokines, meta.cluster = NULL, FCS_FILES$fSOM.CD8pCD69p, "fSOM.CD8Cytokines", write.count.files = TRUE, write.fcs.files = FALSE, fcsfolder.name = NULL)

##counts for CD8 Cytokines
COUNTS_PATHS <- dir(fp, recursive=TRUE, full.names=TRUE, pattern="*fSOM.CD8Cytokines_counts.csv$")

raw.counts <- generate.counts(COUNTS_PATHS)
raw.counts$unique <- sub("_fSOM.*", "", rownames(raw.counts))
write.csv(raw.counts, "ICS_fSOM.CD8cytokines.MetaCluster_Counts.csv")
counts <- as.data.frame(t(apply(raw.counts[, grep("Meta.Cluster", colnames(raw.counts))], 1, function(i) i/sum(i)))) * 100
counts$unique <- sub("_fSOM.*", "", rownames(counts))
ICS.meta <- readRDS("./data_source/ICS_SEB.meta.rds")
counts.meta <- merge(counts, ICS.meta, by = "unique")
meta.explorer(counts.meta)

ggord(prcomp(counts.meta[, grep("Meta.Cluster", colnames(counts.meta))]), counts.meta$ExperimentALIAS, ellipse = FALSE, vec_ext = 0, vec_lab = FALSE, arrow = NULL, vectyp = 'blank')
