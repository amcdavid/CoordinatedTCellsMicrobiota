### animation of normalization/warping effect

## source flow-related libraries and functions
source("./source_R/flow.src.R", chdir = TRUE)

## FCS paths to non-warped CD45RO/CD197 fSOM_CD4pos HD samples
non.warped = list.files("./data_modified/TPHE_warping/targets_HD/scatter.VG.VE.GD.VA.warped/fSOM_CD3.CD4.CD8.LIVEDEAD/fSOM_CD4pos/",
                        recursive = FALSE, full.names = TRUE, pattern = ".fcs")

## FCS paths to warped CD45RO/CD197 fSOM_CD4pos HD samples
warped = list.files("./data_modified/TPHE_warping/targets_HD/scatter.VG.VE.GD.VA.warped/fSOM_CD3.CD4.CD8.LIVEDEAD/fSOM_CD4pos/fSOM_CD4pos_scatter/fSOM_CD4pos_scatter/",
                        recursive = TRUE, full.names = TRUE, pattern = ".fcs")

## select files for animation
frames.original <- non.warped[1]
frames.warped <- warped[1]

frames.original.set <- fset.compensate(fset.trim(read.flowSet(frames.original, transformation = FALSE)))
frames.warped.set <- fset.compensate(read.flowSet(frames.warped, transformation = FALSE))

trans.obj <- readRDS("results_R/FCS.transforms/TPHE/asinh.transforms.TPHE.rds")

frames.original.set <- transform.set(frames.original.set)
frames.warped.set <- transform.set(frames.warped.set)

## pre-define samp index so per-sample cell events match
samp <- sample(1:nrow(exprs(frames.warped.set[[1]])), 40000)

frames.original.set <- fsApply(frames.original.set, function (i) {
  exprs(i) <- exprs(i)[samp, ]
  i
})
frames.warped.set <- fsApply(frames.warped.set, function (i) {
  exprs(i) <- exprs(i)[samp, ]
  i
})

markers <- get.markers(frames.original.set[[1]])

dat.original <- data.frame(exprs(as(frames.original.set, "flowFrame")))
dat.warped <- data.frame(exprs(as(frames.warped.set, "flowFrame")))

colnames(dat.original)[1:length(markers)] <- markers
colnames(dat.warped)[1:length(markers)] <- markers

dat.original$id <- 1
dat.warped$id <- 2

all(dat.original[, 7:(length(dat.original)-2)] == dat.warped[, 7:(length(dat.warped)-2)])

dat.frames <- rbind2(dat.original, dat.warped)

ggplot(dat.original[dat.original$Original == 5, ], aes(x = SSC.A)) + theme_void() + theme(legend.position="none") + geom_histogram(aes(color = factor(Original)), binwidth = 400)
ggplot(dat.warped[dat.warped$Original == 5, ], aes(x = SSC.A)) + theme_void() + theme(legend.position="none") + geom_histogram(aes(color = factor(Original)), binwidth = 400) +
  geom_vline(xintercept = 190000) +
  geom_vline(xintercept = 10000) +
  geom_vline(xintercept = 65000)

ggplot(dat.original, aes(x = FSC.A)) + theme_void() + theme(legend.position="none") + geom_histogram(aes(color = factor(Original)), binwidth = 400) + facet_wrap(~Original, nrow = 5, ncol = 1) +
  geom_vline(xintercept = 115650)
ggplot(dat.warped, aes(x = FSC.A, color = factor(Original))) + theme_void() + theme(legend.position="none") + geom_histogram(binwidth = 400) + facet_wrap(~Original, nrow = 5, ncol = 1)

library(gganimate)
p <- ggplot(dat.frames, aes(x = FSC.A)) + 
  geom_histogram(aes(color = factor(Original)), binwidth = 400) +
  theme_void() +
  theme(legend.position = "none") +
  geom_vline(xintercept = 115650) +
  facet_wrap(~Original, nrow = 5, ncol = 1) +
  transition_states(id, transition_length = 10, state_length = 0, wrap = FALSE) +
  ease_aes('linear')

animate(p)

anim_save("FSCA.warping.gif", path = "./results/plots/")

library(gganimate)
p <- ggplot(dat.frames, aes(x = SSC.A)) + 
  geom_histogram(aes(color = factor(Original)), binwidth = 400) +
  theme_void() +
  theme(legend.position = "none") +
  geom_vline(xintercept = 190000) +
  geom_vline(xintercept = 10000) +
  geom_vline(xintercept = 65000) +
  facet_wrap(~Original, nrow = 5, ncol = 1) +
  transition_states(id, transition_length = 10, state_length = 0, wrap = FALSE) +
  ease_aes('linear')

animate(p)

anim_save("SSCA.warping.gif", path = "./results/plots/")

library(gganimate)
p <- ggplot(dat.frames, aes(x = CD45RO, y = CD197)) + 
  geom_point(size = 0.1) +
  transition_states(id, transition_length = 10, state_length = 0, wrap = FALSE) +
  ease_aes('linear')

animate(p)

anim_save("FSCA.SSCA.warping.gif", path = "./results/plots/")
