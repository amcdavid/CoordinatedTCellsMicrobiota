### animation of normalization/warping effect

## source flow-related libraries and functions
source("src/flow.src.R", chdir = TRUE)

## FCS paths to TPHE files
tphe.original.files = list.files(recursive = TRUE, full.names = TRUE, pattern = "TPHE.fcs")

## NGL063 is 'bad' - staining error CD4/CD31 in same channel; exclude
tphe.original.files.good <- tphe.original.files[-grep("NGL063", tphe.original.files)]

## Let's start with just the healthy donors
tphe.original.files.good.HD <- tphe.original.files.good[grep("HD", tphe.original.files.good)]

## 'Warped' healthy donors
tphe.warped.files.good.HD <- list.files(recursive = TRUE, full.names = TRUE, pattern = "WARPED")

## select files for animation
frames.original <- tphe.original.files.good.HD[seq(1, 22, 5)]
frames.warped <- tphe.warped.files.good.HD[seq(1, 22, 5)]

frames.original.set <- fset.trim(read.flowSet(frames.original, transformation = FALSE))
frames.warped.set <- read.flowSet(frames.warped, transformation = FALSE)

## pre-define samp index so per-sample cell events match

samp <- sample(1:nrow(exprs(frames.original.set[[1]])), 50000)

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
p <- ggplot(dat.frames, aes(x = FSC.A, y = SSC.A)) + 
  geom_point(size = 0.1) +
  transition_states(id, transition_length = 10, state_length = 0, wrap = FALSE) +
  ease_aes('linear')

animate(p)

anim_save("FSCA.SSCA.warping.gif", path = "./results/plots/")
