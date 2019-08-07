library(flowCore); library(FlowSOM); library(pheatmap); library(shiny); library(ggplot2); library(viridis); library(Rtsne); library(gganimate); library(dplyr)

### read in and map a single file
input.files <- list.files("./data_source/ICS2_with.spillover/", recursive = TRUE, full.names = TRUE, pattern = "*SEB.fcs$")
input.files[grep("RPRC0530106", input.files)]

tmp <- read.flowSet(input.files[grep("RPRC0530106", input.files)], transformation = F)

tmp <- fsApply(tmp, function(i) compensate(i, i@description$`$SPILLOVER`))

fsom <- readRDS("./analyses_R/FSOM_analysis/01_fSOM.CD3LIVE.rds")
fsom$markers <- get.markers(tmp[[1]])
fsom$markers[grep("CD14", fsom$markers)] <- "CD14_LIVEDEAD"

tmp <- transform(tmp, fsom$elgcl)

# tmp.inverse <- transform(tmp.trans, inverseLogicleTransform(estimateLogicle(tmp.comp[[1]], colnames(tmp.comp[[1]])[-grep("FSC|SSC|Time|Original", colnames(tmp.comp[[1]]))])))
# 
# tmp.comp[[1]]@exprs[, 7][1] == tmp.inverse[[1]]@exprs[, 7][1]                
# 
# all(round(tmp.comp[[1]]@exprs, digits = 5) == round(tmp.inverse[[1]]@exprs, digits = 5))

fsom.tmp <- NewData(fsom, tmp[[1]])
dat <- fsom.dat(fsom.tmp)
index <- sample(1:nrow(dat), 10000)

tsne.out <- Rtsne(dat[index, ][, fsom.tmp$map$colsUsed], perplexity = 50, verbose = TRUE)
dat.tsne <- cbind(dat[index, ], data.frame(tsne.x = tsne.out$Y[, 1], tsne.y = tsne.out$Y[, 2], frame = 1))

### OR...read in an existing .rds with tsne.x,y parameters
dat.tsne <- readRDS("./analyses_R/ICS_tSNE.Explore/SEB.CD3pos.RDS")$input.data
colnames(dat.tsne)[grep("tsne.x", colnames(dat.tsne), ignore.case = T)] <- "tsne.x"
colnames(dat.tsne)[grep("tsne.y", colnames(dat.tsne), ignore.case = T)] <- "tsne.y"
index <- sample(1:nrow(dat.tsne), 10000)
dat.tsne <- cbind(dat.tsne[index, ], data.frame(frame = 1))
###

ggplot() +
  geom_point(data = dat.tsne,
             aes(x = tsne.x,
                 y = tsne.y,
                 color = factor(MCluster)),
             show.legend = F) +
  theme_void()

## Random dimensions
random_disc_jitter <- function(num_points,
                               disc_radius,
                               random_seed = 42) {
  
  if(!is.null(random_seed)) {
    set.seed(random_seed)
  }
  
  # Random radius positions
  r <- runif(num_points, 0, disc_radius ^ 2)
  # Random angles
  t <- runif(num_points, 0, 2 * pi)
  
  # Convert radius and angles to cartesian coordinates
  data.frame(x = sqrt(r) * cos(t),
             y = sqrt(r) * sin(t))
}

# Get the max x or y value from the tsne dims.
# I'll use this for the radius
max_tsne <- max(abs(c(dat.tsne$tsne.x, dat.tsne$tsne.y)))

dat.tsne <- dat.tsne %>%
  mutate(random_x = random_disc_jitter(num_points = n(),
                                       disc_radius = max_tsne)$x,
         random_y = random_disc_jitter(num_points = n(),
                                       disc_radius = max_tsne)$y)

ggplot() +
  geom_point(data = dat.tsne,
             aes(x = random_x,
                 y = random_y,
                 color = factor(MCluster)),
             show.legend = F) +
  theme_void()

fermat_jitter <- function(num_points, 
                          size, 
                          center_x, 
                          center_y) {
  
  golden_ratio <- (sqrt(5) + 1) / 2
  fibonacci_angle <- 360 / (golden_ratio ^ 2)
  
  ci <- sqrt(size / num_points)
  
  x <- rep(center_x, num_points)
  y <- rep(center_y, num_points)
  
  
  for (m in 1:(num_points - 1)) {
    n <- m - 1
    r <- ci * sqrt(n)
    theta <- fibonacci_angle * (n)
    x[n] <- center_x + r * cos(theta)
    y[n] <- center_y + r * sin(theta)
    
  }
  
  data.frame(x = x, 
             y = y)
  
}

# Some trial and error for this value to find a size that packs the points without getting
# too much overlap
max_size <- 30
# We'll also need the max cluster n to scale the jitter radius
max_n <- max(table(dat.tsne$MCluster))

dat.tsne <- dat.tsne %>%
  group_by(MCluster) %>%
  mutate(centroid_x = mean(tsne.x, trim = 0.5),
         centroid_y = mean(tsne.y, trim = 0.5)) %>%
  mutate(fermat_x = fermat_jitter(num_points = n(),
                                  size = max_size * n() / max_n,
                                  center_x = centroid_x[1],
                                  center_y = centroid_y[1])$x,
         fermat_y = fermat_jitter(num_points = n(),
                                  size = max_size * n() / max_n,
                                  center_x = centroid_x[1],
                                  center_y = centroid_y[1])$y) %>%
  ungroup()

ggplot() +
  geom_point(data = dat.tsne,
             aes(x = fermat_x,
                 y = fermat_y,
                 color = factor(MCluster)),
             show.legend = F) +
  theme_void()


frame1 <- data.frame(x = dat.tsne$random_x,
                     y = dat.tsne$random_y,
                     MCluster = dat.tsne$MCluster,
                     frame = 1)

frame2 <- data.frame(x = dat.tsne$tsne.x,
                     y = dat.tsne$tsne.y,
                     MCluster = dat.tsne$MCluster,
                     frame = 2)

frame3 <- data.frame(x = dat.tsne$fermat_x,
                     y = dat.tsne$fermat_y,
                     MCluster = dat.tsne$MCluster,
                     frame = 3)

dat.tsne.anim <- rbind(frame1, frame2, frame3)

p <- ggplot(dat.tsne.anim, aes(x, y, color = factor(MCluster))) + 
  geom_point(show.legend = F) +
  theme_void() +
  transition_states(frame, transition_length = 10, state_length = 10, wrap = F) 

animate(p, renderer = gifski_renderer(loop = F))

a <- animate(p, renderer = ffmpeg_renderer())
anim_save("random_tsne.cells_fSOM.metaclusters.NOLOOP.mp4", a)
