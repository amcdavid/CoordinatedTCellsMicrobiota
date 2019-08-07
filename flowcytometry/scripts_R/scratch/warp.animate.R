dat.start <- data.frame(exprs(non.warped))
dat.end <- data.frame(exprs(warped))
colnames(dat.start) <- markers
colnames(dat.end) <- markers

index <- sample(1:nrow(dat.start), 20000)

dat.start <- dat.start[index, ]
dat.end <- dat.end[index, ]

dat.start$id <- 1
dat.end$id <- 2

ggplot(dat.start, aes(x = CD31, y = FSC.A)) + geom_point(size = 0.5) + xlim(c(0,7.5))
ggplot(dat.end, aes(x = CD31, y = FSC.A)) + geom_point(size = 0.5)

dat <- rbind(dat.start, dat.end)

library(gganimate)
p <- ggplot(dat, aes(x = CD31, y = FSC.A)) + 
  geom_point(size = 0.5) +
  transition_states(id, transition_length = 10, state_length = 3, wrap = TRUE) +
  ease_aes('linear')

animate(p)

anim_save("CD31.FSCA.warping.gif", path = "./results_R/plots/")
