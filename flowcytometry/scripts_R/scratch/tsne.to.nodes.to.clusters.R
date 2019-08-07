TPHE.CD4p <- readRDS("./results_R/fSOMs/TPHE/TPHE.CD4pos.unified.RDS")

dat <- sample_n(TPHE.CD4p$input.data, 20000)

frame1 <- dat[, c("nodes", "MCluster", "tsne.X", "tsne.Y")]
frame1$id <- 1

dat.nodes <- cbind(aggregate(. ~ nodes, data = dat, FUN = median), node.counts = as.data.frame(table(dat$nodes))$Freq)

frame2 <- merge(frame1[, c("nodes", "MCluster")], dat.nodes[, c("nodes", "tsne.X", "tsne.Y")], by = "nodes")
frame2$id <- 2

frames <- rbind(frame1, frame2)

ggplot(frame1, aes(x = tsne.X, y = tsne.Y, color = factor(MCluster))) + 
  geom_count()
  #geom_point(data = dat.nodes, aes(x = tsne.X, y = tsne.Y, color = "red"))

ggplot(frame2, aes(x = tsne.X, y = tsne.Y, color = factor(MCluster))) + 
  geom_count(show.legend = FALSE) +
  theme_void()


ggplot(frames, aes(x = tsne.X, y = tsne.Y, color = factor(MCluster))) + 
  geom_point(show.legend = FALSE) + 
  transition_states(states = id,
                    transition_length = 3,
                    state_length = 1) +
  theme_void()


for(i in seq(to.list)){
  temp[[i]] <- ls()[to.list][i]
}
