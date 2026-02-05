# ================================================
# 02_variable_clustering.R
# ================================================

library(tidyverse)
library(cluster)
library(pvclust)
library(here)

setwd(here::here())

# Load cleaned data
metrics <- read.csv("../data/cleaned/all_networks.csv") %>%
  select(
    richness, links, connectance, diameter, complexity, distance, basal, top, intermediate, predpreyRatio,
    herbivory, omnivory, cannibal, l_S, GenSD, VulSD, TL, ChLen, ChSD, ChNum, path, LinkSD,
    S1, S2, S4, S5, ρ, centrality, loops, robustness, intervals, MaxSim, Clust
  )

# Scale metrics
metrics_scaled <- scale(metrics)

# Distance matrix and hierarchical clustering
dist_mat <- as.dist(1 - abs(cor(metrics_scaled, use="pairwise.complete.obs")))
hc <- hclust(dist_mat, method="average")
plot(hc, main="Hierarchical Clustering of Network Metrics", cex=0.7)

# Bootstrap for cluster stability
set.seed(42)
pv <- pvclust(metrics_scaled, method.hclust="average", method.dist="correlation", nboot=1000)
plot(pv)
pvrect(pv, alpha=0.95)

# Select representative metrics manually or programmatically
reduced_metrics <- metrics %>% select(connectance, Clust, complexity, robustness)
write.csv(reduced_metrics, 
            "../data/cleaned/reduced_metrics.csv", 
            row.names = FALSE)
