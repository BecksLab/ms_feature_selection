# ================================================
# 04_variable_clustering_focus.R
# ================================================

library(tidyverse)
library(cluster)
library(pvclust)
library(ggdendro)
library(here)

setwd(here::here())

# --- Load cleaned data ---
metrics <- read.csv("data/cleaned/all_networks.csv")

# ================================================
# 1. Scale the metrics
# ================================================
metrics_scaled <- scale(metrics)  # center=TRUE, scale=TRUE

# ================================================
# 2. Compute correlation-based distance
# ================================================
# We use |correlation| to capture similarity regardless of sign
dist_mat <- as.dist(1 - abs(cor(metrics_scaled, use = "pairwise.complete.obs")))

# ================================================
# 3. Hierarchical clustering
# ================================================
hc <- hclust(dist_mat, method = "average")

# Functional categories (as an experiment)
behaviour_metrics <- c("ρ", "complexity", "robustness")
geometry_metrics <- c("connectance", "l_S", "links", "richness", "Clust", "GenSD", "VulSD", "LinkSD", "diameter", "intervals")
path_metrics <- c("ChLen", "ChSD", "ChNum", "path", "S1", "S2", "S4", "S5", "omnivory", "loops", "predpreyRatio", "distance")
node_metrics <- c("basal", "top","intermediate", "herbivory", "cannibal", "TL", "centrality", "MaxSim")

dend_data <- dendro_data(hc) 

tip_names <- label(dend_data) %>%
  mutate(Category = case_when(
    label %in% behaviour_metrics ~ "Behaviour",
    label %in% geometry_metrics ~ "Geometry",
    label %in% path_metrics ~ "Path",
    label %in% node_metrics ~ "Node",
    TRUE ~ "Other"))

dend_plot <- 
  ggplot(segment(dend_data)) +
  geom_segment(aes(x = x, 
                   y = y, 
                   xend = xend, 
                   yend = yend)) +
  geom_text(data = tip_names, 
            aes(x = x, 
                y = -0.02, 
                label = label,
                colour = Category),
            angle = 45, hjust = 1, size = 3) +
  theme_void() +
  labs(y = "1 - |Correlation| (Distance)", 
       x = "") +
  ylim(-0.05, max(segment(dend_data)$y) * 1.1) +
  theme(axis.title.y = element_text(angle = 90))

# Save figure
ggsave("../figures/metric_hclust.png", 
       dend_plot, 
       width = 10, height = 8)

# ================================================
# 4. Bootstrap clusters for stability (optional)
# ================================================
set.seed(66)
pv <- pvclust(metrics_scaled, 
              method.hclust="average", 
              method.dist="correlation", nboot=1000)
plot(pv)
pvrect(pv, alpha=0.95)  # highlights clusters with >= 95% bootstrap support

# Compute correlation-based distance between metrics
dist_mat <- as.dist(1 - abs(cor(metrics_scaled)))

# Compute silhouette widths for k = 2 to 10
sil_width <- numeric(9)
for (k in 2:10) {
  cluster_assign <- cutree(hclust(dist_mat, method = "average"), k)
  sil <- silhouette(cluster_assign, dist_mat)
  sil_width[k-1] <- mean(sil[, 3])
}

plot(2:10, sil_width, type="b", xlab="Number of clusters (k)", ylab="Average silhouette width")
abline(v=which.max(sil_width)+1, col="red", lty=2)
cat("Optimal k suggested by silhouette:", which.max(sil_width)+1, "\n")

# ================================================
# 5. Interpret clusters
# ================================================

# Cut dendrogram into clusters
clusters <- cutree(hc, k = 4)
cluster_df <- data.frame(Metric=names(clusters), Cluster=clusters)

# Summary: metrics per cluster
cluster_summary <- cluster_df %>%
  group_by(Cluster) %>%
  summarise(Metrics = paste(Metric, collapse=", "),
            Count = n())

print(cluster_summary)

# --- Select representative metric from each cluster ---
# Here you might choose:
# - Ecologically meaningful metric
# - Metric with highest mean correlation to others in the cluster

representative_metrics <- cluster_df %>%
  group_by(Cluster) %>%
  summarise(Representative = Metric[1])  # simple: first metric per cluster

print(representative_metrics)

# Save cluster info
write.csv(cluster_df, "../tables/metric_clusters.csv", row.names=FALSE)
write.csv(cluster_summary, "../tables/metric_cluster_summary.csv", row.names=FALSE)
write.csv(representative_metrics, "../tables/reduced_metrics_from_clustering.csv", row.names=FALSE)
