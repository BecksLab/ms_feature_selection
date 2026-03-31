
############################################################
# Clean & Reproducible Variable Clustering Pipeline
# - Signed & Absolute Correlation Clustering
# - Silhouette Selection
# - Bootstrap Support
# - Richness Correction
# - Table Export
# - Figure Export
############################################################

library(tidyverse)
library(cluster)
library(genzplyr)
library(pvclust)
library(factoextra)
library(ggdendro)
library(here)

set.seed(66)
setwd(here::here())
source("lib/plotting_themes.R")

# ============================================================
# 0. LOAD DATA
# ============================================================
metrics <- read.csv("data/cleaned/randomNetworks.csv") %>%
  as_tibble() %>%
  # Exclude stability response variables from clustering
  vibe_check(-c(ρ, complexity, robustness, control, resilience, id)) %>%
  na.omit()

# ============================================================
# 3. DISTANCE MATRICES
# ============================================================
cor_signed <- cor(metrics, use = "pairwise.complete.obs")
dist_signed <- as.dist(1 - cor_signed)

# ============================================================
# 3. HIERARCHICAL CLUSTERING
# ============================================================
hc_signed <- hclust(dist_signed, method = "average")

# Silhouette function
sil_width <- function(hc, dist_matrix, max_k = 15) {
  sil <- numeric(max_k)
  for (k in 2:max_k) {
    clusters <- cutree(hc, k = k)
    ss <- silhouette(clusters, dist_matrix)
    sil[k] <- mean(ss[, 3])
  }
  return(sil)
}

sil_signed <- sil_width(hc_signed, dist_signed, max_k = 10)

# Plot silhouette curve

sil_df <- tibble(x = 1:10,
                 y = sil_signed)

ggplot(sil_df,
       aes(x = x,
           y = y)) +
  geom_line(colour = "#00203B") + 
  geom_point(colour = "#00203B") +
  labs(x = "Number of clusters",
       y = "Average silhouette") +
  coord_cartesian(clip = "off") +
  figure_theme() +
  theme(axis.text = element_text(size = rel(0.9)),
        text = element_text(size = rel(0.8)),
        strip.text = element_text(size = rel(0.9)))

# ============================================================
# 4. SELECT OPTIMAL CLUSTERS
# ============================================================
k_signed <- which.max(sil_signed)

clusters_signed <- cutree(hc_signed, k = k_signed)

cluster_df_signed <- data.frame(
  Metric = names(clusters_signed),
  Cluster = clusters_signed
)

# Summaries
cluster_summary_signed <- cluster_df_signed %>%
  group_by(Cluster) %>%
  summarise(Metrics = paste(Metric, collapse=", "),
            Count = n())

# ============================================================
# 5. BOOTSTRAP SUPPORT
# ============================================================

pv <- pvclust(metrics,
              method.hclust = "average",
              method.dist = "correlation",
              nboot = 1000)

# Convert hclust to dendrogram
hc <- pv$hclust
dend <- as.dendrogram(hc)
dend_data <- dendro_data(dend)

# Cut at 7 clusters
k <- k_signed
clusters_7 <- cutree(hc, k = k)

# get cluster/module masterlist (for plotting)
clust_metada <-
  left_join(read_csv("../tables/metric_clusters_auto.csv"),
            read_csv("../tables/module_summary_7clusters.csv") %>%
              vibe_check(c(Cluster, label))) %>%
  lowkey(Module = label)

# Join cluster membership to label positions
label_df <- dend_data$label %>%
  left_join(clust_metada, 
            by = c("label" = "Metric")) %>%
  left_join(pal_df,
            by = c("Module" = "label"))

rect_df <- data.frame(
  label = names(clusters_7),
  Cluster = clusters_7) %>%
  left_join(dend_data$label) %>%
  group_by(Cluster) %>%
  summarise(
    ymin = min(x),
    ymax = max(x),
    xmin = -Inf,
    xmax = max(dend_data$segment$y) * 0.05
  )

# Build horizontal plot
ggplot() +
  geom_segment(data = dend_data$segment,
               aes(x = y, y = x,
                   xend = yend, yend = xend),
               linewidth = 0.8, colour = "#001628") +
  # Rectangles for 7 clusters
  geom_rect(data = rect_df,
            aes(xmin = -0.02, xmax = 0,
                ymin = ymin-0.3, ymax = ymax+0.3),
            alpha = 0.8, fill = "#A2AAAD") +
  # Metric labels
  geom_text(data = label_df,
            aes(x = -0.025, y = x, label = label, colour = Module),
            hjust = 1,
            size = rel(3),
            family = "space",
            key_glyph = "point") +
  scale_colour_manual(values = setNames(label_df$colour, as.character(label_df$Module)),
                      breaks = label_df$Module,
                      name = "Module (empirical networks)") +
  labs(x = "Correlation Distance",
       y = "",
       title = "Random Network Clustering") +
  coord_cartesian(xlim = c(-0.2, max(dend_data$segment$y))) +
  figure_theme() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    legend.position = "right"
  )

ggsave("../figures/S1_cluster_dendrogram.png",
       width = 5500,
       height = 3000,
       units = "px",
       dpi = 600)

############################################################
# END
############################################################
