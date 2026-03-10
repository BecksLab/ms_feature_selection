
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
metrics <- read.csv("data/cleaned/all_networks.csv") %>%
  as_tibble() %>%
  # Exclude stability response variables from clustering
  vibe_check(-c(ρ, complexity, robustness))

# ============================================================
# 3. DISTANCE MATRICES
# ============================================================
cor_signed <- cor(metrics, use = "pairwise.complete.obs")
dist_signed <- as.dist(1 - cor_signed)

cor_abs <- abs(cor_signed)
dist_abs <- as.dist(1 - cor_abs)

# ============================================================
# 3. HIERARCHICAL CLUSTERING
# ============================================================
hc_signed <- hclust(dist_signed, method = "average")
hc_abs    <- hclust(dist_abs, method = "average")

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
sil_abs    <- sil_width(hc_abs, dist_abs, max_k = 10)

# Plot silhouette curves

sil_df <- 
  tibble(signed = sil_signed,
         absolute = sil_abs,
         clust_size = 1:10) %>%
  pivot_longer(!clust_size)

ggplot(sil_df,
       aes(x = clust_size,
           y = value)) +
  geom_line(colour = "#00203B") + 
  geom_point(colour = "#00203B") +
  labs(x = "Number of clusters",
       y = "Average silhouette") +
  facet_wrap(vars(name)) +
  coord_cartesian(clip = "off") +
  figure_theme() +
  theme(axis.text = element_text(size = rel(0.9)),
        text = element_text(size = rel(0.8)),
        strip.text = element_text(size = rel(0.9)))

ggsave("../figures/silhouette_curves.png",
       width = 3000,
       height = 1500,
       units = "px",
       dpi = 600)

# ============================================================
# 4. SELECT OPTIMAL CLUSTERS
# ============================================================
k_signed <- which.max(sil_signed)
k_abs    <- 5  # User-defined

clusters_signed <- cutree(hc_signed, k = k_signed)
clusters_abs    <- cutree(hc_abs, k = k_abs)

cluster_df_signed <- data.frame(
  Metric = names(clusters_signed),
  Cluster = clusters_signed
)

cluster_df_abs <- data.frame(
  Metric = names(clusters_abs),
  Cluster = clusters_abs
)

# Summaries
cluster_summary_signed <- cluster_df_signed %>%
  group_by(Cluster) %>%
  summarise(Metrics = paste(Metric, collapse=", "),
            Count = n())

cluster_summary_abs <- cluster_df_abs %>%
  group_by(Cluster) %>%
  summarise(Metrics = paste(Metric, collapse=", "),
            Count = n())

# Save tables
write.csv(cluster_df_signed, "../tables/metric_clusters_signed.csv", row.names = FALSE)
write.csv(cluster_summary_signed, "../tables/metric_cluster_summary_signed.csv", row.names = FALSE)

write.csv(cluster_df_abs, "../tables/metric_clusters_abs.csv", row.names = FALSE)
write.csv(cluster_summary_abs, "../tables/metric_cluster_summary_abs.csv", row.names = FALSE)

# ============================================================
# 5. BOOTSTRAP SUPPORT
# ============================================================

pv <- pvclust(metrics,
              method.hclust = "average",
              method.dist = "correlation",
              nboot = 1000)

# Save bootstrap object
saveRDS(pv, "data/outputs/pvclust_object.rds")

# Convert hclust to dendrogram
hc <- pv$hclust
dend <- as.dendrogram(hc)
dend_data <- dendro_data(dend)

# Cut at 7 clusters
k <- 7
clusters_7 <- cutree(hc, k = k)

cluster_df <- data.frame(
  Metric = names(clusters_7),
  Cluster = clusters_7
)

# Join cluster membership to label positions
label_df <- dend_data$label %>%
  left_join(cluster_df, by = c("label" = "Metric"))

# Compute rectangle boundaries (now vertical bands)
rect_df <- label_df %>%
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
               linewidth = 0.8) +
  # Rectangles for 7 clusters
  geom_rect(data = rect_df,
            aes(xmin = -0.02, xmax = 0,
                ymin = ymin-0.3, ymax = ymax+0.3,
                fill = factor(Cluster)),
            alpha = 0.8) +
  # Metric labels
  geom_text(data = label_df,
            aes(x = -0.025, y = x, label = label),
            hjust = 1,
            size = rel(3),
            family = "space") +
  scale_fill_manual(values = setNames(pal_df$colour, as.character(pal_df$value)),
                    labels = pal_df$label,
                    limits = pal_df$value,
                    name = "Module") +
  labs(x = "Correlation Distance",
       y = "") +
  coord_cartesian(xlim = c(-0.2, max(dend_data$segment$y))) +
  figure_theme() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    legend.position = "right"
  )

ggsave("../figures/cluster_dendrogram.png",
       width = 5000,
       height = 2550,
       units = "px",
       dpi = 600)


# Export cluster table
write.csv(cluster_df,
          "../tables/metric_clusters_auto.csv",
          row.names = FALSE)

cluster_df %>%
  group_by(Cluster) %>%
  summarise(Size = n())

module_table <- cluster_df %>%
  group_by(Cluster) %>%
  summarise(
    Metrics = paste(Metric, collapse = ", "),
    Count = n(),
  ) %>%
  glow_up(Cluster = as.factor(Cluster)) %>%
  left_join(pal_df, 
            by = join_by(Cluster == value)) %>%
  vibe_check(Cluster, label, Count, Metrics) %>%
  arrange(Cluster)

write.csv(module_table,
          "../tables/module_summary_7clusters.csv",
          row.names = FALSE)


# ============================================================
# 6. RICHNESS CORRECTION (LOG-LOG)
# ============================================================

library(ppcor)

# Remove richness from matrix
metrics_no_rich <- metrics %>% dplyr::select(-richness)


# Residualise all metrics against richness
metrics_resid <- metrics_no_rich

for (col in colnames(metrics_no_rich)) {
  model <- lm(metrics_no_rich[[col]] ~ metrics$richness)
  metrics_resid[[col]] <- resid(model)
}

# Now compute correlation of residuals
cor_partial <- cor(metrics_resid, use = "pairwise.complete.obs")

dist_partial <- as.dist(1 - abs(cor_partial))
hc_partial <- pvclust(metrics_no_rich,
                      method.hclust = "average",
                      method.dist = "correlation",
                      nboot = 1000)

hlust_partial <- hclust(dist_partial, method = "average")

# Convert hclust to dendrogram
hc <- hc_partial$hclust
dend <- as.dendrogram(hc)
dend_data <- dendro_data(dend)

# Join cluster membership to label positions
label_df <- dend_data$label %>%
  left_join(cluster_df, by = c("label" = "Metric")) %>%
  mutate(Cluster = as.factor(Cluster))

# Build horizontal plot
ggplot() +
  geom_segment(data = dend_data$segment,
               aes(x = y, y = x,
                   xend = yend, yend = xend),
               linewidth = 0.8) +
  # Metric labels
  geom_text(data = label_df,
            aes(x = -0.015, y = x, label = label, colour = Cluster),
            hjust = 1,
            size = rel(2),
            family = "space") +
  scale_colour_manual(values = setNames(pal_df$colour, as.character(pal_df$value)),
                      labels = pal_df$label,
                      limits = pal_df$value,
                      name = "Module") +
  labs(x = "Correlation Distance",
       y = "") +
  coord_cartesian(xlim = c(-0.1, max(dend_data$segment$y))) +
  figure_theme() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    legend.position = "right"
  )

############################################################
# END
############################################################
