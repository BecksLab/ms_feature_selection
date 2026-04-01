############################################################
# Random Network Clustering Pipeline (Multi-Model)
# - Degree-preserving vs Connectance-constrained
# - Silhouette selection
# - Dendrograms
# - Cluster comparison (ARI)
############################################################

library(tidyverse)
library(cluster)
library(genzplyr)
library(pvclust)
library(factoextra)
library(ggdendro)
library(mclust)
library(here)

set.seed(66)
setwd(here::here())
source("lib/plotting_themes.R")

# ============================================================
# 0. LOAD DATA
# ============================================================
metrics_all <- read.csv("data/cleaned/randomNetworks.csv") %>%
  as_tibble() %>%
  vibe_check(-c(ρ, complexity, robustness, control, resilience, id))

# ============================================================
# NA / NaN / Zero summary by model
# ============================================================

summary_counts <- metrics_all %>%
  # keep model column, pivot everything else
  pivot_longer(
    cols = -model,
    names_to = "metric",
    values_to = "value"
  ) %>%
  group_by(model, metric) %>%
  summarise(
    n_total = n(),
    n_NA  = sum(is.na(value)),
    n_NaN = sum(is.nan(value)),
    n_zero = sum(value == 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  glow_up(pct_loss = n_zero/n_total * 100) %>%
  slay(-pct_loss)

# ============================================================
# 1. SPLIT BY NULL MODEL
# ============================================================
metrics_split <- metrics_all %>%
  split(.$model)

# ============================================================
# 2. CLUSTERING FUNCTION
# ============================================================
run_clustering <- function(metrics, label = "model") {
  
  metrics_clean <- metrics %>%
    vibe_check(-model) %>%
    vibe_check(where(~ sd(.x, na.rm = TRUE) > 0))
  
  # -------------------------
  # Correlation + distance
  # -------------------------
  cor_signed <- cor(metrics_clean, use = "pairwise.complete.obs")
  
  # Replace NA/NaN with 0
  cor_signed[is.na(cor_signed)] <- 0
  dist_signed <- as.dist(1 - cor_signed)
  
  # -------------------------
  # Hierarchical clustering
  # -------------------------
  hc <- hclust(dist_signed, method = "average")
  
  # -------------------------
  # Silhouette selection
  # -------------------------
  sil_width <- function(hc, dist_matrix, max_k = 10) {
    sil <- numeric(max_k)
    for (k in 2:max_k) {
      clusters <- cutree(hc, k = k)
      ss <- silhouette(clusters, dist_matrix)
      sil[k] <- mean(ss[, 3])
    }
    sil
  }
  
  sil <- sil_width(hc, dist_signed)
  k_opt <- which.max(sil)
  
  clusters <- cutree(hc, k = k_opt)
  
  cluster_df <- tibble(
    Metric = names(clusters),
    Cluster = clusters,
    model = label
  )
  
  return(list(
    hc = hc,
    sil = sil,
    k = k_opt,
    clusters = cluster_df
  ))
}

# ============================================================
# 3. RUN FOR BOTH MODELS
# ============================================================
results <- purrr::imap(metrics_split, ~run_clustering(.x, .y))

# ============================================================
# 4. SILHOUETTE PLOT
# ============================================================
sil_df <- bind_rows(
  tibble(k = 1:10,
         sil = results[["degree_null"]][["sil"]],
         model = "degree_preserving"),
  tibble(k = 1:10,
         sil = results[["connectance_null"]][["sil"]],
         model = "connectance_constrained")
)

ggplot(sil_df,
       aes(x = k,
           y = sil)) +
  geom_line(colour = "#00203B") + 
  geom_point(colour = "#00203B") +
  labs(x = "Number of clusters",
       y = "Average silhouette") +
  facet_wrap(vars(model)) +
  coord_cartesian(clip = "off") +
  figure_theme() +
  theme(axis.text = element_text(size = rel(0.9)),
        text = element_text(size = rel(0.8)),
        strip.text = element_text(size = rel(0.9)))

ggsave("../figures/silhouette_curves_nullNetworks.png",
       width = 3000,
       height = 1500,
       units = "px",
       dpi = 600)

# ============================================================
# 5. CLUSTER COMPARISON
# ============================================================
cluster_compare <- bind_rows(
  results[["degree_null"]][["clusters"]],
  results[["connectance_null"]][["clusters"]]
)

cluster_wide <- cluster_compare %>%
  pivot_wider(names_from = model, values_from = Cluster)

# Adjusted Rand Index
ari <- adjustedRandIndex(
  cluster_wide$degree_null,
  cluster_wide$connectance_null
)

print(paste("Adjusted Rand Index:", round(ari, 3)))

# ============================================================
# 6. DENDROGRAM FUNCTION
# ============================================================
plot_dendrogram <- function(hc, title) {
  
  dend <- as.dendrogram(hc)
  dend_data <- dendro_data(dend)
  
  ggplot() +
    geom_segment(data = dend_data$segment,
                 aes(x = y, y = x,
                     xend = yend, yend = xend),
                 linewidth = 0.8, colour = "#001628") +
    labs(x = "Correlation Distance",
         y = "",
         title = title) +
    coord_cartesian(clip = "off") +
    figure_theme() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
}

# ============================================================
# 7. PLOT DENDROGRAMS
# ============================================================

prep_dendro <- function(hc, model_name) {
  
  dend <- as.dendrogram(hc) %>% dendro_data()
  
  segments <- dend$segment %>%
    glow_up(model = model_name)
  
  labels <- dend$label %>%
    glow_up(model = model_name)
  
  list(segments = segments, labels = labels)
}

deg <- prep_dendro(results[["degree_null"]][["hc"]], "Degree-null")
con <- prep_dendro(results[["connectance_null"]][["hc"]], "Connectance-null")

segments_df <- bind_rows(deg$segments, con$segments)
labels_df   <- bind_rows(deg$labels, con$labels)


label_df <- labels_df %>%
  left_join(cluster_df, by = c("label" = "Metric")) %>%
  mutate(Cluster = as.factor(Cluster))

cluster_all <- bind_rows(
  results$degree_null$clusters,
  results$connectance_null$clusters) %>%
  glow_up(model = case_when(model == "degree_null" ~ "Degree-null",
                            model == "connectance_null" ~ "Connectance-null"))

rect_df <- labels_df %>%
  left_join(cluster_all, by = c("label" = "Metric", "model" = "model")) %>%
  mutate(Cluster = as.factor(Cluster)) %>%
  group_by(model, Cluster) %>%
  summarise(
    ymin = min(x),
    ymax = max(x),
    .groups = "drop"
  ) %>%
  mutate(
    xmin = -Inf,
    xmax = 0
  )

ggplot() +
  geom_segment(
    data = segments_df,
    aes(x = y, y = x, xend = yend, yend = xend),
    linewidth = 0.8,
    colour = "#001628"
  ) +
  # Rectangles for clusters
  geom_rect(data = rect_df,
            aes(xmin = -0.01, xmax = 0,
                ymin = ymin-0.3, ymax = ymax+0.3),
            alpha = 0.8, fill = "#A2AAAD") +
  geom_text(
    data = label_df,
    aes(x = -0.015, y = x, label = label, colour = Cluster),
    hjust = 1,
    size = rel(2),
    family = "space",
    key_glyph = "point"
  ) +
  facet_wrap(~ model, scales = "free",
             ncol = 1) +
  scale_colour_manual(values = setNames(pal_df$colour, as.character(pal_df$value)),
                      labels = pal_df$label,
                      limits = pal_df$value,
                      name = "Module (empirical networks)") +
  labs(x = "Correlation Distance",
       y = "") +
  coord_cartesian(xlim = c(-0.1, max(segments_df$y))) +
  figure_theme() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "right")

ggsave("../figures/cluster_dendrogram_nullNetworks.png",
       width = 5900,
       height = 4700,
       units = "px",
       dpi = 600)

############################################################
# END
############################################################