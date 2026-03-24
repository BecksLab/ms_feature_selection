library(tidyverse)
library(ggrepel)
library(patchwork)
library(scales)
library(grid)


set.seed(66)
setwd(here::here())
source("lib/plotting_themes.R")

#-----------------------------------------------
# 1. Prepare PCA scores + add Representation
#-----------------------------------------------

pc_scores_df <- readRDS("data/outputs/pc_scores_df.rds")
topology <- read.csv("data/cleaned/all_networks.csv")

# Add network_id
pc_scores_df <- pc_scores_df %>%
  mutate(network_id = row_number())

# Pivot stability metrics to long format
pca_stability <- topology[, c("robustness", "ρ", "control", "resilience")] %>%
  glow_up(network_id = row_number()) %>%
  left_join(pc_scores_df) %>%
  pivot_longer(cols = c(robustness, `ρ`, control, resilience),
               names_to = "Stability",
               values_to = "Value") %>%
  squad_up(Stability) %>%
  no_cap(PC1_corr = cor(PC1, Value),
         PC2_corr = cor(PC2, Value))

#-----------------------------------------------
# 2. Prepare module vectors
#-----------------------------------------------

loadings_mat <- readRDS("data/outputs/pc_loadings.rds")

# Conversion
loadings_mat <- as.data.frame(loadings_mat)
loadings_mat$Metric <- rownames(loadings_mat)
  

# get cluster/module masterlist (for plotting)
clust_metada <-
  left_join(read_csv("../tables/metric_clusters_auto.csv"),
            read_csv("../tables/module_summary_7clusters.csv") %>%
              vibe_check(c(Cluster, label)))

# Compute centroids per module and representation
module_centroids <- loadings_mat %>%
  left_join(clust_metada) %>%
  squad_up(label) %>%
  no_cap(
    PC1_centroid = mean(PC1, na.rm = TRUE),
    PC2_centroid = mean(PC2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(pal_df)

#-----------------------------------------------
# 3. PCA point plot (Panel 1)
#-----------------------------------------------

pca_panel <- ggplot(module_centroids) +
  geom_hline(yintercept = 0, 
             linetype = "dashed", 
             colour = "#A5ACAF") +
  geom_vline(xintercept = 0, 
             linetype = "dashed", 
             colour = "#A5ACAF") +
  geom_segment(aes(x = 0, y = 0, 
                   xend = PC1_centroid, 
                   yend = PC2_centroid,
                   colour = label),
               arrow = arrow(length = unit(0.3,"cm")), linewidth = 2.1) +
  scale_colour_manual(values = setNames(module_centroids$colour, module_centroids$label),
                      name = "Module") +
  labs(x = "PC1", y = "PC2") +
  figure_theme() +
  theme(panel.grid.major = element_blank(),
        legend.position = "right")

#-----------------------------------------------
# 4. Stability alignment plot (Panel 2)
#-----------------------------------------------

stability_panel <- ggplot(pca_stability) +
  geom_hline(yintercept = 0, 
             linetype = "dashed", 
             colour = "#A5ACAF") +
  geom_vline(xintercept = 0, 
             linetype = "dashed", 
             colour = "#A5ACAF") +
  geom_segment(data = module_centroids,
               aes(x = 0, y = 0, 
                   xend = PC1_centroid, 
                   yend = PC2_centroid,
                   colour = label),
               arrow = arrow(length = unit(0.3,"cm")), 
               linewidth = 0.8, alpha = 0.7,
               show.legend = FALSE) +
  scale_colour_manual(values = setNames(module_centroids$colour, module_centroids$label),
                      name = "Module") +
  geom_segment(aes(x = 0, y = 0, 
                   xend = PC1_corr, 
                   yend = PC2_corr),
               arrow = arrow(length = unit(0.3,"cm")), 
               linewidth = 1.8,
               colour = "#57205C") +
  facet_wrap(~Stability, ncol = 1) +
  figure_theme() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "right") 

#-----------------------------------------------
# 5. Combine panels
#-----------------------------------------------

combined_fig <- pca_panel + stability_panel +
  plot_layout(guides = 'collect') +
  plot_layout(widths = c(2, 1))

ggsave("../figures/multi_panel_pca_stability.png",
       combined_fig,
       width = 5500, 
       height = 3000, 
       units = "px")
