library(tidyverse)
library(ggrepel)
library(patchwork)
library(scales)
library(grid)

#-----------------------------------------------
# 1. Prepare PCA scores + add Representation
#-----------------------------------------------

pc_scores_df <- readRDS("data/outputs/pc_scores_df.rds")

# Add network_id
pc_scores_df <- pc_scores_df %>%
  mutate(network_id = row_number())

# Pivot stability metrics to long format
pca_stability <- topology[, c("robustness", "ρ", "complexity", "control")] %>%
  mutate(network_id = row_number()) %>%
  pivot_longer(cols = c(robustness, `ρ`, complexity, control),
               names_to = "Stability",
               values_to = "Value") %>%
  left_join(pc_scores_df, by = "network_id") %>%
  group_by(Stability) %>%
  mutate(Value_scaled = scales::rescale(Value, to = c(0,1))) %>%
  ungroup()

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
               arrow = arrow(length = unit(0.3,"cm")), linewidth = 1.2) +
  scale_colour_manual(values = setNames(module_centroids$colour, module_centroids$label),
                      name = "Module") +
  labs(x = "PC1", y = "PC2", title = "Panel 1: PCA structural space with module loadings") +
  figure_theme() +
  theme(legend.position = "right")

#-----------------------------------------------
# 4. Stability alignment plot (Panel 2)
#-----------------------------------------------

# Compute alignment of modules with stability
stability_vectors <- pca_stability %>%
  squad_up(Stability) %>%
  no_cap(
    PC1_centroid = mean(PC1, na.rm = TRUE),
    PC2_centroid = mean(PC2, na.rm = TRUE),
    .groups = "drop"
  )

stability_panel <- ggplot(stability_vectors) +
  geom_hline(yintercept = 0, 
             linetype = "dashed", 
             colour = "#A5ACAF") +
  geom_vline(xintercept = 0, 
             linetype = "dashed", 
             colour = "#A5ACAF") +
  geom_segment(aes(x = 0, y = 0, 
                   xend = PC1_centroid, 
                   yend = PC2_centroid),
               arrow = arrow(length = unit(0.3,"cm")), 
               linewidth = 1.2,
               colour = "#001628") +
  facet_wrap(~Stability, ncol = 1) +
  labs(x = "Alignment with stability") +
  figure_theme() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "right") 

#-----------------------------------------------
# 5. Combine panels
#-----------------------------------------------

combined_fig <- pca_panel + stability_panel +
  plot_layout(guides = 'collect')

ggsave("../figures/multi_panel_pca_stability.png",
       combined_fig,
       width = 5500, height = 4500, units = "px")
