library(tidyverse)
library(ggrepel)
library(scales)
library(grid)

#---------------------------------------------------
# 1. Prepare PCA scores
#---------------------------------------------------
# pc_scores_df: rows = networks, columns = PC1, PC2, ..., network_id
# topology: network-level stability metrics
# coefficients_all + clust_metada + pal_df: module contribution info

pca_stability <- pc_scores_df %>%
  glow_up(network_id = row_number()) %>%
  left_join(topology[, c("robustness", "ρ", "complexity", "control")] %>%
              glow_up(network_id = row_number()), 
            by = "network_id")  %>%
  pivot_longer(cols = c(robustness, ρ, complexity, control),
               names_to = "Stability",
               values_to = "Value") %>%
  group_by(Stability) %>%
  mutate(
    Value_scaled = scales::rescale(Value, to = c(0,1))  # min-max 0-1 within Stability
  ) %>%
  ungroup()

#---------------------------------------------------
# 2. Compute module vectors for PCA (optional)
#---------------------------------------------------
# We compute the projection of each module onto PC1 and PC2
# using standardized elastic net coefficients.

# import pca object
pca <- readRDS("data/outputs/pca_object.rds")

pca_loadings <- as.data.frame(pca$rotation)
pca_loadings$Metric <- rownames(pca_loadings)

module_vectors_pc <- coefficients_all %>%
  filter(Representation == "PC_Scores") %>%
  mutate(Module = Predictor,      # each PC score is its own “module”
         # Assign unit vector along PC axes
         PC1_vec = ifelse(Module == "PC1", 1, 0),
         PC2_vec = ifelse(Module == "PC2", 1, 0)) %>%
  group_by(Representation, Stability, Module) %>%
  summarise(
    PC1_vec = sum(Coefficient * PC1_vec),
    PC2_vec = sum(Coefficient * PC2_vec),
    .groups = "drop"
  )

module_vectors_others <- coefficients_all %>%
  filter(Representation != "PC_Scores") %>%
  left_join(clust_metada, by = join_by(Predictor == Metric)) %>%
  mutate(Module = if_else(is.na(label), "PCA Axis", label)) %>%
  left_join(pca_loadings, by = c("Predictor" = "Metric")) %>%
  filter(!is.na(PC1) & !is.na(PC2)) %>%
  group_by(Representation, Stability, Module) %>%
  summarise(
    PC1_vec = sum(Coefficient * PC1),
    PC2_vec = sum(Coefficient * PC2),
    .groups = "drop"
  ) %>%
  group_by(Representation, Stability) %>%
  mutate(
    PC1_vec_raw = PC1_vec,
    PC2_vec_raw = PC2_vec,
    vec_len_raw = sqrt(PC1_vec^2 + PC2_vec^2)
  ) %>%
  ungroup()

# Compute a global scaling factor
max_len_global <- max(module_vectors_others$vec_len_raw, na.rm = TRUE)
desired_max <- 15  # length in PCA units
scale_factor_global <- desired_max / max_len_global

module_vectors_all <- module_vectors_others %>%
  mutate(
    PC1_vec_scaled = PC1_vec_raw * scale_factor_global,
    PC2_vec_scaled = PC2_vec_raw * scale_factor_global
  ) %>%
  left_join(pal_df, by = join_by(Module == label)) %>%
  glow_up(label = if_else(is.na(colour), "PCA Axis", Module),
          colour = if_else(is.na(colour), pca_col, colour))

# direction stability increases in 
stability_vectors <- pca_stability %>%
  group_by(Stability) %>%
  summarise(
    PC1_dir = cor(Value, PC1, use = "complete.obs"),
    PC2_dir = cor(Value, PC2, use = "complete.obs"),
    .groups = "drop"
  )

module_vectors <- pca_loadings %>%
  mutate(Module = clusters[Metric]) %>%
  group_by(Module) %>%
  summarise(
    PC1_vec = mean(PC1),
    PC2_vec = mean(PC2),
    .groups = "drop"
  )

stability_module_effect <- stability_vectors %>%
  crossing(module_vectors) %>%
  mutate(
    effect = PC1_dir * PC1_vec +
      PC2_dir * PC2_vec,
    Module = as.factor(Module)) %>%
  group_by(Stability) %>%
  # sort modules by effect within each Stability facet
  mutate(Module_order = reorder_within(Module, -effect, Stability)) 

ggplot(stability_module_effect,
       aes(x = effect,
           y = Module_order,
           fill = Module,
           colour = Module)) +
  geom_segment(aes(xend = 0),
               linewidth = 0.7) +
  geom_vline(xintercept = 0,
             colour = "#001628",
             linewidth = 0.4) +
  geom_point(size = 6,
             shape = 21,
             colour = "white") +
  scale_fill_manual(values = setNames(pal_df$colour, as.character(pal_df$value)),
                    labels = pal_df$label,
                    limits = pal_df$value,
                    name = "Module") +
  scale_colour_manual(values = setNames(pal_df$colour, as.character(pal_df$value)),
                      labels = pal_df$label,
                      limits = pal_df$value,
                      name = "Module") +
  labs(x= "Alignement with stability") +
  facet_wrap(~Stability,
             scales = "free_y") +
  figure_theme() +
  theme(legend.position = "right",
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank())

ggsave("../figures/stability_alignment.png",
       width = 5000,
       height = 2000,
       units = "px")

#---------------------------------------------------
# 3. PCA + Stability plot
#---------------------------------------------------

pca_stab_plot <- ggplot(pca_stability, 
                        aes(x = PC1, y = PC2)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "#A5ACAF") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "#A5ACAF") +
  # networks as points
  geom_point(aes(size = Value_scaled), 
             alpha = 0.7, shape = 21, fill = NA) +
  scale_size_continuous(range = c(3,8)) +
  # optional module arrows
  geom_segment(data = module_vectors_all,
               aes(x = 0, y = 0, 
                   xend = PC1_vec_scaled, 
                   yend = PC2_vec_scaled, 
                   colour = Module),
               arrow = arrow(length = unit(0.2,"cm")),  # slightly bigger arrowhead
               size = 1.2) +
  geom_segment(data = stability_vectors,
               aes(x = 0, y = 0,
                   xend = PC1_dir * 3,
                   yend = PC2_dir * 3),
               arrow = arrow(length = unit(0.3,"cm")),
               colour = "black",
               linewidth = 1.7) +
  scale_colour_manual(values = setNames(module_vectors_all$colour, 
                                        as.character(module_vectors_all$Module)),
                      breaks = module_vectors_all$label[module_vectors_all$label != "PCA Axis"],
                      name = "Module") +
  scale_fill_gradient(low = seattle_anchors[1], high = seattle_anchors[5]) +
  guides(#fill = guide_legend(title = "Stabilty Value (scaled)"), 
    size = guide_legend(title = "Stabilty Value (scaled)")) + 
  facet_grid(Stability ~ Representation) +
  labs(
    x = "PC1",
    y = "PC2",
    title = "Food-web structural variation and stability landscape"
  ) +
  figure_theme() +
  theme(legend.position = "right",
        panel.grid.major = element_blank())

ggsave("../figures/pca_struct_stability.png",
       pca_stab_plot,
       width = 5500,
       height = 4500,
       units = "px")

#---------------------------------------------------
# 4. Notes:
#---------------------------------------------------
# - Points show networks in structural PCA space.
# - Colour/size of points indicates the magnitude of the stability metric (robustness, ρ, complexity).
# - Arrows show contributions of modules to PCA axes, approximating which structural domains drive stability.
# - Facets separate stability metrics and structural representations (Cluster medoids, PCA axes, etc.).