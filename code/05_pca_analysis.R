library(tidyverse)
library(genzplyr)
library(ggalluvial)
library(ggrepel)
library(patchwork)
library(pheatmap)

set.seed(66)
setwd(here::here())
source("lib/plotting_themes.R")

##################################################
# 1. LOAD DATA
##################################################

predictors <- read.csv("data/cleaned/all_networks.csv") %>%
  as_tibble() %>%
  vibe_check(-c(ρ, complexity, robustness, control))

cluster_df <- read.csv("../tables/metric_clusters_auto.csv")

##################################################
# 2. PREPROCESS FOR PCA
##################################################

# Log-transform skewed count variables (optional but recommended)
pca_data <- predictors %>%
  glow_up(across(c(links, intervals, ChNum),
                 ~ log1p(.x)))

##################################################
# 3. RUN PCA
##################################################

pca_all <- prcomp(pca_data,
                  center = TRUE,
                  scale. = TRUE)

loadings <- pca_all$rotation
pc_scores <- pca_all$x

explained_var <- summary(pca_all)$importance[2, ]

# Save pca object
saveRDS(pca_all, "data/outputs/pca_object.rds")

##################################################
# 4. DEFINE MODULES
##################################################

clusters <- cluster_df$Cluster
names(clusters) <- cluster_df$Metric

modules <- unique(clusters)
pcs <- colnames(loadings)

##################################################
# 5. MODULE → PCA ALIGNMENT FUNCTION
##################################################

compute_module_alignment <- function(loadings, clusters) {
  
  modules <- unique(clusters)
  pcs <- colnames(loadings)
  
  out <- matrix(0,
                nrow = length(modules),
                ncol = length(pcs),
                dimnames = list(modules, pcs))
  
  for (m in modules) {
    
    metrics_in_module <- names(clusters[clusters == m])
    
    for (pc in pcs) {
      
      out[m, pc] <- sum(loadings[metrics_in_module, pc]^2,
                        na.rm = TRUE)
    }
  }
  
  return(out)
}

obs_alignment <- compute_module_alignment(loadings, clusters)

##################################################
# 6. PERMUTATION TEST
##################################################

set.seed(66)

n_perm <- 10000

perm_results <- array(0,
                      dim = c(nrow(obs_alignment),
                              ncol(obs_alignment),
                              n_perm))

for (p in seq_len(n_perm)) {
  
  perm_clusters <- sample(clusters)
  names(perm_clusters) <- names(clusters)
  
  perm_alignment <- compute_module_alignment(loadings,
                                             perm_clusters)
  
  perm_results[,,p] <- perm_alignment
}

# ---- Compute p-values ----
p_values <- matrix(NA,
                   nrow = nrow(obs_alignment),
                   ncol = ncol(obs_alignment),
                   dimnames = dimnames(obs_alignment))

for (i in seq_len(nrow(obs_alignment))) {
  for (j in seq_len(ncol(obs_alignment))) {
    
    p_values[i, j] <- mean(
      perm_results[i, j, ] >= obs_alignment[i, j]
    )
  }
}

# ---- Global test ----
obs_T  <- sum(obs_alignment^2)
perm_T <- apply(perm_results, 3,
                function(x) sum(x^2))

p_global <- mean(perm_T >= obs_T)

##################################################
# 7. Z-SCORES
##################################################

z_scores <- (obs_alignment -
               apply(perm_results, c(1,2), mean)) /
  apply(perm_results, c(1,2), sd)

sig_mat <- abs(z_scores) > 1.96

module_summary <- data.frame(
  Module = rownames(z_scores),
  Num_Significant_PCs = rowSums(sig_mat),
  Strongest_PC = apply(abs(z_scores), 1, which.max)
)

module_effect <- data.frame(
  Module = rownames(z_scores),
  Mean_Z = rowMeans(z_scores),
  Max_Z = apply(z_scores, 1, max)
)

# which PCs are module driven
pc_summary <- data.frame(
  PC = colnames(z_scores),
  Num_Significant_Modules = colSums(sig_mat),
  Strongest_Module = apply(abs(z_scores), 2, which.max)
)

##################################################
# 8. VARIANCE EXPLAINED
##################################################

obs_alignment <- matrix(0,
                        nrow = length(unique(clusters)),
                        ncol = ncol(loadings),
                        dimnames = list(
                          unique(clusters),
                          colnames(loadings)
                        ))

for (m in unique(clusters)) {
  
  metrics_in_module <- names(clusters[clusters == m])
  
  for (pc in colnames(loadings)) {
    
    obs_alignment[m, pc] <-
      sum(loadings[metrics_in_module, pc]^2,
          na.rm = TRUE)
  }
}

module_var_frac <- obs_alignment

for (pc in colnames(module_var_frac)) {
  module_var_frac[, pc] <- module_var_frac[, pc] /
    sum(module_var_frac[, pc])
}

module_var_df <- as.data.frame(module_var_frac) %>%
  rownames_to_column("Module") %>%
  pivot_longer(-Module,
               names_to = "PC",
               values_to = "Variance_Fraction")

module_var_df$PC <- factor(module_var_df$PC,
                           levels = unique(module_var_df$PC))

ggplot(data = module_var_df,
       aes(x = PC, y = Variance_Fraction, alluvium = as.character(Module))) +
  geom_alluvium(aes(fill = as.character(Module)),
                alpha = 0.9) +
  scale_fill_manual(values = setNames(pal_df$colour, as.character(pal_df$value)),
                    labels = pal_df$label,
                    limits = pal_df$value,
                    name = "Module") +
  coord_cartesian(clip = "on",
                  expand = FALSE) +
  labs(title = "Variance Explained by Module per PC",
       y = "Fraction of PC Variance",
       fill = "Module",
       x = "") +
  figure_theme() +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("../figures/variance_explained.png",
       width = 7000,
       height = 2750,
       units = "px",
       dpi = 600)

ggplot(module_var_df,
       aes(x = PC,
           y = Module,
           fill = Variance_Fraction)) +
  geom_tile() +
  scale_fill_continuous(
    palette = seattle_abyssal_gen(100),
    name = "Variance Explained"
  ) +
  figure_theme() +
  theme(legend.position = 'right')

ggsave("../figures/variance_heatmap.png",
       width = 7000,
       height = 2750,
       units = "px",
       dpi = 600)

##################################################
# 8. VISUALISE MODULE ALIGNMENT
##################################################

var_prop <- summary(pca_all)$importance[2, ]  # variance explained
cum_var  <- cumsum(var_prop)

threshold <- 0.8  # 80% cumulative variance

pcs_plot <- which(cum_var <= threshold)

# Always include the PC that crosses the threshold
if (max(pcs_plot) < length(cum_var)) {
  pcs_plot <- c(pcs_plot, max(pcs_plot) + 1)
}

pcs_plot <- pcs_plot[pcs_plot <= ncol(loadings)]
pcs_plot_names <- paste0("PC", pcs_plot)

alignment_mat <- obs_alignment[, pcs_plot_names, drop = FALSE] %>%
  as.tibble() %>%
  glow_up(value = as.factor(row_number())) %>%
  left_join(pal_df) %>%
  vibe_check(-value, -colour) %>%
  pivot_longer(-label)

ggplot(alignment_mat,
       aes(x = name,
           y = label,
           fill = value)) +
  geom_tile() +
  scale_fill_continuous(
    palette = seattle_abyssal_gen(100),
    name = "Module Alignment"
  ) +
  labs(x = NULL,
       y = NULL) +
  figure_theme() +
  theme(legend.position = 'right')

ggsave("../figures/module_alignment_heatmap.png",
       width = 6000,
       height = 2000,
       units = "px",
       dpi = 600)

##################################################
# 8.b VISUALISE FULL LOADING STRUCTURE
##################################################

# Extract metric order from dendrogram
metric_order <- cluster_df %>%
  arrange(Cluster) %>%
  pull(Metric)

loadings_mat <- pca_all$rotation[metric_order, ]

# save this for later pca plotting
saveRDS(loadings_mat,
        "data/outputs/pc_loadings.rds")

# Keep only PCs explaining >5% variance
explained_var <- summary(pca_all)$importance[2, ]
impactful_pcs <- which(explained_var > 0.05)

loadings_mat <- loadings_mat[, impactful_pcs, drop = FALSE]

loadings_df <- as.data.frame(loadings_mat) %>%
  rownames_to_column("Metric") %>%
  pivot_longer(
    -Metric,
    names_to = "PC",
    values_to = "Loading"
  ) %>%
  mutate(Module = clusters[Metric])

# Sidebar (Module annotation)
cluster_annotation <- data.frame(
  Metric = metric_order) %>%
  mutate(Module = clusters[Metric]) %>%
  mutate(Metric = factor(Metric, levels = metric_order)) %>%
  arrange(Module)

cluster_bounds <- cluster_annotation %>%
  group_by(Module) %>%
  summarise(
    y = max(as.numeric(Metric) + 0.5)) %>%
  add_row(Module = 0,
          y = 0.5)

p_sidebar <- ggplot(cluster_annotation,
                    aes(x = 1,
                        y = Metric,
                        fill = as.factor(Module))) +
  geom_tile(alpha = 0.8) +
  scale_fill_manual(values = setNames(pal_df$colour, as.character(pal_df$value)),
                    labels = pal_df$label,
                    limits = pal_df$value,
                    name = "Module") +
  theme_void() +
  theme(legend.position = "right",
        text = element_text(family = "space", color = "#001628"))

# Heatmap of loadings
p_loadings <- ggplot(loadings_df,
                     aes(x = PC,
                         y = factor(Metric, levels = metric_order),
                         fill = Loading)) +
  geom_tile() +
  geom_hline(yintercept = cluster_bounds$y,
             colour = "#001628",
             linewidth = 0.7)  +
  coord_cartesian(clip = "on",
                  expand = FALSE) +
  scale_fill_gradientn(
    colors = seattle_div,
    values = scales::rescale(c(-1, 0, 1))
  ) +
  theme_void() +
  theme(axis.text.x = element_text(family = "space", size = 10, color = "#001628"),
        text = element_text(family = "space", color = "#001628"))

p_loadings +
  p_sidebar + 
  plot_layout(widths = c(1, 0.05),
              guides = 'collect')

ggsave("../figures/pca_loadings_heatmap.png",
       width = 6000,
       height = 3000,
       units = "px",
       dpi = 600)

##################################################
# 9. PCA LOADING PLOT
##################################################

loadings_plot <- as.data.frame(loadings[,1:2]) %>%
  rownames_to_column("Metric") %>%
  mutate(Module = as.factor(clusters[Metric]))

var_explained <- summary(pca_all)$importance["Proportion of Variance", ]
pc1_label <- paste0("PC1 (", round(var_explained[1] * 100, 1), "%)")
pc2_label <- paste0("PC2 (", round(var_explained[2] * 100, 1), "%)")

pca_plot <- ggplot(loadings_plot, 
                   aes(PC1, PC2, 
                       colour = Module, 
                       fill = Module, 
                       label = Metric)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "#A5ACAF") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "#A5ACAF") +
  geom_point(size = 3,
             alpha = 0.75,
             shape = 21,
             colour = "white") +
  geom_text_repel(size = rel(4),
                  family = "space",
                  show.legend = FALSE) +
  labs(x = pc1_label,
       y = pc2_label) +
  scale_fill_manual(values = setNames(pal_df$colour, as.character(pal_df$value)),
                    labels = pal_df$label,
                    limits = pal_df$value,
                    name = "Module") +
  scale_colour_manual(values = setNames(pal_df$colour, as.character(pal_df$value)),
                      labels = pal_df$label,
                      limits = pal_df$value,
                      name = "Module") +
  figure_theme() +
  theme(legend.position = 'right',
        panel.grid.major = element_blank())

ggsave("../figures/pca_loadings.png",
       plot = pca_plot,
       width = 6000, 
       height = 4000,
       units = "px",
       dpi = 600)

##################################################
# 9. PCA LOADING PLOT - hulls
##################################################

# Compute convex hulls for each module
hulls <- loadings_plot %>%
  group_by(Module) %>%
  slice(chull(PC1, PC2)) %>%
  ungroup()

pca_hull_plot <- ggplot(loadings_plot, 
                        aes(PC1, PC2, 
                            colour = Module, 
                            fill = Module, 
                            label = Metric)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "#A5ACAF") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "#A5ACAF") +
  geom_polygon(data = hulls, aes(x = PC1, y = PC2, group = Module),
               alpha = 0.2, colour = NA) +
  geom_point(size = 3, alpha = 0.75, shape = 21, colour = "white") +
  geom_text_repel(size = rel(4),
                  family = "space",
                  show.legend = FALSE) +
  labs(x = pc1_label, y = pc2_label) +
  scale_fill_manual(values = setNames(pal_df$colour, as.character(pal_df$value)),
                    labels = pal_df$label,
                    limits = pal_df$value,
                    name = "Module") +
  scale_colour_manual(values = setNames(pal_df$colour, as.character(pal_df$value)),
                      labels = pal_df$label,
                      limits = pal_df$value,
                      name = "Module") +
  figure_theme() +
  theme(legend.position = 'right',
        panel.grid.major = element_blank())

ggsave("../figures/pca_loadings_hulls.png",
       width = 6000, 
       height = 4000, 
       units = "px",
       dpi = 600)

##################################################
# 10. Z-SCORES HEATMAP
##################################################

sig_mask <- (p_values < 0.05) & (abs(z_scores) > 1.96)

sig_mask <- ifelse(sig_mask, "*", "")

# Use Z-scores for visualization (better than raw alignment)
plot_mat <- z_scores %>%
  as.tibble() %>%
  glow_up(value = as.factor(row_number())) %>%
  left_join(pal_df) %>%
  vibe_check(-value, -colour) %>%
  pivot_longer(-label) %>%
  glow_up(number = str_extract(name, "[[:digit:]]")) %>%
  slay() %>%
  glow_up(name = fct_inorder(name))

max_val <- max(abs(plot_mat$value))

# Create annotation of significance
annotation_matrix <- sig_mask %>%
  as.tibble() %>%
  glow_up(value = as.factor(row_number())) %>%
  left_join(pal_df) %>%
  vibe_check(-value, -colour) %>%
  pivot_longer(-label) %>%
  yeet(value == "*")

ggplot(plot_mat) +
  geom_tile(aes(x = name,
                y = label,
                fill = value)) +
  geom_text(data = annotation_matrix,
            aes(x = name,
                y = label),
            colour = "#001628",
            label = "*") +
  scale_fill_gradient2(
    low = seattle_div[1],
    mid = seattle_div[2],
    high = seattle_div[3],
    limits = c(-max_val, max_val),
    name = "Z-score"
  ) +
  labs(x = NULL,
       y = NULL) +
  figure_theme() +
  theme(legend.position = 'right')

ggsave("../figures/heatmap_zscores.png",
       width = 7000,
       height = 2500,
       units = "px",
       dpi = 600)

joint_sig <- (p_values < 0.05) & (abs(z_scores) > 1.96)

sig_summary <- data.frame(
  Module = rownames(z_scores),
  Num_Significant_PCs = rowSums(joint_sig),
  Mean_Z = rowMeans(z_scores),
  Max_Z = apply(z_scores, 1, max)
)


