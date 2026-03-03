# ================================================
# 05_pca_analysis_MODIFIED.R
# ================================================

library(FactoMineR)
library(factoextra)
library(tidyverse)
library(ggrepel)

# Load datasets
predictors <- read.csv("data/cleaned/all_networks.csv") %>%
  as_tibble() %>%
  vibe_check(-c(ρ, complexity, robustness))

# scale
metrics_scaled <- scale(predictors)

# 2. Run the PCA on Predictors only
pca_all <- prcomp(metrics_scaled, center = TRUE, scale. = FALSE)

summary(pca_all)

loadings <- pca_all$rotation
pc_scores <- pca_all$x

# Suppose you use absolute clustering at k = 5
clusters <- cutree(hc_signed, k = 7)

module_pc_strength <- data.frame()

for (m in unique(clusters)) {
  
  metrics_in_module <- names(clusters[clusters == m])
  
  # Subset loadings
  sub_loadings <- loadings[metrics_in_module, , drop = FALSE]
  
  # Compute mean absolute loading per PC
  strength <- colMeans(abs(sub_loadings), na.rm = TRUE)
  
  module_pc_strength <- rbind(
    module_pc_strength,
    data.frame(Module = m, t(strength))
  )
}

module_pc_strength

module_contribution <- data.frame()

for (m in unique(clusters)) {
  
  metrics_in_module <- names(clusters[clusters == m])
  
  # Sum squared loadings
  contrib <- colSums(loadings[metrics_in_module, ]^2, na.rm = TRUE)
  
  module_contribution <- rbind(
    module_contribution,
    data.frame(Module = m, t(contrib))
  )
}

module_contribution

library(pheatmap)

mat <- as.matrix(module_pc_strength[,-1])
rownames(mat) <- paste0("Module_", module_pc_strength$Module)

pheatmap(mat,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Module Alignment with PCA Axes")

library(ape)

# Get variance explained
var_exp <- summary(pca_all)$importance[2, ]  # Proportion variance
cum_var  <- cumsum(var_exp)

# Choose threshold
threshold <- 0.8

pcs_keep <- which(cum_var <= threshold)

# Always include the last one that crosses threshold
if (max(pcs_keep) < length(cum_var)) {
  pcs_keep <- c(pcs_keep, max(pcs_keep) + 1)
}

pcs_keep

elbow <- which(diff(diff(cum_var)) == max(diff(diff(cum_var))))
elbow

# Subset module-PC matrix to impactful PCs
mat_trim <- module_pc_strength[, c("Module", paste0("PC", pcs_keep))]

mat_plot <- as.matrix(mat_trim[,-1])
rownames(mat_plot) <- paste0("Module_", mat_trim$Module)

pheatmap(mat_plot,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         display_numbers = TRUE,
         main = paste0("Module Alignment with Top PCs (",
                       round(threshold*100), "% variance)"))

df_plot <- data.frame(
  PC1 = loadings[,1],
  PC2 = loadings[,2],
  Metric = rownames(loadings),
  Module = factor(clusters[rownames(loadings)])
)

ggplot(df_plot, aes(PC1, PC2, color = Module, label = Metric)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70") +
  geom_segment(
    aes(x = 0, y = 0, xend = PC1, yend = PC2, colour = Module),
    arrow = arrow(length = unit(0.01, "npc")),
    linewidth = 0.9
  ) +
  geom_text_repel(vjust = 1.2, size = 3,
                  family = "space") +
  scale_colour_manual(values = as.vector(kraken_7)) +
  figure_theme()

obs_mat <- as.matrix(module_pc_strength[,-1])

compute_alignment <- function(loadings, clusters) {
  
  modules <- unique(clusters)
  pcs <- colnames(loadings)
  
  out <- matrix(0, nrow = length(modules), ncol = length(pcs))
  rownames(out) <- modules
  colnames(out) <- pcs
  
  for (m in modules) {
    
    metrics_in_module <- names(clusters[clusters == m])
    
    for (pc in pcs) {
      
      out[m, pc] <- sum(loadings[metrics_in_module, pc]^2, na.rm = TRUE)
    }
  }
  
  return(out)
}

obs_alignment <- compute_alignment(loadings, clusters)

set.seed(66)

n_perm <- 1000
perm_results <- array(0,
                      dim = c(nrow(obs_alignment),
                              ncol(obs_alignment),
                              n_perm))

for (p in 1:n_perm) {
  
  # Shuffle module labels
  perm_clusters <- sample(clusters)
  names(perm_clusters) <- names(clusters)
  
  perm_alignment <- compute_alignment(loadings, perm_clusters)
  
  perm_results[,,p] <- perm_alignment
}


p_values <- matrix(1,
                   nrow = nrow(obs_alignment),
                   ncol = ncol(obs_alignment))

rownames(p_values) <- rownames(obs_alignment)
colnames(p_values) <- colnames(obs_alignment)

for (i in 1:nrow(obs_alignment)) {
  for (j in 1:ncol(obs_alignment)) {
    
    perm_dist <- perm_results[i, j, ]
    
    p_values[i, j] <- mean(perm_dist >= obs_alignment[i, j])
  }
}


obs_T <- sum(obs_alignment^2)

perm_T <- apply(perm_results, 3, function(x) sum(x^2))

p_global <- mean(perm_T >= obs_T)

p_global


z_scores <- matrix(0,
                   nrow = nrow(obs_alignment),
                   ncol = ncol(obs_alignment))

rownames(z_scores) <- rownames(obs_alignment)
colnames(z_scores) <- colnames(obs_alignment)

for (i in 1:nrow(obs_alignment)) {
  for (j in 1:ncol(obs_alignment)) {
    
    perm_dist <- perm_results[i, j, ]
    
    z_scores[i, j] <- (obs_alignment[i, j] - mean(perm_dist)) /
      sd(perm_dist)
  }
}

z_mat <- as.matrix(z_scores)

sig_mat <- abs(z_mat) > 1.96

module_summary <- data.frame(
  Module = 1:nrow(z_mat),
  Num_Significant_PCs = rowSums(sig_mat),
  Strongest_PC = apply(abs(z_mat), 1, which.max)
)

module_summary


dominant_module_per_pc <- apply(abs(z_mat), 2, which.max)

dominant_module_per_pc


summary_table <- data.frame(
  Module = 1:nrow(z_mat),
  Dominant_PC = apply(abs(z_mat), 1, which.max),
  Dominant_Z = apply(z_mat, 1, function(x) x[which.max(abs(x))]),
  Num_Sig_PCs = rowSums(sig_mat)
)

summary_table

# Extract metric order from dendrogram
metric_order <- label_df %>%
  arrange(x) %>%
  pull(label)

# Reorder PCA loading matrix
loadings_mat <- pca_all$rotation
loadings_mat <- loadings_mat[metric_order, ]

explained_var <- summary(pca_all)$importance[2, ]
impactful_pcs <- which(explained_var > 0.05)  # >5% variance
loadings_mat <- loadings_mat[, impactful_pcs]

loadings_df <- as.data.frame(loadings_mat)
loadings_df$Metric <- rownames(loadings_df)

loadings_long <- loadings_df %>%
  pivot_longer(-Metric,
               names_to = "PC",
               values_to = "Loading")

# Ensure consistent ordering
loadings_long$Metric <- factor(loadings_long$Metric,
                               levels = metric_order)

cluster_df <- read.csv("../tables/metric_clusters_auto.csv")

cluster_annotation <- cluster_df %>%
  mutate(Metric = factor(Metric, levels = metric_order)) %>%
  arrange(Metric)

cluster_bounds <- cluster_annotation %>%
  group_by(Cluster) %>%
  summarise(
    y = max(as.numeric(Metric) + 0.5)) %>%
  add_row(Cluster = 0,
          y = 0.5)

p_sidebar <- ggplot(cluster_annotation,
                    aes(x = 1,
                        y = Metric,
                        fill = factor(Cluster))) +
  geom_tile() +
  scale_fill_manual(values = cluster_colors,
                    labels = module_names$Module_Name,
                    name = "Module") +
  theme_void() +
  theme(legend.position = "right")


pca_loadings <- ggplot(loadings_long,
                       aes(x = PC,
                           y = Metric,
                           fill = Loading)) +
  geom_tile() +
  geom_hline(yintercept = cluster_bounds$y,
             colour = "#1A1A1A") +
  scale_fill_gradient2(low = "#2C3E50",
                       mid = "white",
                       high = "#C1785A",
                       midpoint = 0) +
  labs(title = "PCA Loadings",
       x = "Principal Components",
       y = "") +
  theme_void() +
  theme(axis.text.x = element_text(family = "space", size = 10, color = "#1A1A1A"))

pca_loadings +
  p_sidebar + 
  plot_layout(widths = c(1, 0.08),
              guides = 'collect')

ggsave("../figures/pca_loadings_heatmap.png",
       width = 6000,
       height = 3000,
       units = "px",
       dpi = 600)
