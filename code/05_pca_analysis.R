library(tidyverse)
library(ggrepel)
library(patchwork)
library(pheatmap)

set.seed(66)

##################################################
# 1. LOAD DATA
##################################################

predictors <- read.csv("data/cleaned/all_networks.csv") %>%
  as_tibble() %>%
  vibe_check(-c(ρ, complexity, robustness))

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

n_perm <- 1000

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
p_values <- apply(perm_results,
                  c(1,2),
                  function(x)
                    mean(x >= obs_alignment))

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

##################################################
# 8. VISUALISE MODEL ALIGNMENT
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

alignment_mat <- obs_alignment[, pcs_plot_names, drop = FALSE]

hm <- pheatmap(
  as.matrix(alignment_mat),
  color = seattle_abyssal_gen(100),
  cluster_rows = TRUE,
  cluster_cols = TRUE)

ggsave("../figures/module_alignment_heatmap.png",
       plot = hm$gtable,
       width = 6000,
       height = 3000,
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
  scale_fill_manual(values = cluster_colors,
                    name = "Module") +
  theme_void() +
  theme(legend.position = "right")

# Heatmap of loadings
p_loadings <- ggplot(loadings_df,
                     aes(x = PC,
                         y = factor(Metric, levels = metric_order),
                         fill = Loading)) +
  geom_tile() +
  geom_hline(yintercept = cluster_bounds$y,
             colour = "#1A1A1A",
             linewidth = 0.4) +
  scale_fill_gradient2(
    low = "#2C3E50",
    mid = "white",
    high = "#C1785A",
    midpoint = 0
  ) +
  theme_void() +
  theme(axis.text.x = element_text(family = "space", size = 10, color = "#1A1A1A"))

p_loadings +
  p_sidebar + 
  plot_layout(widths = c(1, 0.08),
              guides = 'collect')

ggsave("../figures/pca_loadings_heatmap.png",
       width = 6000,
       height = 3000,
       units = "px",
       dpi = 600)

##################################################
# 9. PCA LOADING PLOT (MODULE COLORED)
##################################################

loadings_plot <- as.data.frame(loadings[,1:2]) %>%
  rownames_to_column("Metric") %>%
  mutate(Module = as.factor(clusters[Metric]))

ggplot(loadings_plot, 
       aes(PC1, PC2, 
           color = Module, 
           label = Metric)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "#A5ACAF") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "#A5ACAF") +
  geom_point(size = 2.4,
             alpha = 0.8) +
  geom_text_repel(vjust = 1.2, size = 3,
                  family = "space") +
  scale_colour_manual(values = as.vector(kraken_7)) +
  figure_theme() +
  theme(legend.position = 'right',
        panel.grid.major = element_blank())

ggsave("../figures/pca_loadings.png",
       width = 5000,
       height = 3000,
       units = "px",
       dpi = 600)
