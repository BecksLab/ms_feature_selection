
library(genzplyr)
library(tidyverse)
library(patchwork)
library(reshape2)

set.seed(66)
setwd(here::here())
source("lib/plotting_themes.R")

# ============================================================
# 0. LOAD DATA
# ============================================================
metrics <- read.csv("data/cleaned/all_networks.csv") %>%
  as_tibble() %>%
  # Exclude stability response variables from clustering
  vibe_check(-c(ρ, complexity, robustness, control))

# import pca object
pca <- readRDS("data/outputs/pca_object.rds")

# get cluster/module masterlist
clust_metada <-
  left_join(read_csv("../tables/metric_clusters_auto.csv"),
            read_csv("../tables/module_summary_7clusters.csv") %>%
              vibe_check(c(Cluster, label)))

############################################################
# 1) Extract Cluster Medoids
############################################################

extract_cluster_medoids <- function(metrics_matrix, metric_cluster_df){
  
  medoid_list <- list()
  
  for (cl in unique(metric_cluster_df$Cluster)) {
    
    metrics_in_cluster <- metric_cluster_df %>%
      filter(Cluster == cl) %>%
      pull(Metric)
    
    sub_mat <- metrics_matrix[, metrics_in_cluster, drop = FALSE]
    
    # Correlation-based distance
    cor_mat <- cor(sub_mat, use = "pairwise.complete.obs")
    dist_mat <- 1 - abs(cor_mat)
    
    # Mean distance per metric
    mean_dist <- rowMeans(dist_mat)
    
    medoid_metric <- names(which.min(mean_dist))
    
    medoid_list[[cl]] <- medoid_metric
  }
  
  medoids <- unlist(medoid_list)
  
  return(medoids)
}

cluster_medoids <- extract_cluster_medoids(metrics,
                                           clust_metada)

############################################################
# 2) Extract PC-Dominant Metrics
############################################################

extract_pc_dominant_metrics <- function(pca_obj,
                                        variance_threshold = 0.80){
  
  # Variance explained
  var_explained <- summary(pca_obj)$importance[2, ]
  cum_var <- cumsum(var_explained)
  
  # PCs to retain
  pcs_keep <- which(cum_var <= variance_threshold)
  
  # Ensure we include first PC exceeding threshold
  if(max(pcs_keep) < length(cum_var)){
    pcs_keep <- seq_len(max(pcs_keep) + 1)
  }
  
  loadings <- pca_obj$rotation[, pcs_keep, drop = FALSE]
  
  dominant_metrics <- apply(loadings, 2, function(pc){
    names(which.max(abs(pc)))
  })
  
  return(unique(dominant_metrics))
}

pc_dominant_metrics <- extract_pc_dominant_metrics(pca,
                                                   variance_threshold = 0.80)

############################################################
# 3) Extract PC Scores (≥80% variance)
############################################################

extract_pc_scores <- function(pca_obj,
                              variance_threshold = 0.80){
  
  var_explained <- summary(pca_obj)$importance[2, ]
  cum_var <- cumsum(var_explained)
  
  pcs_keep <- which(cum_var <= variance_threshold)
  
  if(max(pcs_keep) < length(cum_var)){
    pcs_keep <- seq_len(max(pcs_keep) + 1)
  }
  
  pc_scores <- as.data.frame(pca_obj$x[, pcs_keep, drop = FALSE])
  
  return(pc_scores)
}

pc_scores_df <- extract_pc_scores(pca,
                                  variance_threshold = 0.80)

############################################################
# 4) Export outputs
############################################################

write.csv(cluster_medoids, "data/outputs/cluster_medoids.csv",
          row.names = FALSE)
write.csv(pc_dominant_metrics, "data/outputs/pc_dominant_metrics.csv",
          row.names = FALSE)
write.csv(pc_scores_df, "data/outputs/pc_scores_df.csv",
          row.names = FALSE)

############################################################
# 5) Sanity check
############################################################

############################################################
# Dimensional Reduction Characterisation Engine
############################################################

characterise_reduction <- function(metric_matrix,
                                   cluster_set,
                                   pc_dominant_set,
                                   pc_scores_df){
  
  results <- list()
  
  ##########################################################
  # Helper Function: Compute Diagnostics for One Set
  ##########################################################
  
  compute_set_metrics <- function(data_subset, full_matrix){
    
    data_subset <- full_matrix[, data_subset, drop = FALSE]
    data_subset <- scale(data_subset)
    
    # ---- Variance Retained ----
    pca_full <- prcomp(full_matrix, scale. = TRUE)
    total_var <- sum(pca_full$sdev^2)
    
    var_retained <- sum(apply(data_subset, 2, var)) / total_var
    
    # ---- Internal Redundancy ----
    if (ncol(data_subset) > 1){
      cor_mat <- cor(data_subset, use = "pairwise.complete.obs")
      mean_abs_cor <- mean(abs(cor_mat[upper.tri(cor_mat)]))
    } else {
      mean_abs_cor <- NA
    }
    
    # ---- Effective Dimensionality (SVD Spectrum) ----
    svd_obj <- svd(data_subset)
    sing_vals <- svd_obj$d^2 / sum(svd_obj$d^2)
    
    cumulative_var <- cumsum(sing_vals)
    
    eff_dim <- min(which(cumulative_var >= 0.8))
    
    return(list(
      variance_retained = var_retained,
      mean_abs_correlation = mean_abs_cor,
      singular_values = sing_vals,
      cumulative_variance = cumulative_var,
      effective_dimensionality = eff_dim
    ))
  }
  
  ##########################################################
  # 1. Compute Metrics For Each Representation
  ##########################################################
  
  results$cluster <- compute_set_metrics(cluster_set, metric_matrix)
  
  results$pc_dominant <- compute_set_metrics(pc_dominant_set, metric_matrix)
  
  # PC scores already computed externally
  results$pc_scores <- compute_set_metrics(colnames(pc_scores_df),
                                           cbind(metric_matrix, pc_scores_df))
  
  ##########################################################
  # 2. Cross-Representation Similarity
  ##########################################################
  
  compute_rep_similarity <- function(mat1, mat2){
    
    # Standardize each representation independently
    X1 <- scale(mat1)
    X2 <- scale(mat2)
    
    # Compute distance matrices
    d1 <- as.vector(dist(X1))
    d2 <- as.vector(dist(X2))
    
    # Correlate structural geometry
    similarity <- cor(d1, d2, use = "complete.obs")
    
    return(similarity)
  }
  
  
  cluster_mat    <- scale(metric_matrix[, cluster_set, drop = FALSE])
  pcdom_mat      <- scale(metric_matrix[, pc_dominant_set, drop = FALSE])
  pcscore_mat    <- scale(pc_scores_df)
  
  results$similarity_matrix <- matrix(
    c(
      compute_rep_similarity(cluster_mat, pcscore_mat),
      compute_rep_similarity(pcdom_mat, pcscore_mat)
    ),
    nrow = 2,
    byrow = TRUE
  )
  
  rownames(results$similarity_matrix) <- c(
    "Cluster_vs_PCscore",
    "PCdom_vs_PCscore"
  )
  
  ##########################################################
  # 3. Summary Table (For Manuscript)
  ##########################################################
  
  results$summary_table <- data.frame(
    Representation = c("Cluster",
                       "PC_Dominant",
                       "PC_Scores"),
    Num_Predictors = c(length(cluster_set),
                       length(pc_dominant_set),
                       ncol(pc_scores_df)),
    Variance_Retained = c(results$cluster$variance_retained,
                          results$pc_dominant$variance_retained,
                          results$pc_scores$variance_retained),
    Effective_Dim = c(results$cluster$effective_dimensionality,
                      results$pc_dominant$effective_dimensionality,
                      results$pc_scores$effective_dimensionality),
    Mean_Abs_Correlation = c(results$cluster$mean_abs_correlation,
                             results$pc_dominant$mean_abs_correlation,
                             results$pc_scores$mean_abs_correlation)
  )
  
  return(results)
}

# calculate 

dimensional_results <- characterise_reduction(
  metric_matrix = metrics,
  cluster_set = cluster_medoids,
  pc_dominant_set = pc_dominant_metrics,
  pc_scores_df = pc_scores_df
)

############################################################
# Dimensional Reduction Patchwork Figure
############################################################

plot_dimensional_summary <- function(results){
  
  summary_df <- results$summary_table
  
  ########################################
  # Panel A — Variance Retention
  ########################################
  
  p1 <- ggplot(summary_df,
               aes(x = Representation,
                   y = Variance_Retained,
                   fill = Representation)) +
    geom_col(colour = "white") +
    scale_fill_manual(values = secondary_palette) +
    labs(title = "A. Variance Retention",
         y = "Proportion Variance Retained",
         x = "") +
    figure_theme() +
    theme(legend.position = 'none',
          axis.title = element_text(size = rel(0.8), face = "bold"))
  
  ########################################
  # Panel B — Redundancy
  ########################################
  
  p2 <- ggplot(summary_df,
               aes(x = Representation,
                   y = Mean_Abs_Correlation,
                   fill = Representation)) +
    geom_col(colour = "white") +
    scale_fill_manual(values = secondary_palette) +
    labs(title = "B. Internal Redundancy",
         y = "Mean |r|",
         x = "") +
    figure_theme() +
    theme(legend.position = 'none',
          axis.title = element_text(size = rel(0.8), face = "bold"))
  
  ########################################
  # Panel C — Effective Dimensionality
  ########################################
  
  p3 <- ggplot(summary_df,
               aes(x = Representation,
                   y = Effective_Dim,
                   fill = Representation)) +
    geom_col(colour = "white") +
    scale_fill_manual(values = secondary_palette) +
    labs(title = "C. Effective Dimensionality",
         y = "Dimensionality (80% variance)",
         x = "") +
    figure_theme() +
    theme(legend.position = 'none',
          axis.title = element_text(size = rel(0.8), face = "bold"))
  
  ########################################
  # Combine With Patchwork
  ########################################
  
  combined <- (p1 / p2 / p3)
  
  return(combined)
}

dimensional_plot <- plot_dimensional_summary(dimensional_results)

dimensional_plot

ggsave("../figures/fig_dimensional_reduction.png",
       dimensional_plot,
       width = 3000,
       height = 3500,
       units = "px")
