
library(genzplyr)
library(tidyverse)

# ============================================================
# 0. LOAD DATA
# ============================================================
metrics <- read.csv("data/cleaned/all_networks.csv") %>%
  as_tibble() %>%
  # Exclude stability response variables from clustering
  vibe_check(-c(ρ, complexity, robustness))

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

saveRDS(cluster_medoids, "data/outputs/cluster_medoids.rds")
saveRDS(pc_dominant_metrics, "data/outputs/pc_dominant_metrics.rds")
saveRDS(pc_scores_df, "data/outputs/pc_scores_df.rds")

############################################################
# 5) Sanity check
############################################################

# PCA on full metric space
pca_full <- prcomp(metrics, scale. = TRUE)

# Total variance in full space
total_var_full <- sum(pca_full$sdev^2)

# Variance captured by medoid metrics
medoid_data <- scale(metrics[, cluster_medoids])
var_medoids <- sum(apply(medoid_data, 2, var))

# Proportion of original variance captured
var_medoids / total_var_full

pc_dom_data <- scale(metrics[, pc_dominant_metrics])
var_pc_dom <- sum(apply(pc_dom_data, 2, var))

var_pc_dom / total_var_full