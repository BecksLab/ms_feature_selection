############################################################
# Module-Based Prediction of Stability
# Cleaned & Aligned with Structural Modules
############################################################

library(tidyverse)
library(randomForest)
library(pdp)
library(gridExtra)
library(here)
library(genzplyr)

set.seed(66)
setwd(here::here())

############################################################
# 1. Load Data
############################################################

topology <- read.csv("data/cleaned/all_networks.csv")

cluster_medoids <- readRDS("data/outputs/cluster_medoids.rds")
pc_dominant_metrics <- readRDS("data/outputs/pc_dominant_metrics.rds")
pc_scores_df <- readRDS("data/outputs/pc_scores_df.rds")

# get cluster/module masterlist
clust_metada <-
  left_join(read_csv("../tables/metric_clusters_auto.csv"),
            read_csv("../tables/module_summary_7clusters.csv") %>%
              vibe_check(c(Cluster, label)))

############################################################
# Create Cross-Validation Folds
############################################################

create_cv_folds <- function(n, k = 5, seed = 66){
  
  set.seed(seed)
  folds <- sample(rep(1:k, length.out = n))
  return(folds)
}

############################################################
# Cross-Validated RF with Permutation Importance
############################################################

run_rf_cv <- function(response,
                      predictors,
                      data,
                      folds,
                      nperm = 100,
                      ntree = 1000){
  
  df <- data[, c(response, predictors)]
  df <- na.omit(df)
  
  k <- length(unique(folds))
  
  # Storage
  predictions_all <- data.frame()
  fold_perf <- data.frame()
  importance_obs <- list()
  importance_null <- list()
  
  for (f in 1:k){
    
    train_data <- df[folds != f, ]
    test_data  <- df[folds == f, ]
    
    # Fit model
    rf_model <- randomForest(
      as.formula(paste(response, "~ .")),
      data = train_data,
      ntree = ntree,
      importance = TRUE
    )
    
    # Predict
    preds <- predict(rf_model, newdata = test_data)
    
    # Store predictions
    predictions_all <- rbind(
      predictions_all,
      data.frame(
        Fold = f,
        Observed = test_data[[response]],
        Predicted = preds
      )
    )
    
    # Performance metrics
    r2 <- cor(test_data[[response]], preds)^2
    rmse <- sqrt(mean((test_data[[response]] - preds)^2))
    
    fold_perf <- rbind(
      fold_perf,
      data.frame(Fold = f, R2 = r2, RMSE = rmse)
    )
    
    # Observed importance
    obs_imp <- importance(rf_model)[,"%IncMSE"]
    importance_obs[[f]] <- obs_imp
    
    # Permutation null importance
    null_mat <- matrix(NA, nrow = nperm, ncol = length(predictors))
    colnames(null_mat) <- predictors
    
    for (p in 1:nperm){
      
      test_perm <- test_data
      test_perm[[response]] <- sample(test_perm[[response]])
      
      rf_perm <- randomForest(
        as.formula(paste(response, "~ .")),
        data = rbind(train_data, test_perm),
        ntree = ntree,
        importance = TRUE
      )
      
      null_mat[p, ] <- importance(rf_perm)[,"%IncMSE"]
    }
    
    importance_null[[f]] <- null_mat
  }
  
  # Aggregate performance
  mean_R2 <- mean(fold_perf$R2)
  mean_RMSE <- mean(fold_perf$RMSE)
  
  # Aggregate importance
  obs_matrix <- do.call(rbind, importance_obs)
  null_matrix <- do.call(rbind, importance_null)
  
  mean_obs_imp <- colMeans(obs_matrix)
  mean_null_imp <- colMeans(null_matrix)
  
  z_scores <- (mean_obs_imp - mean_null_imp) /
    apply(null_matrix, 2, sd)
  
  p_values <- colMeans(abs(null_matrix) >=
                         matrix(abs(mean_obs_imp),
                                nrow = nrow(null_matrix),
                                ncol = length(mean_obs_imp),
                                byrow = TRUE))
  
  return(list(
    predictions = predictions_all,
    fold_performance = fold_perf,
    mean_performance = data.frame(R2 = mean_R2,
                                  RMSE = mean_RMSE),
    importance_observed = obs_matrix,
    importance_null = null_matrix,
    mean_importance = mean_obs_imp,
    z_scores = z_scores,
    p_values = p_values
  ))
}

############################################################
# Run RF for Multiple Stability Components
############################################################

run_rf_suite <- function(stability_vars,
                         predictor_set,
                         data,
                         folds,
                         nperm = 100){
  
  results <- list()
  
  for (stab in stability_vars){
    
    results[[stab]] <- run_rf_cv(
      response = stab,
      predictors = predictor_set,
      data = data,
      folds = folds,
      nperm = nperm
    )
  }
  
  return(results)
}

############################################################
# Execute for Each Structural Representation
############################################################

# PC scores must be appended to the topology dataframe
topology_pc <- cbind(topology, pc_scores_df)

# create shared folds
folds <- create_cv_folds(nrow(topology), k = 5)

# run models

# Cluster medoids
rf_cluster <- run_rf_suite(
  stability_vars,
  cluster_medoids,
  topology,
  folds,
  nperm = 100
)

# PC dominant metrics
rf_pc_dominant_metrics <- run_rf_suite(
  stability_vars,
  pc_dominant_metrics,
  topology,
  folds,
  nperm = 100
)

# PC scores
rf_pc <- run_rf_suite(
  stability_vars,
  colnames(pc_scores_df),
  topology_pc,
  folds,
  nperm = 100
)

# Save outputs
saveRDS(rf_cluster, "data/outputs/rf_cluster_results.rds")
saveRDS(rf_pc, "data/outputs/rf_pc_results.rds")
saveRDS(rf_pc_dominant_metrics, "data/outputs/rf_pc_dominant_metrics_results.rds")

############################################################
# maybe move this stuff own script
############################################################

############################################################
# Extract Performance Across Representations
############################################################

extract_performance <- function(rf_object, label){
  
  stability_names <- names(rf_object)
  
  perf_list <- lapply(stability_names, function(stab){
    
    fold_df <- rf_object[[stab]]$fold_performance
    
    data.frame(
      Stability = stab,
      Representation = label,
      Fold = fold_df$Fold,
      R2 = fold_df$R2,
      RMSE = fold_df$RMSE
    )
  })
  
  do.call(rbind, perf_list)
}

perf_cluster <- extract_performance(rf_cluster, "Cluster")
perf_pc_dom  <- extract_performance(rf_pc_dominant_metrics, "PC_Dominant")
perf_pc      <- extract_performance(rf_pc, "PC_Scores")

performance_all <- rbind(perf_cluster,
                         perf_pc_dom,
                         perf_pc)

performance_summary <- aggregate(
  cbind(R2, RMSE) ~ Stability + Representation,
  data = performance_all,
  FUN = function(x) c(mean = mean(x), sd = sd(x))
)

compare_models <- function(stability_name){
  
  df <- subset(performance_all,
               Stability == stability_name)
  
  cluster_r2 <- df$R2[df$Representation == "Cluster"]
  pc_r2      <- df$R2[df$Representation == "PC_Scores"]
  
  t.test(cluster_r2, pc_r2, paired = TRUE)
}

compare_models("robustness")
compare_models("ρ")
compare_models("complexity")

############################################################
# Extract Importance Summary
############################################################

extract_importance <- function(rf_object, label){
  
  stability_names <- names(rf_object)
  
  imp_list <- lapply(stability_names, function(stab){
    
    mean_imp <- rf_object[[stab]]$mean_importance
    z_scores <- rf_object[[stab]]$z_scores
    p_vals   <- rf_object[[stab]]$p_values
    
    data.frame(
      Stability = stab,
      Representation = label,
      Predictor = names(mean_imp),
      MeanImportance = mean_imp,
      Z = z_scores,
      P = p_vals
    )
  })
  
  do.call(rbind, imp_list)
}


imp_cluster <- extract_importance(rf_cluster, "Cluster")
imp_pc_dom  <- extract_importance(rf_pc_dominant_metrics, "PC_Dominant")
imp_pc      <- extract_importance(rf_pc, "PC_Scores")

importance_all <- rbind(imp_cluster,
                        imp_pc_dom,
                        imp_pc)

ggplot(performance_all,
       aes(x = Representation,
           y = R2,
           colour = Representation,
           fill = Representation)) +
  geom_jitter(width = 0.1,
              alpha = 0.6,
              size = 3) +
  stat_summary(fun = mean,
               geom = "point",
               size = 5) +
  stat_summary(fun.data = mean_se,
               aes(colour = Representation),
               geom = "errorbar",
               width = 0.2) +
  scale_colour_manual(values = secondary_palette) +
  scale_fill_manual(values = secondary_palette) +
  facet_wrap(~ Stability) +
  labs(y = expression(R^2),
       x = "",
       title = "Predictive Performance Across Structural Representations") +
  figure_theme() +
  theme(legend.position = 'none')

ggsave("../figures/fig_pred_preformance.png",
       width = 5500,
       height = 2500,
       units = "px")

plot_imp_df <- 
  importance_all %>%
  left_join(clust_metada, 
            by = join_by(Predictor == Metric)) %>%
  # get colours
  left_join(pal_df) %>%
  # add manual colour for PCA axes (keep same since they have no visual language)
  glow_up(colour = if_else(is.na(colour), "#50723C", colour),
          label = if_else(is.na(label), "PCA Axis", label))

ggplot(plot_imp_df,
         aes(y = Predictor,
             x = MeanImportance,
             colour = label)) +
  geom_segment(aes(xend = MeanImportance,
                   x = -Inf,
                   y = Predictor),
               linewidth = 1) +
  geom_point(size = 8) +
  scale_color_manual(values = setNames(plot_imp_df$colour, as.character(plot_imp_df$label)),
                     name = "Module") +
  facet_grid(rows = vars(Representation),
             cols = vars(Stability),
             scales = "free_y") +
  labs(x = "Mean Importance") +
  figure_theme() +
  theme(legend.position = 'right')

ggsave("../figures/mean_importance.png",
       width = 5500,
       height = 3000,
       units = "px")

############################################################
# Extract Variable Importance (Simple & Clean)
############################################################

extract_variable_importance <- function(rf_object, label){
  
  stab_names <- names(rf_object)
  
  imp_list <- lapply(stab_names, function(stab){
    
    imp_matrix <- rf_object[[stab]]$importance_observed
    z_vals     <- rf_object[[stab]]$z_scores
    p_vals     <- rf_object[[stab]]$p_values
    
    # Ensure matrix
    imp_matrix <- as.matrix(imp_matrix)
    
    # Aggregate across folds
    mean_imp <- colMeans(imp_matrix, na.rm = TRUE)
    
    data.frame(
      Stability = stab,
      Representation = label,
      Predictor = names(mean_imp),
      MeanImportance = as.numeric(mean_imp),
      Z = as.numeric(z_vals),
      P = as.numeric(p_vals)
    )
  })
  
  do.call(rbind, imp_list)
}

imp_cluster <- extract_variable_importance(rf_cluster, "Cluster")

imp_pc_dom <- extract_variable_importance(
  rf_pc_dominant_metrics,
  "PC_Dominant"
)

imp_pc <- extract_variable_importance(
  rf_pc,
  "PC_Scores"
)

importance_all <- rbind(
  imp_cluster,
  imp_pc_dom,
  imp_pc
)



############################################################
# End
############################################################