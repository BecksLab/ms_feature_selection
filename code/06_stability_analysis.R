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

############################################################
# 2. Define Module Representative Metrics
############################################################

# -----> Replace with your empirically selected metrics
module_reps <- c("S1", "richness", "predpreyRatio", "distance", "herbivory", "MaxSim", "centrality")

# Keep only variables that exist in dataset
module_reps <- module_reps[module_reps %in% colnames(topology)]

############################################################
# 3. Define Stability Components
############################################################

stability_vars <- c("robustness", "ρ", "complexity")

stability_vars <- stability_vars[stability_vars %in% colnames(topology)]

############################################################
# 4. Function to Fit Random Forest Model
############################################################

fit_rf_model <- function(response, predictors, data){
  
  df <- data %>%
    dplyr::select(all_of(c(response, predictors))) %>%
    drop_na()
  
  formula <- as.formula(
    paste(response, "~", paste(predictors, collapse = "+"))
  )
  
  rf <- randomForest(
    formula,
    data = df,
    importance = TRUE,
    ntree = 1000
  )
  
  return(rf)
}

############################################################
# 5. Fit Models for Each Stability Component
############################################################

rf_models <- list()

for (stab in stability_vars) {
  
  rf_models[[stab]] <- fit_rf_model(
    response = stab,
    predictors = module_reps,
    data = topology
  )
}

############################################################
# 6. Extract Variable Importance
############################################################

importance_table <- map2_df(
  rf_models,
  names(rf_models),
  ~ {
    
    imp <- as.data.frame(importance(.x))
    imp <- tibble(
      Module = rownames(imp),
      Importance = imp$`%IncMSE`,
      Stability = .y
    )
    
    return(imp)
  }
)

print(importance_table)

############################################################
# 7. Extract Model Performance (R²)
############################################################

performance <- tibble(
  Stability = names(rf_models),
  R2 = map_dbl(rf_models, ~ max(.x$rsq))
)

print(performance)

############################################################
# 8. Plot Variable Importance Per Stability Component
############################################################

plots <- map(rf_models, ~ {
  
  varImpPlot(.x,
             type = 1,
             main = paste("Drivers of", deparse(substitute(.x))))
})

# Display side-by-side
par(mfrow = c(1, length(plots)))
invisible(plots)
par(mfrow = c(1,1))


############################################################
# 9. Partial Dependence Plots
############################################################

pdp_plots <- list()

for (stab in stability_vars) {
  
  rf_model <- rf_models[[stab]]
  
  # Extract predictor names from the model
  model_vars <- attr(rf_model$terms, "term.labels")
  
  # Reconstruct training data using SAME variables
  train_data <- topology %>%
    dplyr::select(dplyr::all_of(c(stab, model_vars))) %>%
    tidyr::drop_na()
  
  for (pred in intersect(module_reps, model_vars)) {
    
    # Double-check variable actually exists
    if (!(pred %in% colnames(train_data))) next
    
    p <- partial(
      object = rf_model,
      pred.var = pred,
      train = train_data,
      plot = TRUE,
      rug = TRUE,
      grid.resolution = 20,
      main = paste(stab, "~", pred)
    )
    
    pdp_plots[[paste(stab, pred, sep = "_")]] <- p
  }
}


############################################################
# 10. Permutation Test for Importance
############################################################

# Permute stability variable and refit to build null importance

perm_test <- function(response,
                      predictors,
                      data,
                      nperm = 500,
                      alternative = "two.sided") {
  
  # Remove missing values once
  df <- data[, unique(c(response, predictors))]
  df <- na.omit(df)
  
  # ---- Observed Model ----
  obs_rf <- fit_rf_model(response, predictors, df)
  obs_imp <- importance(obs_rf)[, "%IncMSE"]
  
  # Store null importance
  null_mat <- matrix(NA,
                     nrow = nperm,
                     ncol = length(predictors))
  
  colnames(null_mat) <- predictors
  
  # ---- Permutation Loop ----
  for (i in seq_len(nperm)) {
    
    df_perm <- df
    
    # Permute response (break structure–stability link)
    df_perm[[response]] <- sample(df_perm[[response]])
    
    rf_perm <- fit_rf_model(response, predictors, df_perm)
    
    null_mat[i, ] <- importance(rf_perm)[, "%IncMSE"]
  }
  
  # ---- Statistics ----
  
  null_mean <- colMeans(null_mat)
  null_sd   <- apply(null_mat, 2, sd)
  
  z_scores <- (obs_imp - null_mean) / null_sd
  
  # Empirical p-values
  if (alternative == "greater") {
    
    p_vals <- colMeans(null_mat >= matrix(obs_imp,
                                          nrow = nperm,
                                          ncol = length(obs_imp),
                                          byrow = TRUE))
    
  } else if (alternative == "less") {
    
    p_vals <- colMeans(null_mat <= matrix(obs_imp,
                                          nrow = nperm,
                                          ncol = length(obs_imp),
                                          byrow = TRUE))
    
  } else {
    
    # Two-sided
    p_vals <- colMeans(abs(null_mat) >=
                         matrix(abs(obs_imp),
                                nrow = nperm,
                                ncol = length(obs_imp),
                                byrow = TRUE))
  }
  
  # ---- Return Full Object ----
  
  return(list(
    observed_importance = obs_imp,
    null_distribution   = null_mat,
    z_scores            = z_scores,
    p_values            = p_vals
  ))
}

perm_results <- perm_test(
  response   = "robustness",
  predictors = module_reps,
  data       = topology,
  nperm      = 500
)

perm_results$p_values
perm_results$z_scores

boxplot(perm_results$null_distribution)
points(perm_results$observed_importance,
       col = "red",
       pch = 19)

############################################################
# Cross-Validated Permutation Importance
############################################################

cv_perm_importance <- function(response,
                               predictors,
                               data,
                               nfolds = 5,
                               nperm  = 100,
                               seed   = 66) {
  
  set.seed(seed)
  
  df <- data[, unique(c(response, predictors))]
  df <- na.omit(df)
  
  folds <- sample(rep(1:nfolds,
                      length.out = nrow(df)))
  
  # Storage
  fold_results <- list()
  
  for (f in seq_len(nfolds)) {
    
    train_data <- df[folds != f, ]
    test_data  <- df[folds == f, ]
    
    # ----- Train Model on Fold -----
    rf_model <- fit_rf_model(response, predictors, train_data)
    
    obs_imp <- importance(rf_model)[, "%IncMSE"]
    
    # ---- Permutation Null Within Fold ----
    null_mat <- matrix(NA,
                       nrow = nperm,
                       ncol = length(predictors))
    
    colnames(null_mat) <- predictors
    
    for (p in seq_len(nperm)) {
      
      test_perm <- test_data
      test_perm[[response]] <- sample(test_perm[[response]])
      
      rf_perm <- fit_rf_model(response,
                              predictors,
                              rbind(train_data, test_perm))
      
      null_mat[p, ] <- importance(rf_perm)[, "%IncMSE"]
    }
    
    fold_results[[f]] <- list(
      observed = obs_imp,
      null     = null_mat
    )
  }
  
  obs_matrix  <- do.call(rbind,
                         lapply(fold_results, `[[`, "observed"))
  
  null_matrix <- do.call(rbind,
                         lapply(fold_results, `[[`, "null"))
  
  mean_obs <- colMeans(obs_matrix)
  mean_null <- colMeans(null_matrix)
  
  z_scores <- (mean_obs - mean_null) /
    apply(null_matrix, 2, sd)
  
  p_vals <- colMeans(abs(null_matrix) >=
                       matrix(abs(mean_obs),
                              nrow = nrow(null_matrix),
                              ncol = length(mean_obs),
                              byrow = TRUE))
  
  return(list(
    mean_observed = mean_obs,
    mean_null     = mean_null,
    z_scores      = z_scores,
    p_values      = p_vals,
    null_dist     = null_matrix
  ))
}

cv_results <- cv_perm_importance(
  response   = "robustness",
  predictors = module_reps,
  data       = topology,
  nfolds     = 5,
  nperm      = 100
)

cv_results$p_values
cv_results$z_scores

############################################################
# Plot 1: Observed vs Null
############################################################

plot_observed_vs_null <- function(cv_results){
  
  null_mat <- cv_results$null_dist
  
  modules <- colnames(null_mat)
  
  # ---- Build Proper Null DataFrame ----
  null_df <- expand.grid(
    Iteration = seq_len(nrow(null_mat)),
    Module = modules
  )
  
  null_df$Importance <- as.vector(null_mat)
  
  # ---- Observed ----
  obs_df <- tibble(
    Module = names(cv_results$mean_observed),
    Observed = cv_results$mean_observed
  )
  
  # ---- Plot ----
  ggplot(null_df, aes(x = Module, y = Importance)) +
    geom_boxplot(
      data = null_df,
      aes(group = Module),
      fill = "grey80",
      alpha = 0.6
    ) +
    geom_point(
      data = obs_df,
      aes(x = Module, y = Observed),
      color = "red",
      size = 3
    ) +
    theme_minimal() +
    labs(
      title = "Observed vs Null Module Importance",
      y = "% IncMSE",
      x = "Module"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_observed_vs_null(cv_results)

############################################################
# Plot 2: Ranked Module Importance
############################################################

plot_ranked_importance <- function(cv_results){
  
  df <- tibble(
    Module = names(cv_results$mean_observed),
    Importance = cv_results$mean_observed,
    Z = cv_results$z_scores
  ) %>%
    arrange(desc(Importance)) %>%
    mutate(Module = factor(Module,
                           levels = Module))
  
  df$Significant <- cv_results$p_values < 0.05
  
  ggplot(df, aes(x = Module, y = Importance)) +
    geom_col(fill = "steelblue") +
    geom_errorbar(aes(ymin = Importance - abs(Z),
                      ymax = Importance + abs(Z)),
                  width = 0.2) +
    geom_text(aes(label = ifelse(Significant, "*", "")),
              hjust = -0.5) +
    coord_flip() +
    theme_minimal() +
    labs(title = "Ranked Module Importance",
         y = "Cross-Validated %IncMSE",
         x = "Module")
}

plot_ranked_importance(cv_results)

############################################################
# Multi panel figure
############################################################

cv_results_list <-
  list(robustness = cv_perm_importance(response   = "robustness",
                                       predictors = module_reps,
                                       data       = topology,
                                       nfolds     = 5,
                                       nperm      = 100),
       spectral_radius = cv_perm_importance(response   = "ρ",
                                            predictors = module_reps,
                                            data       = topology,
                                            nfolds     = 5,
                                            nperm      = 100),
       complexity = cv_perm_importance(response   = "complexity",
                                       predictors = module_reps,
                                       data       = topology,
                                       nfolds     = 5,
                                       nperm      = 100)
       )

############################################################
# Multi-Panel Summary Figure
############################################################

plot_multi_stability_summary <- function(cv_results_list){
  
  plot_list <- list()
  
  for (stab in names(cv_results_list)) {
    
    res <- cv_results_list[[stab]]
    null_mat <- res$null_dist
    
    modules <- colnames(null_mat)
    
    # ---- Build Null Data ----
    null_df <- expand.grid(
      Iteration = seq_len(nrow(null_mat)),
      Module = modules
    )
    
    null_df$Importance <- as.vector(null_mat)
    
    # ---- Observed ----
    obs_df <- tibble(
      Module = names(res$mean_observed),
      Observed = res$mean_observed
    )
    
    # ---- Plot for This Stability Component ----
    p <- ggplot(null_df, aes(x = Module, y = Importance)) +
      geom_boxplot(fill = "grey80",
                   alpha = 0.6,
                   outlier.shape = NA) +
      geom_point(
        data = obs_df,
        aes(x = Module, y = Observed),
        color = "red",
        size = 3
      ) +
      coord_flip() +
      theme_minimal() +
      labs(
        title = paste("Stability Component:", stab),
        y = "% IncMSE",
        x = NULL
      ) +
      theme(
        axis.text.y = element_text(size = 8),
        plot.title = element_text(size = 10, face = "bold")
      )
    
    plot_list[[stab]] <- p
  }
  
  # ---- Arrange Panels ----
  gridExtra::grid.arrange(
    grobs = plot_list,
    ncol = 1
  )
}

plot_multi_stability_summary(cv_results_list)

############################################################
# End
############################################################