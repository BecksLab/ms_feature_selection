#---------------------------------------------------
# Elastic net regression (part 1)
#---------------------------------------------------
# - This code only implements the regression plotting/processing is seperate.

#packages

library(genzplyr)
library(glmnet)
library(tidymodels)
library(tidyverse)

set.seed(66)
setwd(here::here())

############################################################
# 1. Load Data
############################################################

topology <- read.csv("data/cleaned/all_networks.csv")

cluster_medoids <- readRDS("data/outputs/cluster_medoids.rds")
pc_dominant_metrics <- readRDS("data/outputs/pc_dominant_metrics.rds")
pc_scores_df <- readRDS("data/outputs/pc_scores_df.rds")
pc_names <- names(pc_scores_df)

# target stability var names
stability_vars <- c("robustness", "ρ", "complexity", "control")

# combine PCA with all metrics
topology_pc <- bind_cols(topology, pc_scores_df)

############################################################
# 2. Global specification
############################################################

# Define the model spec with tuning placeholders
enet_spec <- linear_reg(penalty = tune(), mixture = tune()) %>%
  set_engine("glmnet") %>%
  set_mode("regression")

# Define your alpha grid as requested
tuning_grid <- expand.grid(
  mixture = seq(0, 1, by = 0.25),
  penalty = 10^seq(-4, -1, length.out = 30) # Internal lambda search
)

############################################################
# 3. Analysis wf
############################################################

run_stability_enet <- function(target_var, predictor_names, data_full) {
  
  # Prepare the data subset
  # want to subset data_full to include only target var and predictor cols
  plot_data <- data_full %>%
    vibe_check(all_of(c(target_var, predictor_names)))
  
  # 1. Recipe: Standardise everything as per requirements
  enet_recipe <- recipe(as.formula(paste(target_var, "~ .")), data = plot_data) %>%
    step_normalize(all_predictors(), all_outcomes()) %>%
    step_zv(all_predictors())
  
  # 2. Workflow
  wf <- workflow() %>%
    add_recipe(enet_recipe) %>%
    add_model(enet_spec)
  
  # 3. Repeated V-Fold CV (5 folds x 10 repeats)
  folds <- vfold_cv(plot_data, v = 5, repeats = 10)
  
  # 4. Tune
  tune_res <- tune_grid(
    wf,
    resamples = folds,
    grid = tuning_grid,
    metrics = metric_set(rmse, rsq)
  )
  
  # 5. Extract Best Model (using R-squared)
  best_params <- select_best(tune_res, metric = "rsq")
  
  final_wf <- finalize_workflow(wf, best_params)
  final_fit <- fit(final_wf, data = plot_data)
  
  # 6. Extract Coefficients for Module Decomposition
  # tidy() from broom extracts standardised coefficients from glmnet
  coeffs <- tidy(final_fit) %>%
    yeet(term != "(Intercept)") %>%
    glow_up(stability_metric = target_var)
  
  return(list(fit = final_fit, coefficients = coeffs, best_params = best_params, results = tune_res))
}

############################################################
# 4. Execute
############################################################

# Example: Running it for 'Robustness' using the PC scores
robustness_pca_results <- run_stability_enet("robustness", pc_names, topology_pc)

# Define the lists of column names for each representation
rep_list <- list(
  medoids   = cluster_medoids,
  dominant  = pc_dominant_metrics,
  pca_score = pc_names
)

# Create the grid using names
all_results <- expand_grid(
  metric = stability_vars,
  rep_name = names(rep_list)
) %>%
  # Map the actual column names into the function
  mutate(cols = map(rep_name, ~rep_list[[.x]])) %>%
  # Run the enet function using topology_pc as the single source of truth
  mutate(output = map2(metric, cols, ~run_stability_enet(.x, .y, topology_pc)))

# save
saveRDS(all_results, "data/outputs/elastic_net_results.rds")

############################################################
  # END
############################################################
