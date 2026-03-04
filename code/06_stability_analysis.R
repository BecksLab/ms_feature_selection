############################################################
# Module-Based Prediction of Stability
############################################################

library(glmnet)
library(rsample)
library(dplyr)
library(purrr)
library(tidyr)

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

stability_vars <- c("robustness", "ρ", "complexity")

############################################################
# Elastic Net Modeling Suite
############################################################

run_elastic_suite <- function(stability_vars,
                              predictor_names,
                              data_matrix,
                              nfolds = 5,
                              nrepeats = 20){
  
  results_list <- list()
  
  for(stab in stability_vars){
    
    y <- scale(data_matrix[[stab]])[,1]
    X <- scale(data_matrix[, predictor_names])
    
    df_model <- data.frame(y = y, X)
    
    # Create repeated CV splits
    splits <- vfold_cv(df_model, v = nfolds, repeats = nrepeats)
    
    fold_results <- map_dfr(1:nrow(splits), function(i){
      
      split <- splits$splits[[i]]
      train <- analysis(split)
      test  <- assessment(split)
      
      x_train <- as.matrix(vibe_check(train, -y))
      y_train <- train$y
      
      x_test  <- as.matrix(vibe_check(test, -y))
      y_test  <- test$y
      
      # Tune alpha over small grid
      alpha_grid <- seq(0, 1, by = 0.25)
      
      cv_models <- map(alpha_grid, function(a){
        cv.glmnet(x_train, y_train,
                  alpha = a,
                  standardize = FALSE)
      })
      
      cv_errors <- map_dbl(cv_models, ~min(.x$cvm))
      best_idx  <- which.min(cv_errors)
      
      best_model <- cv_models[[best_idx]]
      best_alpha <- alpha_grid[best_idx]
      best_lambda <- best_model$lambda.min
      
      preds <- predict(best_model,
                       newx = x_test,
                       s = "lambda.min")
      
      r2 <- cor(preds, y_test)^2
      
      tibble(
        Stability = stab,
        Fold = i,
        R2 = r2,
        Alpha = best_alpha,
        Lambda = best_lambda
      )
    })
    
    # Fit final model on full dataset with best alpha (mean over folds)
    best_alpha_global <- mean(fold_results$Alpha)
    
    cv_full <- cv.glmnet(as.matrix(X),
                         y,
                         alpha = best_alpha_global,
                         standardize = FALSE)
    
    final_model <- glmnet(as.matrix(X),
                          y,
                          alpha = best_alpha_global,
                          lambda = cv_full$lambda.min,
                          standardize = FALSE)
    
    coefs <- as.matrix(coef(final_model))[-1, , drop=FALSE]
    
    coef_df <- tibble(
      Stability = stab,
      Predictor = rownames(coefs),
      Coefficient = as.numeric(coefs)
    )
    
    results_list[[stab]] <- list(
      performance = fold_results,
      coefficients = coef_df
    )
  }
  
  return(results_list)
}

# run for all representations

elastic_cluster <- run_elastic_suite(
  stability_vars,
  cluster_medoids,
  topology
)

elastic_pc_dom <- run_elastic_suite(
  stability_vars,
  pc_dominant_metrics,
  topology
)

elastic_pc_scores <- run_elastic_suite(
  stability_vars,
  colnames(pc_scores_df),
  topology_pc
)

perf_cluster <- map_dfr(elastic_cluster, "performance") %>%
  mutate(Representation = "Cluster")

perf_pc_dom <- map_dfr(elastic_pc_dom, "performance") %>%
  mutate(Representation = "PC_Dominant")

perf_pc_scores <- map_dfr(elastic_pc_scores, "performance") %>%
  mutate(Representation = "PC_Scores")

performance_all <- bind_rows(
  perf_cluster,
  perf_pc_dom,
  perf_pc_scores
)

coef_cluster <- map_dfr(elastic_cluster, "coefficients") %>%
  mutate(Representation = "Cluster")

coef_pc_dom <- map_dfr(elastic_pc_dom, "coefficients") %>%
  mutate(Representation = "PC_Dominant")

coef_pc_scores <- map_dfr(elastic_pc_scores, "coefficients") %>%
  mutate(Representation = "PC_Scores")

coefficients_all <- bind_rows(
  coef_cluster,
  coef_pc_dom,
  coef_pc_scores
)

# which stability component is most explainable

plot_r2 <- performance_all %>%
  squad_up(Stability, Representation) %>%
  no_cap(mean = mean(R2),
         sd = sd(R2)) %>%
  glow_up(ymax = mean + sd,
          ymin = mean - sd) %>%
  ggplot() +
  geom_point(aes(x = Representation,
                 y = mean,
                 colour = Representation),
             size = 6) +
  geom_errorbar(aes(x = Representation,
                    ymin = ymin,
                    ymax = ymax,
                    colour = Representation), 
                width = 0.2) +
  scale_colour_manual(values = secondary_palette) +
  scale_fill_manual(values = secondary_palette) +
  facet_wrap(~ Stability) +
  labs(y = expression(R^2),
       x = "",
       title = "Variance Explained Across Structural Representations") +
  ylim(0, 1.0) +
  figure_theme() +
  theme(legend.position = 'none')

# Is stability controlled by many small structural 
# effects or a few dominant ones
# 0 - many correlated predictors share weight
# 1 - sparse, few predictors dominate

plot_alpha <- performance_all %>%
  squad_up(Stability, Representation) %>%
  no_cap(mean = mean(Alpha),
         sd = sd(Alpha)) %>%
  ggplot() +
  geom_col(aes(x = Representation,
               y = mean,
               fill = Representation)) +
  scale_fill_manual(values = secondary_palette) +
  facet_wrap(~ Stability) +
  labs(y = "Alpha",
       x = "",
       title = "Stability Control Across Structural Representations") +
  ylim(0, 1.0) +
  figure_theme() +
  theme(legend.position = 'none')

plot_r2 /
  plot_alpha + 
  plot_annotation(tag_levels = 'A')

ggsave("../figures/struct_stability_summ.png",
       width = 5500,
       height = 7000,
       units = "px")

# coefficients

plot_coeff_df <- coefficients_all %>%
  left_join(clust_metada, 
            by = join_by(Predictor == Metric)) %>%
  # get colours
  left_join(pal_df) %>%
  # add manual colour for PCA axes (keep same since they have no visual language)
  glow_up(colour = if_else(is.na(colour), "#50723C", colour),
          label = if_else(is.na(label), "PCA Axis", label))

ggplot(plot_coeff_df,
       aes(x = Coefficient,
           y = Predictor,
           colour = label)) +
  geom_segment(aes(xend = 0),
               linewidth = 1) +
  geom_vline(xintercept = 0,
             colour = "#001628",
             linewidth = 0.6) +
  geom_point(size = 6) +
  scale_color_manual(values = setNames(plot_imp_df$colour, as.character(plot_imp_df$label)),
                     breaks = plot_imp_df$label[plot_imp_df$label != "PCA Axis"],
                     name = "Module") +
  facet_grid(rows = vars(Representation),
             cols = vars(Stability),
             scales = "free_y") +
  figure_theme() +
  theme(legend.position = 'right')

ggsave("../figures/struct_stability_coeff.png",
       width = 6000,
       height = 5000,
       units = "px")

############################################################
# End
############################################################