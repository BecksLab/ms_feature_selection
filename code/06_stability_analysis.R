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
source("lib/plotting_themes.R")

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
             size = 8) +
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
               fill = Representation),
           colour = 'white') +
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
       height = 6000,
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
           fill = label,
           colour = label)) +
  geom_segment(aes(xend = 0),
               linewidth = 1) +
  geom_vline(xintercept = 0,
             colour = "#001628",
             linewidth = 0.6) +
  geom_point(size = 9,
             shape = 21,
             colour = "white") +
  scale_fill_manual(values = setNames(plot_coeff_df$colour, as.character(plot_coeff_df$label)),
                     breaks = plot_coeff_df$label[plot_coeff_df$label != "PCA Axis"],
                     name = "Module") +
  scale_colour_manual(values = setNames(plot_coeff_df$colour, as.character(plot_coeff_df$label)),
                    breaks = plot_coeff_df$label[plot_coeff_df$label != "PCA Axis"],
                    name = "Module") +
  facet_grid(rows = vars(Representation),
             cols = vars(Stability),
             scales = "free_y") +
  figure_theme() +
  theme(legend.position = 'right')

ggsave("../figures/struct_stability_coeff.png",
       width = 5000,
       height = 4000,
       units = "px")

############################################################
# Variance partitioning - metric
############################################################

variance_partition <- coefficients_all %>%
  mutate(VarContribution = Coefficient^2) %>%
  group_by(Stability, Representation) %>%
  mutate(PropVar = VarContribution / sum(VarContribution)) %>%
  left_join(clust_metada, 
            by = join_by(Predictor == Metric)) %>%
  # get colours
  left_join(pal_df) %>%
  # add manual colour for PCA axes (keep same since they have no visual language)
  glow_up(colour = if_else(is.na(colour), "#50723C", colour),
          label = if_else(is.na(label), "PCA Axis", label))


ggplot(variance_partition,
       aes(y = Predictor,
           x = PropVar,
           fill = label,
           colour = label))  +
  geom_segment(aes(xend = 0),
               linewidth = 1) +
  geom_point(size = 9,
             shape = 21,
             colour = "white") +
  facet_grid(rows = vars(Representation),
             cols = vars(Stability),
             scales = "free_y") +
  scale_colour_manual(values = setNames(variance_partition$colour, 
                                      as.character(variance_partition$label)),
                    breaks = variance_partition$label[variance_partition$label != "PCA Axis"],
                    name = "Module") +
  scale_fill_manual(values = setNames(variance_partition$colour, 
                                      as.character(variance_partition$label)),
                    breaks = variance_partition$label[variance_partition$label != "PCA Axis"],
                    name = "Module") +
  labs(x = "Proportion of explained variance",
       y = NULL,) +
  figure_theme()v

ggsave("../figures/struct_stability_variance.png",
       width = 6500,
       height = 4000,
       units = "px")

############################################################
# Variance partitioning - module
############################################################

module_variance <- coefficients_all %>%
  left_join(clust_metada,
            by = join_by(Predictor == Metric)) %>%
  mutate(Module = if_else(is.na(label), "PCA Axis", label),
         VarContribution = Coefficient^2) %>%
  group_by(Stability, Representation, Module) %>%
  summarise(ModuleVar = sum(VarContribution), .groups = "drop") %>%
  group_by(Stability, Representation) %>%
  mutate(PropVar = ModuleVar / sum(ModuleVar)) %>%
  yeet(PropVar != 0) %>%
  # get colours
  left_join(pal_df, 
            by = join_by(Module == label)) %>%
  # add manual colour for PCA axes (keep same since they have no visual language)
  glow_up(colour = if_else(is.na(colour), "#50723C", colour))

ggplot(module_variance,
       aes(x = Representation,
           y = PropVar,
           fill = Module)) +
  geom_col(colour = "white") +
  scale_fill_manual(values = setNames(module_variance$colour, 
                                      as.character(module_variance$Module)),
                    breaks = module_variance$Module[module_variance$Module != "PCA Axis"],
                    name = "Module") +
  facet_wrap(~ Stability) +
  labs(y = "Proportion of explained variance",
       x = "") +
  figure_theme() +
  theme(legend.position = 'right',
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("../figures/struct_stability_variance_module.png",
       width = 6500,
       height = 2500,
       units = "px")

############################################################
# Variance partitioning - scaled
############################################################

model_r2 <- performance_all %>%
  group_by(Stability, Representation) %>%
  summarise(R2 = mean(R2),
            .groups = "drop_last")

module_var <- coefficients_all %>%
  left_join(clust_metada,
            by = join_by(Predictor == Metric)) %>%
  mutate(Module = if_else(is.na(label), "PCA Axis", label),
         VarContribution = Coefficient^2) %>%
  group_by(Stability, Representation, Module) %>%
  summarise(ModuleVar = sum(VarContribution), 
            .groups = "drop") %>%
  group_by(Stability, Representation) %>%
  mutate(PropVar = ModuleVar / sum(ModuleVar))

module_var_scaled <- module_var %>%
  left_join(model_r2,
            by = c("Stability", "Representation")) %>%
  mutate(ScaledVar = PropVar * R2)

ggplot(module_var_scaled,
       aes(x = Representation,
           y = ScaledVar,
           fill = Module)) +
  geom_col(width = 0.7) +
  facet_wrap(~ Stability) +
  labs(
    y = expression("Variance explained (" * R^2 * ")"),
    x = "",
    title = "Structural contributions to ecological stability"
  ) +
  scale_fill_manual(values = pal_df$colour) +
  figure_theme()

############################################################
# Stability - structure PCA link
############################################################

# Join PCA scores with stability
pca_stability <- pc_scores_df %>%
  glow_up(network_id = row_number()) %>%
  left_join(topology[, c("robustness", "ρ", "complexity")]%>%
              glow_up(network_id = row_number()), 
            by = "network_id") 

# Example: colour by robustness, size by complexity
ggplot(pca_stability, aes(x = PC1, y = PC2)) +
  geom_point(aes(colour = robustness, size = complexity), alpha = 0.8) +
  scale_colour_viridis_c(option = "plasma") +
  scale_size_continuous(range = c(3,8)) +
  labs(
    x = "PC1",
    y = "PC2",
    colour = "Robustness",
    size = "Complexity",
    title = "Food-web PCA with stability metrics overlayed"
  )  +
  figure_theme()

stab_corr <- pca_stability %>%
  summarise(
    r_robust_PC1 = cor(robustness, PC1),
    r_robust_PC2 = cor(robustness, PC2),
    r_rho_PC1    = cor(ρ, PC1),
    r_rho_PC2    = cor(ρ, PC2),
    r_comp_PC1   = cor(complexity, PC1),
    r_comp_PC2   = cor(complexity, PC2)
  ) %>% 
  pivot_longer(everything(), names_to = "metric", values_to = "cor") %>%
  separate(metric, into = c("stab","PC"), sep = "_PC") %>%
  pivot_wider(names_from = PC, 
              values_from = cor,
              names_prefix = "PC") %>%
  glow_up(stab = case_when(stab == "r_robust" ~ "robustness",
                           stab == "r_rho" ~ "ρ",
                           .default = "complexity"))

ggplot(pca_stability, aes(x = PC1, y = PC2)) +
  geom_point(alpha = 0.7) +
  geom_segment(data = stab_corr,
               aes(x = 0, y = 0, 
                   xend = PC1*5, 
                   yend = PC2*5, 
                   colour = stab), 
               arrow = arrow(length = unit(0.3,"cm"))) +
  geom_text(data = stab_corr,
            aes(x = PC1*5.5, y = PC2*5.5, label = stab, colour = stab),
            size = 5) +
  scale_colour_manual(values = stability_palette) +
  labs(x = "PC1", y = "PC2", title = "PCA of network structure with stability gradients") +
  figure_theme() +
  theme(legend.position = 'right')



############################################################
# End
############################################################