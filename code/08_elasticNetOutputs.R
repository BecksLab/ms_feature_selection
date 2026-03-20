#---------------------------------------------------
# Elastic net regression (part 2)
#---------------------------------------------------
# Post model processing - import results from 07_elasticNetModel.R
# this avoids the long wait period for runnin code

#packages

library(genzplyr)
library(patchwork)
library(tidyverse)

set.seed(66)
setwd(here::here())
source("lib/plotting_themes.R")

############################################################
# 1. Load Data + extract
############################################################

all_results <- readRDS("data/outputs/elastic_net_results.rds") %>%
  glow_up(rep_name = case_when(rep_name == "medoids" ~ "Cluster Dominant",
                               rep_name == "pca_score" ~ "PC Score",
                               rep_name == "dominant" ~ "PC Dominant"))

# get cluster/module masterlist (for plotting)
clust_metada <-
  left_join(read_csv("../tables/metric_clusters_auto.csv"),
            read_csv("../tables/module_summary_7clusters.csv") %>%
              vibe_check(c(Cluster, label)))

# Extract Global Model Performance
model_performance <- all_results %>%
  glow_up(
    r2_val = map_dbl(output, ~show_best(.x$results, metric = "rsq", n = 1)$mean),
    alpha_val = map_dbl(output, ~.x$best_params$mixture),
    lambda_val = map_dbl(output, ~.x$best_params$penalty)
  ) %>%
  vibe_check(metric, rep_name, r2_val, alpha_val, lambda_val)

print(model_performance)

# partition variance
get_module_variance <- function(enet_output, metadata) {
  cv_r2 <- show_best(enet_output$results, metric = "rsq", n = 1)$mean
  
  enet_output$coefficients %>%
    left_join(metadata, by = c("term" = "Metric")) %>%
    glow_up(
      rel_imp = estimate^2 / sum(estimate^2),
      abs_var = rel_imp * cv_r2
    ) %>%
    group_by(Cluster, label) %>%
    no_cap(module_var = sum(abs_var), .groups = "drop")
}

# Apply and unnest
variance_decomposition <- all_results %>%
  glow_up(mod_vars = map(output, ~get_module_variance(.x, clust_metada))) %>%
  vibe_check(metric, rep_name, mod_vars) %>%
  unnest(mod_vars) %>%
  # add relevant colour codes
  left_join(pal_df) %>%
  glow_up(label = if_else(is.na(colour), "PCA Axis", label),
          colour = if_else(is.na(colour), pca_col, colour))

# Extract all coefficients into a long-format dataframe
all_estimates <- all_results %>%
  glow_up(estimates = map(output, ~ .x$coefficients)) %>%
  vibe_check(metric, rep_name, estimates) %>%
  unnest(estimates) %>%
  # Filter out intercept to focus on structural drivers
  yeet(term != "(Intercept)")

# View the top drivers for each stability metric
all_estimates %>%
  squad_up(metric, rep_name) %>%
  slice_max(order_by = abs(estimate), n = 5)

# Merge with metadata for module-level insight
directed_estimates <- all_estimates %>%
  left_join(clust_metada, by = c("term" = "Metric")) %>%
  glow_up(direction = if_else(estimate > 0, "Positive", "Negative")) %>%
  # add relevant colour codes
  left_join(pal_df) %>%
  glow_up(label = if_else(is.na(colour), "PCA Axis", label),
          colour = if_else(is.na(colour), pca_col, colour))

# Calculate the average effect size per module
module_trends <- directed_estimates %>%
  squad_up(metric, label, direction) %>%
  no_cap(
    mean_estimate = mean(estimate),
    count = n(),
    .groups = "drop"
  )

############################################################
# 5. Combined summary
############################################################

# Identify the 'Dominant' module for each model
dominant_modules <- variance_decomposition %>%
  squad_up(metric, rep_name) %>%
  slice_max(order_by = module_var, n = 1) %>%
  lowkey(top_module = label, top_module_var = module_var) %>%
  vibe_check(metric, rep_name, top_module, top_module_var)

# Identify the 'Top Predictor' (estimate) for each model
top_estimates <- directed_estimates %>%
  squad_up(metric, rep_name) %>%
  slice_max(order_by = abs(estimate), n = 1) %>%
  vibe_check(metric, rep_name, term, estimate) %>%
  lowkey(top_predictor = term, max_estimate = estimate)

# Fold everything into model_performance
final_summary <- model_performance %>%
  left_join(dominant_modules, by = c("metric", "rep_name")) %>%
  left_join(top_estimates, by = c("metric", "rep_name")) %>%
  glow_up(
    # Create an 'Architecture' label based on alpha
    architecture = case_when(
      alpha_val >= 0.75 ~ "Sparse (Lasso-like)",
      alpha_val <= 0.25 ~ "Distributed (Ridge-like)",
      TRUE              ~ "Intermediate (Elastic)"
    )
  )

final_summary %>% 
  slay(desc(r2_val))

############################################################
# Plotting
############################################################

# variance explained
ggplot(variance_decomposition, 
       aes(x = rep_name, 
           y = module_var, 
           fill = label)) +
  geom_col() +
  facet_wrap(~ metric) +
  ylim(0, 1.0) +
  labs(title = "Variance in Stability Explained by Structural Modules",
       y = "Absolute Variance Explained (CV R-squared)",
       x = NULL) +
  scale_fill_manual(values = setNames(variance_decomposition$colour, 
                                      as.character(variance_decomposition$label)),
                    breaks = variance_decomposition$label[variance_decomposition$label != "PCA Axis"],
                    name = "Module") +
  figure_theme()

ggsave("../figures/stability_variance.png",
       width = 6000,
       height = 4000,
       units = "px")


# estimates/coefficients
ggplot(directed_estimates, 
       aes(x = reorder(term, estimate), 
           y = estimate, 
           fill = label, 
           colour = label)) +
  geom_hline(yintercept = 0,
             colour = "#001628",
             linewidth = 0.6) +
  geom_segment(aes(xend = term, yend = 0),
               linewidth = 1) +
  geom_point(size = 8,
             shape = 21,
             colour = "white")+
  coord_flip() +
  facet_grid(rows = vars(rep_name),
             cols = vars(metric),
             scales = "free_y") +
  scale_fill_manual(values = setNames(directed_estimates$colour, 
                                      as.character(directed_estimates$label)),
                    breaks = directed_estimates$label[directed_estimates$label != "PCA Axis"],
                    name = "Module")+
  scale_colour_manual(values = setNames(directed_estimates$colour, 
                                        as.character(directed_estimates$label)),
                      breaks = directed_estimates$label[directed_estimates$label != "PCA Axis"],
                      name = "Module") +
  labs(title = "Structural Drivers of Stability (PCA Scores)",
       subtitle = "Standardised Coefficients (direction and magnitude)",
       x = NULL,
       y = "Standardised Estimate") +
  figure_theme()

ggsave("../figures/stability_estimate.png",
       width = 6000,
       height = 4000,
       units = "px")

# alpha
ggplot(final_summary, 
       aes(x = alpha_val, 
           y = metric, 
           colour = rep_name,
           shape = rep_name)) +
  geom_vline(xintercept = 0.5, 
             linetype = "dashed", 
             alpha = 0.3,
             colour = "#001628") +
  geom_point(size = 8,
             alpha = 0.6) +
  # Add labels to the axes to make the Ridge/Lasso distinction clear
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_colour_manual(values = secondary_palette) +
  labs(
    title = "The Architecture of Structural Control",
    subtitle = "Elastic Net Mixing Parameter (α): 0 = Distributed (Ridge), 1 = Sparse (Lasso)",
    x = expression(paste("Mixing Parameter (", alpha, ")")),
    y = "Stability Metric",
    colour = "Representation",
    shape = "Representation"
  ) +
  figure_theme()

ggsave("../figures/stability_alpha.png",
       width = 6000,
       height = 4000,
       units = "px")

main_plot <- ggplot(final_summary, 
                    aes(x = alpha_val, y = r2_val)) +
  geom_vline(xintercept = 0.5, 
             linetype = "dashed", 
             alpha = 0.9,
             colour = "#001628",
             linewidth = 0.8) +
  geom_point(aes(fill = rep_name, 
                 size = top_module_var),
             shape = 21,
             colour = "white",
             stroke = 0.6) +
  scale_size_continuous(
    range = c(2, 12), 
    breaks = c(0.1, 0.4, 0.7),
    limits = c(0, 0.9),
    labels = c("10%", "40%", "70%")) +
  scale_fill_manual(values = secondary_palette) +
  guides(
    size = guide_legend(
      title = "Module Contribution",
      override.aes = list(
        shape = 21, 
        fill = "#001628",
        colour = "white", 
        size = c(2, 7, 12))),
    fill = guide_legend(
      title = "Representation",
      override.aes = list(size = 6, shape = 21))) +
  facet_wrap(~metric, ncol = 1) +
  labs(
    title = "Predictive Power vs. Control Architecture",
    x = expression(paste("Sparsity (", alpha, ")")),
    y = expression(paste("Variance Explained (CV ", R^2, ")")),
    size = "Module Contribution",
    fill = "Representation"
  ) +
  figure_theme() +
  theme(
    panel.grid.major.x = element_blank(),
    legend.key.height = unit(1.2, "cm"), # Essential for big points
    legend.spacing.y = unit(0.2, "cm")
  )


# Create a tiny plot just for the arrows
arrow_plot <- ggplot() +
  annotate("segment", x = 0.48, xend = 0.05, y = 0, yend = 0, 
           arrow = arrow(length = unit(0.2, "cm")), color = "#001628", linewidth = 0.7) +
  annotate("segment", x = 0.52, xend = 0.95, y = 0, yend = 0, 
           arrow = arrow(length = unit(0.2, "cm")), color = "#001628", linewidth = 0.7) +
  annotate("text", x = 0.25, y = -0.6, label = "Distributed Control", 
           size = 3.5, family = "space", color = "#001628", lineheight = 0.9) +
  annotate("text", x = 0.75, y = -0.6, label = "Sparse Control", 
           size = 3.5, family = "space", color = "#001628", lineheight = 0.9) +
  xlim(0, 1) + 
  ylim(-1.8, 0.5) +
  theme_void() +
  theme(
    plot.margin = margin(t = 10, r = 20, b = 20, l = 20),
    plot.background = element_rect(fill = "transparent", color = NA)
  ) +
  coord_cartesian(clip = "off")

# Combine them
main_plot / arrow_plot + plot_layout(heights = c(10, 1))

ggsave("../figures/stability_power_v_control.png",
       width = 5000,
       height = 4500,
       units = "px")

############################################################
# END
############################################################