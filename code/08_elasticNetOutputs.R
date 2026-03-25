#---------------------------------------------------
# Elastic net regression (part 2)
#---------------------------------------------------
# Post model processing - import results from 07_elasticNetModel.R
# this avoids the long wait period for runnin code

#packages

library(genzplyr)
library(patchwork)
library(tidyverse)
library(tune)

set.seed(66)
setwd(here::here())
source("lib/plotting_themes.R")

############################################################
# 1. Load CSV outputs
############################################################

# Model summary from Julia
model_performance <- read_csv("data/outputs/elasticNet_summary.csv") %>%
  mutate(
    rep_name = case_when(
      rep_name == "medoids"   ~ "Cluster Dominant",
      rep_name == "pca_score" ~ "PC Score",
      rep_name == "dominant"  ~ "PC Dominant",
      rep_name == "complexity"  ~ "Complexity",
      TRUE ~ rep_name
    ),
    alpha_val = alpha,
    r2_val = r2
  )

# Coefficients from Julia
all_estimates <- read_csv("data/outputs/elasticNet_coefficients.csv") %>%
  mutate(
    rep_name = case_when(
      rep_name == "medoids"   ~ "Cluster Dominant",
      rep_name == "pca_score" ~ "PC Score",
      rep_name == "dominant"  ~ "PC Dominant",
      rep_name == "complexity"  ~ "Complexity",
      TRUE ~ rep_name
    )
  )

print(model_performance)

############################################################
# 2. Metadata for module plotting
############################################################

clust_metada <-
  left_join(
    read_csv("../tables/metric_clusters_auto.csv"),
    read_csv("../tables/module_summary_7clusters.csv") %>%
      vibe_check(Cluster, label),
    by = "Cluster"
  )

############################################################
# 4. All coefficients in long format
############################################################

directed_estimates <- all_estimates %>%
  left_join(clust_metada, by = c("term" = "Metric")) %>%
  glow_up(
    direction = if_else(estimate > 0, "Positive", "Negative")
  ) %>%
  left_join(pal_df, by = "label") %>%
  glow_up(
    label = case_when(rep_name == "PC Score" ~ "PCA Axis",
                      rep_name == "Complexity" ~ "Complexity",
                      TRUE ~ label),
    colour = case_when(rep_name == "PC Score" ~ pca_col,
                       rep_name == "Complexity" ~ complexity_col,
                       TRUE ~ colour))

############################################################
# 5. Combined summary
############################################################

final_summary <- model_performance %>%
  mutate(
    architecture = case_when(
      alpha >= 0.75 ~ "Sparse (Lasso-like)",
      alpha <= 0.25 ~ "Distributed (Ridge-like)",
      TRUE              ~ "Intermediate (Elastic)"
    )
  )

final_summary %>%
  arrange(desc(r2))

############################################################
# Plotting
############################################################

# variance explained
ggplot(directed_estimates, 
       aes(x = rep_name, 
           y = abs(variance_explained), 
           fill = label)) +
  geom_col() +
  facet_wrap(~ metric) +
  ylim(0, 1.0) +
  labs(title = "Variance in Stability Explained by Structural Modules",
       y = "Absolute Variance Explained (CV R-squared)",
       x = NULL) +
  scale_fill_manual(values = setNames(variance_decomposition$colour, 
                                      as.character(variance_decomposition$label)),
                    breaks = variance_decomposition$label[!(variance_decomposition$label %in% c("PCA Axis", "Complexity"))],
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
                    breaks = directed_estimates$label[!(directed_estimates$label %in% c("PCA Axis", "Complexity"))],
                    name = "Module")+
  scale_colour_manual(values = setNames(directed_estimates$colour, 
                                        as.character(directed_estimates$label)),
                      breaks = directed_estimates$label[!(directed_estimates$label %in% c("PCA Axis", "Complexity"))],
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


############################################################
# END
############################################################