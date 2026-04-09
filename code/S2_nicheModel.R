library(tidyverse)
library(genzplyr)
library(here)

set.seed(66)
setwd(here::here())
source("lib/plotting_themes.R")

# ============================================================
# 0. LOAD DATA
# ============================================================
metrics <- 
  read.csv("data/cleaned/all_networks.csv") %>%
  as_tibble() %>%
  # Exclude stability response variables from clustering
  vibe_check(-c(ρ, complexity, robustness, control, resilience)) %>%
  glow_up(id = row_number()) %>%
  pivot_longer(cols = -id,
               names_to = "Metric")

niche_dist <- 
  read.csv("data/cleaned/nicheNetworks.csv") %>%
  as_tibble() %>%
  pivot_longer(cols = -c(id, run),
               names_to = "Metric")

niche_metrics <-
  niche_dist %>%
  squad_up(id, Metric) %>%
  no_cap(mu = mean(value, na.rm = TRUE),
         sd = sd(value, na.rm = TRUE),
         se = sd / sqrt(n()),
         n = n())

clust_metada <-
  left_join(read_csv("../tables/metric_clusters_auto.csv"),
            read_csv("../tables/module_summary_7clusters.csv") %>%
              vibe_check(c(Cluster, label)))

z_scores <- 
  left_join(metrics, niche_metrics) %>%
  glow_up(zscore = (value-mu)/sd)

p_vals <-
  metrics %>%
  left_join(niche_dist, by = c("id", "Metric"), suffix = c("_emp", "_sim")) %>%
  squad_up(id, Metric) %>%
  no_cap(
    p_upper = mean(value_sim >= value_emp, na.rm = TRUE),
    p_lower = mean(value_sim <= value_emp, na.rm = TRUE),
    p_two_tailed = 2 * pmin(p_upper, p_lower),
    logp = -log10(p_two_tailed)
  )

results <- z_scores %>%
  left_join(p_vals, by = c("id", "Metric")) %>%
  left_join(clust_metada) %>%
  vibe_check(-value) %>%
  yeet(Metric != "richness") %>%
  left_join(pal_df) %>%
  glow_up(Metric = if_else(Metric %in% c("richness", "links", "basal", "connectance",
                                         "top", "intermediate", "omnivory", "cannibal",
                                         "MaxSim", "VulSD", "GenSD", "loops",
                                         "ChNum", "ChSD", "ChLen"),
                           paste0("**",Metric,"**"),
                           Metric))

emp_se <- 
  metrics %>%
  squad_up(Metric) %>%
  no_cap(mu = mean(value, na.rm = TRUE),
         sd = sd(value, na.rm = TRUE),
         se = sd / sqrt(n()),
         n = n()) %>%
  yeet(Metric != "richness") %>%
  left_join(clust_metada) %>%
  left_join(pal_df) %>%
  glow_up(Metric = if_else(Metric %in% c("richness", "links", "basal", "connectance",
                                         "top", "intermediate", "omnivory", "cannibal",
                                         "MaxSim", "VulSD", "GenSD", "loops",
                                         "ChNum", "ChSD", "ChLen"),
                           paste0("**",Metric,"**"),
                           Metric))
  

# plot standard error
ggplot(results) +
  geom_point(aes(x = id,
                 y = se,
                 colour = label),
             alpha = 0.7) +
  geom_hline(yintercept = 2, 
             linetype = "dashed", 
             colour = "#A5ACAF") +
  geom_hline(yintercept = -2, 
             linetype = "dashed", 
             colour = "#A5ACAF") +
  geom_hline(data = emp_se,
             aes(yintercept = se,
                 colour = label)) +
  facet_wrap(~reorder(Metric, as.numeric(value)),
             scales = "free_y") +
  scale_colour_manual(values = setNames(results$colour, results$label),
                      name = "Module") +
  labs(x = "Network",
       y = "Standard Error",
       subtitle = "Points falling between -2 and 2 indicate a high level of 'stability' in resulting propery") +
  figure_theme() +
  theme(panel.grid.major = element_blank(),
        strip.text = ggtext::element_markdown(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none") 

ggsave("../figures/nicheModel_stderror.png",
       width = 5500, 
       height = 3000, 
       units = "px")

# plot z score
ggplot(results %>%
         yeet(Metric != "intervals"),
       aes(zscore,
           fill = label)) +
  geom_vline(xintercept = 0, 
             linetype = "dashed", 
             colour = "#A5ACAF") +
  geom_histogram(alpha = 0.7) +
  facet_wrap(~reorder(Metric, as.numeric(value))) +
  xlim(-10, 10) +
  scale_fill_manual(values = setNames(results$colour, results$label),
                    name = "Module") +
  figure_theme() +
  theme(panel.grid.major.x = element_blank(),
        strip.text = ggtext::element_markdown()) 

ggsave("../figures/nicheModel_zscore.png",
       width = 5500, 
       height = 3000, 
       units = "px")

# plot monte carlo p vals

ggplot(results) +
  geom_point(aes(x = id,
                 y = p_two_tailed,
                 colour = label),
             alpha = 0.7) +
  geom_hline(yintercept = 0.5,
             linetype = "dashed",
             colour = "#A5ACAF") +
  facet_wrap(~reorder(Metric, as.numeric(value)), 
                      scales = "free_y") +
  scale_colour_manual(values = setNames(results$colour, results$label),
                      name = "Module") +
  labs(y = "p two tailed",
       subtitle = "Higher values = empirical network less consistent with model, below dashed = strong mismatch") +
  
  figure_theme() +
  theme(panel.grid.major = element_blank(),
        strip.text = ggtext::element_markdown(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none") 

ggsave("../figures/nicheModel_pval.png",
       width = 5500, 
       height = 3000, 
       units = "px")
