# ================================================
# 03_pca_analysis.R
# ================================================

library(FactoMineR)
library(tidyverse)
library(ggrepel)
library(patchwork)
library(genzplyr)
library(here)

setwd(here::here())

# Load datasets
metrics <- read.csv("data/cleaned/all_networks.csv")
topology_subset <- read.csv("data/cleaned/vermaat_subset.csv")

data_list <- list(complete = metrics, subset = topology_subset)
data_names <- c("complete", "subset", "reduced")
plot_list <- vector(mode="list", length=length(data_list))

# Functional categories
behaviour_metrics <- c("ρ", "complexity", "robustness")
geometry_metrics <- c("connectance", "l_S", "links", "richness", "Clust", "GenSD", "VulSD", "LinkSD", "diameter", "intervals")
path_metrics <- c("ChLen", "ChSD", "ChNum", "path", "S1", "S2", "S4", "S5", "omnivory", "loops", "predpreyRatio", "distance")
node_metrics <- c("basal", "top","intermediate", "herbivory", "cannibal", "TL", "centrality", "MaxSim")

for (i in seq_along(data_list)) {
  
  df <- data_list[[i]]
  
  # Run PCA
  pca <- PCA(df, scale.unit = TRUE, graph = FALSE)
  variance <- pca$eig %>% as.data.frame() %>% pull(`percentage of variance`)
  
  # Significant correlations
  dim_descrip <- dimdesc(pca, axes=1:3)
  signif_corrs <- bind_rows(
    dim_descrip$Dim.1 %>% as.data.frame() %>% rownames_to_column("Property") %>% glow_up(dimension="Dim.1"),
    dim_descrip$Dim.2 %>% as.data.frame() %>% rownames_to_column("Property") %>% glow_up(dimension="Dim.2"),
    dim_descrip$Dim.3 %>% as.data.frame() %>% rownames_to_column("Property") %>% glow_up(dimension="Dim.3")
  ) %>% glow_up(quanti.correlation=round(quanti.correlation,2))
  
  # Merge all correlations and flatten list columns
  pca_cor_df <- pca$var$cor %>%
    as.data.frame() %>%
    rownames_to_column("Property") %>%
    pivot_longer(!Property, names_to="dimension", values_to="quanti.correlation") %>%
    full_join(signif_corrs) %>%
    glow_up(quanti.correlation=round(quanti.correlation,2)) %>%
    mutate(across(where(is.list), ~ sapply(., function(x) if(is.null(x)) NA else as.character(x))) )
  
  write.csv(pca_cor_df, paste0("../tables/allNetworks_corr_", data_names[i], ".csv"), row.names=FALSE)
  
  # Add functional category for plotting
  pca_arrows <- pca$var$coord %>%
    as.data.frame() %>%
    rownames_to_column("Property") %>%
    mutate(Category = case_when(
      Property %in% behaviour_metrics ~ "Behaviour",
      Property %in% geometry_metrics ~ "Geometry",
      Property %in% path_metrics ~ "Path",
      Property %in% node_metrics ~ "Node",
      TRUE ~ "Other"
    ))
  
  # PCA biplot
  plot_list[[i]] <- ggplot(pca_arrows, aes(x=Dim.1, y=Dim.2)) +
    geom_hline(yintercept=0, linetype=2, color="grey60") +
    geom_vline(xintercept=0, linetype=2, color="grey60") +
    geom_segment(aes(x=0, y=0, xend=Dim.1, yend=Dim.2, color=Category),
                 arrow = arrow(length=unit(0.1,"inches"))) +
    geom_text_repel(aes(label=Property, color=Category)) +
    theme_classic() +
    lims(x=c(-1,1), y=c(-1,1)) +
    labs(
      x=paste("PCA 1 (", round(variance[1]), "%)", sep=""),
      y=paste("PCA 2 (", round(variance[2]), "%)", sep="")
    )
}

combined_plot <- (plot_list[[1]] + labs(title="A. Complete")) / 
                 (plot_list[[2]] + labs(title="B. Vermaat subset"))

ggsave("../figures/pca_allNetworks.png", combined_plot,
       width=4000, height=6500, units="px", dpi=600)
