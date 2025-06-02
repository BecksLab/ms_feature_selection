##Load packages
library(caret)
library(corrplot)
library(FactoMineR)
library(here)
library(psych)
library(ggrepel)
library(patchwork)
library(tidyverse)

setwd(here::here())

# import network summary data
topology <- read.csv("data/vermaat_2009/vermaat_summary.csv") %>%
  rbind(read.csv("data/mangal/mangal_summary.csv")) %>%
  select(!id)

# create subset that only uses the Vermaat vars...
topology_subset <- topology %>%
  select(-c(links, S1, S2, S4, S5, diameter, ρ, centrality, complexity))

# get an idea of the number of NA vals
sapply(topology, function(x) sum(is.na(x))) / nrow(topology) * 100

# same for prop of zeros
colSums(topology==0, na.rm = TRUE)/nrow(topology)*100

#### PCA ####

#' here we will dynamically create and save the outputs of the PCA. The 
#' correlations will be output as a .csv file and the figures will be 
#' saved into a list so that they can be combined into a composite figure

data_list <- list(complete = topology, 
                  subset = topology_subset)
data_names <- c("complete", "subset")

plot_list <- vector(mode = "list", length = 2)

for (i in 1:length(data_list)) {
  
  df = data_list[[i]]
  
  pca <- PCA(df, scale.unit = TRUE, graph = FALSE)
  
  # get variance explained for labels
  variance = pca$eig %>%
    as.data.frame() %>%
    pull(`percentage of variance`)
  
  
  # pull all vars signif correlated with first three dims
  dim_descrip <- dimdesc(pca, axes = 1:3)
  
  signif_corrs <- 
    dim_descrip[["Dim.1"]] %>%
    as.data.frame() %>%
    rownames_to_column(., var = "Property") %>%
    mutate(dimension = "Dim.1") %>%
    rbind(dim_descrip[["Dim.2"]] %>%
            as.data.frame() %>%
            rownames_to_column(., var = "Property") %>%
            mutate(dimension = "Dim.2")) %>%
    rbind(dim_descrip[["Dim.3"]] %>%
            as.data.frame() %>%
            rownames_to_column(., var = "Property") %>%
            mutate(dimension = "Dim.3")) %>%
    mutate(quanti.correlation = round(quanti.correlation, digits = 2))
  
  pca$var$cor %>%
    as.data.frame() %>%
    rownames_to_column(., var = "Property") %>%
    select(Property, Dim.1, Dim.2, Dim.3) %>%
    pivot_longer(!Property,
                 names_to = "dimension",
                 values_to = "quanti.correlation") %>%
    mutate(quanti.correlation = round(quanti.correlation, digits = 2)) %>%
    full_join(signif_corrs) %>%
    mutate(dimension = case_when(dimension == "Dim.1" ~ paste("PCA 1 (", round(variance[1]), "%)", sep = ""),
                                 dimension == "Dim.2" ~ paste("PCA 2 (", round(variance[2]), "%)", sep = ""),
                                 dimension == "Dim.3" ~ paste("PCA 3 (", round(variance[3]), "%)", sep = "")),
           quanti.correlation = case_when(abs(quanti.correlation) > 0.53 & abs(quanti.correlation) <= 0.53 & quanti.p.value <= 0.05 ~ paste("**", quanti.correlation, "**", sep = ""),
                                          abs(quanti.correlation) > 0.66 & quanti.p.value <= 0.01 ~ paste("**", quanti.correlation, "**", sep = ""),
                                          .default = as.character(quanti.correlation))) %>%
    select(!quanti.p.value) %>%
    pivot_wider(names_from = dimension,
                values_from = quanti.correlation) %>%
    write.csv(.,
              paste0("../tables/allNetworks_corr_", data_names[i], ".csv"),
              row.names = FALSE)
  
  # plot pca
  
  # get signif variable for sims one and two - so we can assign colours
  plot_list[[i]] <- 
    pca$var$cor %>%
    as.data.frame() %>%
    rownames_to_column(., var = "Property") %>%
    select(Property, Dim.1, Dim.2, Dim.3) %>%
    pivot_longer(!Property,
                 names_to = "dimension",
                 values_to = "quanti.correlation") %>%
    mutate(quanti.correlation = round(quanti.correlation, digits = 2)) %>%
    full_join(signif_corrs) %>%
    filter(dimension != "Dim.3") %>%
    mutate(colour = case_when(abs(quanti.correlation) > 0.53 & abs(quanti.correlation) <= 0.53 & quanti.p.value <= 0.05 ~ "red",
                              abs(quanti.correlation) > 0.66 & quanti.p.value <= 0.01 ~ "red",
                              .default = "black")) %>%
    select(Property, dimension, colour) %>%
    pivot_wider(names_from = dimension,
                values_from = colour) %>%
    mutate(significant = if_else(Dim.1 == "red" | Dim.2 == "red", "signif", "non-signif")) %>%
    select(Property, significant) %>%
    full_join(.,
              pca[["var"]][["coord"]] %>%
                as.data.frame() %>%
                rownames_to_column(., var = "Property")) %>%
    ggplot(.,
           aes(x = Dim.1,
               y = Dim.2)) +
    geom_hline(yintercept = 0,
               colour = "grey60",
               linetype = 2) +
    geom_vline(xintercept = 0,
               colour = "grey60",
               linetype = 2) +
    geom_point(aes(colour = significant)) +
    geom_text_repel(aes(label = Property)) +
    theme_classic() +
    lims(x= c(-1, 1), y = c(-1, 1)) +
    labs(x = paste("PCA 1 (", round(variance[1]), "%)",sep = ""),
         y = paste("PCA 2 (", round(variance[2]), "%)",sep = ""))
  
}

(plot_list[[1]] + labs(title = "A. Complete")) / 
  (plot_list[[2]] + labs(title = "B. Vermaat subset")) +
  plot_layout(guides = 'collect')

ggsave("../figures/pca_allNetworks.png",
       width = 4000,
       height = 6500,
       units = "px",
       dpi = 600)

