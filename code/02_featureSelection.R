##Load packages
library(caret)
library(corrplot)
library(FactoMineR)
library(here)
library(tidyverse)

setwd(here::here())

# import network summary data
topology_mangal <- read.csv("data/mangal/mangal_summary.csv") %>%
  na.omit()

# 1: look at correlation first

#scale all features
topology_scaled <- scale(topology_mangal[2:ncol(topology_mangal)],
                         center = TRUE, scale = TRUE)

# correlation matrix
cor_matrix <- cor(topology_scaled)

corrplot(cor_matrix, order = "hclust")

# remove highly correlated vars
# select the cutoff 
corr_cutoff <- 0.85
# subset by cutoff
hc = findCorrelation(cor_matrix, cutoff = corr_cutoff)

# subset 
topology_scaled_subset <- topology_scaled[,-hc]

cor_mat_subset <- cor(topology_scaled_subset)

corrplot(cor_mat_subset, order = "hclust")

# 2: look at correlation first

pca <- PCA(topology_scaled)

# find number of dims describing 90 % of variance
dim_num <- pca$eig %>%
  as.data.frame() %>%
  filter(`cumulative percentage of variance` < 90) %>%
  nrow() %>% as.numeric()

dim_descrip <- dimdesc(pca, axes = 1:5)

dim_descrip$Dim.1 %>%
  as.data.frame() %>%
  filter(abs(quanti.correlation) > 0.8) %>%
  rbind(dim_descrip$Dim.2 %>%
          as.data.frame() %>%
          filter(abs(quanti.correlation) > 0.8))

# assign broader categories

grps <- 
  tibble(var = row.names(pca[["var"]][["coord"]])) %>%
  reframe(cat = case_when(var %in% c("S1", "S2", "S5", "S4") ~ "Motif",
                          var %in% c("richness", "links", "connectance", "complexity", "l_S") ~ "Structure",
                          var %in% c("diameter", "trophic_level", "path", "cl_mean", "cl_std", "log_fc", "top", "basal", "distance") ~ "Energy",
                          var %in% c("generality", "vulnerability", "link_SD", "cannibal", "omnivory", "herbivory", "intermediate", "") ~ "Node"))

variance = pca$eig %>%
  as.data.frame() %>%
  pull(`cumulative percentage of variance`)

# try plot as points

ggplot(data = pca[["var"]][["coord"]]) +
  geom_point(aes(x = Dim.1,
                 y = Dim.2,
                 colour = grps$cat)) +
  theme_classic() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  lims(x= c(-1, 1), y = c(-1, 1)) +
  labs(x = paste("PCA 1 ", round(variance[1]), "%",sep = ""),
       y = paste("PCA 2 ", round(variance[2]), "%",sep = ""))
  
