##Load packages
library(caret)
library(corrplot)
library(FactoMineR) 
library(here)
library(tidyverse)

# import network summary data
topology_mangal <- read.csv("data/mangal/mangal_summary.csv") %>%
  na.omit()

# 1: look at correlation first

#scale all features
topology_scaled <- scale(topology_mangal[1:ncol(topology_mangal)],
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

pca <- PCA(topology_mangal[1:ncol(topology_mangal)], scale.unit = T, graph=T)

# find number of dims describing 90 % of variance
dim_num <- pca$eig %>%
  as.data.frame() %>%
  filter(`cumulative percentage of variance` < 91) %>%
  nrow() %>% as.numeric()

dim_descrip <- dimdesc(pca, axes = 1:dim_num)

dim_descrip$Dim.1 %>%
  as.data.frame() %>%
  filter(abs(quanti.correlation) > 0.8) %>%
  rbind(dim_descrip$Dim.2 %>%
          as.data.frame() %>%
          filter(abs(quanti.correlation) > 0.8))
