# ================================================
# 04_stability_analysis.R
# ================================================

library(tidyverse)
library(randomForest)
library(here)

setwd(here::here())

# Load cleaned data
topology <- read.csv("../data/cleaned/all_networks.csv")

# Linear model using PCA scores (from previous PCA)
pca_scores <- read.csv("../tables/allNetworks_corr_reduced.csv") # optional: PCA scores CSV
topology_pca <- cbind(topology, pca_scores)

lm_robust <- lm(robustness ~ Dim.1 + Dim.2 + Dim.3, data=topology_pca)
summary(lm_robust)

# Random Forest on reduced metrics
set.seed(66)
rf_model <- randomForest(robustness ~ connectance + Clust + complexity + MaxSim, 
                         data=topology, importance=TRUE)
varImpPlot(rf_model)
