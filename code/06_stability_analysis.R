# ================================================
# 04_stability_analysis.R
# ================================================

library(tidyverse)
library(randomForest)
library(here)

setwd(here::here())

# Load cleaned data
topology <- read.csv("data/cleaned/all_networks.csv")

# Define your structural anchors
anchors <- c("S1", "richness", "predpreyRatio", "distance", "herbivory", "MaxSim")

# Run the three models
set.seed(42) # For reproducibility

# Model A: Persistence
rf_robust <- randomForest(as.formula(paste("robustness ~", paste(anchors, collapse="+"))), 
                          data=topology, importance=TRUE, ntree=1000)

# Model B: Dampening
rf_rho <- randomForest(as.formula(paste("ρ ~", paste(anchors, collapse="+"))), 
                       data=topology, importance=TRUE, ntree=1000)

# Model C: Organization (SVD Complexity)
rf_comp <- randomForest(as.formula(paste("complexity ~", paste(anchors, collapse="+"))), 
                        data=topology, importance=TRUE, ntree=1000)

# Extract importance for all three models
imp_robust <- as.data.frame(importance(rf_robust))
imp_rho    <- as.data.frame(importance(rf_rho))
imp_comp   <- as.data.frame(importance(rf_comp))

# Combine the %IncMSE column (which is usually the first column) into one table
importance_table <- data.frame(
  Pillar = rownames(imp_robust),
  Robustness_IncMSE = imp_robust$'%IncMSE',
  Rho_IncMSE = imp_rho$'%IncMSE',
  Complexity_IncMSE = imp_comp$'%IncMSE'
)

print(importance_table)

# Create a 3-panel plot for your results
par(mfrow=c(1,3)) 
varImpPlot(rf_robust, type=1, main="Drivers of Persistence")
varImpPlot(rf_rho,    type=1, main="Drivers of Dampening")
varImpPlot(rf_comp,   type=1, main="Drivers of Organization")

# 3. Extract and compare Performance (R-squared)
perf <- data.frame(
  Stability_Type = c("Robustness", "Spectral Radius", "SVD Complexity"),
  R_Squared = c(max(rf_robust$rsq), max(rf_rho$rsq), max(rf_comp$rsq))
)
print(perf)


library(pdp)
p1 <- partial(rf_robust, pred.var = "MaxSim", plot = TRUE)
p2 <- partial(rf_rho,    pred.var = "predpreyRatio", plot = TRUE)
p3 <- partial(rf_comp,   pred.var = "S1", plot = TRUE)

library(gridExtra)
grid.arrange(p1, p2, p3, ncol=3)
