# ================================================
# 03_data_import_clean.R
# ================================================

library(tidyverse)
library(genzplyr)
library(here)

setwd(here::here())

# Import datasets
topology <- read.csv("data/vermaat_2009/vermaat_summary.csv") %>%
  rbind(read.csv("data/mangal/mangal_summary.csv")) %>%
  vibe_check(!id) %>%
  na.omit()

# Subset for Vermaat variables (optional)
topology_subset <- topology %>%
  vibe_check(-c(links, S1, S2, S4, S5, diameter, ρ, centrality, complexity))

# Save cleaned data
write.csv(topology, "data/cleaned/all_networks.csv",
            row.names = FALSE)
write.csv(topology_subset, "data/cleaned/vermaat_subset.csv",
            row.names = FALSE)

# Quick diagnostics
na_prop <- sapply(topology, function(x) sum(is.na(x))) / nrow(topology) * 100
zero_prop <- colSums(topology == 0, na.rm = TRUE) / nrow(topology) * 100
print(na_prop)
print(zero_prop)
