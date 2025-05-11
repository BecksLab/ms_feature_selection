##Load packages
library(caret)
library(corrplot)
library(FactoMineR)
library(here)
library(psych)
library(tidyverse)

setwd(here::here())

# import network summary data
topology_vermaat <- read.csv("data/vermaat_2009/vermaat_summary.csv") %>%
  select(-id) %>%
  select(-c(links, S1, S2, S4, S5, diameter, ρ, centrality, complexity))

# get an idea of the number of NA vals
sapply(topology_vermaat, function(x) sum(is.na(x))) / nrow(topology_vermaat) * 100

# same for prop of zeros
colSums(topology_vermaat==0, na.rm = TRUE)/nrow(topology_vermaat)*100

topology_vermaat_subset <-
  topology_vermaat %>%
  select(-c(loops)) %>%
  na.omit()

# 1: look at correlation first

#scale all features
topology_scaled <- scale(topology_vermaat_subset,
                         center = TRUE, scale = TRUE)

# correlation matrix
cor_matrix <- cor(topology_scaled)

corrplot(cor_matrix, order = "hclust")

# remove highly correlated vars
# select the cutoff 
corr_cutoff <- 0.8
# subset by cutoff
hc = findCorrelation(cor_matrix, cutoff = corr_cutoff, names = TRUE)

# subset 
topology_subset <- select(topology_vermaat_subset,-all_of(hc))
topology_scaled_subset <- scale(topology_subset,
                                center = TRUE, scale = TRUE)

cor_mat_subset <- cor(topology_scaled_subset)

corrplot(cor_mat_subset, order = "hclust")

#KMO Test

KMO(topology_subset)
cortest.bartlett(topology_vermaat_subset)


# Determine number of factors to extract
ev <- eigen(cor(topology_scaled_subset)) # get eigenvalues
ev$values

scree(topology_scaled_subset, pc=FALSE)

fa.parallel(topology_scaled_subset, fa="fa")

Nfacs <- 3  # This is for four factors. You can change this as needed.

fit <- factanal(topology_scaled_subset, Nfacs, rotation = "promax")

ggplot() +
  geom_point(aes(x = fit$loadings[, 1],
                 y = fit$loadings[, 2])) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_text(aes(x = fit$loadings[, 1],
                y = fit$loadings[, 2],
                label = names(fit$loadings[, 1])),
            nudge_y = 0.1) +
  theme_classic() +
  labs(x = "Factor 1",
       y = "Factor 2") +
  lims(x= c(-1.2, 1.2), 
       y = c(-1.2, 1.2))

# 2: look at correlation first

pca <- PCA(topology_vermaat_subset)

# find number of dims describing 90 % of variance
dim_num <- pca$eig %>%
  as.data.frame() %>%
  filter(`cumulative percentage of variance` < 90) %>%
  nrow() %>% as.numeric()

dim_descrip <- dimdesc(pca, axes = 1:5)

# assign broader categories

grps <- 
  tibble(var = row.names(pca[["var"]][["coord"]])) %>%
  reframe(cat = case_when(var %in% c("S1", "S2", "S5", "S4") ~ "Motif",
                          var %in% c("richness", "links", "connectance", "complexity", "l_S", "ρ") ~ "Structure",
                          var %in% c("diameter", "TL", "path", "ChLen", "ChNum", "top", "basal", "distance") ~ "Energy",
                          var %in% c("GenSD", "VulSD", "omnivory", "herbivory", "intermediate", "centrality") ~ "Node"))

variance = pca$eig %>%
  as.data.frame() %>%
  pull(`percentage of variance`)

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
  
