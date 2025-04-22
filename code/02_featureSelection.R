##Load packages
library(caret)
library(corrplot)
library(FactoMineR)
library(here)
library(psych)
library(tidyverse)

setwd(here::here())

# import network summary data
topology_mangal <- read.csv("data/mangal/mangal_summary.csv") %>%
  select(-id)

# get an idea of the number of NA vals
sapply(topology_mangal, function(x) sum(is.na(x))) / nrow(topology_mangal) * 100

# same for prop of zeros
colSums(topology_mangal==0, na.rm = TRUE)/nrow(topology_mangal)*100

topology_mangal_subset <-
  topology_mangal %>%
  select(-c(ChSD, cannibal)) %>%
  na.omit()

# 1: look at correlation first

#scale all features
topology_scaled <- scale(topology_mangal_subset,
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
topology_scaled_subset <- topology_scaled[,-hc] %>%
  as.data.frame()

cor_mat_subset <- cor(topology_scaled_subset)

corrplot(cor_mat_subset, order = "hclust")

#KMO Test

KMO(topology_scaled_subset)
cortest.bartlett(topology_scaled_subset)

# Determine number of factors to extract
ev <- eigen(cor(topology_scaled_subset)) # get eigenvalues
ev$values

scree(topology_scaled_subset, pc=FALSE)

fa.parallel(topology_scaled_subset, fa="fa")

Nfacs <- 4  # This is for four factors. You can change this as needed.

fit <- factanal(topology_scaled_subset, Nfacs, rotation="promax")

# 2: look at correlation first

pca <- PCA(topology_scaled_subset)

# find number of dims describing 90 % of variance
dim_num <- pca$eig %>%
  as.data.frame() %>%
  filter(`cumulative percentage of variance` < 90) %>%
  nrow() %>% as.numeric()

dim_descrip <- dimdesc(pca, axes = 1:5)

dim_descrip$Dim.1 %>%
  as.data.frame() %>%
  filter(abs(quanti.correlation) > 0.75) %>%
  rbind(dim_descrip$Dim.2 %>%
          as.data.frame() %>%
          filter(abs(quanti.correlation) > 0.75)) %>%
  rbind(dim_descrip$Dim.3 %>%
          as.data.frame() %>%
          filter(abs(quanti.correlation) > 0.75)) %>%
  rbind(dim_descrip$Dim.4 %>%
          as.data.frame() %>%
          filter(abs(quanti.correlation) > 0.75)) %>%
  rbind(dim_descrip$Dim.5 %>%
          as.data.frame() %>%
          filter(abs(quanti.correlation) > 0.75))



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
  
