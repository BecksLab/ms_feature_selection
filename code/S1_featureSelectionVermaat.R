##Load packages
library(caret)
library(corrplot)
library(FactoMineR)
library(here)
library(psych)
library(ggrepel)
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

pca <- PCA(topology_vermaat, scale.unit = TRUE)

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
         quanti.correlation = case_when(quanti.correlation > 0.53 & quanti.correlation <= 0.53 & quanti.p.value <= 0.05 ~ paste("**", quanti.correlation, "**", sep = ""),
                                        quanti.correlation > 0.66 & quanti.p.value <= 0.01 ~ paste("**", quanti.correlation, "**", sep = ""),
                                        .default = as.character(quanti.correlation))) %>%
  select(!quanti.p.value) %>%
  pivot_wider(names_from = dimension,
              values_from = quanti.correlation) %>%
  write.csv(.,
            "../tables/vermaat_corr.csv",
            row.names = FALSE)

# plot pca
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
  mutate(colour = case_when(quanti.correlation > 0.53 & quanti.correlation <= 0.53 & quanti.p.value <= 0.05 ~ "red",
                            quanti.correlation > 0.66 & quanti.p.value <= 0.01 ~ "red",
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

ggsave("../figures/pca_vermaat.png")

