library(tidyverse)
library(vegan)
library(psych)
library(FactoMineR)
library(factoextra)

# import network summary data
topology <- read_csv("data/vermaat_2009/vermaat_summary.csv") %>%
  rbind(read_csv("data/mangal/mangal_summary.csv")) %>%
  select(!id) %>%
  na.omit()

tp <- topology |> 
  select(-herbivory)

# use vegan rda
pp <- rda(tp, scale = TRUE, center = TRUE)
biplot(pp)
summary(pp)

# grab correlations
out <- scores(pp, choices = 1:4, 
              display = "species", 
              scaling = "species", correlation = TRUE) |> 
  data.frame()

# examine major axes
out |> 
  select(PC1) |> 
  filter(abs(PC1)>0.6)

out |> 
  select(PC2) |> 
  filter(abs(PC2)>0.6)

out |> 
  select(PC3) |> 
  filter(abs(PC3)>0.6)

out |> 
  select(PC4) |> 
  filter(abs(PC4)>0.6)

##-- factor analysis attemp
tp_scaled <- scale(tp)
fa.parallel(tp_scaled, # perform parallel analysis
            fa = "fa",
            fm = "pa",
            show.legend = TRUE,
            main = "Scree Plot and Parallel Analysis")

efa_promax3 <- fa(tp_scaled, # perform EFA with 3 factors
                  nfactors=3,
                  rotate="promax",
                  fm="pa")
fa.diagram(efa_promax3)

# --factominer

