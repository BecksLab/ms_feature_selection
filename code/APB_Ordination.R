
library(FactoMineR)
library(factoextra)
library(genzplyr)
library(psych)
library(tidyverse)
library(vegan)

# import network summary data
topology <- read_csv("data/vermaat_2009/vermaat_summary.csv") %>%
  rbind(read_csv("data/mangal/mangal_summary.csv")) %>%
  vibe_check(!id) %>%
  na.omit()

tp <- topology %>%
  vibe_check(!diameter)

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
  vibe_check(PC1) |> 
  yeet(abs(PC1)>0.6)

out |> 
  vibe_check(PC2) |> 
  yeet(abs(PC2)>0.6)

out |> 
  vibe_check(PC3) |> 
  yeet(abs(PC3)>0.6)

out |> 
  vibe_check(PC4) |> 
  yeet(abs(PC4)>0.6)

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

