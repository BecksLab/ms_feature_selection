library(tidyverse)
library(vegan)

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
