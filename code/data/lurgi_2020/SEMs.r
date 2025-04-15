## ---------------------------
##
## Script name: SEMs.r
##
## Purpose of script: Run piecewise structural equation models (SEMs) on data describing
## geographical variability of network structure in marine intertidal communities
## against environmental predictor variables.
##
## Author: Dr Miguel Lurgi
## Lecturer in Biosciences (Computational Ecology)
## Computational Ecology Lab - Department of Biosciences
## Swansea University, UK
## 
## and
##
## Centre for Biodiversity Theory and Modelling
## Theoretical and Experimental Ecology Station, CNRS, France
##
## Date Created: 19-12-2019
##
## Copyright (c) Miguel Lurgi, 2019
## Email: miguel.lurgi@swansea.ac.uk
##
## ---------------------------
##
## Notes:
##
## This script is provided as supplementary material for the paper:
## Lurgi et al. (2020) Geographical variation of multiplex ecological networks 
## in marine intertidal communities, Ecology.
##
## ---------------------------

############## This code is for the piecewise structural equation models ################
## To keep things tidy I moved all the SEMs to here so the code is more clear

## Load required libraries
require(vegan)
require(piecewiseSEM)
require(MASS)
require(pastecs)
require(AER)

## After discussions with one of the reviewers (Jarrett Byrnes),
## we decided to incorporate all of the predictor variables into the SEMs.
## Before, we were only using those coming out of the RDA as significant (see manuscript for details)

## The environmental_residuals data frame contains the information of the residuals of the environmental
## variables after regressing them against long term mean SST (see methods in the manuscript)

## The following analyses require 4 data frames that contain the values for the structural properties
## of each of the 3 types of interaction networks plus the residuals of the environmental predictors
## as described in the methods of the paper.

## These data can be obtained by running other scripts that are part of the source code (see README)
## Alternatively, for convenience, they have been stored in the R data file called: data-for-chilean-analysis.rda

## First, we load the data from that file
load('data-for-chilean-analysis.rda')

## These 3 lines make sure that all the data frames are in the right order with respect to each other.
## Make sure you always get oredered sequences from 1 to 17

match(row.names(output_trophic), row.names(output_neg_nti))
match(row.names(output_neg_nti), row.names(output_pos_nti))
match(row.names(output_neg_nti), row.names(environment_residuals))

## we 'standardize' the environmental variables to bring them all to the same scale
X_norm <- decostand(environment_residuals, 'standardize')

############### BEGINNING OF ANALYSES FOR THE NON-TROPHIC POSITIVE INTERACTIONS NETWORKS ###############

## we start the analyses with the data for the non-trophic positive interactions
## values of the network properties for this type of interaction are kept in the data frame
## ouptut_pos_nti_stats, which is loaded when loading the workspace (see comments above)
Y_norm <- output_pos_nti[,4:19]

## To run the models we need both, environmental predictors and network descriptors in the same
## data frame
data_for_sems <- cbind(X_norm, Y_norm)

## this is used later for normalisation
data_sems_sd <- apply(data_for_sems, 2, sd)

## the list of candidate models to test the dependece of each network descriptor to environmental predictors, 
## and species richness (S) is this:

## 1. glm(S ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST, na.action=na.omit, data=data_for_sems, family=quasipoisson)
## 2. glm.nb(L ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems)
## 3. lm(L.S ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S + L, na.action=na.omit, data=data_for_sems)
## 4. lm(C ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S + L.S, na.action=na.omit, data=data_for_sems)
## 5. lm(indeg ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems)
## 6. lm(outdeg ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems)
## 7. lm(sd_gen ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems)
## 8. lm(sd_vul ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems)
## 9. lm(modularity ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems)
## 10. lm(fr_mutualism ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems)
## 11. lm(fr_comensalism ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems)

## After discussions with Jarrett Byrnes in early December 2019, and given his experience in dealing with multiple responses
## variables that might covary, preventing thus a good fit of the SEMs, we decided to perform one separate SEM for each response
## variable / network descriptor

## Thus the first SEM consisted of models 1 and 2 from the list above.

## Before constructing and fitting the actual SEMs we verify that the models used are appropriate:
## Since both S and L are count data they were both modelled using a Generalised linear model (GLM) with the appropriate family

## As indicated by the summary statistics over the S
pastecs::stat.desc(data_for_sems$S)
## this variable suffers from underdispersion, i.e. the variance is considerable smaller than the mean: 13.39 vs. 44.53
## for this reason we fit these data to a quasipoisson distribution, as a way to deal with underdispersion

## Once we fit the model:
model_test <- glm(S ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST, na.action=na.omit, data=data_for_sems, family=quasipoisson(link = 'identity'))
## by looking at the p-value resulting from a chi-square test of model deviance
## we conclude this model is a very good fit to the data, i.e. p-value >> 0.05
1-pchisq(summary(model_test)$deviance, summary(model_test)$df.residual)

## Similarly for L, since we are dealing with count data, we fit a GLM. This time, however, since the data is overdispersed,
## as indicated by the summary statistics:
pastecs::stat.desc(data_for_sems$L)
## var = 168.26 vs. mean = 107.41
## we then use a negative binomial distribution for the data fit

## Once we fit the model:
model_test <- glm.nb(L ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems, control=glm.control(epsilon = 1e-8, maxit = 250, trace = FALSE))

## by looking at the p-value resulting from a chi-square test of model deviance
## we conclude this model is a very good fit to the data, i.e. p-value >> 0.05
1-pchisq(summary(model_test)$deviance, summary(model_test)$df.residual)

## Once single models have been verified we construct and evaluate the SEM
sem_nti_pos <- psem(
  glm(S ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST, na.action=na.omit, data=data_for_sems, family=quasipoisson(link = 'identity')),
  glm.nb(L ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems, control=glm.control(epsilon = 1e-8, maxit = 250, trace = FALSE)),
  data=data_for_sems)

sem_nti_pos <- summary(sem_nti_pos)

## Then we calculate the standardised predictors estimates using the formula: estimate * (sd.x / sd.y)
## since the predictors were standardised (above), we only need to do: estimate / sd.y
for(res in unique(sem_nti_pos$coefficients$Response)){
  sd.y <- data_sems_sd[which(names(data_sems_sd) == res)]
  sem_nti_pos$coefficients[which(sem_nti_pos$coefficients$Response == res),]$Std.Estimate <- round(sem_nti_pos$coefficients[which(sem_nti_pos$coefficients$Response == res),]$Estimate / sd.y, 4)
}

coeffs_table <- sem_nti_pos$coefficients
r_sq <- sem_nti_pos$R2

## Similarly, the remaining SEMs were constructed from models 1 and 3, 1 and 4, 1 and 5... and so on

## Since the remaining variables are all continous we fit linear regressions to these data
## To ensure these models comply with the assumptions of linear regression 
## we validate these using the gvlma function from the gvlma package

## For connectance C
model_test <- lm(C ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems)
gvlma::gvlma(model_test)

sem_nti_pos <- psem(
  glm(S ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST, na.action=na.omit, data=data_for_sems, family=quasipoisson(link = 'identity')),
  lm(C ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems),
  data=data_for_sems)

sem_nti_pos <- summary(sem_nti_pos)

## Then we calculate the standardised predictors estimates using the formula: estimate * (sd.x / sd.y)
## since the predictors were standardised (above), we only need to do: estimate / sd.y
for(res in unique(sem_nti_pos$coefficients$Response)){
  sd.y <- data_sems_sd[which(names(data_sems_sd) == res)]
  sem_nti_pos$coefficients[which(sem_nti_pos$coefficients$Response == res),]$Std.Estimate <- round(sem_nti_pos$coefficients[which(sem_nti_pos$coefficients$Response == res),]$Estimate / sd.y, 4)
}

coeffs_table <- rbind(coeffs_table, sem_nti_pos$coefficients[7:13,])
r_sq <- rbind(r_sq, sem_nti_pos$R2[2,])

## for links per species - L/S
## since the response variable (L/S) is not normally distributed (but skewed) we fit a GLM with Gamma family
model_test <- glm(L.S ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems, family=Gamma)
1-pchisq(summary(model_test)$deviance, summary(model_test)$df.residual)

sem_nti_pos <- psem(
  glm(S ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST, na.action=na.omit, data=data_for_sems, family=quasipoisson(link = 'identity')),
  glm(L.S ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems, family=Gamma))

sem_nti_pos <- summary(sem_nti_pos)

## Then we calculate the standardised predictors estimates using the formula: estimate * (sd.x / sd.y)
## since the predictors were standardised (above), we only need to do: estimate / sd.y
for(res in unique(sem_nti_pos$coefficients$Response)){
  sd.y <- data_sems_sd[which(names(data_sems_sd) == res)]
  sem_nti_pos$coefficients[which(sem_nti_pos$coefficients$Response == res),]$Std.Estimate <- round(sem_nti_pos$coefficients[which(sem_nti_pos$coefficients$Response == res),]$Estimate / sd.y, 4)
}

coeffs_table <- rbind(coeffs_table, sem_nti_pos$coefficients[7:13,])
r_sq <- rbind(r_sq, sem_nti_pos$R2[2,])

## for indegree we also use a Gamma family for the response
model_test <- glm(indeg ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems, family=Gamma)
1-pchisq(summary(model_test)$deviance, summary(model_test)$df.residual)

sem_nti_pos <- psem(
  glm(S ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST, na.action=na.omit, data=data_for_sems, family=quasipoisson(link = 'identity')),
  glm(indeg ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems, family=Gamma))

sem_nti_pos <- summary(sem_nti_pos)

## Then we calculate the standardised predictors estimates using the formula: estimate * (sd.x / sd.y)
## since the predictors were standardised (above), we only need to do: estimate / sd.y
for(res in unique(sem_nti_pos$coefficients$Response)){
  sd.y <- data_sems_sd[which(names(data_sems_sd) == res)]
  sem_nti_pos$coefficients[which(sem_nti_pos$coefficients$Response == res),]$Std.Estimate <- round(sem_nti_pos$coefficients[which(sem_nti_pos$coefficients$Response == res),]$Estimate / sd.y, 4)
}

coeffs_table <- rbind(coeffs_table, sem_nti_pos$coefficients[7:13,])
r_sq <- rbind(r_sq, sem_nti_pos$R2[2,])

## for outdegree we also use a Gamma family for the response
model_test <- glm(outdeg ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems, family=Gamma(link = 'inverse'))
1-pchisq(summary(model_test)$deviance, summary(model_test)$df.residual)

sem_nti_pos <- psem(
  glm(S ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST, na.action=na.omit, data=data_for_sems, family=quasipoisson(link = 'identity')),
  glm(outdeg ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems, family=Gamma))

sem_nti_pos <- summary(sem_nti_pos)

## Then we calculate the standardised predictors estimates using the formula: estimate * (sd.x / sd.y)
## since the predictors were standardised (above), we only need to do: estimate / sd.y
for(res in unique(sem_nti_pos$coefficients$Response)){
  sd.y <- data_sems_sd[which(names(data_sems_sd) == res)]
  sem_nti_pos$coefficients[which(sem_nti_pos$coefficients$Response == res),]$Std.Estimate <- round(sem_nti_pos$coefficients[which(sem_nti_pos$coefficients$Response == res),]$Estimate / sd.y, 4)
}

coeffs_table <- rbind(coeffs_table, sem_nti_pos$coefficients[7:13,])
r_sq <- rbind(r_sq, sem_nti_pos$R2[2,])


## sd_gen
model_test <- lm(sd_gen ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems)
gvlma::gvlma(model_test)

## Because of this, and to be sure about heteroskedasticity, we applied the Breusch-Pagan test, a test commonly used to test against heteroskedasticity
bptest(model_test)
## the Breusch-Pagan test reveals that the potential heteroskedacity in the model detected by the gvlma is not really significant: i.e. p-value > 0.05.
## Because of this we decided to keep the model as it is

sem_nti_pos <- psem(
  glm(S ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST, na.action=na.omit, data=data_for_sems, family=quasipoisson(link = 'identity')),
  lm(sd_gen ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems),
  data=data_for_sems)

sem_nti_pos <- summary(sem_nti_pos)

## Then we calculate the standardised predictors estimates using the formula: estimate * (sd.x / sd.y)
## since the predictors were standardised (above), we only need to do: estimate / sd.y
for(res in unique(sem_nti_pos$coefficients$Response)){
  sd.y <- data_sems_sd[which(names(data_sems_sd) == res)]
  sem_nti_pos$coefficients[which(sem_nti_pos$coefficients$Response == res),]$Std.Estimate <- round(sem_nti_pos$coefficients[which(sem_nti_pos$coefficients$Response == res),]$Estimate / sd.y, 4)
}

coeffs_table <- rbind(coeffs_table, sem_nti_pos$coefficients[7:13,])
r_sq <- rbind(r_sq, sem_nti_pos$R2[2,])

## sd_vul
model_test <- lm(sd_vul ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems)
gvlma::gvlma(model_test)

sem_nti_pos <- psem(
  glm(S ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST, na.action=na.omit, data=data_for_sems, family=quasipoisson(link = 'identity')),
  lm(sd_vul ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems))

sem_nti_pos <- summary(sem_nti_pos)

## Then we calculate the standardised predictors estimates using the formula: estimate * (sd.x / sd.y)
## since the predictors were standardised (above), we only need to do: estimate / sd.y
for(res in unique(sem_nti_pos$coefficients$Response)){
  sd.y <- data_sems_sd[which(names(data_sems_sd) == res)]
  sem_nti_pos$coefficients[which(sem_nti_pos$coefficients$Response == res),]$Std.Estimate <- round(sem_nti_pos$coefficients[which(sem_nti_pos$coefficients$Response == res),]$Estimate / sd.y, 4)
}

coeffs_table <- rbind(coeffs_table, sem_nti_pos$coefficients[7:13,])
r_sq <- rbind(r_sq, sem_nti_pos$R2[2,])

## modularity
model_test <- lm(modularity ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems)
gvlma::gvlma(model_test)

sem_nti_pos <- psem(
  glm(S ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST, na.action=na.omit, data=data_for_sems, family=quasipoisson(link = 'identity')),
  lm(modularity ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems))

sem_nti_pos <- summary(sem_nti_pos)

## Then we calculate the standardised predictors estimates using the formula: estimate * (sd.x / sd.y)
## since the predictors were standardised (above), we only need to do: estimate / sd.y
for(res in unique(sem_nti_pos$coefficients$Response)){
  sd.y <- data_sems_sd[which(names(data_sems_sd) == res)]
  sem_nti_pos$coefficients[which(sem_nti_pos$coefficients$Response == res),]$Std.Estimate <- round(sem_nti_pos$coefficients[which(sem_nti_pos$coefficients$Response == res),]$Estimate / sd.y, 4)
}

coeffs_table <- rbind(coeffs_table, sem_nti_pos$coefficients[7:13,])
r_sq <- rbind(r_sq, sem_nti_pos$R2[2,])

## fraction of mutualism
model_test <- lm(fr_mutualism ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems)
gvlma::gvlma(model_test)

sem_nti_pos <- psem(
  glm(S ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST, na.action=na.omit, data=data_for_sems, family=quasipoisson(link = 'identity')),
  lm(fr_mutualism ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems))

sem_nti_pos <- summary(sem_nti_pos)

## Then we calculate the standardised predictors estimates using the formula: estimate * (sd.x / sd.y)
## since the predictors were standardised (above), we only need to do: estimate / sd.y
for(res in unique(sem_nti_pos$coefficients$Response)){
  sd.y <- data_sems_sd[which(names(data_sems_sd) == res)]
  sem_nti_pos$coefficients[which(sem_nti_pos$coefficients$Response == res),]$Std.Estimate <- round(sem_nti_pos$coefficients[which(sem_nti_pos$coefficients$Response == res),]$Estimate / sd.y, 4)
}

coeffs_table <- rbind(coeffs_table, sem_nti_pos$coefficients[7:13,])
r_sq <- rbind(r_sq, sem_nti_pos$R2[2,])


## fraction of commensalism
model_test <- lm(fr_comensalism ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems)
gvlma::gvlma(model_test)

sem_nti_pos <- psem(
  glm(S ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST, na.action=na.omit, data=data_for_sems, family=quasipoisson(link = 'identity')),
  lm(fr_comensalism ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems))

sem_nti_pos <- summary(sem_nti_pos)

## Then we calculate the standardised predictors estimates using the formula: estimate * (sd.x / sd.y)
## since the predictors were standardised (above), we only need to do: estimate / sd.y
for(res in unique(sem_nti_pos$coefficients$Response)){
  sd.y <- data_sems_sd[which(names(data_sems_sd) == res)]
  sem_nti_pos$coefficients[which(sem_nti_pos$coefficients$Response == res),]$Std.Estimate <- round(sem_nti_pos$coefficients[which(sem_nti_pos$coefficients$Response == res),]$Estimate / sd.y, 4)
}

coeffs_table <- rbind(coeffs_table, sem_nti_pos$coefficients[7:13,])
r_sq <- rbind(r_sq, sem_nti_pos$R2[2,])


coeffs_table$Std.Estimate <- as.numeric(coeffs_table$Std.Estimate)
sd.x <- data_sems_sd[which(names(data_sems_sd) == 'S')]
coeffs_table[which(coeffs_table$Predictor == 'S'),]$Std.Estimate <- round(coeffs_table[which(coeffs_table$Predictor == 'S'),]$Std.Estimate * sd.x, 4)


coeffs_table$ranged_coeffs <- round((coeffs_table[,3] - min(coeffs_table[,3])) / (max(coeffs_table[,3]) - min(coeffs_table[,3])), 4)

a <- abs(coeffs_table[,3])
ranged_abs <- (a - min(a)) / (max(a) - min(a))
ranged_abs[which(coeffs_table[,3] < 0)] <- -ranged_abs[which(coeffs_table[,3] < 0)]

coeffs_table$ranged_abs <- round(ranged_abs, 4)

write.csv(coeffs_table, file = 'sem-output-positive-nti-networks.csv')
write.csv(r_sq, file = 'sem-r-squareds-positive-nti-networks.csv')

############### END OF THE ANALYSES FOR THE POSITIVE NON-TROPHIC NETWORKS ###############

############### BEGINNING OF ANALYSES FOR THE TROPHIC NETWORKS (I.E. FOOD WEBS) ###############

## now, we turn our attention to the trophic networks
## values of the network properties for this type of interaction are kept in the data frame
## ouptut_trophic_stats, which is loaded when loading the workspace (see comments above)
Y_norm <- output_trophic[,4:23]

## To run the models we need both, environmental predictors and network descriptors in the same data frame
data_for_sems <- cbind(X_norm, Y_norm)

## this is used later for normalisation
data_sems_sd <- apply(data_for_sems, 2, sd)

## the list of candidate models to test the dependece of each network descriptor to environmental predictors, 
## and species richness (S) is this:

## 1. glm(S ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST, na.action=na.omit, data=data_for_sems, family=quasipoisson(link = 'identity'))
## 2. glm.nb(L ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems)
## 3. lm(L.S ~ LT.SST + Upwelling + Upwelling.SD + Fr.Annual + Clim + Fr.Days + S, na.action=na.omit, data=data_for_sems)
## 4. lm(C ~ LT.SST + Upwelling + Upwelling.SD + Fr.Annual + Clim + Fr.Days + S, na.action=na.omit, data=data_for_sems)
## 5. lm(indeg ~ LT.SST + Upwelling + Upwelling.SD + Fr.Annual + Clim + Fr.Days + S, na.action=na.omit, data=data_for_sems)
## 6. lm(outdeg ~ LT.SST + Upwelling + Upwelling.SD + Fr.Annual + Clim + Fr.Days + S, na.action=na.omit, data=data_for_sems)
## 7. lm(sd_gen ~ LT.SST + Upwelling + Upwelling.SD + Fr.Annual + Clim + Fr.Days + S, na.action=na.omit, data=data_for_sems)
## 8. lm(sd_vul ~ LT.SST + Upwelling + Upwelling.SD + Fr.Annual + Clim + Fr.Days + S, na.action=na.omit, data=data_for_sems)
## 9. lm(basal ~ LT.SST + Upwelling + Upwelling.SD + Fr.Annual + Clim + Fr.Days + S, na.action=na.omit, data=data_for_sems)
## 10. lm(top ~ LT.SST + Upwelling + Upwelling.SD + Fr.Annual + Clim + Fr.Days + S, na.action=na.omit, data=data_for_sems)
## 11. lm(interm ~ LT.SST + Upwelling + Upwelling.SD + Fr.Annual + Clim + Fr.Days + S, na.action=na.omit, data=data_for_sems)
## 12. lm(modularity ~ LT.SST + Upwelling + Upwelling.SD + Fr.Annual + Clim + Fr.Days + S, na.action=na.omit, data=data_for_sems)
## 13. lm(omniv ~ LT.SST + Upwelling + Upwelling.SD + Fr.Annual + Clim + Fr.Days + S, na.action=na.omit, data=data_for_sems)
## 14. lm(mfcl ~ LT.SST + Upwelling + Upwelling.SD + Fr.Annual + Clim + Fr.Days + S, na.action=na.omit, data=data_for_sems)

## Following the protocol for the analyses adopted for the non-trophic positive interaction networks (above),
## we start with models 1 and 2 with the list above:

## Before constructing and fitting the actual SEMs we verify that the models used are appropriate:
## Since both S and L are count data they were both modelled using a Generalised linear model (GLM) with the appropriate family

## In this case, the summary statistics over the S
pastecs::stat.desc(data_for_sems$S)
## reveal that S is again underdispersed, so we use quasipoisson family

## Once we fit the model:
model_test <- glm(S ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST, na.action=na.omit, data=data_for_sems, family=quasipoisson(link = 'identity'))
## by looking at the p-value resulting from a chi-square test of model deviance
## we conclude this model is a very good fit to the data, i.e. p-value >> 0.05
1-pchisq(summary(model_test)$deviance, summary(model_test)$df.residual)

###### Tests for over / underdispersion, to make sure the model we use is the right one
fit = glm(L ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST, na.action=na.omit, data=data_for_sems, family="poisson") 
fit.overdisp = glm(L ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST, na.action=na.omit, data=data_for_sems, family="quasipoisson") 
summary(fit.overdisp)$dispersion # dispersion coefficient
pchisq(summary(fit.overdisp)$dispersion * fit$df.residual, fit$df.residual, lower = F) # significance of overdispersion

dispersiontest(fit, trafo=1, alternative = 'greater')

## Similarly for L, since we are dealing with count data, we fit a GLM. This time, however, since the data is overdispersed,
## as indicated by the summary statistics:
pastecs::stat.desc(data_for_sems$L)
## var = 11426.49 vs. mean = 520.88
## we then use a negative binomial distribution for the data fit

## Once we fit the model:
model_test <- glm.nb(L ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems, 
                     control=glm.control(epsilon = 1e-14, maxit = 2000, trace = FALSE))

## by looking at the p-value resulting from a chi-square test of model deviance
## we conclude this model is a very good fit to the data, i.e. p-value >> 0.05
1-pchisq(summary(model_test)$deviance, summary(model_test)$df.residual)

## Once single models have been verified we construct and evaluate the SEM
sem_trophic <- psem(
  glm(S ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST, na.action=na.omit, data=data_for_sems, family=quasipoisson(link = 'identity')),
  glm.nb(L ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems, control=glm.control(epsilon = 1e-14, maxit = 250, trace = FALSE)),
  data=data_for_sems)

sem_trophic <- summary(sem_trophic)

## Then we calculate the standardised predictors estimates using the formula: estimate * (sd.x / sd.y)
## since the predictors were standardised (above), we only need to do: estimate / sd.y
for(res in unique(sem_trophic$coefficients$Response)){
  sd.y <- data_sems_sd[which(names(data_sems_sd) == res)]
  sem_trophic$coefficients[which(sem_trophic$coefficients$Response == res),]$Std.Estimate <- round(sem_trophic$coefficients[which(sem_trophic$coefficients$Response == res),]$Estimate / sd.y, 4)
}

coeffs_table <- sem_trophic$coefficients
r_sq <- sem_trophic$R2

## Similarly, the remaining SEMs were constructed from models 1 and 3, 1 and 4, 1 and 5... and so on

## Since the remaining variables are all continous we fit linear regressions to these data
## To ensure these models comply with the assumptions of linear regression 
## we validate these using the gvlma function from the gvlma package

## For connectance C
model_test <- lm(C ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems)
gvlma::gvlma(model_test)

sem_trophic <- psem(
  glm(S ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST, na.action=na.omit, data=data_for_sems, family=quasipoisson(link = 'identity')),
  lm(C ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems),
  data=data_for_sems)

sem_trophic <- summary(sem_trophic)

## Then we calculate the standardised predictors estimates using the formula: estimate * (sd.x / sd.y)
## since the predictors were standardised (above), we only need to do: estimate / sd.y
for(res in unique(sem_trophic$coefficients$Response)){
  sd.y <- data_sems_sd[which(names(data_sems_sd) == res)]
  sem_trophic$coefficients[which(sem_trophic$coefficients$Response == res),]$Std.Estimate <- round(sem_trophic$coefficients[which(sem_trophic$coefficients$Response == res),]$Estimate / sd.y, 4)
}

coeffs_table <- rbind(coeffs_table, sem_trophic$coefficients[7:13,])
r_sq <- rbind(r_sq, sem_trophic$R2[2,])

## for links per species - L/S
model_test <- lm(L.S ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems)
gvlma::gvlma(model_test)

sem_trophic <- psem(
  glm(S ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST, na.action=na.omit, data=data_for_sems, family=quasipoisson(link = 'identity')),
  lm(L.S ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems))

sem_trophic <- summary(sem_trophic)

## Then we calculate the standardised predictors estimates using the formula: estimate * (sd.x / sd.y)
## since the predictors were standardised (above), we only need to do: estimate / sd.y
for(res in unique(sem_trophic$coefficients$Response)){
  sd.y <- data_sems_sd[which(names(data_sems_sd) == res)]
  sem_trophic$coefficients[which(sem_trophic$coefficients$Response == res),]$Std.Estimate <- round(sem_trophic$coefficients[which(sem_trophic$coefficients$Response == res),]$Estimate / sd.y, 4)
}

coeffs_table <- rbind(coeffs_table, sem_trophic$coefficients[7:13,])
r_sq <- rbind(r_sq, sem_trophic$R2[2,])

## for indegree
model_test <- lm(indeg ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems)
gvlma::gvlma(model_test)

sem_trophic <- psem(
  glm(S ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST, na.action=na.omit, data=data_for_sems, family=quasipoisson(link = 'identity')),
  lm(indeg ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems))

sem_trophic <- summary(sem_trophic)

## Then we calculate the standardised predictors estimates using the formula: estimate * (sd.x / sd.y)
## since the predictors were standardised (above), we only need to do: estimate / sd.y
for(res in unique(sem_trophic$coefficients$Response)){
  sd.y <- data_sems_sd[which(names(data_sems_sd) == res)]
  sem_trophic$coefficients[which(sem_trophic$coefficients$Response == res),]$Std.Estimate <- round(sem_trophic$coefficients[which(sem_trophic$coefficients$Response == res),]$Estimate / sd.y, 4)
}

coeffs_table <- rbind(coeffs_table, sem_trophic$coefficients[7:13,])
r_sq <- rbind(r_sq, sem_trophic$R2[2,])

## for outdegree
model_test <- lm(outdeg ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems)
gvlma::gvlma(model_test)

sem_trophic <- psem(
  glm(S ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST, na.action=na.omit, data=data_for_sems, family=quasipoisson(link = 'identity')),
  lm(outdeg ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems))

sem_trophic <- summary(sem_trophic)

## Then we calculate the standardised predictors estimates using the formula: estimate * (sd.x / sd.y)
## since the predictors were standardised (above), we only need to do: estimate / sd.y
for(res in unique(sem_trophic$coefficients$Response)){
  sd.y <- data_sems_sd[which(names(data_sems_sd) == res)]
  sem_trophic$coefficients[which(sem_trophic$coefficients$Response == res),]$Std.Estimate <- round(sem_trophic$coefficients[which(sem_trophic$coefficients$Response == res),]$Estimate / sd.y, 4)
}

coeffs_table <- rbind(coeffs_table, sem_trophic$coefficients[7:13,])
r_sq <- rbind(r_sq, sem_trophic$R2[2,])

## sd_gen
model_test <- lm(sd_gen ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems)
gvlma::gvlma(model_test)

sem_trophic <- psem(
  glm(S ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST, na.action=na.omit, data=data_for_sems, family=quasipoisson(link = 'identity')),
  lm(sd_gen ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems))

sem_trophic <- summary(sem_trophic)

## Then we calculate the standardised predictors estimates using the formula: estimate * (sd.x / sd.y)
## since the predictors were standardised (above), we only need to do: estimate / sd.y
for(res in unique(sem_trophic$coefficients$Response)){
  sd.y <- data_sems_sd[which(names(data_sems_sd) == res)]
  sem_trophic$coefficients[which(sem_trophic$coefficients$Response == res),]$Std.Estimate <- round(sem_trophic$coefficients[which(sem_trophic$coefficients$Response == res),]$Estimate / sd.y, 4)
}

coeffs_table <- rbind(coeffs_table, sem_trophic$coefficients[7:13,])
r_sq <- rbind(r_sq, sem_trophic$R2[2,])

## sd_vul
model_test <- lm(sd_vul ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems)
gvlma::gvlma(model_test)

sem_trophic <- psem(
  glm(S ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST, na.action=na.omit, data=data_for_sems, family=quasipoisson(link = 'identity')),
  lm(sd_vul ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems))

sem_trophic <- summary(sem_trophic)

## Then we calculate the standardised predictors estimates using the formula: estimate * (sd.x / sd.y)
## since the predictors were standardised (above), we only need to do: estimate / sd.y
for(res in unique(sem_trophic$coefficients$Response)){
  sd.y <- data_sems_sd[which(names(data_sems_sd) == res)]
  sem_trophic$coefficients[which(sem_trophic$coefficients$Response == res),]$Std.Estimate <- round(sem_trophic$coefficients[which(sem_trophic$coefficients$Response == res),]$Estimate / sd.y, 4)
}

coeffs_table <- rbind(coeffs_table, sem_trophic$coefficients[7:13,])
r_sq <- rbind(r_sq, sem_trophic$R2[2,])

## basal
model_test <- lm(basal ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems)
gvlma::gvlma(model_test)

sem_trophic <- psem(
  glm(S ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST, na.action=na.omit, data=data_for_sems, family=quasipoisson(link = 'identity')),
  lm(basal ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems))

sem_trophic <- summary(sem_trophic)

## Then we calculate the standardised predictors estimates using the formula: estimate * (sd.x / sd.y)
## since the predictors were standardised (above), we only need to do: estimate / sd.y
for(res in unique(sem_trophic$coefficients$Response)){
  sd.y <- data_sems_sd[which(names(data_sems_sd) == res)]
  sem_trophic$coefficients[which(sem_trophic$coefficients$Response == res),]$Std.Estimate <- round(sem_trophic$coefficients[which(sem_trophic$coefficients$Response == res),]$Estimate / sd.y, 4)
}

coeffs_table <- rbind(coeffs_table, sem_trophic$coefficients[7:13,])
r_sq <- rbind(r_sq, sem_trophic$R2[2,])

## top
model_test <- lm(top ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems)
gvlma::gvlma(model_test)

sem_trophic <- psem(
  glm(S ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST, na.action=na.omit, data=data_for_sems, family=quasipoisson(link = 'identity')),
  lm(top ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems))

sem_trophic <- summary(sem_trophic)

## Then we calculate the standardised predictors estimates using the formula: estimate * (sd.x / sd.y)
## since the predictors were standardised (above), we only need to do: estimate / sd.y
for(res in unique(sem_trophic$coefficients$Response)){
  sd.y <- data_sems_sd[which(names(data_sems_sd) == res)]
  sem_trophic$coefficients[which(sem_trophic$coefficients$Response == res),]$Std.Estimate <- round(sem_trophic$coefficients[which(sem_trophic$coefficients$Response == res),]$Estimate / sd.y, 4)
}

coeffs_table <- rbind(coeffs_table, sem_trophic$coefficients[7:13,])
r_sq <- rbind(r_sq, sem_trophic$R2[2,])

## intermediate
model_test <- lm(interm ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems)
gvlma::gvlma(model_test)

sem_trophic <- psem(
  glm(S ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST, na.action=na.omit, data=data_for_sems, family=quasipoisson(link = 'identity')),
  lm(interm ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems))

sem_trophic <- summary(sem_trophic)

## Then we calculate the standardised predictors estimates using the formula: estimate * (sd.x / sd.y)
## since the predictors were standardised (above), we only need to do: estimate / sd.y
for(res in unique(sem_trophic$coefficients$Response)){
  sd.y <- data_sems_sd[which(names(data_sems_sd) == res)]
  sem_trophic$coefficients[which(sem_trophic$coefficients$Response == res),]$Std.Estimate <- round(sem_trophic$coefficients[which(sem_trophic$coefficients$Response == res),]$Estimate / sd.y, 4)
}

coeffs_table <- rbind(coeffs_table, sem_trophic$coefficients[7:13,])
r_sq <- rbind(r_sq, sem_trophic$R2[2,])

## modularity
model_test <- lm(modularity ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems)
gvlma::gvlma(model_test)

sem_trophic <- psem(
  glm(S ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST, na.action=na.omit, data=data_for_sems, family=quasipoisson(link = 'identity')),
  lm(modularity ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems))

sem_trophic <- summary(sem_trophic)

## Then we calculate the standardised predictors estimates using the formula: estimate * (sd.x / sd.y)
## since the predictors were standardised (above), we only need to do: estimate / sd.y
for(res in unique(sem_trophic$coefficients$Response)){
  sd.y <- data_sems_sd[which(names(data_sems_sd) == res)]
  sem_trophic$coefficients[which(sem_trophic$coefficients$Response == res),]$Std.Estimate <- round(sem_trophic$coefficients[which(sem_trophic$coefficients$Response == res),]$Estimate / sd.y, 4)
}

coeffs_table <- rbind(coeffs_table, sem_trophic$coefficients[7:13,])
r_sq <- rbind(r_sq, sem_trophic$R2[2,])

## omnivory
model_test <- lm(omniv ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems)
gvlma::gvlma(model_test)

sem_trophic <- psem(
  glm(S ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST, na.action=na.omit, data=data_for_sems, family=quasipoisson(link = 'identity')),
  lm(omniv ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems))

sem_trophic <- summary(sem_trophic)

## Then we calculate the standardised predictors estimates using the formula: estimate * (sd.x / sd.y)
## since the predictors were standardised (above), we only need to do: estimate / sd.y
for(res in unique(sem_trophic$coefficients$Response)){
  sd.y <- data_sems_sd[which(names(data_sems_sd) == res)]
  sem_trophic$coefficients[which(sem_trophic$coefficients$Response == res),]$Std.Estimate <- round(sem_trophic$coefficients[which(sem_trophic$coefficients$Response == res),]$Estimate / sd.y, 4)
}

coeffs_table <- rbind(coeffs_table, sem_trophic$coefficients[7:13,])
r_sq <- rbind(r_sq, sem_trophic$R2[2,])

## mean food chain length
model_test <- lm(mfcl ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems)
gvlma::gvlma(model_test)

sem_trophic <- psem(
  glm(S ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST, na.action=na.omit, data=data_for_sems, family=quasipoisson(link = 'identity')),
  lm(mfcl ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems))

sem_trophic <- summary(sem_trophic)

## Then we calculate the standardised predictors estimates using the formula: estimate * (sd.x / sd.y)
## since the predictors were standardised (above), we only need to do: estimate / sd.y
for(res in unique(sem_trophic$coefficients$Response)){
  sd.y <- data_sems_sd[which(names(data_sems_sd) == res)]
  sem_trophic$coefficients[which(sem_trophic$coefficients$Response == res),]$Std.Estimate <- round(sem_trophic$coefficients[which(sem_trophic$coefficients$Response == res),]$Estimate / sd.y, 4)
}

coeffs_table <- rbind(coeffs_table, sem_trophic$coefficients[7:13,])
r_sq <- rbind(r_sq, sem_trophic$R2[2,])


coeffs_table$Std.Estimate <- as.numeric(coeffs_table$Std.Estimate)
sd.x <- data_sems_sd[which(names(data_sems_sd) == 'S')]
coeffs_table[which(coeffs_table$Predictor == 'S'),]$Std.Estimate <- round(coeffs_table[which(coeffs_table$Predictor == 'S'),]$Std.Estimate * sd.x, 4)


coeffs_table$ranged_coeffs <- round((coeffs_table[,3] - min(coeffs_table[,3])) / (max(coeffs_table[,3]) - min(coeffs_table[,3])),4)

a <- abs(coeffs_table[,3])
ranged_abs <- (a - min(a)) / (max(a) - min(a))
ranged_abs[which(coeffs_table[,3] < 0)] <- -ranged_abs[which(coeffs_table[,3] < 0)]

coeffs_table$ranged_abs <- round(ranged_abs,4)

write.csv(coeffs_table, file = 'sem-output-trophic-networks.csv')
write.csv(r_sq, file = 'sem-r-squareds-trophic-networks.csv')

############### END OF THE ANALYSES FOR TROPHIC NETWORKS ###############

############### BEGINNING OF ANALYSES FOR NON-TROPHIC NEGATIVE INTERACTIONS NETWORKS ###############

## we now analyse the patterns of environmental variation for the non-trophic negative interactions networks
## values of the network properties for this type of interaction are kept in the data frame
## ouptut_neg_nti_stats, which is loaded when loading the workspace (see comments above)
Y_norm <- output_neg_nti[,4:19]
Y_norm$sqrt_L <- sqrt(Y_norm$L)

## To run the models we need both, environmental predictors and network descriptors in the same data frame
data_for_sems <- cbind(X_norm, Y_norm)

## this is used later for normalisation
data_sems_sd <- apply(data_for_sems, 2, sd)

## the list of candidate models to test the dependece of each network descriptor to environmental predictors, 
## and species richness (S) is this:

## 1. glm(S ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST, na.action=na.omit, data=data_for_sems, family=quasipoisson(link = 'identity'))
## 2. glm.nb(L ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems)
## 3. lm(L.S ~ LT.SST + Upwelling + Upwelling.SD + Fr.Annual + Clim + Fr.Days + S, na.action=na.omit, data=data_for_sems)
## 4. lm(C ~ LT.SST + Upwelling + Upwelling.SD + Fr.Annual + Clim + Fr.Days + S, na.action=na.omit, data=data_for_sems)
## 5. lm(indeg ~ LT.SST + Upwelling + Upwelling.SD + Fr.Annual + Clim + Fr.Days + S, na.action=na.omit, data=data_for_sems)
## 6. lm(outdeg ~ LT.SST + Upwelling + Upwelling.SD + Fr.Annual + Clim + Fr.Days + S, na.action=na.omit, data=data_for_sems)
## 7. lm(sd_gen ~ LT.SST + Upwelling + Upwelling.SD + Fr.Annual + Clim + Fr.Days + S, na.action=na.omit, data=data_for_sems)
## 8. lm(sd_vul ~ LT.SST + Upwelling + Upwelling.SD + Fr.Annual + Clim + Fr.Days + S, na.action=na.omit, data=data_for_sems)
## 9. lm(modularity ~ LT.SST + Upwelling + Upwelling.SD + Fr.Annual + Clim + Fr.Days + S, na.action=na.omit, data=data_for_sems)
## 10. lm(fr_competition ~ LT.SST + Upwelling + Upwelling.SD + Fr.Annual + Clim + Fr.Days + S, na.action=na.omit, data=data_for_sems)
## 11. lm(fr_amensalism ~ LT.SST + Upwelling + Upwelling.SD + Fr.Annual + Clim + Fr.Days + S, na.action=na.omit, data=data_for_sems)

## Following the protocol for the analyses adopted for the non-trophic positive interaction networks (above),
## we start with models 1 and 2 with the list above:

## Before constructing and fitting the actual SEMs we verify that the models used are appropriate:
## Since both S and L are count data they were both modelled using a Generalised linear model (GLM) with the appropriate family

## In this case, the summary statistics over the S
pastecs::stat.desc(data_for_sems$S)
## reveal that S is again underdispersed, so we use quasipoisson family

## Once we fit the model:
model_test <- glm(S ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST, na.action=na.omit, data=data_for_sems, family=quasipoisson(link = 'identity'))
## by looking at the p-value resulting from a chi-square test of model deviance
## we conclude this model is a very good fit to the data, i.e. p-value >> 0.05
1-pchisq(summary(model_test)$deviance, summary(model_test)$df.residual)


###### Tests for over / underdispersion, to make sure the model we use is the right one
fit = glm(S ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST, na.action=na.omit, data=data_for_sems, family="poisson") 
fit.overdisp = glm(S ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST, na.action=na.omit, data=data_for_sems, family="quasipoisson") 
summary(fit.overdisp)$dispersion # dispersion coefficient
pchisq(summary(fit.overdisp)$dispersion * fit$df.residual, fit$df.residual, lower = F) # significance of overdispersion

dispersiontest(fit, trafo=1, alternative = 'less')

## Similarly for L, since we are dealing with count data, we fit a GLM. This time, however, since the data is overdispersed,
## as indicated by the summary statistics:
pastecs::stat.desc(data_for_sems$L)
## var = 98929.85 vs. mean = 1137.29
## we then use a negative binomial distribution for the data fit

## Once we fit the model:
model_test <- glm.nb(L ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems)

## by looking at the p-value resulting from a chi-square test of model deviance
## we conclude this model is a not a good fit to the data, i.e. p-value < 0.05.
1-pchisq(summary(model_test)$deviance, summary(model_test)$df.residual)

## Hence, we transformed the data for the response variable (L) and apply linear regression, which complies with the assumptions
data_for_sems$sqrt_L <- sqrt(data_for_sems$L)
gvlma::gvlma(lm(sqrt_L ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems))

## Once single models have been verified we construct and evaluate the SEM
sem_nti_neg <- psem(
  glm(S ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST, na.action=na.omit, data=data_for_sems, family=quasipoisson(link='identity')),
  lm(sqrt_L ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems),
  data=data_for_sems)

sem_nti_neg <- summary(sem_nti_neg)

## Then we calculate the standardised predictors estimates using the formula: estimate * (sd.x / sd.y)
## since the predictors were standardised (above), we only need to do: estimate / sd.y
for(res in unique(sem_nti_neg$coefficients$Response)){
  sd.y <- data_sems_sd[which(names(data_sems_sd) == res)]
  sem_nti_neg$coefficients[which(sem_nti_neg$coefficients$Response == res),]$Std.Estimate <- round(sem_nti_neg$coefficients[which(sem_nti_neg$coefficients$Response == res),]$Estimate / sd.y, 4)
}

coeffs_table <- sem_nti_neg$coefficients
r_sq <- sem_nti_neg$R2

## Similarly, the remaining SEMs were constructed from models 1 and 3, 1 and 4, 1 and 5... and so on

## Since the remaining variables are all continous we fit linear regressions to these data
## To ensure these models comply with the assumptions of linear regression 
## we validate these using the gvlma function from the gvlma package

## For connectance C
model_test <- lm(C ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems)
gvlma::gvlma(model_test)

sem_nti_neg <- psem(
  glm(S ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST, na.action=na.omit, data=data_for_sems, family=quasipoisson(link = 'identity')),
  lm(C ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems),
  data=data_for_sems)

sem_nti_neg <- summary(sem_nti_neg)

## Then we calculate the standardised predictors estimates using the formula: estimate * (sd.x / sd.y)
## since the predictors were standardised (above), we only need to do: estimate / sd.y
for(res in unique(sem_nti_neg$coefficients$Response)){
  sd.y <- data_sems_sd[which(names(data_sems_sd) == res)]
  sem_nti_neg$coefficients[which(sem_nti_neg$coefficients$Response == res),]$Std.Estimate <- round(sem_nti_neg$coefficients[which(sem_nti_neg$coefficients$Response == res),]$Estimate / sd.y, 4)
}

coeffs_table <- rbind(coeffs_table, sem_nti_neg$coefficients[7:13,])
r_sq <- rbind(r_sq, sem_nti_neg$R2[2,])

## for links per species - L/S
model_test <- lm(L.S ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems)
gvlma::gvlma(model_test)

sem_nti_neg <- psem(
  glm(S ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST, na.action=na.omit, data=data_for_sems, family=quasipoisson(link = 'identity')),
  lm(L.S ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems))

sem_nti_neg <- summary(sem_nti_neg)

## Then we calculate the standardised predictors estimates using the formula: estimate * (sd.x / sd.y)
## since the predictors were standardised (above), we only need to do: estimate / sd.y
for(res in unique(sem_nti_neg$coefficients$Response)){
  sd.y <- data_sems_sd[which(names(data_sems_sd) == res)]
  sem_nti_neg$coefficients[which(sem_nti_neg$coefficients$Response == res),]$Std.Estimate <- round(sem_nti_neg$coefficients[which(sem_nti_neg$coefficients$Response == res),]$Estimate / sd.y, 4)
}

coeffs_table <- rbind(coeffs_table, sem_nti_neg$coefficients[7:13,])
r_sq <- rbind(r_sq, sem_nti_neg$R2[2,])

## for indegree
model_test <- lm(indeg ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems)
gvlma::gvlma(model_test)

sem_nti_neg <- psem(
  glm(S ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST, na.action=na.omit, data=data_for_sems, family=quasipoisson(link = 'identity')),
  lm(indeg ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems))

sem_nti_neg <- summary(sem_nti_neg)

## Then we calculate the standardised predictors estimates using the formula: estimate * (sd.x / sd.y)
## since the predictors were standardised (above), we only need to do: estimate / sd.y
for(res in unique(sem_nti_neg$coefficients$Response)){
  sd.y <- data_sems_sd[which(names(data_sems_sd) == res)]
  sem_nti_neg$coefficients[which(sem_nti_neg$coefficients$Response == res),]$Std.Estimate <- round(sem_nti_neg$coefficients[which(sem_nti_neg$coefficients$Response == res),]$Estimate / sd.y, 4)
}

coeffs_table <- rbind(coeffs_table, sem_nti_neg$coefficients[7:13,])
r_sq <- rbind(r_sq, sem_nti_neg$R2[2,])

## for outdegree
model_test <- lm(outdeg ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems)
gvlma::gvlma(model_test)

sem_nti_neg <- psem(
  glm(S ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST, na.action=na.omit, data=data_for_sems, family=quasipoisson(link = 'identity')),
  lm(outdeg ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems))

sem_nti_neg <- summary(sem_nti_neg)

## Then we calculate the standardised predictors estimates using the formula: estimate * (sd.x / sd.y)
## since the predictors were standardised (above), we only need to do: estimate / sd.y
for(res in unique(sem_nti_neg$coefficients$Response)){
  sd.y <- data_sems_sd[which(names(data_sems_sd) == res)]
  sem_nti_neg$coefficients[which(sem_nti_neg$coefficients$Response == res),]$Std.Estimate <- round(sem_nti_neg$coefficients[which(sem_nti_neg$coefficients$Response == res),]$Estimate / sd.y, 4)
}

coeffs_table <- rbind(coeffs_table, sem_nti_neg$coefficients[7:13,])
r_sq <- rbind(r_sq, sem_nti_neg$R2[2,])

## sd_gen
model_test <- lm(sd_gen ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems)
gvlma::gvlma(model_test)

sem_nti_neg <- psem(
  glm(S ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST, na.action=na.omit, data=data_for_sems, family=quasipoisson(link = 'identity')),
  lm(sd_gen ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems))

sem_nti_neg <- summary(sem_nti_neg)

## Then we calculate the standardised predictors estimates using the formula: estimate * (sd.x / sd.y)
## since the predictors were standardised (above), we only need to do: estimate / sd.y
for(res in unique(sem_nti_neg$coefficients$Response)){
  sd.y <- data_sems_sd[which(names(data_sems_sd) == res)]
  sem_nti_neg$coefficients[which(sem_nti_neg$coefficients$Response == res),]$Std.Estimate <- round(sem_nti_neg$coefficients[which(sem_nti_neg$coefficients$Response == res),]$Estimate / sd.y, 4)
}

coeffs_table <- rbind(coeffs_table, sem_nti_neg$coefficients[7:13,])
r_sq <- rbind(r_sq, sem_nti_neg$R2[2,])

## sd_vul
model_test <- lm(sd_vul ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems)
gvlma::gvlma(model_test)

sem_nti_neg <- psem(
  glm(S ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST, na.action=na.omit, data=data_for_sems, family=quasipoisson(link = 'identity')),
  lm(sd_vul ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems))

sem_nti_neg <- summary(sem_nti_neg)

## Then we calculate the standardised predictors estimates using the formula: estimate * (sd.x / sd.y)
## since the predictors were standardised (above), we only need to do: estimate / sd.y
for(res in unique(sem_nti_neg$coefficients$Response)){
  sd.y <- data_sems_sd[which(names(data_sems_sd) == res)]
  sem_nti_neg$coefficients[which(sem_nti_neg$coefficients$Response == res),]$Std.Estimate <- round(sem_nti_neg$coefficients[which(sem_nti_neg$coefficients$Response == res),]$Estimate / sd.y, 4)
}

coeffs_table <- rbind(coeffs_table, sem_nti_neg$coefficients[7:13,])
r_sq <- rbind(r_sq, sem_nti_neg$R2[2,])

## modularity
model_test <- lm(modularity ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems)
gvlma::gvlma(model_test)

sem_nti_neg <- psem(
  glm(S ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST, na.action=na.omit, data=data_for_sems, family=quasipoisson(link = 'identity')),
  lm(modularity ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems))

sem_nti_neg <- summary(sem_nti_neg)

## Then we calculate the standardised predictors estimates using the formula: estimate * (sd.x / sd.y)
## since the predictors were standardised (above), we only need to do: estimate / sd.y
for(res in unique(sem_nti_neg$coefficients$Response)){
  sd.y <- data_sems_sd[which(names(data_sems_sd) == res)]
  sem_nti_neg$coefficients[which(sem_nti_neg$coefficients$Response == res),]$Std.Estimate <- round(sem_nti_neg$coefficients[which(sem_nti_neg$coefficients$Response == res),]$Estimate / sd.y, 4)
}

coeffs_table <- rbind(coeffs_table, sem_nti_neg$coefficients[7:13,])
r_sq <- rbind(r_sq, sem_nti_neg$R2[2,])

## fraction of competition
model_test <- lm(fr_competition ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems)
gvlma::gvlma(model_test)

## due to the failing of the link function test on gvlma, we tested assumptions again using
## bootstrapping simulations with the DHARMa package.
## These reveal the model looks fine (i.e. p-value of dispersion of the residuals > 0.05)

## using the Dharma package we can check the residuals of the model:
require(DHARMa)
simulationOutput <- simulateResiduals(fittedModel = model_test, n=2000, plot=F)
testResiduals(simulationOutput)
## they look fine

sem_nti_neg <- psem(
  glm(S ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST, na.action=na.omit, data=data_for_sems, family=quasipoisson(link = 'identity')),
  lm(fr_competition ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems))

sem_nti_neg <- summary(sem_nti_neg)

## Then we calculate the standardised predictors estimates using the formula: estimate * (sd.x / sd.y)
## since the predictors were standardised (above), we only need to do: estimate / sd.y
for(res in unique(sem_nti_neg$coefficients$Response)){
  sd.y <- data_sems_sd[which(names(data_sems_sd) == res)]
  sem_nti_neg$coefficients[which(sem_nti_neg$coefficients$Response == res),]$Std.Estimate <- round(sem_nti_neg$coefficients[which(sem_nti_neg$coefficients$Response == res),]$Estimate / sd.y, 4)
}

coeffs_table <- rbind(coeffs_table, sem_nti_neg$coefficients[7:13,])
r_sq <- rbind(r_sq, sem_nti_neg$R2[2,])

## fraction of amensalism
model_test <- lm(fr_amensalism ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems)
gvlma::gvlma(model_test)
## even though the gvlma checks for the assumptions of the linear regression indicate a potential inappropriate link function
## looking at the plots of the model and testing for this with the DHARMa package suggest that it is acceptable (i.e. p-value of dispersion of the residuals > 0.05)
plot(model_test)

## using the Dharma package we can check the residuals of the model:
require(DHARMa)
simulationOutput <- simulateResiduals(fittedModel = model_test, n=2000, plot=F)
testResiduals(simulationOutput)
## they look fine

sem_nti_neg <- psem(
  glm(S ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST, na.action=na.omit, data=data_for_sems, family=quasipoisson(link = 'identity')),
  lm(fr_amensalism ~ Fr.Days + Upwelling + Upwelling.SD + Fr.Annual + Clim + LT.SST + S, na.action=na.omit, data=data_for_sems))

sem_nti_neg <- summary(sem_nti_neg)

## Then we calculate the standardised predictors estimates using the formula: estimate * (sd.x / sd.y)
## since the predictors were standardised (above), we only need to do: estimate / sd.y
for(res in unique(sem_nti_neg$coefficients$Response)){
  sd.y <- data_sems_sd[which(names(data_sems_sd) == res)]
  sem_nti_neg$coefficients[which(sem_nti_neg$coefficients$Response == res),]$Std.Estimate <- round(sem_nti_neg$coefficients[which(sem_nti_neg$coefficients$Response == res),]$Estimate / sd.y, 4)
}

coeffs_table <- rbind(coeffs_table, sem_nti_neg$coefficients[7:13,])
r_sq <- rbind(r_sq, sem_nti_neg$R2[2,])

coeffs_table$Std.Estimate <- as.numeric(coeffs_table$Std.Estimate)
sd.x <- data_sems_sd[which(names(data_sems_sd) == 'S')]
coeffs_table[which(coeffs_table$Predictor == 'S'),]$Std.Estimate <- round(coeffs_table[which(coeffs_table$Predictor == 'S'),]$Std.Estimate * sd.x, 4)


coeffs_table$ranged_coeffs <- round((coeffs_table[,3] - min(coeffs_table[,3])) / (max(coeffs_table[,3]) - min(coeffs_table[,3])),4)

a <- abs(coeffs_table[,3])
ranged_abs <- (a - min(a)) / (max(a) - min(a))
ranged_abs[which(coeffs_table[,3] < 0)] <- -ranged_abs[which(coeffs_table[,3] < 0)]

coeffs_table$ranged_abs <- round(ranged_abs,4)

write.csv(coeffs_table, file = 'sem-output-negative-nti-networks.csv')
write.csv(r_sq, file = 'sem-r-squareds-negative-nti-networks.csv')


############### END OF THE ANALYSES FOR NON-TROPHIC NEGATIVE INTERACTIONS NETWORKS ###############




