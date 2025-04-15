
## ---------------------------
##
## Script name: rda-procrustes.r
##
## Purpose of script: This script implements the redundancy analysis over the network
## properties using the environmental predictors as constraints to produce a constrained
## ordination of the local networks across the environmental gradient. This analysis
## is replicated over the quantitative (i.e. probabilistic) version of the networks
## in order to compare both patterns (qualitative vs. quantitative) using ordination
## rotation via procrustes and protest analysis.
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

## Load required libraries
require(vegan)
require(RColorBrewer)

## The following analyses require 4 data frames that contain the values for the structural properties
## of each of the 3 types of interaction networks plus the residuals of the environmental predictors
## as described in the methods of the paper.

## These data can be obtained by running other scripts that are part of the source code (see README)
## Alternatively, for convenience, they have been stored in the R data file called: data-for-chilean-analysis.rda

## First, we load the data from that file
load('data-for-chilean-analysis.rda')

## RDA analysis and procrustes comparisons for the trophic networks start here
X <- decostand(environment_residuals, 'standardize')
Y <- output_trophic[c(4:12,14,16,18,19,23)]

redundancy_model <- rda(Y~., data=X, scale=T, na.action = na.omit)
summary(redundancy_model)
anova(redundancy_model, by="terms", permutations = 10000)

## From this ANOVA analysis we see that the predictors significantly influencing this ordination are:
## Fr Days, Fr Annual, and LT SST, so we repeat the RDA with these:

X <- X[c('Fr.Days', 'Fr.Annual', 'LT.SST')]
redundancy_model <- rda(Y~., data=X, scale=T, na.action = na.omit)
summary(redundancy_model)
anova(redundancy_model, by="terms", permutations = 10000)


## This code generates a figure for the redundacy analysis of network properties along the environmental
## gradient. It is useful to visualise the effects of different environmental factors over different
## network properties.
pdf('ordination-trophic-new.pdf', height=7, width = 7)
par(mar=c(4.5,4.5,1,1))
plot(redundancy_model, display=c('sp','cn'), type='none', xlab='RDA 1 (29.34 %)' , ylab='RDA 2 (10.4 %)', cex.lab=1.5, cex.axis=1.5, tck=.03, mgp = c(2.5, .5, 0))
text(redundancy_model, display=c('cn'), cex=1.2, lwd=3)

orditorp(redundancy_model, labels=c('S', 'L', 'L/S', 'C', 'Gen', 'Vul', 'SD Gen', 'SD Vul', 'B', 'T', 'I', 'O', 'Q', 'MFCL'), 
         display = "species", priority=c(1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1), 
         pch=c(1, 1, 1, 1, 1, 16, 17, 1, 1, 1, 1, 1, 1, 1), col='darkorange', cex = 1.5, air = 1)
dev.off()


## To proceed with the comparison between qualitative and quantitative networks we repeat the same
## analysis above but only with the properties for which we have quantitative equivalents

X <- decostand(environment_residuals, 'standardize')
Y <- output_trophic[c("Fraction top level qualitative","Fraction intermediate qualitative",
                  "Fraction basal qualitative", "Link density qualitative", "Connectance qualitative", 
                  "Mean chain length qualitative", "Degree of omnivory qualitative", 
                  "Generality qualitative", "Vulnerability qualitative", 
                  "SD standardised generality qualitative", "SD standardised vulnerability qualitative",
                  "modularity")]

red_analysis_qual <- rda(Y~., data=X, scale=T, na.action = na.omit)
anova(red_analysis_qual, by="terms", permutations = 10000)

## ANOVA identifies the same environmental predictors as strong drivers, so:
X <- X[c('Fr.Days', 'Fr.Annual', 'LT.SST')]

red_analysis_qual <- rda(Y~., data=X, scale=T, na.action = na.omit)

## Now we do the same with the quantitative networks
Y <- output_trophic[c("Fraction top level unweighted","Fraction intermediate unweighted",
                  "Fraction basal unweighted", "Link density unweighted", "Connectance unweighted", 
                  "Mean chain length unweighted", "Degree of omnivory unweighted", 
                  "Generality unweighted", "Vulnerability unweighted", 
                  "SD standardised generality unweighted", "SD standardised vulnerability unweighted",
                  "modularity_quant")]

red_analysis_weigh <- rda(Y~., data=X, scale=T, na.action = na.omit)
anova(red_analysis_weigh, by="terms", permutations = 10000)

trophic.proc <- procrustes(red_analysis_qual, red_analysis_weigh, scores='sites', symmetric=TRUE)

summary(trophic.proc)
plot(trophic.proc)

protest(red_analysis_qual, red_analysis_weigh, scores='sites', permutations = 10000)

## Analysis of trophic networks ends here

## Non-trophic positive interaction networks
X <- decostand(environment_residuals, 'standardize')
Y <- output_pos_nti[c(4:12,18,19)]

redundancy_model <- rda(Y~., data=X, scale=T, na.action = na.omit)
summary(redundancy_model)
anova(redundancy_model, by="terms", permutations = 10000)

## From this ANOVA analysis we see that the only predictors significantly influencing this ordination is:
## Fr Days so we repeat the RDA with this alone:

X <- X[c('Fr.Days')]
redundancy_model <- rda(Y~., data=X, scale=T, na.action = na.omit)
summary(redundancy_model)
anova(redundancy_model, by="terms", permutations = 10000)

## This code generates a figure for the redundacy analysis of network properties along the environmental
## gradient. It is useful to visualise the effects of different environmental factors over different
## network properties.
pdf('ordination-pos-nti-new.pdf', height=7, width = 7)
par(mar=c(4.5,4.5,1,1))
plot(redundancy_model, display=c('sp','cn'), type='none', xlab='RDA 1 (17.16 %)' , ylab='PC 1 (52.33 %)', cex.lab=1.5, cex.axis=1.5, tck=.03, mgp = c(2.5, .5, 0))
text(redundancy_model, display=c('cn'), cex=1.2, lwd=3)

orditorp(redundancy_model, labels=c('S', 'L', 'L/S', 'C', 'Gen', 'Vul', 'SD Gen', 'SD Vul', 'Q', 'Fr Mut', 'Fr Com'),
         display = "species", cex = 1.5, pch = 19, air = 1, col='darkorange')
dev.off()


## To proceed with the comparison between qualitative and quantitative networks we repeat the same
## analysis above but only with the properties for which we have quantitative equivalents

X <- decostand(environment_residuals, 'standardize')
Y <- output_pos_nti[c("Link density qualitative", "Connectance qualitative", 
                      "Generality qualitative","Vulnerability qualitative", 
                      "SD standardised generality qualitative", "SD standardised vulnerability qualitative",
                      "modularity")]


red_analysis_qual <- rda(Y~., data=X, scale=T, na.action = na.omit)
anova(red_analysis_qual, by="terms", permutations = 10000)

## ANOVA identifies the same environmental predictors as strong drivers, so:
X <- X[c('Fr.Days')]

red_analysis_qual <- rda(Y~., data=X, scale=T, na.action = na.omit)

## Now we do the same with the quantitative networks
Y <- output_pos_nti[c("Link density unweighted", "Connectance unweighted", 
                      "Generality unweighted", "Vulnerability unweighted", 
                      "SD standardised generality unweighted", "SD standardised vulnerability unweighted",
                      "modularity_quant")]

red_analysis_weigh <- rda(Y~., data=X, scale=T, na.action = na.omit)
anova(red_analysis_weigh, by="terms", permutations = 10000)

pos.nti.proc <- procrustes(red_analysis_qual, red_analysis_weigh, scores='sites', symmetric=TRUE)
summary(pos.nti.proc)
plot(pos.nti.proc)

protest(red_analysis_qual, red_analysis_weigh, scores='sites', permutations = 10000)

## Analysis of non-trophic positive networks ends here

## Non-trophic negative interaction networks
X <- decostand(environment_residuals, 'standardize')
Y <- output_neg_nti[c(4:12,18,19)]

redundancy_model <- rda(Y~., data=X, scale=T, na.action = na.omit)
summary(redundancy_model)
anova(redundancy_model, by="terms", permutations = 10000)

## From this ANOVA analysis we see that the only predictors significantly influencing this ordination is:
## Fr Days so we repeat the RDA with this alone:

X <- X[c('LT.SST')]
redundancy_model <- rda(Y~., data=X, scale=T, na.action = na.omit)
summary(redundancy_model)
anova(redundancy_model, by="terms", permutations = 10000)

## This code generates a figure for the redundacy analysis of network properties along the environmental
## gradient. It is useful to visualise the effects of different environmental factors over different
## network properties.
pdf('ordination-neg-nti-new.pdf', height=7, width = 7)
par(mar=c(4.5,4.5,1,1))
plot(redundancy_model, display=c('sp','cn'), type='none', xlab='RDA 1 (34.35 %)' , ylab='PC 1 (43.55 %)', cex.lab=1.5, cex.axis=1.5, tck=.03, mgp = c(2.5, .5, 0))
text(redundancy_model, display=c('cn'), cex=1.2,  lwd=3)

orditorp(redundancy_model, labels=c('S', 'L', 'L/S', 'C', 'Gen', 'Vul', 'SD Gen', 'SD Vul', 'Q', 'Fr Comp', 'Fr Amen'),
         display = "species", priority=c(1,1,1,1,1,0,1,0,1,1,1), 
         pch=c(1,1,1,1,1,16,1,17,1,1,1), col='darkorange',
         cex = 1.5, air = .7)
dev.off()


## To proceed with the comparison between qualitative and quantitative networks we repeat the same
## analysis above but only with the properties for which we have quantitative equivalents

X <- decostand(environment_residuals, 'standardize')
Y <- output_neg_nti[c("Link density qualitative", "Connectance qualitative", 
                      "Generality qualitative","Vulnerability qualitative", 
                      "SD standardised generality qualitative", "SD standardised vulnerability qualitative",
                      "modularity")]


red_analysis_qual <- rda(Y~., data=X, scale=T, na.action = na.omit)
anova(red_analysis_qual, by="terms", permutations = 10000)

## ANOVA identifies the same environmental predictors as strong drivers, so:
X <- X[c('Fr.Days', 'Fr.Annual')]

red_analysis_qual <- rda(Y~., data=X, scale=T, na.action = na.omit)

## Now we do the same with the quantitative networks
Y <- output_neg_nti[c("Link density unweighted", "Connectance unweighted", 
                      "Generality unweighted", "Vulnerability unweighted", 
                      "SD standardised generality unweighted", "SD standardised vulnerability unweighted",
                      "modularity_quant")]

red_analysis_weigh <- rda(Y~., data=X, scale=T, na.action = na.omit)
anova(red_analysis_weigh, by="terms", permutations = 10000)

neg.nti.proc <- procrustes(red_analysis_qual, red_analysis_weigh, scores='sites', symmetric=TRUE)
summary(neg.nti.proc)
plot(neg.nti.proc)

protest(red_analysis_qual, red_analysis_weigh, scores='sites', permutations = 10000)

#### End of PROCRUSTES and PROTEST analysis...

