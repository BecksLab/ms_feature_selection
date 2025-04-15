
## ---------------------------
##
## Script name: get-environmental-variables.r
##
## Purpose of script: This script extracts the values of all the environmental variables / predictors
## used to analyse the effects of environmental variability on the structure of multiplex networks
## of ecological interactions in the rocky shore of central Chile.
## It includes an algorithm for calculating upwelling index across different locations along the coast
## from satellite imagery.
##
## See the Methods section of the paper for full details. 
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

## This script reads a series of files from the home directory. Ensure these files are placed in the same
## directory from which this script is executed

## Load the required libraries
require(ncdf4)
require(raster)
require(ggplot2)

## This table contains the code names and coordinates of the sampling sites used in this study
ordered_sites <- read.csv('sites-with-coordinates.csv', row.names = 1)

## This .nc file contains the satellite imagery used to calculate the upwelling as detailed in the methods
## section of the paper

## The image specification is:
## SST (Sea Surface Temperature), Pathfinder Ver 5.2 (L3C), Night, Global, 0.0417°, 1981-2012, Science Quality (8 Day Composite)
## obtained from coastwatch.pfeg.noaa.gov
## This file is too large to be provided as part of the supplementary information of the paper. Download beforehand from
## the indicated source and rename it to match the following line of code and in the 'stack' instruction below.
file <- nc_open('../erdPH2sstn8day_9a89_8086_d7fa.nc')
str(file)

lats <- ncvar_get(file, "latitude")
lon <- ncvar_get(file, "longitude")
time <- as.Date(as.POSIXct(as.numeric(as.character(ncvar_get(file, "time"))),origin="1970-01-01",tz="GMT"))

## For ease of manipulation we use the raster package to access the layers.

stacked_raster <- stack('../erdPH2sstn8day_9a89_8086_d7fa.nc')
plot(stacked_raster[[1]])
plot(stacked_raster[[length(stacked_raster@layers)]])


## The following lines of code implement the algorithm detailed in the methods section of the paper
## that permits the calculation of the upwelling index at the pixel level for each day of the
## specified period of time. Refer to the manuscript for details.

## it was agreed with Sergio that the UI would be calculated on a daily basis
## and then the average would be obtained
months_to_include <- c('01', '02', '09', '10', '11', '12')
buffer_matrix <- matrix(rep(1,64), ncol=8, byrow=TRUE)
buffer_matrix[4,4] <- 0

buffer_on <- matrix(rep(1,9), ncol=3, byrow=TRUE)
buffer_on[2,2] <- 0

buffers_cells <- c()

out_upwellings <- NULL
t_bottom <- 10
for(i in 1:length(time)){
  if(!format(time[i], '%m') %in% months_to_include) next
 
  cur_layer <- stacked_raster[[i]]
  
  for(s in unique(ordered_sites$site)){
    temp_site <- ordered_sites[which(ordered_sites$site == s),]
    t_on <- raster::extract(cur_layer, temp_site[c('fixed_lon2', 'fixed_lat2')])
    if(is.na(t_on) | is.nan(t_on)) next # | (t_on < 10) | (t_on > 20.5)) next
    
    dists_off <- as.data.frame(distanceFromPoints(cur_layer, temp_site[c('fixed_lon2', 'fixed_lat2')]), xy=T)
    dists_off$layer[which(is.nan(dists_off$layer))] <- 0
    
    closest_lat <- min(abs(unique(dists_off$y) - temp_site$fixed_lat2))
    closest_lat <- which(abs(unique(dists_off$y) - temp_site$fixed_lat2) == closest_lat)
    closest_lat <- unique(dists_off$y)[closest_lat]
    
    dists_off <- subset(dists_off, (y == closest_lat))
    
    delta_dists <- abs(dists_off$layer - 350000)
    off_idx <- which(delta_dists == min(delta_dists))
    off_pt <- dists_off[off_idx,]
    
    buffer <- adjacent(cur_layer, cellFromXY(cur_layer, off_pt[c('x','y')]), pairs=F, include=T, directions=buffer_matrix)
    
    buffers_cells <- append(buffers_cells, buffer)
    
    buffer <- raster::extract(cur_layer, buffer)
    t_off <- mean(buffer[!is.nan(buffer)], na.rm=T)
    
    if(is.nan(t_off)) next
    
    upwell_idx <- (t_off - t_on)/(t_off - t_bottom)
    
    if(is.na(upwell_idx) | upwell_idx > 5) next
    
    cur_out <- data.frame(date=time[i], site=s, upwell_idx, t_on)
    if(is.null(out_upwellings)){
      out_upwellings <- cur_out
    }else{
      out_upwellings <- rbind(out_upwellings, cur_out)
    }
  }
  
}

out_upwellings$lat <- ordered_sites[match(out_upwellings$site, ordered_sites$site),]$fixed_lat2
out_upwellings$lon <- ordered_sites[match(out_upwellings$site, ordered_sites$site),]$fixed_lon2

## values smaller than -.5 suggest potential problems with these data points (e.g. miscalculations due to
## extreme values of SST). In our case only one value falls under this. This wouldn't affect our results.
if(length(which((out_upwellings$upwell_idx < -.5))) > 0){
  out_upwellings <- out_upwellings[-which((out_upwellings$upwell_idx < -.5)),]  
}

## This is how the upwelling data looks like
hist(out_upwellings$upwell_idx, breaks=100)

out_upwellings <- out_upwellings[order(out_upwellings$lat),]
out_upwellings$site <- factor(out_upwellings$site, levels=unique(out_upwellings$site))

## These are the upwelling values across latitutde
ggplot(out_upwellings, aes(lat, upwell_idx, colour=site)) + geom_point()

## Here we obtain the average and standard deviation of the upwelling index 
## (see Methods section of the manuscript)
avgs <- aggregate(upwell_idx~site, out_upwellings, FUN=mean)
sds <- aggregate(upwell_idx~site, out_upwellings, FUN=sd)
avgs$lat <- ordered_sites[match(avgs$site, ordered_sites$site),]$fixed_lat2
avgs$lon <- ordered_sites[match(avgs$site, ordered_sites$site),]$fixed_lon2

avgs$sd <- sds$upwell_idx

#### This is to count the fraction of days with temperature below 14 degrees
avgs$fr_days <- 0
for(s in avgs$site){
  temp <- out_upwellings[which(out_upwellings$site == s),]
  avgs[which(avgs$site == s),]$fr_days <- length(which(temp$t_on < 14))/dim(temp)[1]
}


## Now let's see how the muzic index (Tapia et al. 2009) relates to the newly calculated
## upwelling values

## This is the data from Figure 5b in Tapia et al. 2009 
## (Tapia, F.J., et al. Thermal indices of upwelling effects on inner-shelf habitats. 
## Prog. Oceanogr. (2009), doi:10.1016/j.pocean.2009.07.035)
muzic <- read.csv('Data_Fig5b_Tapia_et_al_2009.csv', stringsAsFactors = F)
muzic <- muzic[which(muzic$Name %in% ordered_sites$site),]

avgs$muzic <- NA
avgs$muzic[match(muzic$Name, avgs$site)] <- muzic$CR

## And this plit shows how it relates to the upwelling index in our study sites
## This code produces Supplementary Figure 2 of the paper
pdf('supp-fig2.pdf', height=5, width=6)
ggplot(avgs, aes(muzic, upwell_idx)) + geom_point(colour="#E6AB02") + 
  geom_text(aes(label=site),hjust=0, vjust=0, size=3.5) +
  geom_abline() + theme_bw() + ylab(expression('UI'['Avg'])) + xlab('\nCooling rate (CR)') +
  theme(axis.text.x = element_text(colour='black', size = 12), 
        axis.text.y = element_text(colour='black', size = 12),
        axis.title.x = element_text(colour='black', size = 15), 
        axis.title.y = element_text(colour='black', size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlim(0.25,0.5) + ylim(0.25,0.5)
dev.off()

## The relationship is statistically significant
summary(lm(upwell_idx~muzic, data=avgs, na.action = na.omit))

## Lastly, we assign the upwelling categories to each of our study site
avgs$upwell_cat <- NA
avgs[which(avgs$site == 'PTAL'),]$upwell_cat <- 'STRONG'
avgs[which(avgs$site == 'POSC'),]$upwell_cat <- 'STRONG'
avgs[which(avgs$site == 'CUR'),]$upwell_cat <- 'STRONG'
avgs[which(avgs$site == 'PTL'),]$upwell_cat <- 'STRONG'
avgs[which(avgs$site == 'BUCA'),]$upwell_cat <- 'STRONG'
avgs[which(avgs$site == 'PEL'),]$upwell_cat <- 'STRONG'
avgs[which(avgs$site == 'BUP'),]$upwell_cat <- 'STRONG'

avgs[which(avgs$site == 'QUN'),]$upwell_cat <- 'INTERMEDIATE'
avgs[which(avgs$site == 'ELQ'),]$upwell_cat <- 'INTERMEDIATE'

avgs[which(avgs$site == 'TEM'),]$upwell_cat <- 'WEAK'
avgs[which(avgs$site == 'ARR'),]$upwell_cat <- 'WEAK'
avgs[which(avgs$site == 'GUA'),]$upwell_cat <- 'WEAK'
avgs[which(avgs$site == 'MOLL'),]$upwell_cat <- 'WEAK'
avgs[which(avgs$site == 'MONT'),]$upwell_cat <- 'WEAK'
avgs[which(avgs$site == 'ECIMN'),]$upwell_cat <- 'WEAK'
avgs[which(avgs$site == 'LCRUC'),]$upwell_cat <- 'WEAK'
avgs[which(avgs$site == 'CON'),]$upwell_cat <- 'WEAK'

avgs$upwell_cat <- factor(avgs$upwell_cat, levels=c('WEAK', 'INTERMEDIATE', 'STRONG', NA))

## This code generates Supplementary Figure 3 of the paper.
pdf('supp-fig3.pdf', height=5, width=6)
ggplot(avgs, aes(upwell_cat, upwell_idx)) + geom_point(colour="#E6AB02") + 
  geom_text(aes(label=site),hjust=0, vjust=0, size=3.5) + 
  theme_bw() + ylab(expression('UI'['Avg'])) + xlab('\nUpwelling category') +
  theme(axis.text.x = element_text(colour='black', size = 12), 
        axis.text.y = element_text(colour='black', size = 12),
        axis.title.x = element_text(colour='black', size = 15), 
        axis.title.y = element_text(colour='black', size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dev.off()

## Finally, this is how latitude relates to upwelling index in our dataset
ggplot(avgs, aes(lat, upwell_idx)) + geom_point() + 
  geom_text(aes(label=site),hjust=0, vjust=0, size=2.5)

## Information contained in the avgs table at this point can be used to construct Supplementary Table 1 in the
## manuscript

## Once we have upwelling and its related measures we obtain the environmental variables obtained from AVHRR satellite imagery.
## These environmental variables are: (1) long-term mean SST (LT SST), (2) fraction of the annual cycle not explained by season (Fr Annual),
## and (3) climatology (Clim) - See the Methods section in the paper for a detailed account of these variables.

## These data were kindly provided by Bernardo

avhrr_stats <- read.csv('avhrr_environmental.csv', stringsAsFactors = FALSE)

## we add these variables to our avgs data frame where the others are stored
avgs$fr_annual <- avhrr_stats[match(as.character(avgs$site), avhrr_stats$site),]$FRACTION.ANUAL
avgs$lt_mean <- avhrr_stats[match(as.character(avgs$site), avhrr_stats$site),]$MEAN
avgs$climatology <- avhrr_stats[match(as.character(avgs$site), avhrr_stats$site),]$CLIMATOLOGY.BELOW13

## Once we have all our environmental variables we proceed to obtain their residuals from long-term mean
## This is to ensure there is as little correlation among predictors as possible.
## See the Methods section of the paper for further details.

## Correlation among predictors can be visualised using the following code, which generates Supplementary
## Figure 4 of the paper
env_to_plot <- avgs[c("lt_mean", "upwell_idx", "sd", "fr_days", "fr_annual", "climatology" )]
names(env_to_plot) <- c('LT.SST', 'UI.Avg', 'UI.SD', 'Fr.Days', 'Fr.Annual', 'Clim')
pdf('supp-fig4.pdf', height=5, width=7)
chart.Correlation(env_to_plot, histogram = TRUE, pch=19, cex.axis=10)
dev.off()

environment <- avgs

r <- summary(lm(upwell_idx~lt_mean, data=environment))$residuals
environment_residuals <- data.frame(upw_res = (r))

r <- summary(lm(sd~lt_mean, data=environment))$residuals
environment_residuals <- cbind(environment_residuals, data.frame(sd_res = (r)))

r <- summary(lm(fr_days~lt_mean, data=environment))$residuals
environment_residuals <- cbind(environment_residuals, data.frame(frd_res = (r)))

r <- summary(lm(fr_annual~lt_mean, data=environment))$residuals
environment_residuals <- cbind(environment_residuals, data.frame(fran_res = (r)))

r <- summary(lm(climatology~lt_mean, data=environment))$residuals
environment_residuals <- cbind(environment_residuals, data.frame(clim_res = (r)))

environment_residuals$lt_mean <- environment$lt_mean
row.names(environment_residuals) <- environment$site

names(environment_residuals) <- c('Upwelling', 'Upwelling.SD', 'Fr.Days', 'Fr.Annual', 'Clim', 'LT.SST')

## After processing the environmental variables above we finally obtain the values that are going to be
## for the statistical analyses and ultimately the fit of the SEMs

## These data are stored in 'environment_residuals'
## and we order it in the same order as the ordered sites to ensure match with the network variables
environment_residuals <- environment_residuals[match(ordered_sites$site, row.names(environment_residuals)),]

environment_residuals
