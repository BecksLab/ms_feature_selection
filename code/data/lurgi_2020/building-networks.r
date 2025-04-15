## ---------------------------
##
## Script name: building-networks.r
##
## Purpose of script: This script constructs networks of ecological networks
## at different localities based on data from the rocky shore marine intertidal
## of the coast of central Chile. Additionally, a series of network properties are
## calculated that are then used for further statistical analyses.
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

## Load required libraries
require(vegan)
require(ggrepel)
require(igraph)
require(cheddar)
require(ggplot2)
source('./utils.r')
source('./auxiliary-functions.r')

spatial_comms <- read.csv('chilean-rocky-shores-samples.csv', header=T, stringsAsFactors = F)[-1]

#### let's find out if there are species not present in any site
names(spatial_comms) # this tells us that the species go from column 1 to 112

#### Below, we obtain the equivalences between spatial points and sites names
dict_names_coords <- unique(spatial_comms[c("WESTLONG", "SOUTHLAT", "SITE")])

#### I identified a problem with duplicated site entries due to several different WESTLONG values
#### assigned to the same site
#### The following code corrects that:
#### here we obtain the locations with value with the most number of repetitions for the duplicated
#### locations and we change it on the spatial data frame to have all the records of a given
#### site with the same coordinates
dup_locations <- dict_names_coords[which(duplicated(dict_names_coords$SITE)),]$SITE

### sites that appear more than once are: MONT ECIMN & LCRUC
for(d in dup_locations){
  print(d)
  rank <- table(spatial_comms[which(spatial_comms$SITE == d),]$WESTLONG)
  
  cur_counts <- 0
  for(i in 1:length(rank)){
    if(cur_counts < rank[i]){
      cur_counts <- rank[i];
      cur_value <- as.numeric(names(rank)[i]);
    }
  }
  
  spatial_comms[which(spatial_comms$SITE == d),]$WESTLONG <- cur_value
}

#### now we obtain the equivalences between spatial points and sites names
#### after fixing the naming issue
dict_names_coords <- unique(spatial_comms[c("WESTLONG", "SOUTHLAT", "SITE")])

df <- dict_names_coords
names(df) <- c('lon', 'lat', 'site')
df$lon <- df$lon * (-1)
df$lat <- df$lat * (-1)

###### let's see what's the latitudinal gradient of species richness (if any)
ordered_sites <- df[with(df, order(-lat)), ]
##### as a first approximation we will look at the total number
##### of species reported per site regardless of sampling effort heterogeneity
ordered_sites$richness <- 0

for(s in ordered_sites$site){
  temp <- subset(spatial_comms, (SITE==s))
  temp <- temp[,1:114]
  species_in_site <- names(which(colSums(temp)!=0))
  n_species <- length(species_in_site)
  ordered_sites[which(ordered_sites$site==s),]$richness <- n_species
}

model <- lm(richness~lat, ordered_sites)
summary(model)

plot(ordered_sites$lat, ordered_sites$richness, xlab='latitude', ylab='richness')
lines(ordered_sites$lat, predict.lm(model,interval='prediction',type='resp')[,1], col="red")

## to assess the effects of sampling effort we can have a look at species accumulation curves
## These suggest that the sampling effort is good enough to reach a stable number of species in the communities
## Here is where we generate Supplementary Figure 1
data_for_rarefaction <- NULL
samples <- c()

pdf('sup-fig1.pdf', width=10, height = 7)
plot(0,0, ylim=c(0,80), xlim=c(0, 150), cex=0, xlab='No. of samples', ylab="Species")

c <- 1
cols <- rainbow(length(unique(ordered_sites$site)))

xs <- c()
ys <- c()
for(s in unique(ordered_sites$site)){
  cur_comms <- subset(spatial_comms, (SITE == s & HEIGHT %in% c('L','M')))
  print(paste('number of samples = ', dim(cur_comms)[1]))
  
  if(dim(cur_comms)[1] == 0) next
  
  samples <- append(samples, dim(cur_comms)[1])
  
  temp_comm <- cur_comms[,1:114]
  temp_comm[temp_comm > 0] <- as.integer(1)
  
  if(is.null(data_for_rarefaction)) data_for_rarefaction <- temp_comm
  else data_for_rarefaction <- rbind(data_for_rarefaction, temp_comm)
  
  sp_acc <- specaccum(temp_comm)
  plot(sp_acc,add=T, col=cols[c], ci=0, lwd=2)
  xs <- append(xs, sp_acc$sites[length(sp_acc$sites)] + 5)
  ys <- append(ys, sp_acc$richness[length(sp_acc$richness)])
  c <- c + 1
  
}

text(xs, ys, unique(ordered_sites$site), cex=.5)

dev.off()

###### End of Supplementary Figure 1 #######

### at the regional level...
sp_acc <- specaccum(data_for_rarefaction, col = rainbow(dim(data_for_rarefaction)[1]), cex=.6)
plot(sp_acc)

hist(samples)

## Just for fun, we draw the site locations on a map, to see how it looks like!
chile_sst_2014 <- my_get_sst_chile('2014')

# look at the structure of the result
str(chile_sst_2014)

msst <-  melt_sst(chile_sst_2014)
head(msst)

df <- dict_names_coords
names(df) <- c('lon', 'lat', 'site')
df$lon <- df$lon * (-1)
df$lat <- df$lat * (-1)

ggplot() + 
  geom_raster(data = msst, aes(x = long, y = lat, fill = sst), interpolate = TRUE) +
  scale_fill_gradientn(colours = rev(rainbow(7)), na.value = NA) +
  theme_bw() +
  coord_fixed(1.3) +
  geom_point(data = df, aes(x = lon, y = lat), size = 2, shape = 19) +
  geom_text_repel(aes(x = df$lon, y = df$lat, label = df$site), size=2.5, fontface ='bold', nudge_x = 2)


#### after conversations with Bernardo, it seems also appropriate to remove the site matanzas = MAZ
#### because it is subject to a succesional regime due to sand flooding and hence probably
#### not a good representation of reality
spatial_comms <- spatial_comms[-which(spatial_comms$SITE %in% c('MAZ')),]
ordered_sites <- ordered_sites[-which(ordered_sites$site %in% c('MAZ')),]

#### so, this is the final spatial_comms data frame


############# This code is for obtaining taxonomic information for the species from public repositories
## Execute this code to obtain the information included in the supplementary table 2 of the paper
# require(taxize)
# new_metadata <- read.csv('chilean_metadata_for_paper.csv', stringsAsFactors = F)
# for(sp in new_metadata$Species.name){
#   if(grepl('spp.', sp)){
#     temp_class <- classification(strsplit(sp, ' ')[[1]][1], db='worms')[[1]]
#   }else{
#     temp_class <- classification(sp, db='worms')[[1]]
#   }
# 
#   if(is.na(temp_class)){
#     print(paste('no register for', sp))
#     next
#   }
# 
#   if('Phylum' %in% temp_class$rank){
#     new_metadata[which(new_metadata$Species.name == sp),]$Phylum <- temp_class[which(temp_class$rank == 'Phylum'),]$name
#   }else if('Phylum (Division)' %in% temp_class$rank){
#     new_metadata[which(new_metadata$Species.name == sp),]$Phylum <- temp_class[which(temp_class$rank == 'Phylum (Division)'),]$name
#   }
# 
#   if('Subphylum' %in% temp_class$rank){
#     new_metadata[which(new_metadata$Species.name == sp),]$Subphylum <- temp_class[which(temp_class$rank == 'Subphylum'),]$name
#   }else if('Subphylum (Subdivision)' %in% temp_class$rank){
#     new_metadata[which(new_metadata$Species.name == sp),]$Subphylum <- temp_class[which(temp_class$rank == 'Subphylum (Subdivision)'),]$name
#   }else{
#     new_metadata[which(new_metadata$Species.name == sp),]$Subphylum <- NA
#   }
# 
#   if('Class' %in% temp_class$rank){
#     new_metadata[which(new_metadata$Species.name == sp),]$Class <- temp_class[which(temp_class$rank == 'Class'),]$name
#   }
# }


#################### After the data processing and selection above, network construction and 
#################### analysis begins here!

#### Now the real fun! Let's build networks!!!!
#### these are the species in the database of species and codes provided by Evie and Sergio
species <- read.csv('species-dict.csv', stringsAsFactors = FALSE)

#### these are the species found everywhere that are not recorded in the spatial database
permanent_species <- c('petrolisthes punctatus', 'petrolisthes spinifrons', 'petrolisthes angulosus', 'petrolisthes tuberculatus','petrolisthes tuberculosus', 'gulls', 'cincloides', 'onchidella')

#### This is the network of trophic interactions
trophic_ints <- read.csv('chilean_TI.txt', sep = '\t', stringsAsFactors = F)

#### These are the species in the interaction matrix:
sps_in_matrix <- trophic_ints$X.1
metadata <- read.csv('chilean_metadata.csv', stringsAsFactors = F)

#### These are the species that appear in the metadata with a different name than on the interaction matrix
setdiff(sps_in_matrix, metadata$Species.names)

trophic_ints <- trophic_ints[,-1]
names(trophic_ints)[1] <- 'name'
names(trophic_ints)[2:(dim(trophic_ints)[2])] <- trophic_ints$name

rownames(trophic_ints) <- trophic_ints$name
trophic_ints <- trophic_ints[,-1]

species_in_network <- names(trophic_ints)

## These are the species that are found both, in the network and in our species names dictionary
intersect(unique(species$name), unique(species_in_network))

#we can create a graph from the adjacency matrix of the network
trophic_ints <- as.matrix(trophic_ints)

g_trophic <- graph.adjacency(trophic_ints)
g_trophic <- induced.subgraph(g_trophic, V(g_trophic)[which(degree(g_trophic) != 0)])

## this is how the trophic network looks like
plot.igraph(g_trophic, vertex.label=NA, layout = layout.random, edge.arrow.size=.5, edge.curved=T, vertex.size=7)

#### and explore some properties
S <- vcount(g_trophic)
L <- ecount(g_trophic)
L.S <- L/S
C <- L/((S**2)-S)

############ This code is for assigning the full names to species in order to create the
############ interactions matrix

# V(g)$name <- species[match(V(g)$name, species$name),]$corrected_name
# write.csv(as_adjacency_matrix(g, sparse = FALSE), file = 'trophic-network.csv')

###### now, we create networks along the latitudinal gradients

## here, we remove the species from the database for which we do not have a taxonomic assignment
comms_sp_networks <- spatial_comms[,-which(!(names(spatial_comms)[1:112]) %in% species$code)]

#### species that are on the dataset of the survey samples but that are not on the networks are:
####  "TON_SP", "COLLISELLA_SP", "FIS_SP"           
####  "CALIPTR", "TALIEPUS_DENTATUS"
####  "CRYPTOM", "AHNFELTI", "RHODYM"           
####  "MACR_PRI", "MACR_PYR"

## in this data frame we keep the output
output_trophic <- NULL

## In this loop, all the sites are checked for their samples and the trophic networks at each site are
## built based on the information from the metaweb of trophic interactions (stored in the 'g' graph)
## then, all network properties are calculated and the output written to output_trophic

## we use the ordered sites as a reference to keep the output ordered by latitude
for(site in ordered_sites$site){
  current_comm <- subset(comms_sp_networks, (SITE == site))
  lat <- current_comm$SOUTHLAT[1]
  samples_no <- dim(current_comm)[1]
  
  if(dim(current_comm)[1] == 0) next
  
  abs <- current_comm[,which((names(current_comm)) %in% species$code)]
  samples_height <- dim(abs)[1]
  print(paste('site = ', site, ' - samples ', dim(current_comm)[1]))
  
  sp_local_network <- names(abs)[which(colSums(abs) != 0)]
  sp_local_network <- as.character(species[which(species$code %in% sp_local_network),]$name)
  current_g <- induced.subgraph(g_trophic, sp_local_network)
  
  #### now we remove the nodes that have 0 degree
  current_g <- induced.subgraph(current_g, V(current_g)[which(degree(current_g) != 0)])
  
  #### in case the network is empty
  if(length(V(current_g)) == 0) next
  
  S <- vcount(current_g)
  L <- ecount(current_g)
  L.S <- L/S
  C <- L/((S**2)-S)
  
  M <- as.matrix(get.adjacency(current_g))
  
  indeg <- sum(colSums(M))/sum(colSums(M) != 0)
  outdeg <-  sum(rowSums(M))/sum(rowSums(M) != 0)
  
  sd_gen <- SDGenerality(M)
  sd_vul <- SDVulnerability(M)
  
  basal <- FractionOfBasal(M)
  n_basal <- NumberOfBasal(M)
  
  top <- FractionOfTop(M)
  n_top <- NumberOfTop(M)
  
  interm <- FractionOfIntermediate(M)
  n_interm <- NumberOfIntermediate(M)
  
  mod <- tryCatch({
    cluster_louvain(as.undirected(current_g))
  }, warning = function(w){
    NA
  }, error = function(e){
    NA
  }, finally = {
    
  })
  
  if(is.na(mod)){
    modularity <- -1
    n_modules <- 0
  }else{
    modularity <- max(mod$modularity)
    n_modules <- max(mod$membership)
    
    V(current_g)$module <- mod$membership
    
    ## uncomment this line if you want to save the networks to the home folder
    # write.graph(current_g, file=paste0('./networks/trophic-net-',site,'.graphml'), format = 'graphml')
  }
  
  mfcl <- tryCatch({
    MeanFoodChainLength(M)
  }, warning = function(w) {
    NA
  }, error = function(e) {
    NA
  }, finally = {
    
  })
  
  ###### this bit was added to get the quantified version of the network and its corresponding properties
  g_edges <- ends(current_g, E(current_g))
  comm_size <- dim(abs)[1]
  edges_weights <- c()
  for(e in 1:dim(g_edges)[1]){
    prey <- as.character(species[which(species$name == g_edges[e,][1]),]$code)
    predator <- as.character(species[which(species$name == g_edges[e,][2]),]$code)
    co_occur <- (length(which(rowSums(abs[prey]) != 0 )) / comm_size) * (length(which(rowSums(abs[predator])!= 0 )) / comm_size)
    edges_weights <- append(edges_weights, co_occur)
  }
  current_g <- set.edge.attribute(current_g, 'interaction.probability', value=edges_weights)
  ###### end of network links quantification
  
  undirec_g <- as.undirected(current_g, edge.attr.comb = 'first')
  mod_quant <- tryCatch({
    cluster_louvain(undirec_g, weights=E(undirec_g)$interaction.probability)
  }, warning = function(w){
    NA
  }, error = function(e){
    NA
  }, finally = {
    
  })
  
  if(is.na(mod_quant)){
    modularity_quant <- -1
    n_modules_quant <- 0
  }else{
    modularity_quant <- max(mod_quant$modularity)
    n_modules_quant <- max(mod_quant$membership)
  }
  
  ###### now that we have quantified interactions we can use cheddar to calculate Bersier's properties
  M_Q <- get.adjacency(current_g, attr='interaction.probability', sparse=FALSE)
  cheddar_network <- cheddar::Community(nodes=data.frame(node=rownames(M_Q)), trophic.links=PredationMatrixToLinks(M_Q, link.property = 'interaction.probability'), properties=list(title=site));
  
  omniv <- cheddar::Omnivory(cheddar_network)
  quants <- t(QuantitativeDescriptors(cheddar_network, 'interaction.probability'))
  
  cur_out <- data.frame(site, lat, samples_no, S, L, L.S, C, indeg, outdeg, sd_gen, sd_vul, basal, n_basal, top, n_top, interm, n_interm, omniv, modularity, n_modules, modularity_quant, n_modules_quant, mfcl)
  
  ##### here we add the results of the quantitative analysis to the output data frame
  temp <- quants[1,]
  names(temp) <- paste0(names(temp),' qualitative')
  cur_out <- cbind(cur_out, t(temp))
  temp <- quants[2,]
  names(temp) <- paste0(names(temp),' unweighted')
  cur_out <- cbind(cur_out, t(temp))
  temp <- quants[3,]
  names(temp) <- paste0(names(temp),' weighted')
  cur_out <- cbind(cur_out, t(temp))
  
  ###### that's it
  
  if(is.null(output_trophic)) output_trophic <- cur_out
  else output_trophic <- rbind(output_trophic, cur_out)
}

output_trophic$lat <- -1*output_trophic$lat
row.names(output_trophic) <- as.character(output_trophic$site)


## we perform a similar procedure to construct and analyse the non-trophic negative interactions networks
## this time using the information from the NTI negative matrix

nt_neg_ints <- read.csv('chilean_NTIneg.txt', sep = '\t', stringsAsFactors = F)

nt_neg_ints <- nt_neg_ints[,-1]
names(nt_neg_ints)[1] <- 'name'
names(nt_neg_ints)[2:(dim(nt_neg_ints)[2])] <- nt_neg_ints$name

rownames(nt_neg_ints) <- nt_neg_ints$name
nt_neg_ints <- nt_neg_ints[,-1]

species_in_network <- names(nt_neg_ints)

setdiff(unique(species_in_network), unique(species$name))

#we can create a graph from the adjacency matrix of the network
nt_neg_ints  <- as.matrix(nt_neg_ints)
g_nt_neg <- graph.adjacency(nt_neg_ints)

### we remove the nodes that do not have any interactions
g_nt_neg <- induced.subgraph(g_nt_neg, V(g_nt_neg)[which(degree(g_nt_neg) != 0)])

S_neg <- vcount(g_nt_neg)
L_neg <- ecount(g_nt_neg)

par(mar=c(0,0,0,0))
plot.igraph(g_nt_neg, vertex.label=NA, layout = layout.random, edge.arrow.size=.5, edge.curved=T, vertex.size=7)

output_neg_nti <- NULL
for(site in ordered_sites$site){
  current_comm <- subset(comms_sp_networks, SITE == site)
  lat <- current_comm$SOUTHLAT[1]
  samples_no <- dim(current_comm)[1]
  
  abs <- current_comm[,which((names(current_comm)) %in% species$code)]
  print(paste('site = ', site, ' - samples ', samples_no))
  
  sp_local_network <- names(abs)[which(colSums(abs) != 0)]
  sp_local_network <- as.character(species[which(species$code %in% sp_local_network),]$name)
  
  sp_local_network <- intersect(V(g_nt_neg)$name, sp_local_network)
  current_g <- induced.subgraph(g_nt_neg, sp_local_network)
  
  #### now we remove the nodes that have 0 degree
  current_g <- induced.subgraph(current_g, V(current_g)[which(degree(current_g) != 0)])
  
  S <- vcount(current_g)
  L <- ecount(current_g)
  L.S <- L/S
  C <- L/((S**2)-S)
  
  ## this is the local network in adjacency matrix format
  M <- as.matrix(get.adjacency(current_g))
  
  ## some further network properties
  indeg <- sum(colSums(M))/sum(colSums(M) != 0)
  outdeg <-  sum(rowSums(M))/sum(rowSums(M) != 0)
  sd_gen <- SDGenerality(M)
  sd_vul <- SDVulnerability(M)
  
  mod <- tryCatch({
    cluster_louvain(as.undirected(current_g))
  }, warning = function(w){
    NA
  }, error = function(e){
    NA
  }, finally = {
    
  })
  
  if(is.na(mod)){
    modularity <- -1
    n_modules <- 0
  }else{
    modularity <- max(mod$modularity)
    n_modules <- max(mod$membership)
    
    V(current_g)$module <- mod$membership
    
  }
  
  ## uncomment this line if you want to save the local networks to the hard drive : in graphml format
  ## write.graph(current_g, file=paste0('./networks/nti-negative-net-',site,'.graphml'), format = 'graphml')

  ## here we calculate the number of competitive vs amensalistic interactions  
  competition <- 0
  amensalism <- 0
  adj_mat <- get.adjacency(current_g, sparse=F)
  simplified_g <- current_g
  
  for(sp1 in 1:vcount(current_g)){
    for(sp2 in 1:vcount(current_g)){
      if(sp1 <= sp2) next
      if(adj_mat[sp1,sp2] != 0 & adj_mat[sp2,sp1] != 0){
        competition <- competition + 1
        simplified_g <- simplified_g - E(simplified_g)[sp1 %->% sp2]
      }else if(adj_mat[sp1,sp2] == 0 & adj_mat[sp2,sp1] != 0) amensalism <- amensalism + 1
      else if(adj_mat[sp1,sp2] != 0 & adj_mat[sp2,sp1] == 0) amensalism <- amensalism + 1
      
    }
  }
  
  ###### this bit was added to get the quantified version of the network and its corresponding properties
  g_edges <- ends(current_g, E(current_g))
  comm_size <- dim(abs)[1]
  edges_weights <- c()
  for(e in 1:dim(g_edges)[1]){
    prey <- as.character(species[which(species$name == g_edges[e,][1]),]$code)
    predator <- as.character(species[which(species$name == g_edges[e,][2]),]$code)
    co_occur <- (length(which(rowSums(abs[prey]) != 0 )) / comm_size) * (length(which(rowSums(abs[predator])!= 0 )) / comm_size)
    edges_weights <- append(edges_weights, co_occur)
  }
  current_g <- set.edge.attribute(current_g, 'interaction.probability', value=edges_weights)
  ###### end of network links quantification
  
  undirec_g <- as.undirected(current_g, edge.attr.comb = 'first')
  mod_quant <- tryCatch({
    cluster_louvain(undirec_g, weights=E(undirec_g)$interaction.probability)
  }, warning = function(w){
    NA
  }, error = function(e){
    NA
  }, finally = {
    
  })
  
  if(is.na(mod_quant)){
    modularity_quant <- -1
    n_modules_quant <- 0
  }else{
    modularity_quant <- max(mod_quant$modularity)
    n_modules_quant <- max(mod_quant$membership)
  }
  
  ###### now that we have quantified interactions we can use cheddar to calculate Bersier's properties
  M_Q <- get.adjacency(current_g,  attr='interaction.probability', sparse=FALSE)
  cheddar_network <- cheddar::Community(nodes=data.frame(node=rownames(M_Q)), trophic.links=PredationMatrixToLinks(M_Q, link.property = 'interaction.probability'), properties=list(title=site));
  quants <- t(QuantitativeDescriptorsNoChain(cheddar_network, 'interaction.probability'))
  
  cur_out <- data.frame(site, lat, samples_no, S, L, L.S, C, indeg, outdeg, sd_gen, sd_vul, modularity, n_modules, modularity_quant, n_modules_quant, competition, amensalism, fr_competition=(competition/(L-competition)), fr_amensalism=(amensalism/(L-competition)))
  ##### here we add the results of the quantitative analysis to the output data frame
  temp <- quants[1,]
  names(temp) <- paste0(names(temp),' qualitative')
  cur_out <- cbind(cur_out, t(temp))
  temp <- quants[2,]
  names(temp) <- paste0(names(temp),' unweighted')
  cur_out <- cbind(cur_out, t(temp))
  temp <- quants[3,]
  names(temp) <- paste0(names(temp),' weighted')
  cur_out <- cbind(cur_out, t(temp))
  
  ###### that's it
  
  ##### let's try to get quantified versions of competition and amensalism
  competition <- 0
  amensalism <- 0
  adj_mat <- M_Q
  
  for(sp1 in 1:vcount(current_g)){
    for(sp2 in 1:vcount(current_g)){
      if(sp1 <= sp2) next
      
      if(adj_mat[sp1,sp2] != 0 & adj_mat[sp2,sp1] != 0){
        competition <- competition + (adj_mat[sp1,sp2] + adj_mat[sp2,sp1])
      }else if(adj_mat[sp1,sp2] == 0 & adj_mat[sp2,sp1] != 0) amensalism <- amensalism + adj_mat[sp2,sp1]
      else if(adj_mat[sp1,sp2] != 0 & adj_mat[sp2,sp1] == 0) amensalism <- amensalism + adj_mat[sp1,sp2]
      
    }
  }
  
  cur_out <- cbind(cur_out, data.frame(fr_competition_quant=(competition/(sum(adj_mat))), fr_amensalism_quant=(amensalism/(sum(adj_mat)))))
  
  if(is.null(output_neg_nti)) output_neg_nti <- cur_out
  else output_neg_nti <- rbind(output_neg_nti, cur_out)
}

output_neg_nti$lat <- -1*output_neg_nti$lat
row.names(output_neg_nti) <- as.character(output_neg_nti$site)

ggplot(output_neg_nti, aes(lat, fr_competition)) + 
  geom_point() + theme_bw() + geom_smooth(method='lm', se=F) +
  xlab('latitude')

############ This code is for assigning the full names to species in order to create the
############ interactions matrix and save it to a csv file
# V(g_nt_neg)$name <- species[match(V(g_nt_neg)$name, species$name),]$corrected_name
# write.csv(as_adjacency_matrix(g_nt_neg, sparse = FALSE), file = 'non-trophic-negative-network.csv')

## Lastly, we build and analyse the non-trophic positive interactions networks along our latitudinal gradient
nt_pos_ints <- read.csv('chilean_NTIpos.txt', sep = '\t', stringsAsFactors = F)
nt_pos_ints <- nt_pos_ints[,-1]
names(nt_pos_ints)[1] <- 'name'
names(nt_pos_ints)[2:(dim(nt_pos_ints)[2])] <- nt_pos_ints$name

rownames(nt_pos_ints) <- nt_pos_ints$name
nt_pos_ints <- nt_pos_ints[,-1]

species_in_network <- names(nt_pos_ints)
setdiff(unique(species_in_network), unique(species$name))

#we can create a graph from the adjacency matrix of the network
nt_pos_ints  <- as.matrix(nt_pos_ints)
g_nt_pos <- graph.adjacency(nt_pos_ints)
V(g_nt_pos)$code <- as.character(species$code[match(V(g_nt_pos)$name, species$name)])

### we remove the nodes that do not have any interactions
g_nt_pos <- induced.subgraph(g_nt_pos, V(g_nt_pos)[which(degree(g_nt_pos) != 0)])
S_pos <- vcount(g_nt_pos)
L_pos <- ecount(g_nt_pos)

par(mar=c(0,0,0,0))
plot.igraph(g_nt_pos, vertex.label=NA, layout = layout.random, edge.arrow.size=.5, edge.curved=T, vertex.size=7)

## this is the data frame where we store the output
output_pos_nti <- NULL
for(site in ordered_sites$site){
  current_comm <- subset(comms_sp_networks, SITE == site)
  lat <- current_comm$SOUTHLAT[1]
  temp_code <- current_comm$TEMP_CODE[1]
  samples_no <- dim(current_comm)[1]
  
  abs <- current_comm[,which((names(current_comm)) %in% species$code)]
  print(paste('site = ', site, ' - samples ', dim(current_comm)[1]))
  
  sp_local_network <- names(abs)[which(colSums(abs) != 0)]
  sp_local_network <- as.character(species[which(species$code %in% sp_local_network),]$name)
  sp_local_network <- intersect(V(g_nt_pos)$name, sp_local_network)
  current_g <- induced.subgraph(g_nt_pos, sp_local_network)
  
  #### now we remove the nodes that have 0 degree
  current_g <- induced.subgraph(current_g, V(current_g)[which(degree(current_g) != 0)])
  
  S <- vcount(current_g)
  L <- ecount(current_g)
  L.S <- L/S
  C <- L/((S**2)-S)
  
  # this is the adjacency matrix
  M <- as.matrix(get.adjacency(current_g))
  
  indeg <- sum(colSums(M))/sum(colSums(M) != 0)
  outdeg <-  sum(rowSums(M))/sum(rowSums(M) != 0)
  
  sd_gen <- SDGenerality(M)
  sd_vul <- SDVulnerability(M)
  
  mod <- tryCatch({
    cluster_louvain(as.undirected(current_g))
  }, warning = function(w){
    NA
  }, error = function(e){
    NA
  }, finally = {
    
  })
  
  if(is.na(mod)){
    modularity <- -1
    n_modules <- 0
  }else{
    modularity <- max(mod$modularity)
    n_modules <- max(mod$membership)
    V(current_g)$module <- mod$membership
  }
  
  # uncomment the line below if you want to save the local networks to disk : in graphml format
  #write.graph(current_g, file=paste0('./networks/nti-positive-net-',site,'.graphml'), format = 'graphml')
  
  mutualism <- 0
  comensalism <- 0
  adj_mat <- get.adjacency(current_g, sparse=F)
  for(sp1 in 1:vcount(current_g)){
    for(sp2 in 1:vcount(current_g)){
      if(sp1 <= sp2) next
      if(adj_mat[sp1,sp2] != 0 & adj_mat[sp2,sp1] != 0) mutualism <- mutualism + 1
      else if(adj_mat[sp1,sp2] == 0 & adj_mat[sp2,sp1] != 0) comensalism <- comensalism + 1
      else if(adj_mat[sp1,sp2] != 0 & adj_mat[sp2,sp1] == 0) comensalism <- comensalism + 1
      
    }
  }
  
  ###### this bit was added to get the quantified version of the network and its corresponding properties
  g_edges <- ends(current_g, E(current_g))
  comm_size <- dim(abs)[1]
  edges_weights <- c()
  for(e in 1:dim(g_edges)[1]){
    prey <- as.character(species[which(species$name == g_edges[e,][1]),]$code)
    predator <- as.character(species[which(species$name == g_edges[e,][2]),]$code)
    co_occur <- (length(which(rowSums(abs[prey]) != 0 )) / comm_size) * (length(which(rowSums(abs[predator])!= 0 )) / comm_size)
    edges_weights <- append(edges_weights, co_occur)
  }
  current_g <- set.edge.attribute(current_g, 'interaction.probability', value=edges_weights)
  ###### end of network links quantification
  
  undirec_g <- as.undirected(current_g, edge.attr.comb = 'first')
  mod_quant <- tryCatch({
    cluster_louvain(undirec_g, weights=E(undirec_g)$interaction.probability)
  }, warning = function(w){
    NA
  }, error = function(e){
    NA
  }, finally = {
    
  })
  
  if(is.na(mod_quant)){
    modularity_quant <- -1
    n_modules_quant <- 0
  }else{
    modularity_quant <- max(mod_quant$modularity)
    n_modules_quant <- max(mod_quant$membership)
  }
  
  ###### now that we have quantified interactions we can use cheddar to calculate Bersier's properties
  M_Q <- get.adjacency(current_g, attr='interaction.probability', sparse=FALSE)
  cheddar_network <- cheddar::Community(nodes=data.frame(node=rownames(M_Q)), trophic.links=PredationMatrixToLinks(M_Q, link.property = 'interaction.probability'), properties=list(title=site));
  quants <- t(QuantitativeDescriptorsNoChain(cheddar_network, 'interaction.probability'))
  
  cur_out <- data.frame(site, lat, samples_no, S, L, L.S, C, indeg, outdeg, sd_gen, sd_vul, modularity, n_modules, modularity_quant, n_modules_quant, mutualism, comensalism, fr_mutualism=mutualism/(L-mutualism), fr_comensalism=comensalism/(L-mutualism))
  
  ##### here we add the results of the quantitative analysis to the output data frame
  temp <- quants[1,]
  names(temp) <- paste0(names(temp),' qualitative')
  cur_out <- cbind(cur_out, t(temp))
  temp <- quants[2,]
  names(temp) <- paste0(names(temp),' unweighted')
  cur_out <- cbind(cur_out, t(temp))
  temp <- quants[3,]
  names(temp) <- paste0(names(temp),' weighted')
  cur_out <- cbind(cur_out, t(temp))
  
  ###### that's it
  if(is.null(output_pos_nti)) output_pos_nti <- cur_out
  else output_pos_nti <- rbind(output_pos_nti, cur_out)
}


output_pos_nti$lat <- -1*output_pos_nti$lat
row.names(output_pos_nti) <- as.character(output_pos_nti$site)

ggplot(output_pos_nti, aes(lat, fr_comensalism)) + 
  geom_point() + theme_bw() + geom_smooth(method='lm', se=F) + 
  xlab('latitude')

############ This code is for assigning the full names to species in order to create the
############ interactions matrix
# V(g_nt_pos)$name <- species[match(V(g_nt_pos)$name, species$name),]$corrected_name
# write.csv(as_adjacency_matrix(g_nt_pos, sparse = FALSE), file = 'non-trophic-positive-network.csv')

#### Now that we have the networks for all the different interaction types we build the networks
#### comprising all interaction types. Interestingly, this procedure does not produce multiple edges
#### i.e. there are no two species with the same directionality in their interaction across the three networks

## first, we build the whole interaction network by merging the three networks above: trophic, non-trophic
## negative, and non-trophic positive

g_whole <- igraph::union(g_trophic, g_nt_neg, g_nt_pos)

## this is the data frame where we store the output
output_all_ints <- NULL
for(site in ordered_sites$site){
  current_comm <- subset(comms_sp_networks, SITE == site)
  lat <- current_comm$SOUTHLAT[1]
  temp_code <- current_comm$TEMP_CODE[1]
  samples_no <- dim(current_comm)[1]
  
  abs <- current_comm[,which((names(current_comm)) %in% species$code)]
  print(paste('site = ', site, ' - samples ', dim(current_comm)[1]))
  
  sp_local_network <- names(abs)[which(colSums(abs) != 0)]
  sp_local_network <- as.character(species[which(species$code %in% sp_local_network),]$name)
  sp_local_network <- intersect(V(g_whole)$name, sp_local_network)
  current_g <- induced.subgraph(g_whole, sp_local_network)
  
  #### now we remove the nodes that have 0 degree
  current_g <- induced.subgraph(current_g, V(current_g)[which(degree(current_g) != 0)])
  
  S <- vcount(current_g)
  L <- ecount(current_g)
  L.S <- L/S
  C <- L/((S**2)-S)
  
  # this is the adjacency matrix
  M <- as.matrix(get.adjacency(current_g))
  
  indeg <- sum(colSums(M))/sum(colSums(M) != 0)
  outdeg <-  sum(rowSums(M))/sum(rowSums(M) != 0)
  
  mod <- tryCatch({
    cluster_louvain(as.undirected(current_g))
  }, warning = function(w){
    NA
  }, error = function(e){
    NA
  }, finally = {
    
  })
  
  if(is.na(mod)){
    modularity <- -1
    n_modules <- 0
  }else{
    modularity <- max(mod$modularity)
    n_modules <- max(mod$membership)
    V(current_g)$module <- mod$membership
  }
  
  average_path <- average.path.length(current_g)
  
  ###### this bit was added to get the quantified version of the network and its corresponding properties
  g_edges <- ends(current_g, E(current_g))
  comm_size <- dim(abs)[1]
  edges_weights <- c()
  for(e in 1:dim(g_edges)[1]){
    prey <- as.character(species[which(species$name == g_edges[e,][1]),]$code)
    predator <- as.character(species[which(species$name == g_edges[e,][2]),]$code)
    co_occur <- (length(which(rowSums(abs[prey]) != 0 )) / comm_size) * (length(which(rowSums(abs[predator])!= 0 )) / comm_size)
    edges_weights <- append(edges_weights, co_occur)
  }
  current_g <- set.edge.attribute(current_g, 'interaction.probability', value=edges_weights)
  ###### end of network links quantification
  
  undirec_g <- as.undirected(current_g, edge.attr.comb = 'first')
  mod_quant <- tryCatch({
    cluster_louvain(undirec_g, weights=E(undirec_g)$interaction.probability)
  }, warning = function(w){
    NA
  }, error = function(e){
    NA
  }, finally = {
    
  })
  
  if(is.na(mod_quant)){
    modularity_quant <- -1
    n_modules_quant <- 0
  }else{
    modularity_quant <- max(mod_quant$modularity)
    n_modules_quant <- max(mod_quant$membership)
  }
  
  ###### now that we have quantified interactions we can use cheddar to calculate Bersier's properties
  M_Q <- get.adjacency(current_g, attr='interaction.probability', sparse=FALSE)
  cheddar_network <- cheddar::Community(nodes=data.frame(node=rownames(M_Q)), trophic.links=PredationMatrixToLinks(M_Q, link.property = 'interaction.probability'), properties=list(title=site));
  quants <- t(QuantitativeDescriptorsNoChain(cheddar_network, 'interaction.probability'))
  
  cur_out <- data.frame(site, lat, samples_no, S, L, L.S, C, indeg, outdeg, modularity, n_modules, modularity_quant, n_modules_quant, average_path)
  
  ##### here we add the results of the quantitative analysis to the output data frame
  temp <- quants[1,]
  names(temp) <- paste0(names(temp),' qualitative')
  cur_out <- cbind(cur_out, t(temp))
  temp <- quants[2,]
  names(temp) <- paste0(names(temp),' unweighted')
  cur_out <- cbind(cur_out, t(temp))
  temp <- quants[3,]
  names(temp) <- paste0(names(temp),' weighted')
  cur_out <- cbind(cur_out, t(temp))
  
  ###### that's it
  if(is.null(output_all_ints)) output_all_ints <- cur_out
  else output_all_ints <- rbind(output_all_ints, cur_out)
}


output_all_ints$lat <- -1*output_all_ints$lat
row.names(output_all_ints) <- as.character(output_all_ints$site)

