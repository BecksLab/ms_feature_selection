
## ---------------------------
##
## Script name: auxiliary-functions.r
##
## Purpose of script: This script provides a series of auxiliary functions
## used in different parts of the other scripts.
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


## function for downloading and formatting satellite image map for Sea Surface Temperature off the central coast of Chile
my_get_sst_chile <- function(year) {
  require(stringr)
  require(ncdf4)
  # base url with YEAR meant to be replaced
  turl <- "http://coastwatch.pfeg.noaa.gov/erddap/griddap/erdMBsstd8day.nc?sst[(YEAR-12-12T00:00:00Z)][(0.0)][(-40.0):(-25)][(280.0):(295.0)]"
  
  # the URL with YEAR replaced with a value
  the_url <- str_replace(turl, "YEAR", year)
  
  # the filename to save the downloaded data in
  the_file <- paste("sst_", year, ".nc", sep="")
  
  # if the file isn't here, download it
  if(!file.exists(the_file)) {
    download.file(the_url, the_file, mode='wb')
  } else {
    message(paste("Using existing file", the_file, collapse = " "))
  }
  
  # now, grab stuff out of the netcdf file and return it in a list
  # called ret
  sstFile <- nc_open(the_file)
  ret <- list()
  ret$lats <- ncvar_get(sstFile, "latitude")
  ret$lon <- ncvar_get(sstFile, "longitude") - 360 # we need them as negative values
  ret$time <- ncvar_get(sstFile, "time")
  ret$sst <- ncvar_get(sstFile, "sst")
  ret$date <- paste("12-12-", year, sep = "")
  nc_close(sstFile)
  
  ret
}

## melting data function for the map
melt_sst <- function(L) {
  require(reshape2)
  dimnames(L$sst) <- list(long = L$lon, lat = L$lats)
  ret <- melt(L$sst, value.name = "sst")
  cbind(date = L$date, ret, degF = ret$sst * 9/5 + 32)
}


############ The next two functions were copied and modified from the cheddar package
############ to tailor their output to our needs

NodeQuantitativeDescriptorsNoChain <- function (community, weight) {
  if (!is.Community(community)) 
    stop("Not a Community")
  # .RequireTrophicLinks(community)
  vlog2v <- function(v) v * log2(v)
  b <- PredationMatrix(community, weight)
  bIn <- colSums(b)
  bOut <- rowSums(b)
  HN <- -rowSums(vlog2v(t(b)/bIn), na.rm = TRUE)
  HP <- -rowSums(vlog2v(b/bOut), na.rm = TRUE)
  nN <- 2^HN
  nN[0 == bIn] <- 0
  nP <- 2^HP
  nP[0 == bOut] <- 0
  d.prime <- nN/(nN + nP)
  d <- bIn * nN/(bIn * nN + bOut * nP)
  # t <- LongestTrophicLevel(community)
  res <- ResourcesByNode(community)
  # o.prime <- -1 + 2^sapply(1:NumberOfNodes(community), function(k) {
  #   r <- res[[k]]
  #   n <- table(t[r])
  #   return(-sum(vlog2v(n/sum(n))))
  # })
  # o <- -1 + 2^sapply(1:NumberOfNodes(community), function(k) {
  #   r <- res[[k]]
  #   bt <- tapply(b[r, k], t[r], sum)
  #   return(-sum(vlog2v(bt/sum(bt))))
  # })
  g.prime <- nN * NumberOfNodes(community)/sum(nN)
  g <- bIn * nN * NumberOfNodes(community)/sum(bIn * nN)
  v.prime <- nP * NumberOfNodes(community)/sum(nP)
  v <- bOut * nP * NumberOfNodes(community)/sum(bOut * nP)
  res <- cbind(NResources = cheddar::NumberOfResources(community), 
               NConsumers = cheddar::NumberOfConsumers(community), bIn, bOut, 
               nN, nP, d.prime, d, g.prime, g, v.prime, v)
  rownames(res) <- unname(NP(community, "node"))
  return(res)
}

QuantitativeDescriptorsNoChain <- function (community, weight, top.level.threshold = 0.99) {
  vlog2v <- function(v) v * log2(v)
  np <- NodeQuantitativeDescriptorsNoChain(community, weight)
  b <- PredationMatrix(community, weight)
  sumb <- sum(b)
  HN <- -rowSums(vlog2v(t(b)/np[, "bIn"]), na.rm = TRUE)
  HP <- -rowSums(vlog2v(b/np[, "bOut"]), na.rm = TRUE)
  fracT.q.prime <- mean(np[, "d.prime"] >= top.level.threshold)
  fracI.q.prime <- mean(0 < np[, "d.prime"] & np[, "d.prime"] < top.level.threshold)
  fracB.q.prime <- mean(0 == np[, "d.prime"])
  fracT.q <- mean(np[, "d"] >= top.level.threshold)
  fracI.q <- mean(0 < np[, "d"] & np[, "d"] < top.level.threshold)
  fracB.q <- mean(0 == np[, "d"])
  NP.q.prime <- 2^(-sum(vlog2v(np[, "nP"]/sum(np[, "nP"])), na.rm = TRUE))/2^(-sum(vlog2v(np[, "nN"]/sum(np[, "nN"])), na.rm = TRUE))
  NP.q <- 2^(-sum(vlog2v((np[, "bOut"] * np[, "nP"])/sum(np[, "bOut"] * np[, "nP"])), na.rm = TRUE))/2^(-sum(vlog2v((np[, "bIn"] * np[, "nN"])/sum(np[, "bIn"] * np[, "nN"])), na.rm = TRUE))
  LD.q.prime <- (sum(np[, "nP"]) + sum(np[, "nN"]))/(2 * NumberOfNodes(community))
  LD.q <- (sum(np[, "bOut"] * np[, "nP"]/sumb, na.rm = TRUE) + sum(np[, "bIn"] * np[, "nN"]/sumb, na.rm = TRUE))/2
  C.q.prime <- LD.q.prime/NumberOfNodes(community)
  C.q <- LD.q/NumberOfNodes(community)
  Phi <- sum(HP * np[, "bOut"]/sumb) + sum(HN * np[, "bIn"]/sumb)
  m <- 2^(Phi/2)
  PhiAB <- function(A, B) {
    A <- A(community)
    B <- B(community)
    bOut <- rowSums(b)[A]
    bIn <- colSums(b)[B]
    return(sum((bOut/sumb) * -vlog2v(b[A, B]/bOut), na.rm = TRUE) + 
             sum((bIn/sumb) * -vlog2v(t(b[A, B])/bIn), na.rm = TRUE))
  }
  fracTI.q <- PhiAB(IntermediateNodes, TopLevelNodes)/Phi
  fracTB.q <- PhiAB(BasalNodes, TopLevelNodes)/Phi
  fracII.q <- PhiAB(IntermediateNodes, IntermediateNodes)/Phi
  fracIB.q <- PhiAB(BasalNodes, IntermediateNodes)/Phi
  PhiPrime <- sum(HP/NumberOfNodes(community)) + sum(HN/NumberOfNodes(community))
  PhiABPrime <- function(A, B) {
    A <- A(community)
    B <- B(community)
    bOut <- rowSums(b)[A]
    bIn <- colSums(b)[B]
    s <- NumberOfNodes(community)
    return(sum((1/s) * -vlog2v(b[A, B]/bOut), na.rm = TRUE) + 
             sum((1/s) * -vlog2v(t(b[A, B])/bIn), na.rm = TRUE))
  }
  fracTI.q.prime <- PhiABPrime(IntermediateNodes, TopLevelNodes)/PhiPrime
  fracTB.q.prime <- PhiABPrime(BasalNodes, TopLevelNodes)/PhiPrime
  fracII.q.prime <- PhiABPrime(IntermediateNodes, IntermediateNodes)/PhiPrime
  fracIB.q.prime <- PhiABPrime(BasalNodes, IntermediateNodes)/PhiPrime
  tlps <- TLPS(community, link.properties = weight)
  
  nT <- length(TopLevelNodes(community))
  nI <- length(IntermediateNodes(community))
  nB <- length(BasalNodes(community))
  G.q.prime <- sum(np[, "nN"])/(nT + nI)
  G.q <- sum(np[, "nN"] * np[, "bIn"]/sumb, na.rm = TRUE)
  V.q.prime <- sum(np[, "nP"])/(nI + nB)
  V.q <- sum(np[, "nP"] * np[, "bOut"]/sumb, na.rm = TRUE)
  tlps <- TLPS(community, node.properties = c("IsTopLevelNode", "IsIntermediateNode", "IsBasalNode"))
  fracTI <- with(tlps, sum(resource.IsIntermediateNode & consumer.IsTopLevelNode))/nrow(tlps)
  fracTB <- with(tlps, sum(resource.IsBasalNode & consumer.IsTopLevelNode))/nrow(tlps)
  fracII <- with(tlps, sum(resource.IsIntermediateNode & consumer.IsIntermediateNode))/nrow(tlps)
  fracIB <- with(tlps, sum(resource.IsBasalNode & consumer.IsIntermediateNode))/nrow(tlps)
  Qualitative <- c(FractionTopLevelNodes(community), FractionIntermediateNodes(community), 
                   FractionBasalNodes(community), 
                   sum(cheddar::NumberOfConsumers(community) > 0)/sum(cheddar::NumberOfResources(community) > 0), 
                   LinkageDensity(community), 
                   DirectedConnectance(community), fracTI, fracTB, fracII, 
                   fracIB,
                   mean(cheddar::TrophicGenerality(community)[cheddar::NumberOfResources(community) > 0]), 
                   mean(cheddar::TrophicVulnerability(community)[cheddar::NumberOfConsumers(community) > 0]), 
                   sd(cheddar::NormalisedTrophicGenerality(community)), sd(cheddar::NormalisedTrophicVulnerability(community)))
  Unweighted <- c(fracT.q.prime, fracI.q.prime, fracB.q.prime, 
                  NP.q.prime, LD.q.prime, C.q.prime, fracTI.q.prime, fracTB.q.prime, 
                  fracII.q.prime, fracIB.q.prime, 
                  G.q.prime, V.q.prime, sd(np[, "g.prime"]), sd(np[, "v.prime"]))
  Weighted <- c(fracT.q, fracI.q, fracB.q, NP.q, LD.q, C.q, 
                fracTI.q, fracTB.q, fracII.q, fracIB.q,
                G.q, V.q, sd(np[, "g"]), sd(np[, "v"]))
  res <- cbind(Qualitative, Unweighted, Weighted)
  rownames(res) <- c("Fraction top level", "Fraction intermediate", 
                     "Fraction basal", "Ratio resources:consumers", "Link density", 
                     "Connectance", "Fraction links top:intermediate", "Fraction links top:basal", 
                     "Fraction links intermediate:intermediate", "Fraction links intermediate:basal", 
                     "Generality", "Vulnerability", 
                     "SD standardised generality", "SD standardised vulnerability")
  return(res)
}

##### End of modified cheddar functions

