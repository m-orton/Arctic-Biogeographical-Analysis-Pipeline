###############
# Chironomid Biogeography Pipeline 
# (Using Subarctic Shapefile Filter)

# Authored by Matthew G. Orton and Sally J. Adamowicz

# Credit to Torbjorn Ekrem and Elisabeth Ster for private Chironomidae data and some of the public
# Chironomidae data from BOLD and for helping us on the function of the code and plots.

##############
# Packages

# install.packages("foreach")
library(foreach)
# install.packages("ape")
library(ape)
# read_tsv function.
# install.packages("readr")
library(readr)
# source("https://bioconductor.org/biocLite.R")
# biocLite("Biostrings")
# biocLite("muscle")
# biocLite("DECIPHER")
library(DECIPHER)
library(Biostrings) 
library(muscle)
# install.packages("plotly")
library(plotly)
# install.packages ("ggplot2") 
require(ggplot2)
# install.packages("raster")
library(raster)
# install.packages("rgdal")
library(rgdal)
# install.packages("rgeos")
library(rgeos)
# install.packages("vegan")
library(vegan)
# install.packages("tidyr")
library(tidyr)
# install.packages("dplyr")
library(dplyr)
# install.packages("data.table")
library(data.table)
# install.packages("vegan")
library(vegan)
# install.packages("stringr")
library(stringr)

##############
# Parsing from BOLD

# Public records for each of the three regions - for most recent data from BOLD:
# dfNearctic <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Chironomidae&geo=Alaska|Canada&format=tsv")
# dfGreenland <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Chironomidae&geo=Greenland&format=tsv")
# dfPalearctic <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Chironomidae&geo=Norway|Denmark|Iceland|Sweden|Finland&format=tsv")

# To read in datasets from initial download on Oct 2017:
dfGreenland <- read.csv("dfGreenland_Oct17.csv")
dfGreenland <- dfGreenland[,2:81]
# Must unzip Nearctic csv file first
dfNearctic <- read.csv("dfNearctic_Oct17.csv")
dfNearctic <- dfNearctic[,2:81]
dfPalearctic <- read.csv("dfPalearctic_Oct17.csv")
dfPalearctic <- dfPalearctic[,2:81]

# Note that BOLD makes the distinction between Greenland and Denmark and distinguishes them as separate countries
# (even though they are not) so the datasets between Denmark and Greenland are nonoverlapping

##############
# Record filtering - lat/lon coordinates, BIN, presence of a sequence, COI-5P

# Labeling regions
dfNearctic$globalRegion <- "Nearctic"
dfGreenland$globalRegion <- "Greenland"
dfPalearctic$globalRegion <- "Palearctic"

# Combine dataframes with regional identifiers
dfChironomid <- rbind(dfNearctic, dfGreenland, dfPalearctic)

# Filter for presence of BIN assignment (grep by colon since all BIN identifiers have this):
containBin <- grep( "[:]", dfChironomid$bin_uri)
dfChironomidAll <- dfChironomid[containBin,]

# Filter out BINs without sequence data since we need sequence data for determining outlier sequences:
containNucleotides <- grep( "[ACGT]", dfChironomid$nucleotides)
dfChironomid <- dfChironomid[containNucleotides,]

# Filter out BINs without coordinate data:
containLatLon <- grep( "[0-9]", dfChironomid$lat)
dfChironomid <- dfChironomid[containLatLon,]

# Filter according to COI-5P
containCOI <- grep( "^CO", dfChironomid$markercode)
dfChironomid <- dfChironomid[containCOI,]

# Can use this command to check to make sure all markers are COI-5P
unique(dfChironomid$markercode)

# Conversion to numeric for lat and lon values:
latNum <- with(dfChironomid, as.numeric(as.character(lat))) 
dfChironomid$latNum <- latNum
lonNum <- with(dfChironomid, as.numeric(as.character(lon))) 
dfChironomid$lonNum <- lonNum

##############
# Reading of Private Data and Combining with Public Data

# Filtering is redone for private records due to differences in column names between public and private

# Private dataset used from April 24th
dfPrivateData <- read_csv("Private_Chironomid_Data_ModifiedSingleSheet.csv")

# Read in the sequence data for the private sequence dataset
privateSeqs <- readDNAStringSet("PrivateSequenceData.fas")
privateSeqNames <- as.character(privateSeqs@ranges@NAMES)
privateSeqNames2 <- strsplit(privateSeqNames, "[|]")
privateSeqNames2 <-  sapply( privateSeqNames2, "[", 1 )
nucleotides <- unname(as.character(privateSeqs))
dfPrivateData2 <- as.data.frame(nucleotides)
dfPrivateData2$Process_id <- privateSeqNames2
dfPrivateData <- merge(dfPrivateData, dfPrivateData2, by.x = "Process ID", by.y = "Process_id")

# Make all columns char type for all dataframes (makes df manipulations easier)
dfPrivateData <- data.frame(lapply(dfPrivateData, as.character), stringsAsFactors=FALSE)

# Filter by BIN assignment (grep by colon since all BIN identifiers have this):
containBin2 <- grep( "[:]", dfPrivateData$BIN)
dfPrivateData <- dfPrivateData[containBin2,]

# Filter out BINs without coordinate data:
containLatLon2 <- grep( "[0-9]", dfPrivateData$Lat)
dfPrivateData <- dfPrivateData[containLatLon2,]

# Check the records for which countries are included
unique(dfPrivateData$Country.Ocean)

# Only Canada, Greenland, Iceland and Norway were found from Private data

# Assigning regions: Greenland, Nearctic or Palearctic
for(i in seq(from = 1, to = nrow(dfPrivateData), by = 1)) {
  if(dfPrivateData$Country.Ocean[i] == "Canada") {
    dfPrivateData$globalRegion[i] <- "Nearctic"
  } else if(dfPrivateData$Country.Ocean[i] == "Greenland"){
    dfPrivateData$globalRegion[i] <- "Greenland"
  } else if(dfPrivateData$Country.Ocean[i] == "Iceland"){
    dfPrivateData$globalRegion[i] <- "Palearctic"
  } else if(dfPrivateData$Country.Ocean[i] == "Norway"){
    dfPrivateData$globalRegion[i] <- "Palearctic"
  }
}

# Check if there is any intersection between dfPrivate data and dfChironomid for sample ID
sampleIdIntersect <- intersect(dfPrivateData$Sample.ID, dfChironomid$sampleid)

# Process ID
processIdIntersect <- intersect(dfPrivateData$Process.ID, dfChironomid$processid)

# Same number of elements in both - 2418, will subset dfChironomid for the duplicate records
# by sampleID
duplicateSubset <- which(dfChironomid$sampleid %in% sampleIdIntersect)
dfChironomid <- dfChironomid[-duplicateSubset,]

# Conversion to numeric coordinates for private data:
latNum <- with(dfPrivateData, as.numeric(as.character(Lat))) 
dfPrivateData$latNum <- latNum
lonNum <- with(dfPrivateData, as.numeric(as.character(Lon))) 
dfPrivateData$lonNum <- lonNum

# Combining private data with dfChironomid and excluding unecessary columns for the analysis
# and mapping
colnames(dfPrivateData)[1] <- "processid"
colnames(dfPrivateData)[5] <- "subfamily_name"
colnames(dfPrivateData)[8] <- "species_name"
colnames(dfPrivateData)[25] <- "collectors"
colnames(dfPrivateData)[27] <- "country"
colnames(dfPrivateData)[28] <- "province_state"
colnames(dfPrivateData)[29] <- "region"
colnames(dfPrivateData)[30] <- "sector"
colnames(dfPrivateData)[31] <- "exactsite"
colnames(dfPrivateData)[48] <- "bin_uri"

# Ensuring same column headings for both private and public data
dfPrivateData <- (dfPrivateData[,c("processid","globalRegion","bin_uri","species_name","subfamily_name","latNum","lonNum",
                                   "country","province_state","region","sector","exactsite","collectors","nucleotides")])

dfChironomid <- (dfChironomid[,c("processid","globalRegion","bin_uri","species_name","subfamily_name","latNum","lonNum",
                                 "country","province_state","region","sector","exactsite","collectors","nucleotides")])

# Combine both dataframes together
dfChironomidAll <- rbind(dfPrivateData, dfChironomid)

# *** Upon checking with Elisabeth, the only BIN that we should eliminate is ACZ1013
# as it was an outlier BINs
binCheck <- which(dfChironomidAll$bin_uri == "BOLD:ACZ1013")
dfChironomidAll <- dfChironomidAll[-binCheck,]

##############
# Subarctic filtering according to subarctic shapefile

# Subarctic filtering of all data including private data
# Read in the subarctic shapefile from the CAFF website:
# http://geo.abds.is/geonetwork/srv/eng/catalog.search#/metadata/2ad7a7cb-2ad7-4517-a26e-7878ef134239
rawShapefile <- shapefile("Arctic_Zones")
subarcticZone <- rawShapefile[3]

# Loading required package sp to make a spatial points dataframe from our coordinate data
library(sp)
xy <- data.frame(ID = dfChironomidAll$bin_uri, X = dfChironomidAll$lonNum, Y = dfChironomidAll$latNum)
coordinates(xy) <- c("X", "Y")
proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")
# Have to convert to epsg:3408 - its a northern map view with UTM coordinates
convertUTM <- spTransform(xy, CRS("+init=epsg:3408"))

# Transform our shapefile
subarcticZone <- spTransform(subarcticZone, CRSobj = "+init=epsg:3408")

# Projection of each spatial dataframe
projection(convertUTM) 
projection(subarcticZone) 

# Find which points overlap, commands may take a while:
pointOverlap <- sp::over(convertUTM, subarcticZone, fn = NULL)
pointOverlaprgeos <- rgeos::gIntersection(convertUTM, subarcticZone)

# Plot points in polygon, points within the polygon plot green
plot(subarcticZone)
plot(convertUTM, pch = 19, cex = 1, add = TRUE)
plot(pointOverlaprgeos, pch = 19, cex = 0.5, col = 'green', add = TRUE)
box()

# Append pointOverlap to dfChironomidAll
dfChironomidAll$shapeLen <- pointOverlap$Shape_Leng
withinPoly <- grep( "[0-9]", dfChironomidAll$shapeLen)
# Filter by subarctic zone!
dfChironomidAll <- dfChironomidAll[withinPoly,]

##############
# Selecting One Sequence per BIN 
# (for SuperBIN clustering only)

# New dataframe for this section only
dfChironomidFilter <- dfChironomidAll

# Make a list of all BINs
binList <- lapply(unique(dfChironomidFilter$bin_uri), function(x) 
  dfChironomidFilter[dfChironomidFilter$bin_uri == x,])
  
# Also need to find the number of unique BINs in the binlist 
binNumber<- unique(dfChironomidFilter$bin_uri)
binNumber <- length(binNumber)
  
# Extract record id from each BIN
binRecordId <- foreach(i=1:binNumber) %do% unique(binList[[i]]$processid)

# Count sequence length per BIN
seqLengthBIN <- foreach(i=1:binNumber) %do% nchar(binList[[i]]$nucleotides)

# Name seqLengthBIN with the record ids
for (i in seq(from = 1, to = binNumber, by = 1)) {
  names(seqLengthBIN[[i]]) <- binRecordId[[i]]
}
# Which sequence is closest to 658 bp
binSelect <- foreach(i=1:binNumber) %do% 
  which(abs(seqLengthBIN[[i]]-658)==min(abs(seqLengthBIN[[i]]-658)))
# If multiple sequences then pick first
binSelect <- foreach(i=1:binNumber) %do% head(seqLengthBIN[[i]], 1)

# Unlist, keep name and retype
binSelect <- unlist(binSelect)
binSelect <- names(binSelect)
binSelect <- as.character(binSelect)

# Subset dfChironomidFilter by processids selected as representatives
dfSingleSeq <- subset(dfChironomidFilter, processid %in% binSelect)
dfSingleSeq$seqLength <- nchar(dfSingleSeq$nucleotides)
dfSingleSeq <- dfSingleSeq[order(dfSingleSeq[,'processid'],-dfSingleSeq[,'seqLength']),]
dfSingleSeq <- dfSingleSeq[!duplicated(dfSingleSeq$processid),]

##############
# Single-linkage Clustering using the Decipher Package

# Alignment step
# dnaStringSet <- DNAStringSet(dfSingleSeq$nucleotides)
# alignment <- muscle(dnaStringSet, maxiters = 2, diags = TRUE)
# dnaStringSet2 <- DNAStringSet(alignment)

# Name the stringset with record ids
# bin_uri <- dfSingleSeq$bin_uri
# names(dnaStringSet2) <- bin_uri

# Write out to fasta
# fileName <- paste("ChironomidAlignmentApril24_Subarctic.fas")
# writeXStringSet(dnaStringSet2, file=fileName, format = "fasta", width = 658)

# Alignment using subarctic filtered chironomid data
alignmentSubarctic <- readDNAStringSet("ChironomidAlignmentApril24_Subarctic.fas")

# DNAbin format
dnaBin <- as.DNAbin(alignmentSubarctic)
# Distance matrix using TN93 before clustering 
distanceMatrix <- dist.dna(dnaBin, model = "TN93", as.matrix = TRUE, 
  pairwise.deletion = TRUE)

# Clustering according to 4% divergence threshold
clustSingle4 <- IdClusters(distanceMatrix,
                          method = "single",
                          cutoff= 0.04,
                          showPlot = TRUE,
                          type = "both",
                          processors = 2,
                          verbose = TRUE)

# Number of unique clusters for 4%
length(unique(clustSingle4$cluster))

# 4.5% divergence threshold
clustSingle45 <- IdClusters(distanceMatrix,
                           method = "single",
                           cutoff= 0.045,
                           showPlot = TRUE,
                           type = "both",
                           processors = 2,
                           verbose = TRUE)

# Number of unique clusters for 4.5%
length(unique(clustSingle45$cluster))

# 5% divergence threshold
clustSingle5 <- IdClusters(distanceMatrix,
                            method = "single",
                            cutoff= 0.05,
                            showPlot = TRUE,
                            type = "both",
                            processors = 2,
                            verbose = TRUE)

# Number of unique clusters for 5%
length(unique(clustSingle5$cluster))

# Renaming column for 4.5% cluster
clustSingle45 <- setDT(clustSingle45, keep.rownames = TRUE)[]
colnames(clustSingle45)[2] <- "cluster_45"

# Merge clusters to original Chironomid dataset (pre-alignment filter)
# Every record of every BIN will now have a defined "superBIN"!
dfChironomidSBIN <- merge(dfChironomidAll, clustSingle45, by.x ="bin_uri", by.y ="rn")

############
# Taxonomy Curation of Greenland (Elisabeth's revisions to Greenland species)
# For species level analyses

# Read in Elizabeths csv file for Greenland (modified to incorporate her revisions) -
# Certain species removed that were misclassified
dfSpeciesEdit <- read_csv("Greenland records chironomid (1).csv")
colnames(dfSpeciesEdit)[1] <- "col_1"
splitSpecies <- foreach(i=1:nrow(dfSpeciesEdit)) %do% strsplit(dfSpeciesEdit$col_1[i], ",")
# Convert to dataframe format
dfSplitSpecies <- do.call("rbind", lapply(splitSpecies, "[[", 1))
dfSplitSpecies <- as.data.frame(dfSplitSpecies)
dfSplitSpecies$species_names <- as.character(dfSplitSpecies$V3)

############
# Divide into separate regions (BINs) for downstream analysis

# If subsetting for private data only:
# privateSubset <- intersect(dfChironomidAll$processid, dfPrivateData$processid)
# dfChironomidAll <- subset(dfChironomidAll, processid %in% privateSubset)
# dfChironomidSBIN <- subset(dfChironomidSBIN, processid %in% privateSubset)

containGreenland <- which(dfChironomidAll$globalRegion=="Greenland")
dfGreenlandBIN <- dfChironomidAll[containGreenland,]

# Extract out bad species names for Greenland before using for species level analyses
# dfGreenlandBIN <- subset(dfGreenlandBIN, species_name %in% dfSplitSpecies$V3)

containNearctic <- which(dfChironomidAll$globalRegion=="Nearctic")
dfNearcticBIN <- dfChironomidAll[containNearctic,]

containPalearctic <-  which(dfChironomidAll$globalRegion=="Palearctic")
dfPalearcticBIN <- dfChironomidAll[containPalearctic,]

# Divide into separate regions (SBINs)
containGreenlandSBIN <- which(dfChironomidSBIN$globalRegion=="Greenland")
dfGreenlandSBIN <- dfChironomidSBIN[containGreenlandSBIN,]

containNearcticSBIN <- which(dfChironomidSBIN$globalRegion=="Nearctic")
dfNearcticSBIN <- dfChironomidSBIN[containNearcticSBIN,]

containPalearcticSBIN <-  which(dfChironomidSBIN$globalRegion=="Palearctic")
dfPalearcticSBIN <- dfChironomidSBIN[containPalearcticSBIN,]

#############
# Dplyr/Tidyr and Vegan Analyses - BIN, SBIN and Species level analyses

# First separate bin_uri and global region from other columns
dfNSubset_BIN <- (dfNearcticBIN[,c("globalRegion","bin_uri")])
dfPSubset_BIN <- (dfPalearcticBIN[,c("globalRegion","bin_uri")])
dfGSubset_BIN <- (dfGreenlandBIN[,c("globalRegion","bin_uri")])

# Or separate by SBIN - 4.5%
dfNSubset_SBIN <- (dfNearcticSBIN[,c("globalRegion","cluster_45")])
dfPSubset_SBIN <- (dfPalearcticSBIN[,c("globalRegion","cluster_45")])
dfGSubset_SBIN <- (dfGreenlandSBIN[,c("globalRegion","cluster_45")])

# Separate by species 
dfNSubset_Sp <- (dfNearcticSBIN[,c("globalRegion","species_name")])
dfPSubset_Sp <- (dfPalearcticSBIN[,c("globalRegion","species_name")])
dfGSubset_Sp <- (dfGreenlandSBIN[,c("globalRegion","species_name")])

# Group by BIN
nearcticGroup_BIN <- group_by(dfNSubset_BIN, bin_uri)
palearcticGroup_BIN <- group_by(dfPSubset_BIN, bin_uri)
greenlandGroup_BIN <- group_by(dfGSubset_BIN, bin_uri)

# Group by SBIN
nearcticGroup_SBIN <- group_by(dfNSubset_SBIN, cluster_45)
palearcticGroup_SBIN <- group_by(dfPSubset_SBIN, cluster_45)
greenlandGroup_SBIN <- group_by(dfGSubset_SBIN, cluster_45)

# Group by species
nearcticGroup_Sp <- group_by(dfNSubset_Sp, species_name)
palearcticGroup_Sp <- group_by(dfPSubset_Sp, species_name)
greenlandGroup_Sp <- group_by(dfGSubset_Sp, species_name)

# BIN counts per region
countsN_BIN <- summarize(nearcticGroup_BIN, count = n())
countsP_BIN <- summarize(palearcticGroup_BIN, count = n())
countsG_BIN <- summarize(greenlandGroup_BIN, count = n())

# SBIN counts per region
countsN_SBIN <- summarize(nearcticGroup_SBIN, count = n())
countsP_SBIN <- summarize(palearcticGroup_SBIN, count = n())
countsG_SBIN <- summarize(greenlandGroup_SBIN, count = n())

# Species counts per region
countsN_Sp <- summarize(nearcticGroup_Sp, count = n())
countsP_Sp <- summarize(palearcticGroup_Sp, count = n())
countsG_Sp <- summarize(greenlandGroup_Sp, count = n())

# Assign regions again 
# BIN
for (i in 1:nrow(countsN_BIN)){
  countsN_BIN$region[i] <- "Nearctic"
}
for (i in 1:nrow(countsP_BIN)){
  countsP_BIN$region[i] <- "Palearctic"
}
for (i in 1:nrow(countsG_BIN)){
  countsG_BIN$region[i] <- "Greenland"
}
# SBINs
for (i in 1:nrow(countsN_SBIN)){
  countsN_SBIN$region[i] <- "Nearctic"
}
for (i in 1:nrow(countsP_SBIN)){
  countsP_SBIN$region[i] <- "Palearctic"
}
for (i in 1:nrow(countsG_SBIN)){
  countsG_SBIN$region[i] <- "Greenland"
}
# Species
for (i in 1:nrow(countsN_Sp)){
  countsN_Sp$region[i] <- "Nearctic"
}
for (i in 1:nrow(countsP_Sp)){
  countsP_Sp$region[i] <- "Palearctic"
}
for (i in 1:nrow(countsG_Sp)){
  countsG_Sp$region[i] <- "Greenland"
}

# Combine together again - now its in the right format for spread function
countsAll_BIN <- rbind(countsN_BIN, countsP_BIN, countsG_BIN)
countsAll_SBIN <- rbind(countsN_SBIN, countsP_SBIN, countsG_SBIN)
countsAll_Sp <- rbind(countsN_Sp, countsP_Sp, countsG_Sp)

# First converting to the right format using tidyr
counts_spread_BIN <- spread(countsAll_BIN, key = bin_uri, value = count)
counts_spread_SBIN <- spread(countsAll_SBIN, key = cluster_45, value = count)
counts_spread_Sp <- spread(countsAll_Sp, key = species_name, value = count)

# If NA in a cell - assign a 0
counts_spread_BIN[is.na(counts_spread_BIN)] <- 0
counts_spread_SBIN[is.na(counts_spread_SBIN)] <- 0
counts_spread_Sp[is.na(counts_spread_Sp)] <- 0

# Make the region column the rowname
counts_spread1_BIN <- counts_spread_BIN[,-1]
row.names(counts_spread1_BIN) <- counts_spread_BIN$region

counts_spread1_SBIN <- counts_spread_SBIN[,-1]
row.names(counts_spread1_SBIN) <- counts_spread_SBIN$region

counts_spread1_Sp <- counts_spread_Sp[,-1]
row.names(counts_spread1_Sp) <- counts_spread_Sp$region

# Dissimilarity measures using chao for BIN, SBIN and species
chaoBIN <- vegdist(counts_spread1_BIN, method="chao")
chaoBIN

chaoSBIN <- vegdist(counts_spread1_SBIN, method="chao")
chaoSBIN

chaoSp <- vegdist(counts_spread1_Sp, method="chao")
chaoSp

################
# Mapping with plotly (BINs)

# New dataframe for filtered by subarctic
dfChironomidMap <- ChironomidAll
# Not filtred by subarctic
dfChironomidNonArctic <- dfChironomidFilter

containGreenland1 <- which(dfChironomidNonArctic$globalRegion=="Greenland")
dfGreenlandBIN1 <- dfChironomidNonArctic[containGreenland1,]

containNearctic1 <- which(dfChironomidNonArctic$globalRegion=="Nearctic")
dfNearcticBIN1 <- dfChironomidNonArctic[containNearctic1,]

containPalearctic1 <-  which(dfChironomidNonArctic$globalRegion=="Palearctic")
dfPalearcticBIN1 <- dfChironomidNonArctic[containPalearctic1,]

containGreenland2 <- which(dfChironomidMap$globalRegion=="Greenland")
dfGreenlandBIN2 <- dfChironomidMap[containGreenland2,]

containNearctic2 <- which(dfChironomidMap$globalRegion=="Nearctic")
dfNearcticBIN2 <- dfChironomidMap[containNearctic2,]

containPalearctic2 <-  which(dfChironomidMap$globalRegion=="Palearctic")
dfPalearcticBIN2 <- dfChironomidMap[containPalearctic2,]

# Using only region and lat/lon values
dfGMap1 <- (dfGreenlandBIN1[,c("globalRegion","latNum","lonNum")])
dfNMap1 <- (dfNearcticBIN1[,c("globalRegion","latNum","lonNum")])
dfPMap1 <- (dfPalearcticBIN1[,c("globalRegion","latNum","lonNum")])
dfGMap2 <- (dfGreenlandBIN2[,c("globalRegion","latNum","lonNum")])
dfNMap2 <- (dfNearcticBIN2[,c("globalRegion","latNum","lonNum")])
dfPMap2 <- (dfPalearcticBIN2[,c("globalRegion","latNum","lonNum")])
dfGMap3 <- rbind(dfGMap1, dfGMap2)

# delete duplicates - otherwise will overload plotly
dfGMap3 <- dfGMap3[!(duplicated(dfGMap3[c("latNum","lonNum")])), ]
dfNMap1 <- dfNMap1[!(duplicated(dfNMap1[c("latNum","lonNum")])), ]
dfNMap2 <- dfNMap2[!(duplicated(dfNMap2[c("latNum","lonNum")])), ]
dfPMap1 <- dfPMap1[!(duplicated(dfPMap1[c("latNum","lonNum")])), ]
dfPMap2 <- dfPMap2[!(duplicated(dfPMap2[c("latNum","lonNum")])), ]

# Export csv's for import into plotly for further formatting of the map
# on the plotly server:
write.csv(dfGMap3, file = "GMap3.csv")
write.csv(dfNMap1, file = "NMap1.csv")
write.csv(dfPMap1, file = "PMap1.csv")
write.csv(dfNMap2, file = "NMap2.csv")
write.csv(dfPMap2, file = "PMap2.csv")

##############
# Venn Diagram Calculations

# Venn Diagram of BINs

# Counts for each overlap region
GPN_BIN <- length(intersect(intersect(dfGreenlandBIN$bin_uri, dfPalearcticBIN$bin_uri), dfNearcticBIN$bin_uri)) 
GN_BIN <- length(intersect(dfGreenlandBIN$bin_uri, dfNearcticBIN$bin_uri)) - GPN_BIN
GP_BIN <- length(intersect(dfGreenlandBIN$bin_uri, dfPalearcticBIN$bin_uri)) - GPN_BIN
NP_BIN <- length(intersect(dfNearcticBIN$bin_uri, dfPalearcticBIN$bin_uri)) - GPN_BIN

# Counts for each circle
G_BIN <- length(unique(dfGreenlandBIN$bin_uri)) - (GN_BIN + GP_BIN + GPN_BIN)
N_BIN <- length(unique(dfNearcticBIN$bin_uri)) - (GN_BIN + NP_BIN + GPN_BIN)
P_BIN <- length(unique(dfPalearcticBIN$bin_uri)) - (GP_BIN + NP_BIN + GPN_BIN)

# Venn Diagram of Species

# Counts for each overlap region for species
GPN_Sp <- length(intersect(intersect(dfGreenlandBIN$species_name, dfPalearcticBIN$species_name), dfNearcticBIN$species_name))
GN_Sp <- length(intersect(dfGreenlandBIN$species_name, dfNearcticBIN$species_name)) - GPN_Sp
GP_Sp <- length(intersect(dfGreenlandBIN$species_name, dfPalearcticBIN$species_name)) - GPN_Sp
NP_Sp <- length(intersect(dfNearcticBIN$species_name, dfPalearcticBIN$species_name)) - GPN_Sp

# Counts for each circle for species
G_Sp <- length((unique(dfGreenlandBIN$species_name))) - (GN_Sp + GP_Sp + GPN_Sp)
N_Sp <- length(unique(dfNearcticBIN$species_name)) - (GN_Sp + NP_Sp + GPN_Sp)
P_Sp <- length(unique(dfPalearcticBIN$species_name)) - (GP_Sp + NP_Sp + GPN_Sp)

# Venn Diagram of SBINs at 4.5%

# Counts for each overlap region
GPN_SBIN <- length(intersect(intersect(dfGreenlandSBIN$cluster_45, dfPalearcticSBIN$cluster_45), dfNearcticSBIN$cluster_45))
GN_SBIN <- length(intersect(dfGreenlandSBIN$cluster_45, dfNearcticSBIN$cluster_45)) - GPN_SBIN
GP_SBIN <- length(intersect(dfGreenlandSBIN$cluster_45, dfPalearcticSBIN$cluster_45)) - GPN_SBIN
NP_SBIN <- length(intersect(dfNearcticSBIN$cluster_45, dfPalearcticSBIN$cluster_45)) - GPN_SBIN

# Counts for each circle
G_SBIN <- length((unique(dfGreenlandSBIN$cluster_45))) - (GN_SBIN + GP_SBIN + GPN_SBIN)
N_SBIN <- length(unique(dfNearcticSBIN$cluster_45)) - (GN_SBIN + NP_SBIN + GPN_SBIN)
P_SBIN <- length(unique(dfPalearcticSBIN$cluster_45)) - (GP_SBIN + NP_SBIN + GPN_SBIN)

# Using these counts in this shiny app that makes Venn diagrams:
# http://jolars.co/eulerr/
