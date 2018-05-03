###############
# Chironomid Biogeography Pipeline 
# (Not using Subarctic Shapefile Filter - All NOR and CAN records included)

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


##############
# Parsing from BOLD

# Public records for each of the three regions
dfNearctic <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Chironomidae&geo=Alaska|Canada&format=tsv")
dfGreenland <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Chironomidae&geo=Greenland&format=tsv")
dfPalearctic <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Chironomidae&geo=Norway|Denmark|Iceland|Sweden|Finland&format=tsv")

# Note that BOLD makes the distinction between Greenland and Denmark and distinguishes them as separate countries
# (even though they are not) so the datasets between Denmark and Greenland are nonoverlapping

##############
# Record filtering - lat/lon coordinates, BIN, presence of a sequence, COI-5P

# Labeling regions
dfNearctic$globalRegion <- "Nearctic"
dfGreenland$globalRegion <- "Greenland"
dfPalearctic$globalRegion <- "Palearctic"

#Combine dataframes with regional identifiers
dfChironomid <- rbind(dfNearctic, dfGreenland, dfPalearctic)

# Filter by BIN assignment (grep by colon since all BIN identifiers have this):
containBin <- grep( "[:]", dfChironomid$bin_uri)
dfChironomid <- dfChironomid[containBin,]

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

# Check if there is any intersection between dfPrivate data and dfChrionomid for sample ID
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
dfSingleSeq = dfSingleSeq[order(dfSingleSeq[,'processid'],-dfSingleSeq[,'seqLength']),]
dfSingleSeq = dfSingleSeq[!duplicated(dfSingleSeq$processid),]

##############
# Single-linkage Clustering using the Decipher Package

# Alignment step
dnaStringSet <- DNAStringSet(dfSingleSeq$nucleotides)
alignment <- muscle(dnaStringSet, maxiters = 2, diags = TRUE, gapopen = -3000)
dnaStringSet2 <- DNAStringSet(alignment)

# Name the stringset with record ids
bin_uri <- dfSingleSeq$bin_uri
names(dnaStringSet2) <- bin_uri

# Write out to fasta
fileName <- paste("ChironomidAlignment_AllCanNor.fas")
writeXStringSet(dnaStringSet2, file=fileName, format = "fasta", width = 658)

# DNAbin format
dnaBin <- as.DNAbin(dnaStringSet2)
# Distance matrix using TN93 before clustering 
distanceMatrix <- dist.dna(dnaBin, model = "TN93", as.matrix = TRUE, 
                           pairwise.deletion = TRUE)

# Clustering according to 4% divergence threshold
clustSingle4 <- IdClusters(distanceMatrix,
                           method = "single",
                           cutoff= 0.04,
                           showPlot = TRUE,
                           type = "clusters",
                           processors = 2,
                           verbose = TRUE)

# Number of unique clusters for 4%
length(unique(clustSingle4$cluster))

# 4.5% divergence threshold
clustSingle45 <- IdClusters(distanceMatrix,
                            method = "single",
                            cutoff= 0.045,
                            showPlot = TRUE,
                            type = "clusters",
                            processors = 2,
                            verbose = TRUE)

# Number of unique clusters for 4.5%
length(unique(clustSingle45$cluster))

# 5% divergence threshold
clustSingle5 <- IdClusters(distanceMatrix,
                           method = "single",
                           cutoff= 0.05,
                           showPlot = TRUE,
                           type = "clusters",
                           processors = 2,
                           verbose = TRUE)

# Number of unique clusters for 5%
length(unique(clustSingle5$cluster))

# Renaming column for 4% cluster
clustSingle4 <- setDT(clustSingle4, keep.rownames = TRUE)[]
colnames(clustSingle4)[2] <- "cluster_4"

# Merge clusters to original Chironomid dataset (pre-alignment filter)
# Every record of every BIN will now have a defined "superBIN"!
dfChironomidSBIN <- merge(dfChironomidAll, clustSingle4, by.x ="bin_uri", by.y ="rn")

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

# Extract out bad species names for Greenland before using for species analysis
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

# Or separate by SBIN - 4%
dfNSubset_SBIN <- (dfNearcticSBIN[,c("globalRegion","cluster_4")])
dfPSubset_SBIN <- (dfPalearcticSBIN[,c("globalRegion","cluster_4")])
dfGSubset_SBIN <- (dfGreenlandSBIN[,c("globalRegion","cluster_4")])

# Separate by species 
dfNSubset_Sp <- (dfNearcticSBIN[,c("globalRegion","species_name")])
dfPSubset_Sp <- (dfPalearcticSBIN[,c("globalRegion","species_name")])
dfGSubset_Sp <- (dfGreenlandSBIN[,c("globalRegion","species_name")])

# Group by BIN
nearcticGroup_BIN <- group_by(dfNSubset_BIN, bin_uri)
palearcticGroup_BIN <- group_by(dfPSubset_BIN, bin_uri)
greenlandGroup_BIN <- group_by(dfGSubset_BIN, bin_uri)

# Group by SBIN
nearcticGroup_SBIN <- group_by(dfNSubset_SBIN, cluster_4)
palearcticGroup_SBIN <- group_by(dfPSubset_SBIN, cluster_4)
greenlandGroup_SBIN <- group_by(dfGSubset_SBIN, cluster_4)

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
counts_spread_SBIN <- spread(countsAll_SBIN, key = cluster_4, value = count)
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

#############
# Accumulation Curve Analysis (for BINs - all CAN and NOR)
# 100 replicates per region for smooth curving

# Checking the number of records in dfGreenland
nGreen <-length(dfGreenland$bin_uri)
nGreen

# Removing records from dfGreenland that don't contain a BIN.
# Creating a new df, to be used for the accumulation curve analysis of Greenland.
dfAccGreenland <- dfChironomidAll[containGreenland,]
dfAccGreenland$bin_uri
nGreen <-length(dfAccGreenland$bin_uri)
nGreen
length(unique(dfAccGreenland$bin_uri))

# Creating an array to hold the BIN counts based upon random draws from the records
meanBINGreenland <- array (NA,dim=c(length(dfAccGreenland$bin_uri))-1)

# Reducing dfAccGreenland to just the BIN column
dfAccGreenland1 <- dfAccGreenland$bin_uri
dfAccGreenland1

# A resampling analysis was performed to assess how BINs accumulate in Greenland
# as individuals are sampled. A steep accumulation curve would indicate BINs remain
# to be collected. A curve that levels off would indicate the sampling is approaching
# completeness, for a given sampling method.

# Testing and building up components for building the rarefaction curve below.
sample(dfAccGreenland1, size=10)
unique(sample(dfAccGreenland1, size=10))
length(unique(sample(dfAccGreenland1, size=10)))

# Testing that this bit is doing what we want:
greenlandBIN <-replicate(10, {
  length(unique(sample(dfAccGreenland1, size=10)))
})
greenlandBIN

# 100 replicates
for (i in 1:(nGreen-1))
{greenlandBIN <- replicate(100, {
  length(unique(sample(dfAccGreenland1, size=i)))
})
{meanBINGreenland[i] <- mean(greenlandBIN, na.rm=TRUE)
}
}

# Repeating analysis for Nearctic 

# checking the number of records in dfNearctic
nNearctic <-length(dfNearctic$bin_uri)
nNearctic

# Subsetting to get records in Nearctic
dfAccNearctic <- dfChironomidAll[containNearctic,]
dfAccNearctic$bin_uri
nNearctic <-length(dfAccNearctic$bin_uri)
nNearctic
length(unique(dfAccNearctic$bin_uri))

# Creating an array to hold the BIN counts based upon random draws from the records
meanBINNearctic <- array (NA,dim=c(length(dfAccNearctic$bin_uri))-1)

# Reducing dfAccNearctic to just the BIN column
dfAccNearctic1 <- dfAccNearctic$bin_uri
dfAccNearctic1

# Testing and building up components for building the rarefaction curve below.
sample(dfAccNearctic1, size=10)
unique(sample(dfAccNearctic1, size=10))
length(unique(sample(dfAccNearctic1, size=10)))

# Testing that this bit is doing what we want:
nearcticBIN <- replicate(10, {
  length(unique(sample(dfAccNearctic1, size=10)))
})
nearcticBIN

# Testing replicates for Nearctic - 100 replicates
for (i in 1:(nNearctic-1))
{nearcticBIN <- replicate(100, {
  length(unique(sample(dfAccNearctic1, size=i)))
})
{meanBINNearctic[i] <- mean(nearcticBIN, na.rm=TRUE)
}
}

# Palearctic

# checking the number of records in dfPalearctic
nPalearctic <-length(dfPalearctic$bin_uri)
nPalearctic

# Subsetting to get records in Palearctic
dfAccPalearctic <- dfChironomidAll[containPalearctic,]
dfAccPalearctic$bin_uri
nPalearctic <-length(dfAccPalearctic$bin_uri)
nPalearctic
length(unique(dfAccPalearctic$bin_uri))

# Creating an array to hold the BIN counts based upon random draws from the records
meanBINPalearctic <- array (NA,dim=c(length(dfAccPalearctic$bin_uri))-1)

# Reducing dfAccNearctic to just the BIN column
dfAccPalearctic1 <- dfAccPalearctic$bin_uri
dfAccPalearctic1

# Testing and building up components for building the rarefaction curve below.
sample(dfAccPalearctic1, size=10)
unique(sample(dfAccPalearctic1, size=10))
length(unique(sample(dfAccPalearctic1, size=10)))

# Testing that this bit is doing what we want:
palearcticBIN<-replicate(10, {
  length(unique(sample(dfAccPalearctic1, size=10)))
})
palearcticBIN

# Testing replicates for Palearctic - 100 replicates
for (i in 1:(nPalearctic-1))
{palearcticBIN <- replicate(100, {
  length(unique(sample(dfAccPalearctic1, size=i)))
})
{meanBINPalearctic[i] <- mean(palearcticBIN, na.rm=TRUE)
}
}

# Plotly visualization of meanBIN's for each region on the same plot

# First making a dataframe of meanBIN's per region

# Naming each list of means
for (i in 1:(length(meanBINGreenland))){
  names(meanBINGreenland)[i] <- dfAccGreenland$globalRegion[i]
}

for (i in 1:(length(meanBINPalearctic))){
  names(meanBINPalearctic)[i] <- dfAccPalearctic$globalRegion[i]
}

for (i in 1:(length(meanBINNearctic))){
  names(meanBINNearctic)[i] <- dfAccNearctic$globalRegion[i]
}

# Combining means together
combinedRegions <- append(meanBINPalearctic, meanBINGreenland)
combinedRegions <- append(combinedRegions, meanBINNearctic)                          

# Make a dataframe with means before plotting
dfRegion <- data.frame(combinedRegions)
dfRegion$region <- names(combinedRegions)

# Removed the title since it will be included in the figure legend

# X and Y axis 
f <- list(
  family = "open-sans",
  size = 18,
  color = "#7f7f7f"
)

f2 <- list(
  family = "open-sans",
  size = 13,
  color = "#7f7f7f"
)

x <- list(
  title = "Number of Specimens Barcoded",
  titlefont = f,
  tickfont = f2
)

y <- list(
  title = "Number of Barcode Index Numbers (BIN)",
  titlefont = f,
  tickfont = f2
)

# Legend
l <- list(
  font = list(
    family = "open-sans",
    size = 16,
    color = "#000"),
  bgcolor = "#E2E2E2",
  bordercolor = "#FFFFFF",
  borderwidth = 2)

# Making a plot for the region and storing in a variable
pRegion <- plot_ly(data = dfRegion, y = dfRegion$combinedRegions, color = dfRegion$region, type = "scatter", mode = "markers") %>%
  layout(xaxis = x, yaxis = y, legend = l)

pRegion

# Uploading to plotly for further formatting on their web server
# api_create(pRegion, filename = "AccCurve100Rep")
# Sys.setenv("plotly_username"="") 
# Sys.setenv("plotly_api_key"="")

#############
# Greenland East/West Division Dissimilarity Measure - BINs only (not used in first submission)

# Dividing Greenland into 2 regions for one set of dissimilarity measures - East and West - Dividing by -30 lon
# This will divide between Zackenberg Research Station on the east and all points on west
greenlandEast <- which(dfAccGreenland$lonNum>-30)
dfGreenlandEast  <- dfAccGreenland[greenlandEast,]
dfGEastSubset <- (dfGreenlandEast[,c("globalRegion","bin_uri")])

dfGreenlandWest <- dfAccGreenland[-greenlandEast,]
dfGWestSubset <- (dfGreenlandWest[,c("globalRegion","bin_uri")])

greenEast <- group_by(dfGEastSubset, bin_uri)
greenWest <- group_by(dfGWestSubset, bin_uri)

countsGE <- summarize(greenEast, count = n())
countsGW <- summarize(greenWest, count = n())

for (i in 1:nrow(countsGE)){
  countsGE$region[i] <- "GreenlandEast"
}
for (i in 1:nrow(countsGW)){
  countsGW$region[i] <- "GreenlandWest"
}

countsAllGDivide <- rbind(countsN_BIN, countsP_BIN, countsGE, countsGW)

countsAllGDivide_spread <- spread(countsAllGDivide, key = bin_uri, value = count)

countsAllGDivide_spread[is.na(countsAllGDivide_spread)] <- 0

countsGDivide_spread1 <- countsAllGDivide_spread[,-1]
row.names(countsGDivide_spread1) <- countsAllGDivide_spread$region

# Dissimilarity measure using chao when dividing Greenland
chao_GDivide <- vegdist(countsGDivide_spread1, method="chao")
chao_GDivide

##############
# Venn Diagram Calculation

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

# Venn Diagram of SBINs at 4%

# Counts for each overlap region
GPN_SBIN <- length(intersect(intersect(dfGreenlandSBIN$cluster_4, dfPalearcticSBIN$cluster_4), dfNearcticSBIN$cluster_4))
GN_SBIN <- length(intersect(dfGreenlandSBIN$cluster_4, dfNearcticSBIN$cluster_4)) - GPN_SBIN
GP_SBIN <- length(intersect(dfGreenlandSBIN$cluster_4, dfPalearcticSBIN$cluster_4)) - GPN_SBIN
NP_SBIN <- length(intersect(dfNearcticSBIN$cluster_4, dfPalearcticSBIN$cluster_4)) - GPN_SBIN

# Counts for each circle
G_SBIN <- length((unique(dfGreenlandSBIN$cluster_4))) - (GN_SBIN + GP_SBIN + GPN_SBIN)
N_SBIN <- length(unique(dfNearcticSBIN$cluster_4)) - (GN_SBIN + NP_SBIN + GPN_SBIN)
P_SBIN <- length(unique(dfPalearcticSBIN$cluster_4)) - (GP_SBIN + NP_SBIN + GPN_SBIN)

# Using these counts in this shiny app that makes Venn diagrams:
# http://jolars.co/eulerr/
