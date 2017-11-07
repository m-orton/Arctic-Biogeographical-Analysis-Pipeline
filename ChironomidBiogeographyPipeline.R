###############
# Chironomidae Biogeography Pipeline

# Authored by Matthew G. Orton and Sally J. Adamowicz

# Credit to Torbjorn Ekrem and Elisabeth Ster for private Chironomidae data and some of the public
# Chironomidae data from BOLD and for helping us on the function of the code and plots.

# ***Note: to get all components of code to work, ensure that all Arctic zone files and the Private data csv are within the
# same current working R directory.***

##############
# Packages
# install.packages("foreach")
library(foreach)
# For genetic distance determination using the TN93 model, the ape package is
# required.
# install.packages("ape")
library(ape)
# The readr package will speed up parsing of the tsv file using the 
# read_tsv function.
# install.packages("readr")
library(readr)
# For sequence alignments we need the biostrings (DNAStringSet function) 
# and muscle libraries, as follows:
# source("https://bioconductor.org/biocLite.R")
# biocLite("Biostrings")
# biocLite("muscle")
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

##############
# Parsing from BOLD

# Public records for each of the three regions, for now im
# keeping all columns so I dont get rid of any important info by mistake
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

# Filter according to COI-5P - are we just restricting to COI-5P? If not, we will have to perform
# multiple distance matrix calaculations in the outlier section
containCOI <- grep( "^CO", dfChironomid$markercode)
dfChironomid <- dfChironomid[containCOI,]

#Can use this command to check to make sure all markers are COI-5P
unique(dfChironomid$markercode)

# Conversion to numeric:
latNum <- with(dfChironomid, as.numeric(as.character(lat))) 
dfChironomid$latNum <- latNum
lonNum <- with(dfChironomid, as.numeric(as.character(lon))) 
dfChironomid$lonNum <- lonNum

##############
# Select One Sequence per BIN - selecting sequence closest to 658 bp (can be changed)

# Make a list of BINs
binList <- lapply(unique(dfChironomid$bin_uri), function(x) 
  dfChironomid[dfChironomid$bin_uri == x,])
  
# Also need to find the number of unique BINs in dfBinList. 
binNumber<- unique(dfChironomid$bin_uri)
binNumber <- length(binNumber)
  
# Extract record id from each BIN.
binRecordId <- foreach(i=1:binNumber) %do% (binList[[i]]$recordID)

# Select one sequence per BIN with 658 bp
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
binSelect <- as.numeric(binSelect)

# Subset dfChironomid
dfSingleSeq <- subset(dfChironomid, recordID %in% binSelect)

# Remove the Bin list since it is a large variable and no longer needed
rm(binList)
gc()

##############
# Identification of outlier points based on genetic
# distance, based on spread, will eliminate records which are clearly not
# belonging to Chironomidae if by mistake a record was misindentified

# To do this doing an alignment and a distance matrix 
# (just using the TN93 model)

dnaStringSet <- DNAStringSet(dfSingleSeq$nucleotides)
# For testing just set to maxiters=2 and diags=TRUE
alignment <- muscle(dnaStringSet, maxiters = 2, diags = TRUE)
# Back to a string set
dnaStringSet2 <- DNAStringSet(alignment)

# Name the stringset with record ids
bin_uri <- dfSingleSeq$bin_uri
names(dnaStringSet2) <- bin_uri

# DNAbin format
dnaBin <- as.DNAbin(dnaStringSet2)
# Distance matrix
distanceMatrix <- dist.dna(dnaBin, model = "TN93", as.matrix = TRUE, 
  pairwise.deletion = TRUE)

# Compute mean per row (each row corresponding to a spearate BIN), this will give 
# the average distance per BIN
rowmeans <- rowMeans(distanceMatrix, na.rm = FALSE, dims = 1)

# Make a new dataframe for the rowmeans to check for outliers
dfOutlierCheck <- data.frame(rowmeans)
dfOutlierCheck$bin_uri <- names(rowmeans)

# Title for the Boxplot
boxTitle <- "Boxplot of Average Pairwise Distance per BIN"

# Plotting a boxplot of the distances here to find the outliers
# Doing in plotly to identify BINs easily by hovering over plot
plot_ly(y = dfOutlierCheck$rowmeans, type = "box") %>%
  layout(title = paste0(boxTitle))

# Summary of quantiles
summary(dfOutlierCheck$rowmeans)

# Still need to figure out how to identify the outlier point by BIN to eliminate them
# since there are a few you can see past 0.4

# Using the IQR to detect outliers
lowerQuantile <- quantile(dfOutlierCheck$rowmeans)[2]
lowerQuantile
upperQuantile <- quantile(dfOutlierCheck$rowmeans)[4]
upperQuantile
iqr <- upperQuantile - lowerQuantile
upperThreshold <- (iqr * 1.5) + upperQuantile
upperThreshold

# Identify BINs with no relatives within "typical" range of genetic
# divergence (i.e. all of their genetic distances fall outside 1.5 x IQR
# upper threshold.)

# Checking number of values above upperThreshold
length(dfOutlierCheck$rowmeans[dfOutlierCheck$rowmeans>upperThreshold])

# Subsetting by values to get BINs with values above this threshold
outlier1.5 <- which(dfOutlierCheck$rowmeans>upperThreshold)
# Making a separate dataframe with these values and BINs 
dfOutlier1.5 <- dfOutlierCheck[outlier1.5,]

# Checking outliers above the iqr*3 rule
upperThreshold2 <- (iqr * 3) + upperQuantile
upperThreshold2
length(dfOutlierCheck$rowmeans[dfOutlierCheck$rowmeans>upperThreshold2])

# Subsetting by values to get BINs with values above this threshold
outlier3 <- which(dfOutlierCheck$rowmeans>upperThreshold2)
# Making a separate dataframe with these values and BINs 
dfOutlier3 <- dfOutlierCheck[outlier3,]

# ***Upon using avergae pairwise distance per BIN I am now getting very few numbers of 
# outliers, 19 at 1.5 and 4 at 3, each should have its own dataframe with
# BINs listed that we can show to Elisabeth to see what she thinks

##############
# Reading in and filtering of private data for use in mapping and analyses, can assume we can 
# skip the outlier check with the private data?

# Filtering is redone for private records due to differences in column names between public and private
# Note: am currently getting 12 parsing errors in the records
dfPrivateData <- read_csv("Private_Chironomid_Data_ModifiedSingleSheet.csv")

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

# Only Canada, Greenland, Iceland and Norway were found

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
colnames(dfPrivateData)[5] <- "subfamily_name"
colnames(dfPrivateData)[8] <- "species_name"
colnames(dfPrivateData)[25] <- "collectors"
colnames(dfPrivateData)[27] <- "country"
colnames(dfPrivateData)[28] <- "province_state"
colnames(dfPrivateData)[29] <- "region"
colnames(dfPrivateData)[30] <- "sector"
colnames(dfPrivateData)[31] <- "exactsite"
colnames(dfPrivateData)[48] <- "bin_uri"

dfPrivateData <- (dfPrivateData[,c("globalRegion","bin_uri","species_name","subfamily_name","latNum","lonNum",
                                   "country","province_state","region","sector","exactsite","collectors")])

dfChironomid <- (dfChironomid[,c("globalRegion","bin_uri","species_name","subfamily_name","latNum","lonNum",
                                 "country","province_state","region","sector","exactsite","collectors")])

# Combine both dataframes together
dfChironomidAll <- rbind(dfPrivateData, dfChironomid)

# *** Upon checking with Elisabeth, the only BIN that we should eliminate is ACZ1013:
binCheck <- which(dfChironomidAll$bin_uri == "BOLD:ACZ1013")
dfChironomidAll <- dfChironomidAll[-binCheck,]

##############
# Subarctic filtering according to subarctic shapefile

# Subarctic filtering of all data including private data
# Read in the subarctic shapefile
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

length(unique(dfChironomidAll$species_name))
length(unique(dfChironomidAll$bin_uri))

#############
# Accumulation Curve Analysis

# Need to filter according to shapefile

# Checking the number of records in dfGreenland
nGreen <-length(dfGreenland$bin_uri)
nGreen

# Removing records from dfGreenland that don't contain a BIN.
# Creating a new df, to be used for the accumulation curve analysis of Greenland.
containGreenland <- which(dfChironomidAll$globalRegion=="Greenland")
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

# For testing large datasets, I suggest to change the number after replicate to a very small number,
# such as 10. For the final version, we will want a large number to create smooth lines, such as 
# 1,000. For Greenland, I ran 100 replicates in a couple of minutes on my laptop. ***Update, setting to 10 
# for testing****
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
containNearctic <- which(dfChironomidAll$globalRegion=="Nearctic")
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
nearcticBIN <- replicate(100, {
  length(unique(sample(dfAccNearctic1, size=10)))
})
nearcticBIN

# Testing replicates for Nearctic, ***Update, setting to 5 for now
# for testing because Nearctic is so large****
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
containPalearctic <-  which(dfChironomidAll$globalRegion=="Palearctic")
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

# Testing replicates for Nearctic, ***Update, setting to 10 
# for testing****
for (i in 1:(nPalearctic-1))
{palearcticBIN <- replicate(100, {
  length(unique(sample(dfAccPalearctic1, size=i)))
})
{meanBINPalearctic[i] <- mean(palearcticBIN, na.rm=TRUE)
}
}

# Plotly visualization of meanBIN's for each region on the same plot

# I chose plotly because you can hover over points on the graph to see the precise point
# value

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

api_create(pRegion, filename = "AccCurve100Rep")
Sys.setenv("plotly_username"="Matt14") 
Sys.setenv("plotly_api_key"="W0jEdwXDUtAAVsmY9nRJ")


# Truncated Plot of Accumulation according to region, truncated by smallest region Palearctic
x2 <- list(
  title = "No. of Specimens Barcoded",
  titlefont = f,
  tickfont = f2,
  range = c(0, 1600)
)

pRegion2 <- plot_ly(data = dfRegion, y = dfRegion$combinedRegions, color = dfRegion$region, type = "scatter", mode = "markers") %>%
  layout(xaxis = x2, yaxis = y, legend = l)

pRegion2

#############
# Dplyr/Tidyr and Vegan

# First separate bin_uri and global region from other columns
dfNSubset <- (dfAccNearctic[,c("globalRegion","bin_uri")])
dfPSubset <- (dfAccPalearctic[,c("globalRegion","bin_uri")])
dfGSubset <- (dfAccGreenland[,c("globalRegion","bin_uri")])

# Dividing Greenland into 2 regions for one set of dissimilarity measures - East and West - Dividing by -30 lon
# This will divide between Zackenberg Research Station on the east and all points on west
greenlandEast <- which(dfAccGreenland$lonNum>-30)
dfGreenlandEast  <- dfAccGreenland[greenlandEast,]
dfGEastSubset <- (dfGreenlandEast[,c("globalRegion","bin_uri")])
dfGreenlandWest <- dfAccGreenland[-greenlandEast,]
dfGWestSubset <- (dfGreenlandWest[,c("globalRegion","bin_uri")])

# Group by BIN
nearcticGroup <- group_by(dfNSubset, bin_uri)
palearcticGroup <- group_by(dfPSubset, bin_uri)
greenlandGroup <- group_by(dfGSubset, bin_uri)
greenEast <- group_by(dfGEastSubset, bin_uri)
greenWest <- group_by(dfGWestSubset, bin_uri)

# BIN counts per region
countsN <- summarize(nearcticGroup, count = n())
countsP <- summarize(palearcticGroup, count = n())
countsG <- summarize(greenlandGroup, count = n())
countsGE <- summarize(greenEast, count = n())
countsGW <- summarize(greenWest, count = n())

# Assign regions again - ***gives a warning but does work***
for (i in 1:nrow(countsN)){
 countsN$region[i] <- "Nearctic"
}
for (i in 1:nrow(countsP)){
  countsP$region[i] <- "Palearctic"
}
for (i in 1:nrow(countsG)){
  countsG$region[i] <- "Greenland"
}
for (i in 1:nrow(countsGE)){
  countsGE$region[i] <- "GreenlandEast"
}
for (i in 1:nrow(countsGW)){
  countsGW$region[i] <- "GreenlandWest"
}

# Combine together again - now its in the right format for spread function
countsAll <- rbind(countsN, countsP, countsG)
countsAllGDivide <- rbind(countsN, countsP, countsGE, countsGW)

# First converting to the right format using tidyr
counts_spread <- spread(countsAll, key = bin_uri, value = count)
countsAllGDivide_spread <- spread(countsAllGDivide, key = bin_uri, value = count)

# If NA in a cell - assign a 0
counts_spread[is.na(counts_spread)] <- 0
countsAllGDivide_spread[is.na(countsAllGDivide_spread)] <- 0

# Make the region column the rowname
counts_spread1 <- counts_spread[,-1]
row.names(counts_spread1) <- counts_spread$region

countsGDivide_spread1 <- countsAllGDivide_spread[,-1]
row.names(countsGDivide_spread1) <- countsAllGDivide_spread$region


# Dissimilarity measure using chao
chao <- vegdist(counts_spread1, method="chao")
chao

# Dissimilarity measure using cao
cao <- vegdist(counts_spread1, method="cao")
cao

# Dissimilarity measure using chao when dividing Greenland
chao <- vegdist(countsGDivide_spread1, method="chao")
chao

# Dissimilarity measure using cao when dividing Greenland
cao <- vegdist(countsGDivide_spread1, method="cao")
cao


###############
# Plot on plotly after Vegan analysis 

# Also allows visualization of where points may be located

# # *** After generating the map, you will need to click on the icon in the 
# viewer (bottom right corner) that says "show in new window" 
# (little box with arrow beside the refresh icon). 
# Unfortunately, this does not show the actual map directly in Rstudio.
# The map will appear in a web browser window, though you don't have to be 
# online to do this.***

mapLayout <- list(
  showland = TRUE,
  showlakes = FALSE,
  showcountries = FALSE,
  showocean = FALSE,
  countrywidth = 0.5,
  landcolor = toRGB("grey90"),
  lakecolor = toRGB("white"),
  oceancolor = toRGB("white"),
  resolution = list(type = '1000'),
  projection = list(type = 'transverse mercator'),
  lonaxis = list(
    #range = c(30, -160),
    showgrid = TRUE,
    gridcolor = toRGB("gray40"),
    gridwidth = 0.5
  ),
  lataxis = list(
    #range = c(50, 90),
    showgrid = TRUE,
    gridcolor = toRGB("gray40"),
    gridwidth = 0.5
  )
)

# New dataframe column with data for hovering over points on the map.
# Can add more columns to hover if you want more detail on the map for each
# pairing.
dfChironomidAll$hover <- 
  paste("BIN:",dfChironomidAll$bin_uri,
        "Subfamily:",dfChironomidAll$subfamily_name,
        "Country:",dfChironomidAll$country,
        "Province/State:",dfChironomidAll$province_state,
        "Region:",dfChironomidAll$region,
        "Sector:",dfChironomidAll$sector,
        "Site:",dfChironomidAll$exactsite,
        "Collectors:",dfChironomidAll$collectors,
        sep = "<br>")

# Got rid of the title since it will be included in the figure legend

# This command will ensure the pairing results dataframe can be read by plotly.
attach(dfChironomidMap)
# This command will show a map organized by region (each color to a different region)
# of all Arctic Chironomid records
plot2 <- plot_ly(dfChironomidMap, lat = dfChironomidMap$latNum, lon = dfChironomidMap$lonNum, 
        text = hover, color = globalRegion,
        mode = "markers", type = 'scattergeo') %>%
  layout(geo = mapLayout)
  
# Export csv's for import into plotly - this will be each region points
dfGMap <- (dfAccGreenland[,c("globalRegion","latNum","lonNum")])
dfNMap <- (dfAccNearctic[,c("globalRegion","latNum","lonNum")])
dfPMap <- (dfAccPalearctic[,c("globalRegion","latNum","lonNum")])

write.csv(dfGMap, file = "GMap.csv")
write.csv(dfNMap, file = "NMap.csv")
write.csv(dfPMap, file = "PMap.csv")

# Then importing in the points for trace lines:
# Find the mode for Greenland:
modePoint1 <- table(dfAccGreenland$latNum)
# 74.467,-20.567  - most number of occurences at this point - 4344 occurences - Zackenberg Research station

# Taking Churchill for Nearctic - 58.7690, -94.1600

modePoint2 <- table(dfAccPalearctic$latNum)
# 78.071, 13.793 - 74 occurences - Nordenskioldland

##############
# Venn Diagram of BINs

# Counts for each overlap region
ABC <- length(intersect(intersect(dfAccGreenland$bin_uri, dfAccPalearctic$bin_uri), dfAccNearctic$bin_uri)) 
AB <- length(intersect(dfAccGreenland$bin_uri, dfAccNearctic$bin_uri)) - ABC
AC <- length(intersect(dfAccGreenland$bin_uri, dfAccPalearctic$bin_uri)) - ABC
BC <- length(intersect(dfAccNearctic$bin_uri, dfAccPalearctic$bin_uri)) - ABC

# Counts for each circle
A <- length(unique(dfAccGreenland$bin_uri)) - (ABC + AB + AC)
B <- length(unique(dfAccNearctic$bin_uri)) - (AB + BC + ABC)
C <- length(unique(dfAccPalearctic$bin_uri)) - (AC + BC + ABC)

# Should now equal 1520 unique BINs

# Species

# Counts for each overlap region for species
ABC <- length(intersect(intersect(dfAccGreenland$species_name, dfAccPalearctic$species_name), dfAccNearctic$species_name))
AB <- length(intersect(dfAccGreenland$species_name, dfAccNearctic$species_name)) - ABC
AC <- length(intersect(dfAccGreenland$species_name, dfAccPalearctic$species_name)) - ABC
BC <- length(intersect(dfAccNearctic$species_name, dfAccPalearctic$species_name)) - ABC

# Counts for each circle for species
A <- length(unique(dfAccGreenland$species_name)) - (ABC + AB + AC)
B <- length(unique(dfAccNearctic$species_name)) - (AB + BC + ABC)
C <- length(unique(dfAccPalearctic$species_name)) - (AC + BC + ABC)

# Should now equal 572 unique species

# Using these counts in this shiny app that makes Venn diagrams:
# http://jolars.co/eulerr/
