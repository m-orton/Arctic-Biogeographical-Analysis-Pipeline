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
# Subarctic Filtering - filter any points below subartic line using subartic shapefile

rawShapefile <- shapefile("Arctic_Zones")
subarcticZone <- rawShapefile[3]

# Can show a plot of the subarctic zone using the shapefile and geom_polygon
p <- ggplot()+geom_polygon(data=subarcticZone, aes(x=long, y=lat, group=group),
                      fill="cadetblue", color="grey")#+
  #geom_point(data = dfAccPalearctic, aes(x = dfAccPalearctic$lonNum, y = dfChironomid$latNum), color = "black", size = 1) 
subarcticPlot <- ggplotly(p)

# Note that the plot scale is not in lat/lon but in UTM (Universal Transmercator units)
subarcticPlot


# For now just filtering above 55 degrees lat to make the dataset more manageable
# ***Only doing this for testing purposes***
latCheck <- which(dfChironomid$latNum>=55)
dfChironomid <- dfChironomid[latCheck,]

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

# *** Upon checking with Elisabeth, the only BIN that we should eliminate is ACZ1013:
binCheck <- which(dfChironomid$bin_uri == "BOLD:ACZ1013")
dfChironomid <- dfChironomid[-binCheck,]

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

# Conversion to numeric coordinates:
latNum <- with(dfPrivateData, as.numeric(as.character(Lat))) 
dfPrivateData$latNum <- latNum
lonNum <- with(dfPrivateData, as.numeric(as.character(Lon))) 
dfPrivateData$lonNum <- lonNum

# Filtering again according to the shapefile


# Combining private data with dfChironomid and excluding unecessary columns for the analysis
# and mapping
colnames(dfPrivateData)[4] <- "subfamily_name"
colnames(dfPrivateData)[24] <- "collectors"
colnames(dfPrivateData)[26] <- "country"
colnames(dfPrivateData)[27] <- "province_state"
colnames(dfPrivateData)[28] <- "region"
colnames(dfPrivateData)[29] <- "sector"
colnames(dfPrivateData)[30] <- "exactsite"
colnames(dfPrivateData)[47] <- "bin_uri"

dfPrivateData <- (dfPrivateData[,c("globalRegion","bin_uri","subfamily_name","latNum","lonNum",
                                   "country","province_state","region","sector","exactsite","collectors")])

dfChironomid <- (dfChironomid[,c("globalRegion","bin_uri","subfamily_name","latNum","lonNum",
                                 "country","province_state","region","sector","exactsite","collectors")])

dfChironomidAll <- rbind(dfPrivateData, dfChironomid)

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
{greenlandBIN <- replicate(5, {
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
dfAccNearctic1 <- dfNearctic$bin_uri
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

# Testing replicates for Nearctic, ***Update, setting to 5 for now
# for testing because Nearctic is so large****
for (i in 1:(nNearctic-1))
{nearcticBIN <- replicate(5, {
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
{palearcticBIN <- replicate(5, {
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

# Plot title
pRegionTitle = "Accumuluation Curves for Chironomidae Subarctic Regions - 5 replicates"

# Making a plot for the region and storing in a variable
pRegion <- plot_ly(data = dfRegion, y = dfRegion$combinedRegions, color = region, type = "scatter", mode = "markers") %>%
  layout(title = paste0(pRegionTitle))
  
plot_ly(data = dfRegion, y = dfRegion$combinedRegions, color = region, type = "scatter", mode = "markers") %>%
  layout(title = paste0(pRegionTitle))

# Truncated Plot of Accumulation according to region
pRegionTitle2 = "Accumuluation Curves for Chironomidae Subarctic Regions - 5 replicates, Truncated"

pRegion2 <- plot_ly(data = dfRegion, y = dfRegion$combinedRegions, color = region, type = "scatter", mode = "markers") %>%
  layout(title = paste0(pRegionTitle2), xaxis = list(range = c(0, 20000)))

plot_ly(data = dfRegion, y = dfRegion$combinedRegions, color = region, type = "scatter", mode = "markers") %>%
  layout(title = paste0(pRegionTitle2), xaxis = list(range = c(0, 20000)))

# mean BIN per major subfamily on same plot

# First finding the subfamilies
unique(dfChironomidAll$subfamily_name)


#############
# Vegan

# possibly useful package to consider. However, the data format is different from ours.
# So, we would need to consider that before we pursue vegan for our further analyses.
# install.packages("vegan")
# library(vegan)



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
  showlakes = TRUE,
  showcountries = TRUE,
  showocean = TRUE,
  countrywidth = 0.5,
  landcolor = toRGB("grey90"),
  lakecolor = toRGB("white"),
  oceancolor = toRGB("white"),
  resolution = list(type = '50'),
  projection = list(type = 'transverse mercator'),
  lonaxis = list(
    range = list(range = c(-157, 23)),
    showgrid = TRUE,
    gridcolor = toRGB("gray40"),
    gridwidth = 0.5
  ),
  lataxis = list(
    range = list(range = c(55, 90)),
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

# Title for the map
mapTitle <- paste("Subarctic Zone Chironomidae Records per Region")

# This command will ensure the pairing results dataframe can be read by plotly.
attach(dfChironomidAll)
# This command will show a map organized by region (each color to a different region)
# of all Arctic Chironomid records
plot_ly(dfChironomidAll, lat = dfChironomidAll$latNum, lon = dfChironomidAll$lonNum, 
        text = hover, color = globalRegion,
        mode = "markers", type = 'scattergeo') %>%
  layout(title = paste0(mapTitle) , geo = mapLayout)
