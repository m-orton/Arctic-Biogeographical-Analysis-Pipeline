###############
# Identification of outlier points based on genetic
# distance, based on spread, will eliminate records which are clearly not
# belonging to Chironomidae if by mistake a record was misindentified

# Authored by Matthew G. Orton

# Packages used
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

# After selection of sequence representatives per BIN

dnaStringSet <- DNAStringSet(dfSingleSeq$nucleotides)
# For testing just set to maxiters=2 and diags=TRUE
alignment <- muscle(dnaStringSet, maxiters = 2, diags = TRUE)
# Back to a string set
dnaStringSet2 <- DNAStringSet(alignment)

# Name the stringset with record ids
bin_uri <- dfSingleSeq$bin_uri
names(dnaStringSet2) <- bin_uri

# DNAbin format
dnaBin <- as.DNAbin(alignment)
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
