# R scripts for the Analysis of Arctic Community Similarity using DNA Barcoding Data

## Prerequisite Files

Please see files in the Chironomid Datasets and Arctic Shapefile folders.

For convenience, the links to the files have been provided here as well.

To run either of the biogeography analysis pipelines with the Chironomid datasets, please ensure you have the following files in your 
current working R directory:

### Public BOLD data from October/2017 download:

[Greenland](ChironomidDatasets/dfGreenland_Oct17.csv)

dfNearctic_Oct17.csv ([Nearctic](ChironomidDatasets/dfNearctic_Oct17.7z) must be unzipped first with 7zip software)

[Palearctic](ChironomidDatasets/dfPalearctic_Oct17.csv)

### Alternatively current BOLD data can be used via a direct download from the BOLD API by running the commands:

```
dfNearctic <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Chironomidae&geo=Alaska|Canada&format=tsv")
dfGreenland <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Chironomidae&geo=Greenland&format=tsv")
dfPalearctic <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combinedtaxon=Chironomidae&geo=Norway|
                          Denmark|Iceland|Sweden|Finland&format=tsv")
```

### For using the datasets provided by Torbjørn Ekrem and Elisabeth Stur (datasets are now publicly available on BOLD):

[Taxonomic Data](ChironomidDatasets/Private_Chironomid_Data_ModifiedSingleSheet.csv)

[Sequence Data](ChironomidDatasets/PrivateSequenceData.fas)

[Curated Greenland Species Data](ChironomidDatasets/Greenland_records_chironomid.csv) 

### Multiple sequence alignments for single-linkage clustering with a 4, 4.5 or 5% cutoff threshold:

[Full Dataset Alignment](ChironomidDatasets/ChironomidAlignmentMay1_AllCanNor_2nd.fas)

[Subarctic Filtered Alignment](ChironomidDatasets/ChironomidAlignmentApril24_Subarctic.fas)


### To run the [Subarctic Shapefile Filtered Pipeline](ChironomidBiogeographySubarcticFilter.R) please ensure the following shapefiles are also included in the current working R directory:

[.dbf](ArcticShapefiles/Arctic_Zones.dbf)

[.prj](ArcticShapefiles/Arctic_Zones.prj)

[.sbn](ArcticShapefiles/Arctic_Zones.sbn)

[.sbx](ArcticShapefiles/Arctic_Zones.sbx)

[.shp](ArcticShapefiles/Arctic_Zones.shp)

[.shp.xml](ArcticShapefiles/Arctic_Zones.shp.xml)

[.shx](ArcticShapefiles/Arctic_Zones.shx)

## Installation

Please ensure the following packages are installed in R/RStudio:

```
install.packages("foreach")
install.packages("ape")
install.packages("readr")
source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings")
biocLite("muscle")
biocLite("DECIPHER")
install.packages("plotly")
install.packages ("ggplot2") 
install.packages("raster")
install.packages("rgdal")
install.packages("rgeos")
install.packages("vegan")
install.packages("tidyr")
install.packages("dplyr")
install.packages("data.table")
install.packages("vegan")
```

## Authors of Pipeline

Matthew Orton

Dr. Sally Adamowicz

## Acknowledgments
Dr. Torbjørn Ekrem and Dr. Elisabeth Stur for kindly providing their Chironomid datasets.


