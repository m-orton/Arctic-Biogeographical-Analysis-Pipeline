# R scripts for the Analysis of Arctic Community Similarity using DNA Barcoding Data

### Prerequisite Files

Please see files in the Chironomid Datasets and Arctic Shapefile folders.

For convenience, the links to the files have been provided here.

To run either of the biogeography analysis pipelines with the Chironomid datasets, please ensure you have the following files in your 
current working R directory:

Public BOLD data from October/2017 download:

[Greenland](ChironomidDatasets/dfGreenland_Oct17.csv)

dfNearctic_Oct17.csv ([Nearctic](ChironomidDatasets/dfNearctic_Oct17.7z) must be unzipped first with 7zip software)

[Palearctic](ChironomidDatasets/dfPalearctic_Oct17.csv)

For using the datasets provided by Torbjørn Ekrem and Elisabeth Stur (datasets are now publicly available on BOLD):

[Taxonomic Data](ChironomidDatasets/Private_Chironomid_Data_ModifiedSingleSheet.csv)

[Sequence_Data](ChironomidDatasets/PrivateSequenceData.fas)

[Curated Greenland Species Data](ChironomidDatasets/Greenland records chironomid (1).csv) 

Multiple sequence alignments for single-linkage clustering with a 4, 4.5 or 5% cutoff threshold:

[Full Datasets Alignment](ChironomidDatasets/ChironomidAlignmentMay1_AllCanNor_2nd.fas)

[Subacrctic Filtered Alignment](ChironomidDatasets/ChironomidAlignmentApril24_Subarctic.fas)

To run the second pipeline: [Subarctic Pipeline](ChironomidBiogeographySubarcticFilter.R) please ensure the following shapefiles are also included in the current working R directory:

[.dbf file](ArcticShapefiles/Arctic_Zones.dbf)

[.prj file](ArcticShapefiles/Arctic_Zones.prj)

[.sbn file](ArcticShapefiles/Arctic_Zones.sbn)

[.sbx file](ArcticShapefiles/Arctic_Zones.sbx)

[.shp file](ArcticShapefiles/Arctic_Zones.shp)

[.shp.xml file](ArcticShapefiles/Arctic_Zones.shp.xml)

[.shx file](ArcticShapefiles/Arctic_Zones.shx)

### Installation

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

## Acknowledgments
Dr. Torbjørn Ekrem and Dr. Elisabeth Stur for kindly providing their Chironomid datasets.
Dr. Sally Adamowicz for her work on the analysis pipelines.

