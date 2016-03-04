# Basic analysis of DNA methylation microarray data
In this tutorial we are going to be analysing some publicly available data from the Gene Expression Omnibus (http://www.ncbi.nlm.nih.gov/geo/). In particular we'll be focusing on a study that contains prostate cancer samples from benign and tumour tissues (GSE47915, http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47915). 

We'll firstly download the raw data which contains the IDAT files using the Bioconductor package GEOquery.

You can load a R package into the R environment using ```library(name-of-package)```. For example let's load the ```GEOquery``` package
```
library(GEOquery)
```
We can now download the supplementary data from GSE47915 which contains all the raw data that we need. This may take some time depending on your internet connection
```
getGEOSuppFiles('GSE47915')
GSE47915 <- getGEO(GEO = 'GSE47915')
``` 

The getGEOSuppFiles function creates a folder titled GSE47915 in your working directory which contains GSE47915_RAW.tar which we shall extract in the next step. The getGEO function imports information regarding the study. We can extract the metadata from this object using the command below.
```
pd <- pData(phenoData(GSE47915[[1]]))
```
pd is a data frame that contains the meta data of the study. You should be able to see the data in your Environment. Since it is a small data frame you can click on it to view it or use the command below to have a look at the first 6 lines.
```
head(pd)
```
To extract the IDAT files from GSE47915.tar we must untar the file using the command below.
```
untar(tarfile = 'GSE47915/GSE47915_RAW.tar')
```
This has created 8 Grn.idat.gz and 8 Red.idat.gz files which contains the unmethylated and methylated channels for each of the 8 samples. We must gunzip each file such that we can read the files into the R environment.
```
idat_files <- list.files(pattern = 'idat.gz')
for(i in 1:length(idat_files)){
  gunzip(filename = idat_files[i], destname = gsub("[.]gz$", "", idat_files[i]))
}
```
The idat files have now been gunzip and now we can load the minfi package to import the data. The aim here is to get the methylation level of each CpG site for each sample into a single value. In this tutorial we'll be using the beta value.
```
library(minfi)
library(IlluminaHumanMethylation450kmanifest)
RGset <- read.450k.exp('~/Tutorial/')
MSet <- preprocessRaw(RGset)
BetaMatrix <- getBeta(MSet, type = 'Illumina')
```
The BetaMatrix object is large data frame which contains the non-normalised beta values for the 485,512 probes and 8 samples. Beta values range from 0-1 where 0 and 1 indicates that the CpG site is completely unmethylated and methylated respectively.

Before we can start doing some differential methylation analyses we need to normalise the data and check for any outliers.

The Illumina Infinium HumanMethylation450 BeadChip uses two different probe designs which can create some technical variation. In this tutorial we'll be using a normalisation method referred to as BMIQ (Beta MIxture Quantile dilation) which adjusts the data for the two different probe designs followed by quantile normalisation. This method is implemented in the ChAMP package.
```
library(ChAMP)
Norm_Beta_450k <- champ.norm(beta=BetaMatrix, rgSet=FALSE, pd = pd, mset = FALSE,
                             sampleSheet = FALSE, plotBMIQ = TRUE, 
                             resultsDir = "~/Tutorial/", filterXY = FALSE)
```

The champ.norm function can take many arguments. If you are interested in each of these arguments you can run ```?champ.norm``` in the console of RStudio which will display information on the function in the help window.

