# Basic analysis of DNA methylation microarray data
In this tutorial we are going to be analysing some publicly available data from the Gene Expression Omnibus (http://www.ncbi.nlm.nih.gov/geo/). In particular will be focusing on a study that contains prostate cancer samples from benign and tumour tissues (GSE47915, http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47915). 

We'll firstly download the raw data which contains the IDAT files using the Bioconductor package GEOquery.

You can download GEOquery using biocLite
```
biocLite("GEOquery")
```
Next load the package into the R environment 
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





