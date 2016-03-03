# Basic analysis of DNA methylation microarray data
In this tutorial we are going to be analysing some publicly available data from the Gene Expression Omnibus (http://www.ncbi.nlm.nih.gov/geo/). In particular will be focusing on a study that contains prostate cancer samples from benign and tumour tissues (GSE47915, http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47915). 

We'll firstly download the raw data which contains the IDAT files using the Bioconductor package GEOquery.

You can download GEOquery using biocLite
```
biocLite("GEOquery")
```

