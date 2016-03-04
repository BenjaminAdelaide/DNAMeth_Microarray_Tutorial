# Basic analysis of DNA methylation microarray data
In this tutorial we are going to be analysing some publicly available data from the Gene Expression Omnibus (http://www.ncbi.nlm.nih.gov/geo/). In particular we'll be focusing on a study that contains prostate cancer samples from benign and tumour tissues (GSE47915, http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47915). 

We'll firstly download the raw data which contains the IDAT files using the Bioconductor package GEOquery.

You can load a R package into the R environment using ```library(name-of-package)```. For example let's load the ```GEOquery``` package.
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
BetaMatrix_norm <- champ.norm(beta=BetaMatrix, rgSet=FALSE, pd = pd, mset = FALSE,
                             sampleSheet = FALSE, plotBMIQ = TRUE, 
                             resultsDir = "~/Tutorial/", filterXY = FALSE)
datMeth <- BetaMatrix_norm$beta
```

The champ.norm function can take many arguments. If you are interested in each of these arguments you can run ```?champ.norm``` in the console of RStudio which will display information on the function in the help window. The datMeth object is a matrix which contains the normalised beta values.

The DNA methylation data can be checked for potential outliers by analysing a multidimensional scaling (MDS) plot of the 1000 most variable positions or probes. This can also be used as an informative method to see if any sample groups cluster separately based on their DNA methylation. 
```
pdf('MDS_plot.pdf', width=8, height=8)
mdsPlot(dat=datMeth, sampGroups=pd$characteristics_ch1.2)
dev.off()
```
The above command will output a pdf into your working directory which should look like the image below.
![](https://cloud.githubusercontent.com/assets/10754973/13515017/ce3e681e-e1fd-11e5-9e57-b6d8604e3d95.png)

As you can see from the plot it would suggest that the benign and tumour prostate samples cluster separately. However, we are working with a small sample size so it is difficult to draw any valid conclusions.

For differential methylation analysis we will be using ```limma``` which is commonly used in differential expression analyses in gene expression microarrays. 

We'll firstly create model matrix which will contain the grouping information for each sample. However, let's make sure the order of samples in datMeth matches the order in pd.

Use the ```substr``` function to trim letters off the column names such that they'll match the sample IDs in pd. The use ```match``` to reorder pd to match the order of the samples in datMeth.
```
colnames(datMeth) <- substr(x=colnames(datMeth), start=0, stop=10)
pd <- pd[match(colnames(datMeth), row.names(pd)),]
all(colnames(datMeth) == row.names(pd))
```
If the last line of code returns ```TRUE``` then  the sample order matches the order between datMeth and pd.

We can now create the model matrix. In this analysis we'll be comparing benign to tumour samples which is provided in the characteristics_ch1.2 of pd which can be viewed by typing ```pd$characteristics_ch1.2$```.
```
group <- factor(pd$characteristics_ch1.2,
                levels=c("tissue: benign prostate tissues","tissue: prostate cancer tumor"))
design <- model.matrix(~group)
```
Using limma we can now produce a list of differentially methylated CpG sites.
```
fit.reduced <- lmFit(datMeth,design)
fit.reduced <- eBayes(fit.reduced)
top <- topTable(fit.reduced, number=50)
```
The top object now contains the top 50 most significantly differentially methylated CpG sites. To view more or less change the number is the ```topTable``` function.

We can now make a heatmap to illustrate the differences in methylation of the top 50 CpG sites between the two groups. 
```
datMatrix <- datMeth[match(row.names(top), row.names(datMeth)),]
colnames(datMatrix) <- ifelse(test=pd$characteristics_ch1.2 == 'tissue: benign prostate tissues', 
                              yes='Benign', no='Tumour')
pdf('Heatmap.pdf', width=8, height=8)
heatmap(t(datMatrix))
dev.off()
```
![](https://cloud.githubusercontent.com/assets/10754973/13515613/d3b14d16-e202-11e5-835a-4e8a6aeba0d7.png)
The heatmap function also performs unsupervised clustering of the probes and samples.
