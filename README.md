# DNA methylation microarray tutorial
A tutorial for the UofA Bioinformatics Hub and EpicSA
## Basic R script to analyse Illumina Infinium HumanMethylation450 BeadChip data
Created by: Benjamin Mayne (benjamin.mayne@adelaide.edu.au)

### Requirements
This tutorial will be conducted using R, which is a free software used in statistical computing. In order to complete the tutorial you must download the latest version of R from here (https://cran.r-project.org/). Make sure to download the correct version depending on whether you are using a Windows, MAC or a UNIX system.

In addition, also download the latest version of RStudio which can be found here (https://www.rstudio.com/home/). RStudio is a more user friendly interactive version of R. Once you have installed both R and R studio open up RStudio.

R works using R packages which can be installed and loaded into the R environment. In this tutorial we'll be using a type of R packages called Bioconductor packages (http://bioconductor.org/). Bioconductor packages are open-sourced software packages developed by bioinformaticians. To install Bioconductor packages you must firstly install Bioconductor into R. This can be done by typing the command below into R.
```
source("http://bioconductor.org/biocLite.R")
biocLite()
```
For more information on how to install Bioconductor refer to the following link (http://bioconductor.org/install/).

Before starting the tutorial you'll need to download these Bioconductor packages. Bioconductor packages can be install using the ```biocLite()``` command. The following Bioconductor packages that are needed are listed below.
```
biocLite('GEOquery')
biocLite('minfi')
biocLite('IlluminaHumanMethylation450kmanifest')
biocLite('ChAMP')
biocLite('limma')
biocLite('IlluminaHumanMethylation450k.db')
```
