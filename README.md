# DNA methylation microarray tutorial
A tutorial for the UofA Bioinformatics Hub
## Basic R script to analyse Infinium HumanMethylation450 BeadChip data
Created by: Benjamin Mayne (benjamin.mayne@adelaide.edu.au)

### Requirements
This tutorial will be conducted in RStudio (https://www.rstudio.com/home/) and requires R (https://cran.r-project.org/) and several Bioconductor packages (http://bioconductor.org/).

Before starting the tutorial you'll need to download these Bioconductor packages. To install biocLite and for more information on how to install Bioconductor packages visit the website (https://www.bioconductor.org/install/)
```
biocLite('GEOquery')
biocLite('minfi')
biocLite('IlluminaHumanMethylation450kmanifest')
biocLite('ChAMP')
biocLite('limma')
biocLite('IlluminaHumanMethylation450k.db')
```
