# DNA methylation microarray tutorial
A tutorial for the UofA Bioinformatics Hub
## Basic R script to analyse Illumina Infinium HumanMethylation450 BeadChip data
Created by: Benjamin Mayne (benjamin.mayne@adelaide.edu.au)

### Requirements
This tutorial will be conducted in RStudio (https://www.rstudio.com/home/) which requires R (https://cran.r-project.org/) and several Bioconductor packages (http://bioconductor.org/).

Before starting the tutorial you'll need to download these Bioconductor packages. For more information on how to use biocLite and installing Bioconductor packages visit the Bioconductor website (https://www.bioconductor.org/install/)
```
biocLite('GEOquery')
biocLite('minfi')
biocLite('IlluminaHumanMethylation450kmanifest')
biocLite('ChAMP')
biocLite('limma')
biocLite('IlluminaHumanMethylation450k.db')
```
