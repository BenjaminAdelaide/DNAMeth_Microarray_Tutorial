# Install the required Bioconductor packages

source("http://bioconductor.org/biocLite.R")
biocLite('GEOquery')
biocLite('minfi')
biocLite('IlluminaHumanMethylation450kmanifest')
biocLite('ChAMP')
biocLite('limma')
biocLite('IlluminaHumanMethylation450k.db')

# Create a directory called Tutorial and move into it
dir.create('Tutorial')
setwd('Tutorial/')

# Load GEOquery
library(GEOquery)

# Download the data set
getGEOSuppFiles('GSE47915')
GSE47915 <- getGEO(GEO = 'GSE47915')

# Extract the meta data
pd <- pData(phenoData(GSE47915[[1]]))

# Take a look at the meta data
head(pd)

# Extract the IDAT files
untar(tarfile = 'GSE47915/GSE47915_RAW.tar')

# Gunzip the IDAT files
idat_files <- list.files(pattern = 'idat.gz')
for(i in 1:length(idat_files)){
    gunzip(filename = idat_files[i], destname = gsub("[.]gz$", "", idat_files[i]))
}

# Preprocess and obtain the beta values
library(minfi)
library(IlluminaHumanMethylation450kmanifest)
RGset <- read.450k.exp('~/Tutorial/')
MSet <- preprocessRaw(RGset)
BetaMatrix <- getBeta(MSet, type = 'Illumina')

# Normalisation
library(ChAMP)
BetaMatrix_norm <- champ.norm(beta=BetaMatrix, rgSet=FALSE, pd = pd, mset = FALSE,
                              sampleSheet = FALSE, plotBMIQ = TRUE, 
                              resultsDir = "~/Tutorial/", filterXY = FALSE)
datMeth <- BetaMatrix_norm$beta

# MDS plot
pdf('MDS_plot.pdf', width=8, height=8)
mdsPlot(dat=datMeth, sampGroups=pd$characteristics_ch1.2)
dev.off()

# Make sure the samples match up
colnames(datMeth) <- substr(x=colnames(datMeth), start=0, stop=10)
pd <- pd[match(colnames(datMeth), row.names(pd)),]
all(colnames(datMeth) == row.names(pd))

# Setup the design matrix
group <- factor(pd$characteristics_ch1.2,
                levels=c("tissue: benign prostate tissues","tissue: prostate cancer tumor"))
design <- model.matrix(~group)

# Differential methylation
fit.reduced <- lmFit(datMeth,design)
fit.reduced <- eBayes(fit.reduced)
top <- topTable(fit.reduced, number=50)

# Heatmap
datMatrix <- datMeth[match(row.names(top), row.names(datMeth)),]
colnames(datMatrix) <- ifelse(test=pd$characteristics_ch1.2 == 'tissue: benign prostate tissues', 
                              yes='Benign', no='Tumour')
pdf('Heatmap.pdf', width=8, height=8)
heatmap(t(datMatrix))
dev.off()

# Annotation
library(IlluminaHumanMethylation450k.db)

# Extract the chromosome data

chromosomes <- as.list(IlluminaHumanMethylation450kCHR)
chromosomes <- data.frame(unlist(chromosomes))
colnames(chromosomes) <- 'CHR'

# Extract the co-ordinates

location <- as.list(IlluminaHumanMethylation450kCPGCOORDINATE)
location <- data.frame(unlist(location))
colnames(location) <- 'Co_ordinate'
# Extract the genes

symbol <- as.list(IlluminaHumanMethylation450kSYMBOL)
symbol <- data.frame(unlist(symbol))
colnames(symbol) <- 'ID'

# Match the probes with CpG site locations

top$Chr <- chromosomes$CHR[match(row.names(top), row.names(chromosomes))]
top$Position <- location$Co_ordinate[match(row.names(top), row.names(location))]
top$Gene <- symbol$ID[match(row.names(top), row.names(symbol))]

# Export the results
write.csv(top, file = 'Diff_meth.csv')
