dir.create('Tutorial')
setwd('Tutorial/')

library(GEOquery)

getGEOSuppFiles('GSE47915')
GSE47915 <- getGEO(GEO = 'GSE47915')
pd <- pData(phenoData(GSE47915[[1]]))
head(pd)

untar(tarfile = 'GSE47915/GSE47915_RAW.tar')
idat_files <- list.files(pattern = 'idat.gz')

for(i in 1:length(idat_files)){
  gunzip(filename = idat_files[i], destname = gsub("[.]gz$", "", idat_files[i]))
}

library(minfi)
library(IlluminaHumanMethylation450kmanifest)

RGset <- read.450k.exp('~/Tutorial/')
MSet <- preprocessRaw(RGset)
BetaMatrix <- getBeta(MSet, type = 'Illumina')

library(ChAMP)
BetaMatrix_norm <- champ.norm(beta=BetaMatrix, rgSet=FALSE, pd = pd, mset = FALSE,
                             sampleSheet = FALSE, plotBMIQ = TRUE, 
                             resultsDir = "~/Tutorial/", filterXY = FALSE)

datMeth <- BetaMatrix_norm$beta

pdf('MDS_plot.pdf', width=8, height=8)
mdsPlot(dat=datMeth, sampGroups=pd$characteristics_ch1.2)
dev.off()

library(limma)

colnames(datMeth) <- substr(x=colnames(datMeth), start=0, stop=10)
pd <- pd[match(colnames(datMeth), row.names(pd)),]
all(colnames(datMeth) == row.names(pd))

group <- factor(pd$characteristics_ch1.2,
                levels=c("tissue: benign prostate tissues","tissue: prostate cancer tumor"))
design <- model.matrix(~group)
fit.reduced <- lmFit(datMeth,design)
fit.reduced <- eBayes(fit.reduced)
top <- topTable(fit.reduced, number=50)

datMatrix <- datMeth[match(row.names(top), row.names(datMeth)),]
colnames(datMatrix) <- ifelse(test=pd$characteristics_ch1.2 == 'tissue: benign prostate tissues', 
                              yes='Benign', no='Tumour')
pdf('Heatmap.pdf', width=8, height=8)
heatmap(t(datMatrix))
dev.off()

library(IlluminaHumanMethylation450k.db)
  
chromosomes <- as.list(IlluminaHumanMethylation450kCHR)
chromosomes <- data.frame(unlist(chromosomes))
colnames(chromosomes) <- 'CHR'

location <- as.list(IlluminaHumanMethylation450kCPGCOORDINATE)
location <- data.frame(unlist(location))
colnames(location) <- 'Co_ordinate'

symbol <- as.list(IlluminaHumanMethylation450kSYMBOL)
symbol <- data.frame(unlist(symbol))
colnames(symbol) <- 'ID'

top$Chr <- chromosomes$CHR[match(row.names(top), row.names(chromosomes))]
top$Position <- location$Co_ordinate[match(row.names(top), row.names(location))]
top$Gene <- symbol$ID[match(row.names(top), row.names(symbol))]

write.csv(top, file = 'Diff_meth.csv')

