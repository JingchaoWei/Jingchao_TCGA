library(RTCGAToolbox)
valid_aliases <- getFirehoseDatasets()
valid_aliases
stddate <- getFirehoseRunningDates()
stddate
gisticDate <- getFirehoseAnalyzeDates()
gisticDate
brcaData <- getFirehoseData (dataset="BRCA", runDate="20140416",
                             gistic2Date="20140115", clinic=TRUE, RNAseqGene=TRUE, mRNAArray=TRUE,
                             Mutation=TRUE)

data(RTCGASample)
RTCGASample
getData(RTCGASample, "RNASeqGene")
par(mfrow=c(2,2))
diffGeneExprs <- getDiffExpressedGenes(dataObject=RTCGASample, DrawPlots=TRUE,
                                       adj.method="BH", adj.pval=0.05, raw.pval=0.05, logFC=2, hmTopUpN=10,
                                       hmTopDownN=100)
diffGeneExprs
showResults(diffGeneExprs[[1]])
toptableOut <- showResults(diffGeneExprs[[1]])
corrGECN <- getCNGECorrelation(dataObject=RTCGASample, adj.method="BH",
                               adj.pval=0.9, raw.pval=0.05)
corrGECN
showResults(corrGECN[[1]])
corRes <- showResults(corrGECN[[1]])
mutFrq <- getMutationRate(dataObject=RTCGASample)
head(mutFrq[order(mutFrq[, 2], decreasing=TRUE), ])
clinicData <- getData(RTCGASample,"clinical")
head(clinicData)
clinicData <- clinicData[, 3:5]
clinicData[is.na(clinicData[, 3]), 3] <- clinicData[is.na(clinicData[, 3]), 2]
survData <- data.frame(Samples=rownames(clinicData),
                       Time=as.numeric(clinicData[, 3]), Censor=as.numeric(clinicData[, 1]))
getSurvival(dataObject=RTCGASample, geneSymbols=c("FCGBP"), sampleTimeCensor=survData)





