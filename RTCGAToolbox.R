library(RTCGAToolbox)
valid_aliases <- getFirehoseDatasets()
valid_aliases
stddate <- getFirehoseRunningDates()
stddate
gisticDate <- getFirehoseAnalyzeDates()
gisticDate
PRAD_data <- getFirehoseData(dataset="PRAD", runDate="20151101",forceDownload = F,
                             clinical=TRUE, RNASeq2GeneNorm = T,Mutation = T)
data(RTCGASample)
RTCGASample
getData(RTCGASample, "RNASeqGene")
