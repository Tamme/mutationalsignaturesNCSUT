#setwd("~/Desktop/RFN/")
setwd("C://Users/Tamme/Desktop/Yli/Magistritöö/github/mutationalsignaturesNCSUT/")
library("SomaticSignatures")
library(ggplot2)
library(R.matlab)
library(lsa)
load(file = "mutationTypes.rda") #variabe mutationTypes
load(file = "dataMmRownamesRightOrder.rda")
source("plot_functions.R")

results = matrix(0, nrow = 3, ncol=15)
rawData = readMat("RFN/data/119_genomes_96_subs_data.mat")
rawData = readMat("506large96SubsData.mat")
nrOfSamples = 119
#dataOrig =  matrix(unlist(rawData[1]), ncol = nrOfSamples, byrow = FALSE)
dataOrig = rawData$originalGenomes
for (sig in 1:15) {
  sig = 5
  raw_data = readMat(paste0("RFN/models/", nrOfSamples, "_breast_C_L1_0_", sig, "_100_out.mat"))
  W= raw_data$W
  H= raw_data$H
  R =raw_data$R
  diff = dataOrig - W %*% H
  frobenius = sqrt(sum(diag(t(diff) %*% diff)))
  print(frobenius)
  #results[3, sig] = frobenius
}

data = W
rownames(data) = dataMmRownamesRightOrder
data = data[mutationTypes,]
plot = plotMutationalSignatures(data, 0.25)
plot

cor(W) 
cosine(W)
max(cosine(W) - diag(dim(W)[2]))


#Snippet for calculating sparsities
rawData = readMat("data/21_genomes_96_subs_data.mat")
rawData = readMat("data/119_genomes_96_subs_data.mat")
rawData = readMat("data/506_genomes_96_subs_data.mat")
nrOfSamples = 506
nrOfSignatures = 7
dataOrig = rawData$originalGenomes

ctr = 1
sparsitySequence = seq(0, 3.5, 0.05)
resultArray = array(data = 0, dim = c(dim(dataOrig)[1], nrOfSignatures, length(sparsitySequence)), dimnames = NULL)
sparsityResults = matrix(0, nrow= length(sparsitySequence), ncol=3)
for (sparsity in sparsitySequence) {
  sparsity
  raw_data = readMat(paste0("models/", nrOfSamples, "_breast_C_L1_S_", sprintf("%1.1f", sparsity), "_", nrOfSignatures, "_50_out.mat"))
  W= raw_data$W
  H= raw_data$H
  R =raw_data$R
  
  diff = dataOrig - R
  frobenius = sqrt(sum(diag(t(diff) %*% diff)))
  spar = sum(W==0)/(dim(W)[1]*dim(W)[2])
  spar
  print(paste(frobenius, sum(W==0)/(4*96)))
  sparsityResults[ctr, 1] = round(frobenius, 2)
  sparsityResults[ctr, 2] = spar
  sparsityResults[ctr, 3] = sparsity
  resultArray[, , ctr] = W
  ctr = ctr + 1
}
zeroSparsity = sparsityResults[1,]
#frobenius over 15% off
sparsityResults[sparsityResults[,1]>1.5*zeroSparsity[1]] = c(0, 0, 0)
#order by resulting sparsity decreasing
idxes = order(sparsityResults[,2], decreasing = TRUE)
#idxes = c(1,2)
for (index in idxes) {
  if (sparsityResults[index,1] != 0) {
    print(index)
    data = resultArray[ , , index]
    rownames(data) = dataMmRownamesRightOrder
    data = data[mutationTypes,]
    
    res = new("MutationalSignatures",
              signatures = data,
              samples = as.matrix(1:13),
              fitted = as.matrix(0),
              observed = as.matrix(0),
              nSignatures = dim(data)[2]
    )
    #cor(W) 
    print(max(cor(data) - diag(dim(data)[2])))
    #cosine(W)
    print(max(cosine(data) - diag(dim(data)[2])))
    print(sparsityResults[index,])
    #x11()
    paste(sparsityResults[index,])
    Sys.sleep(0)
    print(plotSignatures(res, normalize = TRUE, percent=TRUE) + ggtitle(paste(sparsityResults[index,1], sparsityResults[index,2], sparsityResults[index,3])) +  coord_cartesian(ylim = c(0, 25)))
    readline()
  }
}
sparsityResults = sparsityResults[order(sparsityResults[,2], decreasing = TRUE),]
sparsityResults
writeMat("RFN_1121breast_5sigs_sparsity.mat", sparsityResults = sparsityResults)

#Normalized proportions
library(reshape2)
library(RColorBrewer)
h = matrix(unlist(H), ncol = dim(H)[2], byrow = FALSE)
H = t(H)
normalizedH = H/rowSums(H)
rowSums(normalizedH)

mydf2 = as.data.frame(normalizedH)
mydf2$Samples = row.names(mydf2)
mydf2.molten <- melt(mydf2, value.name="Identity", variable.name="Signature")
colourCount = dim(H)[2]
getPalette = colorRampPalette(brewer.pal(5, "Set3"))
p = ggplot( data=mydf2.molten, aes(x = Samples, y = value, fill=variable))
p +  geom_bar(stat = "identity") + scale_fill_manual(values = getPalette(colourCount)) + ggtitle(paste0("RFN C impl ", dim(normalizedH)[2], " signatures proportions H"))
