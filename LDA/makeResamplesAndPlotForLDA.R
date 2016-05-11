#R Version R-3.2.2
#script for creating data sets in correct format for lda-c

library(R.matlab)
#set global folder
#in linux
setwd("~/...")
#in windows
setwd("C:/Users/...")
library("SomaticSignatures")
library(ggplot2)
load(file = "mutationTypes.rda") #variabe mutationTypes
load(file = "dataMmRownamesRightOrder.rda")
source("plot_functions.R")

rawData = readMat("LDA/data/21_genomes_96_subs_data.mat")
rawData = readMat("LDA/data/119_genomes_96_subs_data.mat")
rawData = readMat("LDA/data/506_genomes_96_subs_data.mat")
dataOrig = rawData$originalGenomes
rownames(dataOrig) = dataMmRownamesRightOrder
dataOrig = dataOrig[mutationTypes,]
#Turn into topic modelling input format
nrOfSamples = dim(dataOrig)[2]
nrOfTypes = dim(dataOrig)[1]

typeDictionary <- vector(mode="list", length=nrOfTypes)
names(typeDictionary) <- rownames(dataOrig)
for (i in 1:nrOfTypes) {
  typeDictionary[[i]] = i
}

nrOfSampledSets = 50
for (ctr in 1:nrOfSampledSets) {
  data_mm_boot = matrix(nrow=nrOfTypes, ncol=nrOfSamples, 0)
  for (sampl in 1:nrOfSamples) {
    prob = dataOrig[,sampl]/sum(dataOrig[,sampl])
    n = sum(dataOrig[,sampl])
    names(prob) = NULL
    data_mm_boot[,sampl] = rowSums(rmultinom(n, 1, prob))
  }
  rownames(data_mm_boot) = rownames(dataOrig)
  colnames(data_mm_boot) = colnames(dataOrig)

  data = data_mm_boot
  
  filename = paste0("data/lda_original_506_breast", ctr, ".txt")
  file.create(filename)
  print(ctr)
  for (i in 1:nrOfSamples) {
    nonZeroIndexes = which(data[,i] != 0)
    
    row = paste(colSums(data!=0)[i], paste(typeDictionary[rownames(data[nonZeroIndexes,])], data[nonZeroIndexes,i], sep=":", collapse=" "))
    if (length(nonZeroIndexes) == 1) {
      row = paste(colSums(data!=0)[i], paste(typeDictionary[names(nonZeroIndexes)], data[nonZeroIndexes,i], sep=":", collapse=" "))
    }
    write(row, file = filename, ncolumns = 1, append = TRUE, sep = "\n")
  }
}

#output 

#retrieve H and W from topic modelling lagorithms
result = matrix(0, ncol = 12, nrow = 1)
for (sig in 2:12) {
  sig=4
  dataW <- read.table(paste0("LDA/lda-c/21breast_rs_", sig, "sigs_full_moreconv/final.beta.prob.txt"), header=FALSE,sep=" ")
  dataH <- read.table(paste0("LDA/lda-c/21breast_rs_", sig, "sigs_full_moreconv/final.gamma"), header=FALSE,sep=" ")
  
  H = as.matrix(t(dataH))
  dataW = dataW[,-1]
  W = t(dataW)
  #now W and H are in correct format with shapes featuresXpatterns and patternsXsamples
  R = W %*% H
  diff = R - dataOrig
  frobenius = sqrt(sum(diag(t(diff) %*% diff)))
  print(frobenius)
  result[1, sig] = frobenius
  rownames(W) = mutationTypes
  data = W
}
plot = plotMutationalSignatures(data, 0.3)
plot

