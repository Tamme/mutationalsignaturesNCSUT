setwd("~/Desktop/RFN/")
setwd("C://Users/Tamme/Desktop/Yli/Magistritöö/github/mutationalsignaturesNCSUT/")
library("SomaticSignatures")
library(ggplot2)
library(R.matlab)
library(lsa)
load(file = "mutationTypes.rda") #variabe mutationTypes
load(file = "dataMmRownamesRightOrder.rda")
source("plot_functions.R")

results = matrix(0, nrow = 3, ncol=15)

nrOfRealMutations = 2000
rawData = readMat(paste0("RFN/data/synthetic_5sigs_96mut_2000tot_", nrOfRealMutations, "samp.mat"))
dataOrigW =  rawData$sampledW
dataOrig = rawData$originalGenomes

for (sig in 3:8) {
  sig = 5
  raw_data = readMat(paste0("RFN/models/synthetic_", nrOfRealMutations, "_C_L1_0_", sig, "_100_out.mat"))
  W= raw_data$W
  H= raw_data$H
  R =raw_data$R
  diff = dataOrig - W %*% H
  diff=t(dataOrigW - data)
  frobenius = sqrt(sum(diag(t(diff) %*% diff)))
  print(frobenius)
  sum(data)
  data = t(t(data)/colSums(data))
}

data = W
rownames(data) = mutationTypes
plot = plotMutationalSignatures(data, 0.25)
plot

#check correlations and other things
cor(W) 
max(cor(W) - diag(dim(W)[2]))
cosine(W)
max(cosine(W) - diag(dim(W)[2]))