#created by Lauri Tammeveski

setwd("C:/...")
library(ggplot2)
library(R.matlab)
library(reshape2)
library(RColorBrewer)
load(file = "mutationTypes.rda") #variabe mutationTypes
load(file = "dataMmRownamesRightOrder.rda")
source("plot_functions.R")

nrOfMutations = 96
nrOfPatterns = 5
nrOfSamples = 500
#making signatures
patterns = matrix(0, nrow = nrOfPatterns, ncol = nrOfMutations)
for (i in 1:(nrOfPatterns)) {
  #a handcrafted way of making overlapping regions to make well distinguishable signatures
  interactionLength = floor((nrOfMutations / (nrOfPatterns)))
  interactionLength2 = floor(1.3 * (nrOfMutations / (nrOfPatterns)))
  p = rep(0, nrOfMutations)
  if (i != nrOfPatterns) {
    r = abs(rnorm(interactionLength2, mean = i/interactionLength, sd = 0.1))
    p[(((i-1)*interactionLength)+1):((((i-1)*interactionLength)+interactionLength2))] = r #no overlapping  
  } else {
    r = abs(rnorm(interactionLength, mean = i/interactionLength, sd = 0.1))
    p[(((i-1)*interactionLength)+1):(((i-1)*interactionLength)+interactionLength)] = r #no overlapping  
    
  }
  
  #add general small noise everywhere
  noise = runif(nrOfMutations, 0, 0.03) #decrease if more mutations
  #add some noise to 5 sampled mutation types
  p[sample(1:96, 5)] = abs(rnorm(5, sd=0.5))
  p = p + noise
  p = p /sum(p)
  patterns[i,] = p
}
patterns=t(patterns)
rownames(patterns) = mutationTypes[1:nrOfMutations]
plot = plotMutationalSignatures(patterns, 0.25)
plot

#making proportions
proportions = matrix(0, nrow = nrOfSamples, ncol = nrOfPatterns)
for (i in 1:nrOfSamples) {
  r = abs(rnorm(nrOfPatterns, mean = 0, sd = 100))
  rProb = r / sum(r)
  proportions[i,] = rProb
}
#plot some of them for viewing
mydf2 = as.data.frame(proportions[1:180, ])
mydf2$Samples = row.names(mydf2)
mydf2.molten <- melt(mydf2, value.name="Identity", variable.name="Signature")
colorCount = dim(proportions)[2]
getPalette = colorRampPalette(brewer.pal(9, "Set3"))
p = ggplot( data=mydf2.molten, aes(x = Samples, y = value, fill=variable))
p +  geom_bar(stat = "identity") + scale_fill_brewer(palette = "Set3")  + ggtitle("Synthetic data sample proportions")


countMatrix = matrix(0, nrow=nrOfMutations, ncol=nrOfSamples)
countMatrixForPatterns = matrix(0, nrow=nrOfMutations, ncol=nrOfPatterns)
totalCountPerSample = 2000
for (samp in 1:nrOfSamples) {
  x = sample(1:nrOfPatterns, totalCountPerSample, replace=TRUE, prob=proportions[samp, ])
  addCounter = unlist(lapply(x, function(z) sample(1:nrOfMutations, 1, prob=patterns[,z])))
  for (a in addCounter) {
    countMatrix[a, samp] = countMatrix[a, samp] + 1
  }
  for (idx in 1:length(addCounter)) {
    countMatrixForPatterns[addCounter[idx], x[idx]] = countMatrixForPatterns[addCounter[idx], x[idx]] + 1
  }
}
countMatrix2 = countMatrix
counts = t(t(countMatrixForPatterns) / colSums(countMatrixForPatterns))

#counts
#rownames(counts) = mutationTypes
#plot = plotMutationalSignatures(counts, 0.25)
#plot
#colSums(counts)

#head(t(patterns))
#head(counts)
#head(data)
diff = patterns - counts
frobenius = sqrt(sum(diag(t(diff) %*% diff)))
frobenius

Y = countMatrix2
W = patterns
X = proportions
writeMat(paste0("Autoencoder/data/EasierHsynthetic_5sigs_96mut_500samples_2000tot_", totalCountPerSample, "samp.mat"), X=X, originalGenomes=Y, W=W, sampledW=counts, sampledWtotal=countMatrixForPatterns)