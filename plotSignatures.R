#Instructions
#Weights should be with shape (features, patterns) and H should be (patterns, samples)

setwd("C:/Users/.../mutationalsignaturesNCSUT/")
library("SomaticSignatures")
library(ggplot2)
library(R.matlab)
library(reshape2)
library(RColorBrewer)
library(ggplot2)
library(lsa)
source("plot_functions.R")
load(file = "mutationTypes.rda") #variabe mutationTypes
load(file = "dataMmRownamesRightOrder.rda")

rawData = readMat("data/21_genomes_96_subs_data.mat")
dataOrig21 =  rawData$originalGenomes
rawData = readMat("data/119_genomes_96_subs_data.mat")
dataOrig119 =  rawData$originalGenomes
rawData = readMat("data/506_genomes_96_subs_data.mat")
dataOrig506 =  rawData$originalGenomes

if (isAE == TRUE) {
  nrOfSamples=506
  input = readMat(paste0("Autoencoder/data/", nrOfSamples, "_genomes_96_subs_data.mat"))
  for (sig in 1:15) {
    sig=7
    raw_data = readMat(paste0("Autoencoder/models/AE_", nrOfSamples, "_sigs_", sig, ".mat"))
    W= raw_data$W
    H= raw_data$H
    H = H * colSums(W)
    W = t(t(W) / colSums(W))
    diff = dataOrig506 - W %*% H
    frobenius = sqrt(sum(diag(t(diff) %*% diff)))
    print(frobenius)
    data = W
    rownames(data) = dataMmRownamesRightOrder
    data = data[mutationTypes,]
    plot = plotMutationalSignatures(data, 0.25)
    plot
    #ggsave(filename = paste0("../final_plots/21_ae_4sigs.png"), plot = plot, path = NULL, scale = 1.5, width = 350, height = 175, units = "mm", dpi = 300)
  }
}
if (isNMF == TRUE) {
  raw_data = readMat("../../path_to_models.mat")
  W= raw_data$processes
  H= raw_data$exposures
  data = W
  diff = dataOrig506 - W %*% H
  frobenius = sqrt(sum(diag(t(diff) %*% diff)))
  print(frobenius)
  rownames(data) = dataMmRownamesRightOrder
  data = data[mutationTypes,]
  plot = plotMutationalSignatures(data, 0.25)
  plot
  #ggsave(filename = paste0("../final_plots/21_breast_nmf_4sigs_final.png"), plot = last_plot(), path = NULL, scale = 1.5, width = 350, height = 175, units = "mm", dpi = 300)
}
if (isLDA == TRUE) {
  dataW <- read.table("LDA/lda-c/119breast_rs_6sigs_full_moreconv/final.beta.prob.txt", header=FALSE,sep=" ")
  dataH <- read.table(paste0("LDA/lda-c/119breast_rs_6sigs_full_moreconv/final.gamma"), header=FALSE,sep=" ")
  H = as.matrix(t(dataH))
  dataW = dataW[,-1]
  W = t(dataW)
  
  rownames(data) = mutationTypes
  plot = plotMutationalSignatures(data, 0.25)
  plot
  #ggsave(filename = paste0("../final_plots/506_lda_8sigs.png"), plot = plot, path = NULL, scale = 1.5, width = 350, height = 175, units = "mm", dpi = 300)
  
  diff = dataOrig119 - W %*% H
  frobenius = sqrt(sum(diag(t(diff) %*% diff)))
  print(frobenius)
}
if (isRFN == TRUE) {
  raw_data_rfn <- readMat("RFN/models/21_breast_C_L1_0_4_100_out.mat")
  raw_data_rfn <- readMat("RFN/models/119_breast_C_L1_0_5_100_out.mat")
  raw_data_rfn <- readMat("RFN/models/506_breast_C_L1_0_6_100_out.mat")
  
  W= raw_data_rfn$W
  H= raw_data_rfn$H
  H = H * colSums(W)
  W = t(t(W)/colSums(W))
  data = W

  rownames(data) = dataMmRownamesRightOrder
  data = data[mutationTypes,]
  plot = plotMutationalSignatures(data, 0.25)
  plot
  #ggsave(filename = paste0("../final_plots/119_rfn_5sigs.png"), plot = plot, path = NULL, scale = 1.5, width = 350, height = 175, units = "mm", dpi = 300)
  
  diff = dataOrig119 - Wrfn %*% H
  frobenius = sqrt(sum(diag(t(diff) %*% diff)))
  print(frobenius)
  
  cor(W)  
  cosine(W)
  print(max(cosine(W) - diag(sig)))
  print(max(cor(W) - diag(sig)))
}
if (isHLDA == TRUE) {
  str = "26"
  data <- read.table(paste0("LDA/hlda-c/test", str, ".distr"), header=FALSE,sep=" ")
  maxLevel = max(data[,1])
  for (level in 0:maxLevel) {
    levelData = subset(data, data[,1] == level)
    levelData = levelData[,c(-1, -2, -3)]
    levelData = t(levelData)
    rownames(levelData) = mutationTypes
    plot = plotMutationalSignatures(levelData, 0.25)
  
    ggsave(filename = paste0("LDA/hlda-c/plots/21_breast_hlda_2", str, "_levelplot_", level, ".jpg"), plot = plot,
           path = NULL, scale = 1, width = 300, height = 150, units = "mm", dpi = 200)
  }
}

#################### H PLOT
#Normalized proportions from RFN
h = t(H)
normalizedH = h/rowSums(h)
mydf2 = as.data.frame(normalizedH)
mydf2$Samples = row.names(mydf2)
mydf2.molten <- melt(mydf2, value.name="Identity", variable.name="Signature")
colourCount = dim(h)[2]
colorPalette <- c("forestgreen",  "yellowgreen", "greenyellow", "gray60","snow")
p = ggplot( data=mydf2.molten, aes(x = Samples, y = value, fill=variable))
p = p + scale_fill_manual(values = alpha(colorPalette, 0.8)) + theme(text = element_text(size = 15)) + xlab("Samples") + ylab("Proportion") 
p = p  +  geom_bar(stat = "identity") + ggtitle(paste0("RFN C impl ", dim(normalizedH)[2], " signatures proportions H"))
p
#ggsave(filename = paste0("../final_plots/119_rfn_5sigs_proportions.png"), plot = p, path = NULL, scale = 1.5, width = 350, height = 175, units = "mm", dpi = 300)