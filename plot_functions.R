library(reshape)
library(ggplot2)

meltSignatures <- function(x, vars = c("motif", "signature")) {
  
  w_df = melt(x, varnames = vars)
  w_df
  w_df$alteration = sub("([ACGTN])([ACGTN]) .+", "\\1>\\2", w_df$motif)
  w_df$context = sub("[ACGTN][ACGTN] (.+)", "\\1", w_df$motif)
  levels = unique(w_df$signature)
  labels = signatureLabels(length(levels))
  if(all(levels %in% labels))
    levels = labels
  #print(levels)
  #la = c(6, 3, 16, 2, 66, 1, 8)
  w_df$signature = factor(w_df$signature, levels)#, labels=la)
  return(w_df)
}

signatureLabels <- function(n) {
  labels = sprintf("wwwS%d", 1:n)
  return(labels)
}


theme_ss <- function() {
  
  t = theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          axis.text.y = element_text(hjust = 0.5))
  
  return(t)
}

theme_small_axis <- function(x = TRUE, y = TRUE, size = 6, family = "mono") {
  ## decrease the x/y-axis label size
  template = element_text(size = size, family = family)
  t = theme_ss()
  if(x)
    t = t + theme(axis.text.x = template)
  if(y)
    t = t + theme(axis.text.y = template)
  
  return(t)
}


plotMutationalSignatures <- function(data, limit) {
  data = t(t(data) / rowSums(t(data)))
  w_df = meltSignatures(data)
  #wstd_df = meltSignatures(Wstd)
  #w_df$errorMax = w_df$value + wstd_df$value
  #w_df$errorMin = w_df$value - wstd_df$value
  #w_df$errorMin[w_df$errorMin < 0] = 0
  #w_df[, c(2,3,6)] <- sapply(w_df[, c(2,3,6)], as.numeric) #make needed cols numeric for subtraction
  #stabs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7)
  #stabilities = raw_data$stability
  #stabilities = stability
  sig_labbeler = function(value){
    return(paste0("Signature ", value, "\nStability ", round(stabilities[as.numeric(value)], 3)))
  }
  colorPalette <- c("dodgerblue", "gray20", "tomato", "orange1", "olivedrab2", "orchid1")
  p = ggplot(w_df, aes(x = context, y = value, fill = alteration))
  p = p + geom_bar(stat = "identity", position = "identity", width = 0.8)
  p = p + facet_grid(signature ~ alteration)
  p = p + theme_ss()# + theme_small_axis()
  p = p + theme(legend.position = "none") + theme(strip.background = element_rect(fill = 'gray88'), strip.text.x = element_text(size = 15))
  p = p + scale_fill_manual(values = alpha(colorPalette, 0.7))
  p = p + xlab("Motif") + ylab("Contribution") + theme(text = element_text(size=16), axis.text.y = element_text(size=15), axis.text.x = element_text(size=14))
  p = p + coord_cartesian(ylim = c(0, limit))  + scale_y_continuous(labels = scales::percent)#  +
  return(p)
}