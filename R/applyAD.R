#!/usr/bin/env Rscript
library(ggplot2)
library(plyr)
library(stringr)


accuracy = function (tp, tn, fp, fn){
  return ((tp + tn)/(tp + fp + tn +fn))
}

precision = function (tp, fp){
  return (tp/(tp + fp))
}

recall = function (tp, fn){
  return (tp/(tp + fn))
}

specificity = function (tn, fp){
  sp = tn/(tn + fp)
  if(is.na(sp)){
    sp = 0
  }
  return (sp)
}

sensibility = function (tp, fn){
  se = tp/(tp + fn)
  if(is.na(se)){
    se = 0
  }
  return (se)
}

MCC = function (tp, tn, fp, fn){
  numerator = tp*tn-fp*fn
  A1 = as.double(tp+fp)
  A2 = tp+fn
  A3 = tn+fp
  A4 = tn+fn
  denumerator = A1 * A2 *A3 *A4
  mcc = numerator / sqrt(denumerator)
  if(is.na(mcc)){
    mcc = 0
  }
  return (mcc)
}


includeInADD = function(d, d_AD, cutoff_zscore){

  l_inAD = rownames(d_AD)[which(d_AD$Zscore < cutoff_zscore)]

  l_ID_inAD = intersect(rownames(d), l_inAD)

  # define light d
  d = d[l_ID_inAD,]

  d[which(d[,"Real"] == 0),"Real"] = NA
  dact = na.omit(d)
  dinact = d[which(is.na(d[,3])),]



  # SUMMARY PREDICTION
  h = c("NB active", "NB inact", "TP", "FN", "TN", "FP", "Mact prob", "SDact prob", "Minact prob", "SDinact prob", "Acc", "Sp", "Se", "MCC")

  nbact = dim(dact)[1]
  nbinact = dim(dinact)[1]
  TP = length(which(dact[,1] >= 0.5))
  FN = length(which(dact[,1] < 0.5))
  Mact = mean(dact[,1])
  SDact = sd(dact[,1])

  TN = length(which(dinact[,1] < 0.5))
  FP = length(which(dinact[,1] >= 0.5))
  Minact = mean(dinact[,1])
  SDinact = sd(dinact[,1])

  Acc = accuracy(TP, TN, FP, FN)
  Sp = specificity(TN, FP)
  Se = sensibility(TP, FN)
  mcc = MCC(TP, TN, FP, FN)


  cval = c(nbact, nbinact, TP, FN, TN, FP, Mact,SDact, Minact, SDinact, Acc, Sp, Se, mcc)
  names(cval) = h 


  return(cval)

}

################
#     MAIN     #
################

args <- commandArgs(TRUE)
p_pred = args[1]
p_AD = args[2]
pr_out = args[3]


#p_pred = "../../results/AID_588834/predict_model_RF/Merge_pred.csv"
#p_AD = "../../results/AID_588834/AD/AD_zscore.csv"
#pr_out = "../../results/AID_588834/predict_AD//"

# open prediction
d = read.csv(p_pred, sep = "\t")
if (dim(d)[2] == 1){
  d = read.csv(p_pred, sep = ",")  
}
rownames(d) = d[,1]
d = d[,-1]


# open AD
d_AD = read.csv(p_AD, sep = ",")
rownames(d_AD) = d_AD[,1]

cval = includeInADD(d, d_AD, 0.75)
write.csv(cval, file = paste(pr_out, "sum_Pred_075_AD.csv", sep = ""))

cval = includeInADD(d, d_AD, 1)
write.csv(cval, file = paste(pr_out, "sum_Pred_1_AD.csv", sep = ""))

cval = includeInADD(d, d_AD, 2)
write.csv(cval, file = paste(pr_out, "sum_Pred_2_AD.csv", sep = ""))

cval = includeInADD(d, d_AD, 4)
write.csv(cval, file = paste(pr_out, "sum_Pred_4_AD.csv", sep = ""))

cval = includeInADD(d, d_AD, 10)
write.csv(cval, file = paste(pr_out, "sum_Pred_10_AD.csv", sep = ""))