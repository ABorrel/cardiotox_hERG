#!/usr/bin/env Rscript
library (randomForest)
library(pls)
library(nnet)
library(ggplot2)
source("../../../../development/QSAR-QSPR/performance.R") # need to change that

calR2 = function(dreal, dpredict){
  
  dreal = as.double(as.character(dreal))
  dpredict = as.double(as.character(dpredict))
  
  #print(dreal)
  #print(dpredict)

  #print("Nb val in perf:")
  #print(length(dreal))
  
  #print("Nb val predict:")
  #print(dim(dperf))

  M = mean(dreal)
  SCEy = 0.0  
  SCEtot = 0.0
  for (i in seq(1, length(dreal))){
    #print (i)
    SCEy = SCEy + (dreal[i] - dpredict[i])**2
    SCEtot = SCEtot + (dreal[i] - M)**2
  }
  
  #print(SCEy)
  #print(SCEtot)

  r2 = 1 - (SCEy/SCEtot)
  #print(r2)

  print(cor(dreal, dpredict))

  return (as.double(r2))
}




################
#     MAIN     #
################

args <- commandArgs(TRUE)
p_desc = args[1]
p_aff = args[2]
p_AD = args[3]
p_model = args[4]
ML = args[5]
pr_out = args[6]



#p_desc = "c://Users/aborr/research/ILS/HERG/results/AID588834_filtered/desc_global.csv"
#p_aff = "c://Users/aborr/research/ILS/HERG/results/AID588834_filtered/aff_formated.csv"
#p_model = "c://Users/aborr/research/ILS/HERG/results/NCAST_CHEMBL/QSARReg_NCAST_CHEMBL__0.9-90-5-10-0.15-0/2/RFreg/model.RData"
#pr_out = "c://Users/aborr/research/ILS/HERG/results/AID588834_filtered/predict_model_NCAST_CHEMBL_active_reg/PCR/"
#p_AD = "c://Users/aborr/research/ILS/HERG/results/AID588834_filtered/AD/AD_Test_zscore.csv"
#ML = "RFreg"

# open file
d_aff = read.csv(p_aff, sep = ",", row.names=1)
d_aff = na.omit(d_aff)
d_aff = d_aff[-which(d_aff$Aff == 0),] # remove not active chemicals
din = read.csv(p_desc, sep = "\t", header = TRUE, row.names = 1)
l_chem = intersect(rownames(d_aff), rownames(din))
din = din[l_chem,]
din = din[-c(1),]
d_aff = d_aff[l_chem,]

# color by AD
d_ad = read.csv(p_AD, sep = ",", row.names = 1)
d_ad = d_ad[l_chem,]
AD = d_ad$Zscore
#d_aff = d_aff[which(d_ad$Zscore < 2),]
#din = din[which(d_ad$Zscore < 2),]



# prediction
load(p_model)
print(ls())
model = outmodel$model

print(attributes(model))

if(ML == "PLSreg" || ML == "PCRreg"){
  lpred = predict(model, din, ncomp = model$ncomp)
  print(lpred)
}else if(ML == "NNreg") {
  din = na.omit(din)
  lpred = predict(model, din)
  lpred = lpred[,1]
}else{
  lpred = predict(model, din)
}

dout = cbind(lpred, abs(d_aff$LogAC50))

colnames(dout) = c("Pred", "Real")
dout = as.data.frame(dout)
dout$Pred = as.double(as.character(dout$Pred))
dout$Real = as.double(as.character(dout$Real))
dout = cbind(dout, AD)
rownames(dout) = rownames(d_aff)
dout = na.omit(dout)

#dout = dout[which(dout$AD < 1), ]

write.csv(dout, paste(pr_out, "predicted.csv", sep=""))


# quality prediction
valr2 = calR2(dout$Real, dout$Pred)
corval = cor(dout$Real, dout$Pred)
RMSEP = vrmsep(dout$Real, dout$Pred)
MAEval = MAE(dout$Real, dout$Pred)
R02val = R02(dout$Real, dout$Pred)

l_perf = c(valr2, R02val, corval, RMSEP, MAEval)
names(l_perf) = c("R2", "R02", "r", "RMSEP", "MAE")
write.csv(l_perf, paste(pr_out, "perf_predicted.csv", sep=""))

# plot predictted VS real
theme_set(theme_grey())
p = ggplot(dout, aes(Real, Pred))+
  geom_point(size=1.5, colour="black", shape=21) + 
  geom_text(x=4.3, y=7.2, label = paste("R2=",round(valr2,2), sep = ""), size = 8)+
  labs(x = "Experimental pIC50", y = "Predicted pIC50") +
  theme(axis.text.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.text.x = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.x =  element_text(size = 25, hjust = 0.5, vjust =0.1))+
  xlim (c(4, 7.5)) +
  geom_segment(aes(x = 4, y = 4, xend = 7.5, yend = 7.5), linetype=2, size = 0.5) + 
  ylim (c(4, 7.5))
ggsave(paste(pr_out, "PerfRFregpoint.png", sep = ""), width = 7,height = 7, dpi = 300)


# with label - CV10
p = ggplot(dout, aes(Real, Pred, label=rownames(dout)))+
  geom_point(size=1.5, colour="black", shape=21) + 
  geom_text(x=4.3, y=7.2, label = paste("R2=",round(valr2,2), sep = ""), size = 8)+
  labs(x = "pAff", y = "Predicted pAff") +
  geom_text(size = 2.6, aes(label= rownames(dout)), parse = TRUE, color="black", nudge_y = 0.06) + 
  labs(x = "Experimental pIC50", y = "Predicted pIC50") + 
  theme(axis.text.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.text.x = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.x =  element_text(size = 25, hjust = 0.5, vjust =0.1))+
  xlim (c(4, 7.5))  +
  geom_segment(aes(x = 4, y = 4, xend = 7.5, yend = 7.5), linetype=2, size = 0.5) + 
  ylim (c(4, 7.5))  
#print(p)
ggsave(paste(pr_out, "PerfRFregname.png", sep = ""), width = 8,height = 8, dpi = 300)




theme_set(theme_grey())
p = ggplot(dout, aes(Real, Pred, color = AD))+
  geom_point(size=1.5, shape=19) +
  scale_color_gradient(low="blue", high="red") + 
  geom_text(x=4.3, y=7.2, label = paste("R2=",round(valr2,2), sep = ""), size = 8)+
  labs(x = "Experimental pIC50", y = "Predicted pIC50") +
  theme(axis.text.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.text.x = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.x =  element_text(size = 25, hjust = 0.5, vjust =0.1))+
  xlim (c(4, 7.5)) +
  geom_segment(aes(x = 4, y = 4, xend = 7.5, yend = 7.5), linetype=2, size = 0.5) + 
  ylim (c(4, 7.5))
ggsave(paste(pr_out, "PerfRFregpoint_AD.png", sep = ""), width = 8, height = 7, dpi = 300)


