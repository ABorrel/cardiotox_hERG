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



################
#     MAIN     #
################

args <- commandArgs(TRUE)
pprob = args[1]


#pprob = "../../results/AID_588834/predict_model_RF/Merge_pred.csv"




d = read.csv(pprob, sep = "\t")
if (dim(d)[2] == 1){
  d = read.csv(pprob, sep = ",")  
}
rownames(d) = d[,1]
d = d[,-1]

d[which(d[,"Real"] == 0),"Real"] = NA
dact = na.omit(d)
dinact = d[which(is.na(d[,3])),]


# define name short
pprob = str_sub(pprob, 1,-5)

# Use geom_pointrange
ggplot(dact, aes(x=Real, y=Mpred)) + 
  geom_pointrange(aes(ymin=Mpred-SDpred, ymax=Mpred+SDpred))

ggsave(paste(pprob, "_act.png", sep = ""),  width = 8, height = 8, dpi = 300, bg="transparent")



lact = rep("act", dim(d)[1])
lact[which(is.na(d[,3]))] = "inact"

dhist = d
dhist[,3] = lact

mu <- ddply(dhist, "Real", summarise, grp.mean=mean(Mpred))



d = na.omit(d)
print(d)

ggplot(d, aes(x=Mpred, color=Real, fill=Real)) +
  geom_histogram(aes(y=..density..), position="identity", alpha=0.15)+
  geom_density(alpha=0.3)+
  #geom_vline(data=mu, aes(xintercept=grp.mean, color=Real),
  #           linetype=c("dashed", "solid"))+
  #xlim(0,120)+
  theme(text = element_text(size=19))+
  #scale_fill_manual(values=c("#eb9999", "#290000"), labels = c("active", "inactive"))+
  #scale_color_manual(values=c("#eb9999", "#290000", "#d21919", "#a40000"), labels = c("hek293 cell based", "hek2$
  labs(title="",x="Prob", y = "Density")

ggsave(paste(pprob, "_hist.png", sep = ""),  width = 8, height = 7, dpi = 300, bg="transparent")




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



write.csv(cval, file = paste(pprob, "_sum.csv", sep = ""))
