#!/usr/bin/env Rscript
library (randomForest)
library (MASS)
library(rpart)
library(rpart.plot)
library(e1071)
#library(neuralnet)
library(nnet)
library(clusterGeneration)
library(stringr)
library(reshape2)



################
#     MAIN     #
################
args <- commandArgs(TRUE)
pdesc = args[1]
pmodel = args[2]
ML = args[3]
pr_out = args[4]


#pdesc = "c://Users/aborr/research/ILS/HERG/results/AID588834_filtered/desc_global.csv"
#pmodel = "c://Users/aborr/research/ILS/HERG/results/NCAST_CHEMBL/QSAR_NCAST_CHEMBL__0.9-90-5-10-0.15-0/noSampling/2/SVMclass_linear/model.RData"
#ML = "SVMclass"
#pr_out = "c://Users/aborr/research/ILS/HERG/results/AID588834_filtered/predict_model_NCAST_CHEMBL_classif_nosampling/SVM-linear/"



load(pmodel)
model = outmodel$model
din = read.csv(pdesc, sep = "\t", header = TRUE)
rownames(din) = din[,1]
din = din[,-1]


if(ML == "RFclass"){
  lpred = predict(model, din)
  dout = cbind(rownames(din), lpred)
}else if(ML == "SVMclass"){
  scale_val = model$x.scale
  scale_center = scale_val$'scaled:center'
  scale_scale = scale_val$'scaled:scale'
  din = din[,names(scale_center)]
  #din = scale(din, scale=scale_scale, center=scale_center)
  lpred = predict(model,din, decision.values = TRUE)
  lpred = as.double(lpred)
  lpred[which(lpred == 1 )] = 0
  lpred[which(lpred == 2 )] = 1
  dout = cbind(rownames(din), as.vector(lpred))
}else if(ML == "LDAclass"){
  lpred = predict(model, din)
  dout = cbind(rownames(din), as.vector(lpred$class))
# to update  
}else if (ML == "CARTclass"){
  lpred = predict(model, din)
  dout = cbind(lpred[,2], din$Aff)
}


# name exist
xmodel = strsplit(pmodel, "/")
xmodel = xmodel[[1]]
xmodel = xmodel[length(xmodel)]
xmodel = strsplit(xmodel, "[.]")
xmodel = xmodel[[1]]
xmodel = xmodel[1]

colnames(dout) = c("ID", "Pred")
rownames(dout) = rownames(din)
write.csv(dout, paste(pr_out, "pred_", xmodel, sep = ""))

