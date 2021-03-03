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



#pdesc = "../../results/ChEMBL/desc_1D2D.csv"
#pmodel = "../../results/QSARs_models/RF/2.RData"
#ML = "RFclass"
#pr_out = "../../results/ChEMBL_predict/predict_model_RF/"



load(pmodel)
model = outmodel$model
din = read.csv(pdesc, sep = "\t", header = TRUE)
rownames(din) = din[,1]


if(ML == "RFclass"){
  lpred = predict(model, din)
  dout = cbind(rownames(din), lpred)
}else if(ML == "SVMclass"){
  lpred = predict(model, din)
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

