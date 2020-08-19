#!/usr/bin/env Rscript
library(Toolbox)


computeTableSignif = function(ddesc, dAC50, prresult){
  
  lact = rownames(dAC50[which(!is.na(dAC50[,2])),])
  linact = rownames(dAC50[which(is.na(dAC50[,2])),])
    
  dout = data.frame()
  for(j in seq(2, dim(ddesc)[2])){
    print(j)
    dact = ddesc[lact,j]
    dinact = ddesc[linact,j]
    typeTest = conditionTtest(dact, dinact)
    if (typeTest == 0){
      pval = comparisonTest(dact, dinact, "non-parametric")
    }else{
      pval = comparisonTest(dact, dinact, "parametric")
    }
    if(!is.na(pval)){
      dout[j,1] = length(lact)
      dout[j,2] = length(linact)
      dout[j,3] = round(mean(dact, na.rm = TRUE),2)
      dout[j,4] = round(mean(dinact, na.rm = TRUE),2)
      dout[j,5] = pval
      dout[j,6] = signifPvalue(pval)  
    }else{
      dout[j,1] = length(lact)
      dout[j,2] = length(linact)
      dout[j,3] = round(mean(dact, na.rm = TRUE),2)
      dout[j,4] = round(mean(dinact, na.rm = TRUE),2)
      dout[j,5] = "NA"
      dout[j,6] = "NA"
    }
  }      
  rownames(dout) = colnames(ddesc)
  colnames(dout) = c("Nbact", "Nbinact", "Mact", "Minact", "pval", "Signif")
  orderPval = order(dout[,5],decreasing=F)
  dout = dout[orderPval,]
  pfilout = paste(presult, colnames(dAC50)[i], ".csv", sep = "")
  write.csv(dout, pfilout)
}




mergeCol = function(ddesc, lmerge, delNA = 1){
  
  dtemp = ddesc[,lmerge]
  if(delNA == 1){
    dout = rowMeans(dtemp, na.rm = TRUE)
  }else{
    dout = rowMeans(dtemp, na.rm = FALSE)  
  }
  
  dout[which(dout == "NaN")] = NA
  
  return (dout)
}



################
#     MAIN     #
################

args <- commandArgs(TRUE)
pdesc = args[1]
pAC50 = args[2]
pr_out = args[3]

#pdesc = "./../../results/DESC/desc_1D2D.csv"
#pAC50 = "./../../data/AC50_7403.txt"
#presult = "./../../results/SignifDesc/Rdkit_"

# AC50
dAC50 = read.csv(pAC50, sep=",", header = TRUE)
rownames(dAC50) = dAC50[,1]

# descriptor
ddesc = read.csv(pdesc, sep = "\t", header = TRUE)
if(dim(ddesc == 0)){
  ddesc = read.csv(pdesc, sep = ",", header = TRUE)
}

rownames(ddesc) = ddesc[,1]
ddesc = ddesc[,-1]

# for any conditions
computeTableSignif(ddesc, dAC50, prresult)
