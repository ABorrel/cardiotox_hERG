#!/usr/bin/env Rscript
library(Toolbox)


computeTableSignif = function(ddesc, dAC50, pr_out){
  
  lact = rownames(dAC50[which(!is.na(dAC50[,2])),])
  linact = rownames(dAC50[which(is.na(dAC50[,2])),])
    
  dout = data.frame()
  for(j in seq(1, dim(ddesc)[2])){
    print(j)
    print(colnames(ddesc)[j])
    dact = ddesc[lact,j]
    dact = na.omit(dact)
    dinact = ddesc[linact,j]
    dinact = na.omit(dinact)
    typeTest = conditionTtest(dact, dinact)
    if (typeTest == 0){
      pval = comparisonTest(dact, dinact, "non-parametric")
    }else if (typeTest == 1){
      pval = comparisonTest(dact, dinact, "parametric")
    }else{
      pval = NA
    }
    pval = as.double(as.character(pval))
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
      dout[j,5] = NA
      dout[j,6] = NA
    }
  }      
  rownames(dout) = colnames(ddesc)
  colnames(dout) = c("Nbact", "Nbinact", "Mact", "Minact", "pval", "Signif")
  orderPval = order(dout[,5],decreasing=F)
  dout = dout[orderPval,]
  pfilout = paste(pr_out, "desc_signif.csv", sep = "")
  write.csv(dout, pfilout)
}





################
#     MAIN     #
################

args <- commandArgs(TRUE)
pdesc = args[1]
pAC50 = args[2]
pr_out = args[3]

#pdesc = "../../results/DESC/desc_OPERA.csv"
#pAC50 = "./../../data/AC50_7403.txt"
#pr_out = "./../../results/SignifDesc/OPERA_"

# AC50
dAC50 = read.csv(pAC50, sep=",", header = TRUE)
if(dim(dAC50)[2] == 1){
  dAC50 = read.csv(pAC50, sep = "\t", header = TRUE)
}
rownames(dAC50) = dAC50[,1]

# descriptor
ddesc = read.csv(pdesc, sep = "\t", header = TRUE)
if(dim(ddesc)[2] == 1){
  ddesc = read.csv(pdesc, sep = ",", header = TRUE)
}

rownames(ddesc) = ddesc[,1]
ddesc = ddesc[,-1]

# reduce desc for OPERA
if (!is.integer0(grep("OPERA", pr_out))){
  l_opera_pred = NULL
  for(opera_desc in colnames(ddesc)){
    if(!is.integer0(grep("_pred", opera_desc, fixed = TRUE))){
      if(is.integer0(grep("_predRange", opera_desc, fixed = TRUE)) && is.integer0(grep("pKa_b_pred", opera_desc, fixed = TRUE)) && is.integer0(grep("pKa_a_pred", opera_desc, fixed = TRUE)) && is.integer0(grep("CATMoS", opera_desc, fixed = TRUE)) && is.integer0(grep("CERAPP", opera_desc, fixed = TRUE))&& is.integer0(grep("CoMPARA", opera_desc, fixed = TRUE))){
        l_opera_pred = append(l_opera_pred, opera_desc)
      }
    }
  }
  ddesc = ddesc[,l_opera_pred]
}

# for any conditions
computeTableSignif(ddesc, dAC50, pr_out)
