#!/usr/bin/env Rscript
library(Toolbox)


################
#     MAIN     #
################

args <- commandArgs(TRUE)
p_desc_rdkit = args[1]
p_desc_OPERA = args[2]
pr_out = args[3]

#p_desc_rdkit = "../../results/DESC/desc_1D2D.csv"
#p_desc_OPERA = "../../results/DESC/desc_OPERA.csv"
#pr_out = "../../results/QSAR/"


##################
# MERGE dataset  #
##################
d_rdkit = read.csv(p_desc_rdkit, sep = "\t")
rownames(d_rdkit) = d_rdkit[,1]
d_opera = read.csv(p_desc_OPERA, sep = ",", row.names = 1)
if(dim(d_opera)[2] == 0){
  d_opera = read.csv(p_desc_OPERA, sep = "\t", row.names = 1)
}


l_chem = intersect(rownames(d_rdkit), rownames(d_opera))
d_rdkit = d_rdkit[l_chem,]
d_opera = d_opera[l_chem,]

# remove SMILES
if("SMILES" %in% colnames(d_rdkit) == TRUE){
  d_rdkit = d_rdkit[,-which(colnames(d_rdkit)=="SMILES")]
}

if("SMILES" %in% colnames(d_opera) == TRUE){
  d_opera = d_opera[,-which(colnames(d_opera)=="SMILES")]
}

l_opera_pred = NULL
for(opera_desc in colnames(d_opera)){
  if(!is.integer0(grep("_pred", opera_desc, fixed = TRUE))){
    if(is.integer0(grep("_predRange", opera_desc, fixed = TRUE)) && is.integer0(grep("pKa_b_pred", opera_desc, fixed = TRUE)) && is.integer0(grep("pKa_a_pred", opera_desc, fixed = TRUE)) && is.integer0(grep("CATMoS", opera_desc, fixed = TRUE)) && is.integer0(grep("CERAPP", opera_desc, fixed = TRUE))&& is.integer0(grep("CoMPARA", opera_desc, fixed = TRUE))){
      l_opera_pred = append(l_opera_pred, opera_desc)
    }
  }
}

d_opera = d_opera[,l_opera_pred]
d_global = cbind(d_rdkit, d_opera)

p_global = paste(pr_out, "desc_global.csv", sep = "")
write.table(d_global, p_global, sep = "\t", row.names = FALSE)

