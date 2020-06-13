#!/usr/bin/env Rscript
library(ape)
library(phangorn)
library(ggtree)
library(factoextra)



dendogramCircle = function(ddes, daff, prout){
  
  #calibrate affinity for color
  daff = as.data.frame(daff)
  
  
  matTrans1 <- scale(ddes)
  d <- dist(matTrans1, method = "euc")
  tupgma2 <- upgma(d, method="ward.D2")
  legend.title = "AC50"
  
  
  # all chemicals
  pfilout = paste(prout, "HClust_dendo.png", sep = "")
  t4 <- ggtree(tupgma2, layout="circular", size=1)
  t4 <- t4 %<+% daff + geom_text(aes(color=Aff, label=label, angle=angle), hjust=-0.5, size=1) +
    geom_tippoint(aes(color=Aff), alpha=0.75, size=0.5)+
    scale_color_continuous(low='red', high='lightgreen') +
    theme(legend.position="right")+
    theme(plot.margin = unit(c(0,0,0,0), "cm")) +
    geom_treescale(x = 1, y = 1, width = NULL, offset = NULL,
                   color = "white", linesize = 1E-100, fontsize = 1E-100)
  ggsave(pfilout, dpi=300, height = 10, width = 14)
  
  
  # only active chemicals
  daff_act = na.omit(daff)
  matTrans1_act = matTrans1[rownames(daff_act),]
  d <- dist(matTrans1_act, method = "euc")
  tupgma2 <- upgma(d, method="ward.D2")
  legend.title = "AC50"
  
  pfilout = paste(prout, "HClust_dendo_active.png", sep = "")
  t4 <- ggtree(tupgma2, layout="circular", size=1)
  t4 <- t4 %<+% daff + geom_text(aes(color=Aff, label=label, angle=angle), hjust=-0.5, size=1) +
    geom_tippoint(aes(color=Aff), alpha=0.75, size=0.5)+
    scale_color_continuous(low='red', high='lightgreen') +
    theme(legend.position="right")+
    theme(plot.margin = unit(c(0,0,0,0), "cm")) +
    geom_treescale(x = 1, y = 1, width = NULL, offset = NULL,
                   color = "white", linesize = 1E-100, fontsize = 1E-100)
  
  ggsave(pfilout, dpi=300, height = 10, width = 14)
  
}



################
#     MAIN     #
################

args <- commandArgs(TRUE)
p_desc = args[1]
p_AC50 = args[2]
pr_out = args[3]

#p_desc = "../../results/Cleaned_Data/desc1D2D_cleaned.csv"
#p_AC50 = "../../results/Cleaned_Data/AC50_cleaned.csv"
#pr_out = "../../results/HClust/"


# open files
ddesc = read.csv(p_desc, sep = ",", header = TRUE)
rownames(ddesc) = ddesc[,1]
ddesc = ddesc[,-1]

dAC50 = read.csv(p_AC50, sep = ",", header = TRUE)
rownames(dAC50) = dAC50[,1]
dAC50 = dAC50[,-1]


dendogramCircle(ddesc, dAC50, pr_out)
  

