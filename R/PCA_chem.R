#!/usr/bin/env Rscript
library(ggplot2)
library(Toolbox)



PCAplot = function (din, vcolor, prout){
  coord = generatePCAcoords(din)
  data_plot = coord[[1]]
  var_cap = coord[[2]]
  cp = coord[[3]]
  
  col.desc = "black"
  
  
  png (paste (prout, "PCA_text.png", sep = ""), 1700, 1500)
  factor = factorACP (data_plot, cp)
  
  color_arrow = col.desc[rownames(cp)]
  par(mar=c(8,8,8,8))
  plot(data_plot[,1],data_plot[,2], pch=20, main = paste (var_cap[1],var_cap[2], sep = "_" ), xlab = paste("PC1: ", signif (var_cap[1], 4), "%", sep = ""), ylab = paste("PC2: ", signif (var_cap[2], 4), "%", sep = ""), cex.lab = 4, cex.main = 4, cex.axis = 1.75, cex = 4, type = "n")
  text (data_plot[,1],data_plot[,2], label = rownames (din), cex = 1.2)
  abline(h=0,v=0)
  warnings ()
  dev.off()
  

  # PCA color  
  png (paste (prout, "PCA_color.png", sep = ""), 1700, 1500)
  factor = factorACP (data_plot, cp)
  color_arrow =col.desc
  par(mar=c(8,8,8,8))
  plot(data_plot[,1],data_plot[,2], pch=20, col = vcolor, xlab = paste("PC1: ", round (var_cap[1], 2), "%", sep = ""), ylab = paste("PC2: ", round(var_cap[2], 2), "%", sep = ""), cex.lab = 4, cex.main = 4, cex.axis = 1.75, cex = 4)
  points(data_plot[which(vcolor != "gray90"),1], data_plot[which(vcolor != "gray90"),2], pch=20, col = vcolor[which(vcolor != "gray90")], cex = 4)
  abline(h=0,v=0)
  warnings ()
  dev.off()
  
  
  png (paste (prout, "PCA_descriptor.png", sep = ""), 1700, 1500)
  par(mar=c(8,8,8,8))
  plot(data_plot[,1],data_plot[,2], xlab = paste("CP1: ", signif (var_cap[1], 4), "%", sep = ""), ylab = paste("CP2: ", signif (var_cap[2], 4), "%", sep = ""), pch=20, cex.lab = 4, cex.main = 4, cex.axis = 1.75, cex = 4, type = "n")
  #points (data_plot[,1][length (color_point):dim(data_plot)[1]],data_plot[,2][length (color_point):dim(data_plot)[1]], pch=17, cex = 4, col = color_point2)
  abline(h=0,v=0)
  arrows (0,0,cp[,1]*factor,cp[,2]*factor, col = color_arrow, lwd = 4 )
  text (cp[,1]*factor, cp[,2]*factor,rownames(cp), col = color_arrow, cex = 2.5)
  dev.off()
  
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
#pr_out = "../../results/PCA/"


ddesc = read.csv(p_desc, sep = ",", header = TRUE)
rownames(ddesc) = ddesc[,1]
ddesc = ddesc[,-1]

dAC50 = read.csv(p_AC50, sep = ",", header = TRUE)

rownames(dAC50) = dAC50[,1]
dAC50 = dAC50[,-1]

lID = intersect(rownames(dAC50), rownames(ddesc))
print(lID)

ddesc = ddesc[lID,]
dAC50 = dAC50[lID,]

vcolor = rep("gray90", dim(ddesc)[1])
vcolor[which(!is.na(dAC50$Aff))] = "blue"

PCAplot(ddesc, vcolor, pr_out)