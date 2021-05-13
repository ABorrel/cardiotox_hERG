#!/usr/bin/env Rscript
library(ggplot2)
library(factoextra)



generatePCAcoords = function(din){
  
  dinScale = scale(din)
  data.cor=cor(dinScale)
  data.eigen=eigen(data.cor)
  lambda = data.eigen$values
  var_cap = lambda/sum(lambda)*100
  cp = data.eigen$vectors
  rownames (cp) = colnames (dinScale)
  colnames (cp) = colnames (dinScale)
  data_plot = as.matrix(dinScale)%*%cp
  
  return(list(data_plot, var_cap, cp))
}



addTrans <- function(color,trans)
{
  # This function adds transparancy to a color.
  # Define transparancy with an integer between 0 and 255
  # 0 being fully transparant and 255 being fully visable
  # Works with either color and trans a vector of equal length,
  # or one of the two of length 1.
  
  if (length(color)!=length(trans)&!any(c(length(color),length(trans))==1)) stop("Vector lengths not correct")
  if (length(color)==1 & length(trans)>1) color <- rep(color,length(trans))
  if (length(trans)==1 & length(color)>1) trans <- rep(trans,length(color))
  
  num2hex <- function(x)
  {
    hex <- unlist(strsplit("0123456789ABCDEF",split=""))
    return(paste(hex[(x-x%%16)/16+1],hex[x%%16+1],sep=""))
  }
  rgb <- rbind(col2rgb(color),trans)
  res <- paste("#",apply(apply(rgb,2,num2hex),2,paste,collapse=""),sep="")
  return(res)
}


################
#     MAIN     #
################

args <- commandArgs(TRUE)
pdesc1 = args[1]
pAff1 = args[2]
pdesc2 = args[3]
pAff2 = args[4]
prout = args[5]

#pdesc1 = "../../results/Cleaned_Data/desc1D2D_cleaned.csv"
#pAff1 = "../../results/Cleaned_Data/AC50_cleaned.csv"
#pdesc2 = "../../results/AID_588834/desc_cleaned.csv"
#pAff2 = "../../results/AID_588834/aff_cleaned.csv"
#prout = "../../results/AID_588834/PCA_vs/"


ddesc1 = read.csv(pdesc1, sep = ",", row.names = 1)
ddesc2 = read.csv(pdesc2, sep = "\t", row.names = 1)



# remove SMILES
#ddesc2 = ddesc2[,-1]

daff1 = read.csv(pAff1, header = TRUE, row.names = 1)
daff2 = read.csv(pAff2, header = TRUE, row.names = 1)


l_chem1 = intersect(rownames(ddesc1), rownames(daff1))
l_chem2 = intersect(rownames(ddesc2), rownames(daff2))


# organize same order
daff1 = daff1[l_chem1,]
daff2 = daff2[l_chem2,]

ddesc1 = ddesc1[l_chem1,]
ddesc2 = ddesc2[l_chem2,]

res.pca <- prcomp(ddesc1, scale = TRUE)
var_cap = generatePCAcoords(ddesc1)[[2]]

add_coord = predict(res.pca, newdata = ddesc2)

vcolor = rep(addTrans("#595959", 60), dim(ddesc1)[1])
vcolor[which(!is.na(daff1$Aff))] = "#191997"

vcolor_add = rep("#90EE90", dim(ddesc2)[1])
vcolor_add[union(which(daff2$Aff == 0), which(is.na(daff2$Aff)))] = "#2c7825"
vshape = rep(18, dim(ddesc2)[1])
#vshape[union(which(daff2$Aff == 0), which(is.na(daff2$Aff)))] = 8

png(paste(prout, "PCA_color.png", sep = ""), 1700, 1500)
par(mar=c(8,8,8,8))
plot(rbind(add_coord[,1],res.pca$x[,1]) ,rbind(add_coord[,2],res.pca$x[,2]) , pch=19, col = "white", xlab = paste("PC1: ", round (var_cap[1], 2), "%", sep = ""), ylab = paste("PC2: ", round(var_cap[2], 2), "%", sep = ""), cex.lab = 4, cex.main = 4, cex.axis = 1.75, cex = 2.5)
points(res.pca$x[,1], res.pca$x[,2], pch=19, col = vcolor, cex = 2.5)
#points(res.pca$x[which(vcolor != "darkgray"),1], res.pca$x[which(vcolor != "darkgray"),2], pch=19, col = vcolor[which(vcolor != "darkgray")], cex = 2.5)
points(res.pca$x[which(vcolor != "#5959593C"),1], res.pca$x[which(vcolor != "#5959593C"),2], pch=19, col = vcolor[which(vcolor != "#5959593C")], cex = 2.5)
points(add_coord[,1], add_coord[,2], col = vcolor_add, pch = vshape, cex = 3)
abline(h=0,v=0)
warnings ()
legend(min(res.pca$x[,1]), max(res.pca$x[,2]), legend=c("Train set - active", "Train set - inactive", "Test set - active", "Test set - inactive"),
       col=c("#191997", "#595959", "#2c7825", "#90EE90"),pch = 19, cex=2)

dev.off()

#points(data_plot[which(vcolor != "gray90"),1], data_plot[which(vcolor != "gray90"),2], pch=20, col = vcolor[which(vcolor != "gray90")], cex = 4)
#plot(res.pca$x[,1], res.pca$x[,2], xlim = c(-40, 20), ylim = c(-20, 20))

#points(add_coord[,1], add_coord[,2], col = "blue")
#points(res.pca$x[,1], res.pca$x[,2], col = "red")