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


computeDistFromCenter = function(coord_centroid, res.pca, name_plot=""){
  
  # plot distribution distance centroid for the all set
  l_dist = NULL
  i = 1
  imax = dim(res.pca)[1]
  
  while(i <= imax){
    x = rbind(coord_centroid, res.pca[i,])
    d = dist(x)
    l_dist = append(l_dist, d)
    i = i + 1 
  }
  
  l_n = rownames(res.pca)
  d_dist = cbind(l_n, l_dist)
  d_dist = cbind(l_n, d_dist)
  colnames(d_dist) = c("ID", "ID-2", "dist")
  rownames(d_dist) = l_n
  d_dist = as.data.frame(d_dist)
  d_dist$dist = as.double(as.character(d_dist$dist))
  
  if(name_plot != ""){
    ggplot(d_dist, aes(x=dist)) +
      geom_histogram(aes(y=..density..), position="identity", alpha=0.20, color="blue", fill="blue")+
      geom_density(alpha=0.3, fill="blue")+
      theme(text = element_text(size=19))+
      scale_color_manual(values=c("#cde2ff")) + 
      labs(title="",x="Distance to centroid", y = "Density")
    ggsave(paste(pr_out, "hist_", name_plot, "_dist.png", sep = ""),  width = 8, height = 7, dpi = 300, bg="transparent")
  }
  return(d_dist)
}



computeZscore = function(d_dist, M, SD, name_plot=""){
  
  #z-score = (x-μ)/σ
  #x is a raw score to be standardized; 
  #μ is the mean of the population; 
  #σ is the standard deviation of the population.
  
  l_dist = d_dist$dist
  l_n = rownames(d_dist)
  l_zscores = (l_dist - M)/SD
  l_zscores = abs(l_zscores)
  d_Zscore = cbind(l_n, l_zscores)
  d_Zscore = cbind(l_n, d_Zscore)
  colnames(d_Zscore) = c("ID", "ID-2", "Zscore")
  rownames(d_Zscore) = l_n
  d_Zscore = as.data.frame(d_Zscore)
  d_Zscore$Zscore = as.double(as.character(d_Zscore$Zscore))
  
  if(name_plot != ""){
    ggplot(d_Zscore, aes(x=Zscore)) +
      geom_histogram(aes(y=..density..), position="identity", alpha=0.20, color="blue", fill="blue")+
      geom_density(alpha=0.3, fill="blue") +
      theme(text = element_text(size=19)) +
      scale_color_manual(values=c("#cde2ff")) + 
      labs(title="",x="Z-scores to centroid", y = "Density")
    ggsave(paste(pr_out, "hist_", name_plot, "_Zscore.png", sep = ""),  width = 8, height = 7, dpi = 300, bg="transparent")
  }
  return(d_Zscore)
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
p_desc_model = args[1]
p_desc_test = args[2]
pr_out = args[3]


#p_desc_model = "C:\\Users\\Aborrel\\research\\ILS\\HERG\\results\\Cleaned_Data\\desc1D2D_cleaned.csv"
#p_desc_test = "C:\\Users\\Aborrel\\research\\ILS\\HERG\\results\\DrugSet\\DESC\\desc_1D2D.csv"
#pr_out = "C:\\Users\\Aborrel\\research\\ILS\\HERG\\results\\DrugSet\\AD\\"


ddesc1 = read.csv(p_desc_model, sep = ",", row.names = 1)
ddesc2 = read.csv(p_desc_test, sep = "\t", row.names = 1)

# remove SMILES
ddesc2 = ddesc2[,-1]

res.pca <- prcomp(ddesc1, scale = TRUE)
var_cap = generatePCAcoords(ddesc1)[[2]]


# average each dimension
coord_centroid = apply(res.pca$x, 2, "mean")

# for training set
d_dist = computeDistFromCenter(coord_centroid, res.pca$x, "trainning")
M = mean(d_dist$dist)
SD = sd(d_dist$dist)

# transform in z-score
d_zscore_train = computeZscore(d_dist, M, SD, "trainning")


#######################
# for the external set

add_coord = predict(res.pca, newdata = ddesc2)
d_dist_test = computeDistFromCenter(coord_centroid, add_coord, "test")
d_zscore_test = computeZscore(d_dist_test, M, SD, "test")

# define a in/out AD => threshold at 1
l_AD = d_zscore_test$Zscore
l_AD[which(d_zscore_test$Zscore<=1)] = 1
l_AD[which(d_zscore_test$Zscore >= 2)] = 0
l_AD[which(d_zscore_test$Zscore > 1 & d_zscore_test$Zscore < 2)] = 0.5

d_zscore_test = d_zscore_test[,-1]
d_zscore_test = cbind(d_zscore_test, l_AD)
d_zscore_test = d_zscore_test[,-1]
colnames(d_zscore_test)[dim(d_zscore_test)[2]] = "AD"

write.csv(d_zscore_test, file=paste(pr_out, "AD_zscore.csv", sep=""))


