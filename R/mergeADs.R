#!/usr/bin/env Rscript
library(ggplot2)
library(factoextra)
library(Toolbox)# personal library with toolbox functions




histZscore = function(d_Zscore, p_out){
  
  ggplot(d_Zscore, aes(x=Zscore)) +
    geom_histogram(aes(y=..density..), position="identity", alpha=0.20, color="blue", fill="blue")+
    geom_density(alpha=0.3, fill="blue") +
    theme(text = element_text(size=19)) +
    scale_color_manual(values=c("#cde2ff")) + 
    xlim(0, 7)+
    labs(title="",x="Z-scores", y = "Density")
  ggsave(p_out,  width = 6, height = 7, dpi = 300, bg="transparent")
  
}


compteSummary = function(d_Zscore, p_sum){
  
  l_sum = c(min(d_train_AD$Zscore), max(d_train_AD$Zscore), mean(d_train_AD$Zscore), sd(d_train_AD$Zscore), median(d_train_AD$Zscore))
  names(l_sum) = c("Min Zscores", 'Max Zscores', "M Zscores", "SD Zscores", "Median Zscores")
  write.csv(l_sum, p_sum)
  
}


################
#     MAIN     #
################

args <- commandArgs(TRUE)
p_train_AD = args[1]
p_test_AD = args[2]
p_desc = args[3]
pr_out = args[4]


p_train_AD = "./../../results/QSAR/Merge_AD/Zscores_train.csv"
p_test_AD = "./../../results/QSAR/Merge_AD/Zscores_test.csv"
p_desc = "./../../results/Cleaned_Data/desc1D2D_cleaned.csv"

pr_out = "./../../results/QSAR/Merge_AD/"


# desc 
d_desc = read.csv(p_desc, sep = ",", row.names = 1)


# train
d_train_AD = read.csv(p_train_AD, sep = "\t")
rownames(d_train_AD) = d_train_AD[,1]
d_desc_train = d_desc[rownames(d_train_AD),]
histZscore(d_train_AD, paste(pr_out, "train_Zscore.png", sep = ""))
compteSummary(d_train_AD, paste(pr_out, "sum_train_Zscore.csv", sep = ""))


res.pca <- prcomp(d_desc_train, scale = TRUE)
var_cap = generatePCAcoords(d_desc_train)[[2]]


# test
d_test_AD = read.csv(p_test_AD, sep = "\t")
rownames(d_test_AD) = d_test_AD[,1]
d_desc_test = d_desc[rownames(d_test_AD),]
histZscore(d_test_AD, paste(pr_out, "test_Zscore.png", sep = ""))
compteSummary(d_test_AD, paste(pr_out, "sum_test_Zscore.csv", sep = ""))

add_coord = predict(res.pca, newdata = d_desc_test)


# plot PCA
vcolor = rep(addTrans("#3336FF", 80), dim(d_desc_train)[1])
vcolor_add = rep(addTrans("#36D529", 85), dim(d_desc_test)[1])


png(paste(pr_out, "PCA_color.png", sep = ""), 1700, 1500)
par(mar=c(8,8,8,8))
plot(rbind(add_coord[,1],res.pca$x[,1]) ,rbind(add_coord[,2],res.pca$x[,2]) , pch=19, col = "white", xlab = paste("CP1: ", round (var_cap[1], 2), "%", sep = ""), ylab = paste("CP2: ", round(var_cap[2], 2), "%", sep = ""), cex.lab = 4, cex.main = 4, cex.axis = 1.75, cex = 2.5)
points(res.pca$x[,1], res.pca$x[,2], pch=19, col = vcolor, cex = 2.5)
#points(res.pca$x[which(vcolor != "darkgray"),1], res.pca$x[which(vcolor != "darkgray"),2], pch=19, col = vcolor[which(vcolor != "darkgray")], cex = 2.5)
points(add_coord[,1], add_coord[,2], col = vcolor_add, pch = 19, cex = 2.5)
abline(h=0,v=0)
warnings ()
dev.off()


