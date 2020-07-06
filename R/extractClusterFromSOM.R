#!/usr/bin/env Rscript
library(kohonen)
library(ggplot2)


#################
# Color reverse #
#################

applySOM = function(som_model, d_AC50, pr_out){
  
  #write cluster
  dclust = som_model$unit.classif
  names(dclust) = rownames(som_model$data[[1]])
  
  xdim = som_model$grid$xdim
  ydim = som_model$grid$ydim
  
  # remove inactive
  d_AC50 = na.omit(d_AC50)
  dclust = dclust[rownames(d_AC50)]# take only active chemical
  dclust = cbind(names(dclust), dclust)
  colnames(dclust) = c("CASRN", "Cluster")
  
  lAct = rep(0, xdim*ydim)
  names(lAct) = seq(1, xdim*ydim)
  
  ltable = table(dclust[,2])
  lAct[names(ltable)] = ltable
  
  ltabinit = table(som_model$unit.classif)
  linitial = rep(0, xdim*ydim)
  names(linitial) = seq(1, xdim*ydim)
  linitial[names(ltabinit)] = ltabinit
  
  lprob = lAct / linitial
  
  
  write.csv(dclust, paste(pr_out, "SOM_Clusters_act", sep = ""))
  
  # count of active
  svg(paste(pr_out, "SOM_count_act.svg", sep = ""))
  plot(som_model, type = "property", property=lAct, palette.name=coolBlueHotRed, main = "", dpi=300, height = 20, width = 20, bg = "transparent")
  dev.off()
  
  # count of active calibrate
  lAct[which(lAct == max(lAct))] = max(table(som_model$unit.classif))# have to calibrate based on the max of the original SOM
  svg(paste(pr_out, "SOM_count_act_calibrate.svg", sep = ""))
  plot(som_model, type = "property", property=lAct, palette.name=coolBlueHotRed, main = "", dpi=300, height = 20, width = 20, bg = "transparent")
  dev.off()
  
  # plot with proba
  svg(paste(pr_out, "SOM_prob_act.svg", sep = ""))
  plot(som_model, type = "property", property=lprob, palette.name=coolBlueHotRed, main = "Prob active", dpi=300, height = 20, width = 20, bg = "transparent")
  dev.off()
  
  write.csv(lprob, paste(pr_out, "SOM_Clusters_prob", sep = ""))
  
}



colors <- function(n, alpha = 1) {
  rev(heat.colors(n, alpha))
}

coolBlueHotRed <- function(n, alpha = 1) {rainbow(n, end=4/6, alpha=alpha)[n:1]}


################
#     MAIN     #
################

args <- commandArgs(TRUE)
p_desc = args[1]
p_SOMmodel = args[2]
pr_out = args[3]

#p_desc = "../../results/SOM/class_from_interferences/Drug/desc.csv"
#p_SOMmodel = "./../../results/SOM/SOM_model.RData"
#pr_out = "../../results/SOM/class_from_interferences/Drug/"


SOMmodel = load(p_SOMmodel)
d_desc = read.csv(p_desc, sep = "\t", row.names = 1)
d_desc = d_desc[,-1]
d_desc = as.matrix(d_desc)


applySOM(som_model, d_desc, pr_out)


