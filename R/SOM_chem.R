#!/usr/bin/env Rscript
require(kohonen)
library(ggplot2)



#################
# Color reverse #
#################

colors <- function(n, alpha = 1) {
  rev(heat.colors(n, alpha))
}

coolBlueHotRed <- function(n, alpha = 1) {rainbow(n, end=4/6, alpha=alpha)[n:1]}

#######
# SOM #
#######

generateSOM = function(ddesc, xdim, ydim, pr_out){
  
  ddesc = as.matrix(ddesc)
  som_grid <- somgrid(xdim=xdim, ydim=ydim, topo="hexagonal")
  som_model <- som(ddesc, 
                   grid=som_grid, 
                   rlen=100, 
                   alpha=c(0.05,0.01), 
                   keep.data = TRUE) 
  
  save(som_model, file = paste(pr_out, "SOM_model.RData", sep = ""))
  
  
  # draw SOM
  svg(paste(pr_out, "SOM_model_count.svg", sep = ""))
  plot(som_model, type = "count", palette.name=coolBlueHotRed, main = "", dpi=300, height = 20, width = 20, bg = "transparent")
  dev.off()
  
  # cluster
  dclust = som_model$unit.classif
  names(dclust) = rownames(som_model$data[[1]])
  dclust = dclust[rownames(ddesc)]# take only active chemical
  dclust = cbind(names(dclust), dclust)
  write.table(dclust, paste(pr_out, "SOM_Clusters", sep = ""), sep = ",", row.names = FALSE)
  
  return(som_model)
}



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
  
  # enrichment with a fisher score
  nbInact = length(which(is.na(dAC50)))
  nbAct = dim(dAC50)[1] - nbInact
  lpval = NULL
  nbclust = ydim*xdim
  for(clust in seq(1,nbclust)){
    dclust = dAC50[which(som_model$unit.classif == clust),]
    ninactClust = length(which(is.na(dclust[,2])))
    nactClust = dim(dclust)[1] - ninactClust
    FisherMat = matrix(c(nbAct, nbInact, nactClust, ninactClust),
                        nrow = 2,
                        dimnames = list(Pop = c("Act", "Inact"),
                                        Clust = c("Act", "Inact")))
    p = fisher.test(FisherMat, alternative = "greater")
    pval = p$p.value
    lpval = append(lpval, p$p.value)
  }
  svg(paste(pr_out, "SOM_enrich_act.svg", sep = ""))
  plot(som_model, type = "property", property=log10(lpval), palette.name=coolBlueHotRed, main = "Enrichment log(pvalue(Fisher test))", dpi=300, height = 20, width = 20)
  dev.off()
}




################
#     MAIN     #
################

args <- commandArgs(TRUE)
p_desc = args[1]
p_AC50 = args[2]
pr_out = args[3]
size = as.integer(args[4])


p_desc = "../../results/Cleaned_Data/desc1D2D_cleaned.csv"
p_AC50 = "../../results/Cleaned_Data/AC50_cleaned.csv"
pr_out = "../../results/SOM/"
size = 15


# open files
ddesc = read.csv(p_desc, sep = ",", header = TRUE)
rownames(ddesc) = ddesc[,1]
ddesc = ddesc[,-1]

dAC50 = read.csv(p_AC50, sep = ",", header = TRUE)
rownames(dAC50) = dAC50[,1]
dAC50 = dAC50[,-1]


# define model
p_model = paste(pr_out , "SOM_model.RData", sep = "")
if(!file.exists(p_model)){
  som_model = generateSOM(ddesc, size, size, pr_out)
}else{
  load(p_model)
}


# put active in the SOM
applySOM(som_model, dAC50, pr_out)




