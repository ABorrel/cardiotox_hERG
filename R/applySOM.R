#!/usr/bin/env Rscript
library(kohonen)
library(ggplot2)


#################
# Color reverse #
#################

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

#p_desc = "./../../results/ChEMBL/desc_1D2D.csv"
#p_SOMmodel = "./../../results/SOM/SOM_model.RData"
#pr_out = "./../../results/ChEMBL_predict/SOM/"



SOMmodel = load(p_SOMmodel)
d_built = som_model$data[[1]]

d_desc = read.csv(p_desc,sep = "\t", row.names = 1)
d_desc = d_desc[,-1]
d_desc = as.matrix(d_desc)
d_desc = d_desc[, colnames(d_built)]

som_map = kohonen::map(som_model, d_desc)
l_count = table(som_map$unit.classif)

l_plot = rep(0,som_model$grid$xdim*som_model$grid$ydim)
names(l_plot) = seq(1, som_model$grid$xdim*som_model$grid$ydim)
l_plot[names(l_count)] = l_count

d_clust = som_map$unit.classif
names(d_clust) = rownames(d_desc)
write.csv(d_clust, paste(pr_out, "SOM_Clusters", sep = ""))

svg(paste(pr_out, "SOM_applied.svg", sep = ""))
plot(som_model, type = "property", property=l_plot, palette.name=coolBlueHotRed, main = "", dpi=300, height = 20, width = 20, bg = "transparent")
dev.off()
