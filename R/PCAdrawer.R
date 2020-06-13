#!/usr/bin/env Rscript
source ("graph.R")
library(fastICA)
library(cowplot)
library(plotly)
library(dplyr)


##################
# PCA functions  #
##################

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

factorACP = function (coor_point, vector_arrows){
  
  factor = 1
  orgin = vector_arrows
  while (max (vector_arrows[,1]) < max (coor_point[,1]) && max (vector_arrows[,2]) < max (coor_point[,2]) && min (vector_arrows[,1]) > min (coor_point[,1]) && min (vector_arrows[,2]) > min (coor_point[,2]) ){
    factor = factor + 1
    vector_arrows[,1] = vector_arrows[,1] + orgin[,1]
    vector_arrows[,2] = vector_arrows[,2] + orgin[,2]
    
  }
  return (factor-1)
  
}

generateIPAcoords = function(din, path_result){
  
  a = fastICA(din, 3, alg.typ = "parallel", fun = "logcosh", alpha = 1, method = "R", row.norm = FALSE, maxit = 200, tol = 0.0001, verbose = TRUE)
  print (a)
  gifGeneration(paste(path_result, "ICA3D", sep = ""), a$S)
}

PCAcombined3plans = function(d1D, d2D, d3D, pfilout){
  
  lcoord1D = generatePCAcoords(d1D)
  lcoord2D = generatePCAcoords(d2D)
  lcoord3D = generatePCAcoords(d3D)
  
  coordSpace = cbind(lcoord1D[[1]][,1], lcoord2D[[1]][,1])
  coordSpace = cbind(coordSpace, lcoord3D[[1]][,1])
  
  print(paste(lcoord1D[[2]], lcoord2D[[2]], lcoord3D[[2]], sep = "%  "))
  
  gifGeneration(paste(pfilout, "PCA3D", sep = ""), coordSpace)
}

PCAcombined2plans = function(d1D, d2D, pfilout){
  
  lcoord1D = generatePCAcoords(d1D)
  lcoord3D = generatePCAcoords(d2D)
  
  coordSpace = cbind(lcoord1D[[1]][,1], lcoord1D[[1]][,2])
  coordSpace = cbind(coordSpace, lcoord3D[[1]][,1])
  
  #print(paste(lcoord1D[[2]], lcoord2D[[2]], lcoord3D[[2]], sep = "%  "))
  
  gifGeneration(paste(pfilout, "PCA3D2plan", sep = ""), coordSpace)
  
  
}

PCA3D = function(din, path_result){
  
  dinScale = scale(din)
  
  data.cor=cor(dinScale)
  data.eigen=eigen(data.cor)
  lambda = data.eigen$values
  var_cap = lambda/sum(lambda)*100
  cp = data.eigen$vectors
  rownames (cp) = colnames (dinScale)
  colnames (cp) = colnames (dinScale)
  data_plot = as.matrix(dinScale)%*%cp
  
  #col.desc = colorDesc(colnames(din))
  
  col.desc = "black"
  
  gifGeneration(paste(path_result, "PCA3D", sep = ""), data_plot)
  
}

# PCA interactive html

plotPCA <- function(
  my_pca,
  colorby = NULL,
  xaxs = "PC1", yaxs = "PC2",
  title = NULL,
  palette = function(x) rainbow(x, s = 0.6),
  continuous_color = FALSE,
  ...
) {
  
  # obtaining % of variance explained by each axis
  varxaxs <- (summary(my_pca)$importance["Proportion of Variance", xaxs] * 100) %>% round(digits = 1)
  varyaxs <- (summary(my_pca)$importance["Proportion of Variance", yaxs] * 100) %>% round(digits = 1)
  
  myplot <- tibble(
    myxaxs = my_pca$x[, xaxs],
    myyaxs = my_pca$x[, yaxs],
    texts = rownames(my_pca$x),
    colors = if (is.null(colorby)) {NA} else {colorby}
  ) %>%
    ggplot(aes(x = myxaxs, y = myyaxs, color = colors, text = texts)) +
    geom_point() +
    labs(x = paste0(xaxs, " (", varxaxs, "%)"), y = paste0(yaxs, " (", varyaxs, "%)"), title = title, color = NULL) +
    background_grid()
  
  if (continuous_color) {
    myplot <- myplot + scale_color_gradientn(colors = palette(64))
  } else {
    myplot <- myplot + scale_color_manual(values = palette(n_distinct(colorby)))
  }
  
  ggplotly(myplot, tooltip = "text", ...)
  
}

# generating data for the example
metadata <- tibble(
  group = rep(c("group1", "group2", "group3"), each = 50),
  names = paste("sample", 1:150),
  latent_variable = c(
    rnorm(100) + 4,
    rnorm(50)
  )
)
#my_data <- matrix(rnorm(150 * 15), nrow = 150)
#my_data[1:100, 1:10] <- my_data[1:50, 1:10] + 3 # PC1 axis differences
#my_data[1:50, 11:15] <- my_data[1:50, 11:15] + 2 # PC2 differences 
#rownames(my_data) <- metadata$names

# basic example (click on legend to hide/display specific groups):
#plotPCA(prcomp(my_data), colorby = metadata$group)

# changing axis and title
#plotPCA(prcomp(my_data), colorby = metadata$group, xaxs = "PC2", yaxs = "PC3", title = "PC1 not shown")

# changing color palette:
#plotPCA(prcomp(my_data), colorby = metadata$group, palette = viridis::viridis)

# continuous colors:
#plotPCA(prcomp(my_data), colorby = metadata$latent_variable, palette = viridis::viridis, continuous_color = TRUE)

# saving as html:
#p <- plotPCA(prcomp(my_data), colorby = metadata$group)
#htmlwidgets::saveWidget(p, "my_great_pca.html", selfcontained = TRUE)
