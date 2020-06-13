#!/usr/bin/env Rscript

# By Alexandre BORREL
# 02-2016 => update 08-2018

library(MASS)
require(plotrix)
require(lattice)
library(scatterplot3d)
library(ggplot2)

source("PCAdrawer.R")
source("elimcor_sansY.R")


####################
# Opener functions #
####################

# separate data with class group
separeData = function (data, descriptor_class){
  
  data1 = data [which(data[,descriptor_class] == 0),]
  data2 = data [which(data[,descriptor_class] == 1),]
  
  return (list (data1, data2))
}

# most used
openData = function (pfilin, valcor, prout, NbmaxNA=10){
  desc = read.csv (pfilin, header = TRUE, sep = "\t", stringsAsFactors = TRUE)
  
  if(dim(desc)[2] ==1){
    # case of fucking ,
    desc = read.csv (pfilin, header = TRUE, sep = ",", stringsAsFactors = TRUE)
  }
  
  if(length(which(duplicated(desc[,1]))) == 0){
    rownames(desc) = desc[,1]
    desc = desc[,-1]  
  }
  
  
  # deleted col with NA
  lcoldel = NULL
  print (dim(desc))
  for (icol in seq(1, dim(desc)[2])){
    #print (desc[,icol])
    if (sum(is.na(as.vector(desc[,icol]))) > NbmaxNA){
      lcoldel = append(lcoldel, icol)
    }
  }
  
  print("=== colnames deleted ===")
  print(colnames(desc)[lcoldel])
  
  if(is.null(lcoldel)){
    desc = na.omit(desc)
  }else{
    desc = desc[,-lcoldel]	  
    desc = na.omit(desc)  
  }
  
  
  
  # dell when sd = 0
  sd_desc = apply (desc[,1:(dim(desc)[2])], 2, sd)
  
  #print (sd_desc)
  #print ("--------")
  sd_0 = which (sd_desc == 0.0)
  
  #print ("------------")
  #print (mode(sd_0))
  #print (length (sd_0))
  #print ("------------")
  if (length(sd_0) != 0){
    #print (as.factor (sd_0))
    #desc = desc[,-sd_0]
    desc=subset(desc,select=-sd_0)
  }
  if (valcor != 0){
    out_elimcor = elimcor_sansY (desc, valcor)
    descriptor = out_elimcor$possetap
    
    MDSElimcor (desc, out_elimcor, paste (prout, "MDSDesc_", valcor, sep = ""), "corr")
    descriptor = colnames (desc) [descriptor]
    desc = desc[,descriptor]
    #print (dim(desc))
  }
  
  # again with SD null
  sd_desc =apply (desc[,1:(dim(desc)[2])], 2, sd)
  sd_0 = which (sd_desc == 0)
  if (length(sd_0) != 0){
    desc=subset(desc,select=-sd_0)
  }
  return (list((desc),colnames (desc)))
}

MDSElimcor = function (din, out_elimcor, path_file, type_dist){
  
  groupe_elimcor = out_elimcor$groupes
  descriptor_selected = out_elimcor$possetap
  
  #print(descriptor_selected)
  name_descriptor = colnames(din)[descriptor_selected]
  din = din[,descriptor_selected]
  
  din = na.omit (din)
  if (type_dist == "corr"){
    MC = cor(din)
    dist1 = abs(1-MC)  
  }
  
  color_desc = "black" #colorDesc(name_descriptor)
  #print (color_desc)
  
  png (paste (path_file, "_", type_dist, ".png", sep = ""), 3500, 3500, res = 180, bg = "white")
  par( mar=c(6,6,6,6))
  fit <- cmdscale(as.dist(dist1),eig=TRUE, k=2)
  #c = as.vector (col.desc[,colnames (data)])
  plot (fit$points[,1], fit$points[,2], main="MDS Descriptor", xlab = "DIM 1", ylab = "DIM 2", cex.lab = 2.75, cex.axis = 2, cex.main = 2.75, type = "n", xlim = c(-1.3, 1.3), ylim = c(-1, 1))
  #text (fit$points[,1], fit$points[,2]+0.05, labels = num_desc,  cex = 1.8, col = color_desc, font = 12)
  #text (fit$points[descriptor_selected,1], fit$points[descriptor_selected,2]+0.02, labels = "*",  cex = 4, col = color_desc)
  text (fit$points[,1], fit$points[,2], labels = name_descriptor,  cex = 2.6, col = color_desc)
  
  dev.off()
}

openDataVexcluded = function (pfilin, valcor, prout, vexclude, NbmaxNA=100){
  
  desc = read.csv (pfilin, header = TRUE, sep = "\t", stringsAsFactors = TRUE)
  print(dim(desc))
  
  # remove chemical with only NA
  desc = delete.na(desc, dim(desc)[2]-20)
  print(dim(desc))
  
  # remove col not well computed
  desc = t(delete.na(t(desc), NbmaxNA))
  print(dim(desc))
  #print(desc[2,])
  
  # deleted line with NA
  #rownames (desc) = seq (1, dim(desc)[1])
  desc = na.omit(desc)
  #print (dim(desc))
  cexclude = desc[,vexclude]
  desc = desc[,-vexclude]
  
  print(dim(desc))
  
  # dell when sd = 0
  sd_desc = apply (desc[,1:(dim(desc)[2])], 2, sd)
  
  #print (sd_desc)
  #print ("--------")
  sd_0 = which (sd_desc == 0)
  
  #print (sd_0)
  
  #print ("------------")
  #print (mode(sd_0))
  #print (length (sd_0))
  #print ("------------")
  if (length(sd_0) != 0){
    #print (as.factor (sd_0))
    #desc = desc[,-sd_0]
    desc=subset(desc,select=-sd_0)
    #cexclude = subset(cexclude,select=-sd_0)
    #print(dim(desc_new))
  }
  desc = apply(desc,2,as.double)
  
  print(dim(desc))
  print(valcor)
  if (valcor != 0){
    out_elimcor = elimcor_sansY (desc, valcor)
    descriptor = out_elimcor$possetap
    
    #MDSElimcor (desc, out_elimcor, paste (prout, "MDSDesc_", valcor, sep = ""), "corr")
    descriptor = colnames (desc) [descriptor]
    desc = desc[,descriptor]
  }
  
  # again with SD null
  sd_desc = apply (desc[,1:(dim(desc)[2])], 2, sd)
  sd_0 = which (sd_desc == 0)
  if (length(sd_0) != 0){
    desc=subset(desc,select=-sd_0)
  }
  
  
  desc = cbind(cexclude, desc)
  return (list((desc),colnames (desc)))
}

delSDNull = function(desc){
  
  sd_desc = apply (desc[,1:(dim(desc)[2])], 2, sd)
  
  #print (sd_desc)
  #print ("--------")
  sd_0 = which (sd_desc == 0)
  
  #print (sd_0)
  
  #print ("------------")
  #print (mode(sd_0))
  #print (length (sd_0))
  #print ("------------")
  if (length(sd_0) != 0){
    #print (as.factor (sd_0))
    #desc = desc[,-sd_0]
    desc=subset(desc,select=-sd_0)
    #cexclude = subset(cexclude,select=-sd_0)
    #print(dim(desc_new))
  }
  return(desc)
}

reduceNBdesc = function(desc, nbCV){ # reduce number of descriptor -nbCV percent less than number of chemicals
  
  # control number of descriptors VS number of data
  nbdesc = dim(desc)[2]
  nbchemical = dim(desc)[1]
  valcor = 0.9
  dtemp = desc
  
  while((((100-nbCV)/100*nbchemical)-10) <= nbdesc){ # less descriptors than CV 10 fold
    ldesc = colnames(dtemp)[elimcor_sansY (dtemp, valcor)$possetap]
    dtemp = dtemp[,ldesc]
    nbdesc = dim(dtemp)[2]
    print(paste("Limit nb desc regresssion: ", nbdesc, sep = ""))
    print(paste("Corval: ", valcor, sep = ""))
    valcor = valcor-0.01
  }
  return(colnames(dtemp))
  
}



##################
# base functions #
##################

delete.na <- function(DF, n=0) {
  DF[rowSums(is.na(DF)) <= n,]
}

is.integer0 <- function(x)
{
  is.integer(x) && length(x) == 0L
}

eucdist = function(x1, x2, y1, y2){
  distout = sqrt((x1-x2)^2 + (y1-y2)^2)
  return (distout)
}


###################################
# base functions for data manager #
###################################


changeProbaList = function (list_element){
  Y = list_element
  Y2 = rep(1,length(Y))
  Y2[which(Y<0.5)]=0
  return (Y2)
}

# del colum
delCol = function (data, list_name_col){
  for (name_col in list_name_col){
    data = data[,-which(colnames(data)==name_col)]
  }
  return (data)
}


delnohomogeniousdistribution = function(din, cutoff = 80){
  
  countMax = dim(din)[1]*cutoff/100  
  
  i = 1
  imax = dim(din)[2]
  while(i <= imax){
    #print (i)
    #print (as.vector(din[,i]))
    qt = hist(din[,i], breaks = 10, plot = FALSE)$counts
    
    for (qtc in qt){
      if (qtc >= countMax){
        din = din[,-i]
        imax = imax - 1
        i = i - 1
        break()
      }
    }
    i = i + 1
  }
  return(din)
}


mergeIdenticRow = function(din){
  dout = NULL
  for(n in unique(din[,1])){
    dtemp = din[which(din[,1]==n),]
    if(length(which(is.na(dtemp[,2])== dim(dtemp)[1]))){
      dout = rbind(dout, dtemp[1,])
    }else{
      dtemp[1,2] = mean(dtemp[-which(is.na(dtemp[,2])),2])
      dout = rbind(dout, dtemp[1,])
    }
  }
  return(dout)
}



#####################
# Sampler functions #
#####################

sampligDataMulticlassCV = function(din, nbfold, nameCol){
  
  dclass = din[,nameCol]
  lclass = unique(dclass)
  
  lfold = list()
  
  
  lNbElement = NULL
  for(class in lclass){
    nbelem = length(dclass[which(dclass == class)])
    lNbElement = append(lNbElement, round(nbelem /nbfold))
  }
  names(lNbElement) = lclass
  
  lrandom = sample(dim(din)[1])
  
  
  
  for(i in lrandom){
    Score = din[i,nameCol]
    ifold = 1
    while (ifold <= nbfold){
      if(length(lfold) < ifold){
        lfold[[ifold]] = din[i,]
        ifold = 11
      }else{
        nbinfold = length(which(lfold[[ifold]][,nameCol] == Score))
        if(ifold == 10){
          lfold[[ifold]] = rbind(lfold[[ifold]], din[i,])  
          ifold = 11
        }else if(nbinfold < lNbElement[as.character(Score)]){
          lfold[[ifold]] = rbind(lfold[[ifold]], din[i,])  
          ifold = 11
        }else{
          ifold = ifold + 1
        }
      }
    }
  }
  return(lfold) 
}



sampligDataMulticlass = function(din, fract, nameCol){
  
  dclass = din[,nameCol]
  lclass = unique(dclass)
  
  lNbElement = NULL
  for(class in lclass){
    nbelem = length(dclass[which(dclass == class)])
    lNbElement = append(lNbElement, round(nbelem * fract))
  }
  names(lNbElement) = lclass
  
  lrandom = sample(dim(din)[1])
  
  dtrain = NULL
  dtest = NULL
  for(i in lrandom){
    Score = din[i,nameCol]
    nbintest = length(which(dtest[,nameCol] == Score))
    if(nbintest < lNbElement[as.character(Score)]){
      dtest = rbind(dtest, din[i,])
    }else{
      dtrain = rbind(dtrain, din[i,])
    }
    
  }
  return(list(dtrain, dtest)) 
}



samplingDataNgroupClass = function (t_din, i_nb_group, s_nameclass){
  
  # divise two classes
  v_class = as.factor(t_din[,s_nameclass])
  t_dc0 = t_din [which(v_class == 0),]
  t_dc1 = t_din [which(v_class == 1),]
  
  
  # sample data
  v_sampledc0 = sample (dim (t_dc0)[1])
  v_sampledc1 = sample (dim (t_dc1)[1])
  
  # ind limit
  i_limitc0 = as.integer (dim(t_dc0)[1] / i_nb_group)
  i_limitc1 = as.integer (dim(t_dc1)[1] / i_nb_group)
  
  #print (i_limitc0)
  #print (i_limitc1)
  
  output = list ()
  for (g in 1:i_nb_group){
    #print (g)
    # start selct 1
    if (g == 1 ){
      t_group = rbind (t_dc0[v_sampledc0[1:i_limitc0],], t_dc1[v_sampledc1[1:i_limitc1],])
    }
    # last end to number of coulumn
    else if (g == i_nb_group){
      #print ("inf")
      #print (i_limitc0 * (g-1) + 1)
      #print (i_limitc1 * (g-1) + 1)
      #print ("sup")
      #print (length (v_sampledc0))
      #print (length (v_sampledc1))
      #print ("**IC**")
      #print ((i_limitc0 * (g-1) + 1):(length (v_sampledc0)))
      #print ((i_limitc1 * (g-1) + 1):(length (v_sampledc1)))
      
      
      t_group = rbind (t_dc0[v_sampledc0[((i_limitc0 * (g-1) + 1):(length (v_sampledc0)))],], t_dc1[(v_sampledc1[(i_limitc1 * (g-1) + 1):(length (v_sampledc1))]),])
    }
    else{
      t_group = rbind (t_dc0[(v_sampledc0[(i_limitc0 * (g-1) + 1):(i_limitc0 * g)]),], t_dc1[(v_sampledc1[(i_limitc1 * (g-1) + 1):(i_limitc1 * g)]),])
    }
    # append list
    output[[g]] = t_group
  }
  
  return (output)
}

sampligDataFractionCluster = function(t_din, fract, pcluster){
  
  dclust = read.csv(pcluster, sep=",", header = TRUE)
  print (t_din)
  lclust = unique(dclust[,2])
  
  dtrain = NULL
  dtest = NULL
  for (clust in lclust){
    lcpIDclust = dclust[which(dclust[,2] == clust),1]
    descclust = t_din[lcpIDclust,]
    
    # sample data
    v_sample = sample (dim (descclust)[1])
    
    # ind limit
    i_limitc = round((dim (descclust)[1]) * fract)
    
    dtrain = rbind(dtrain, descclust[v_sample[(i_limitc + 1):(length (v_sample))],])
    dtest = rbind(dtest, descclust[v_sample[1:i_limitc],])
  }
  # sample data
  return (list(dtrain, dtest))
  
}

samplingDataFraction = function (t_din, fract){
  
  # sample data
  v_sample = sample (dim (t_din)[1])
  
  # ind limit
  i_limitc = round((dim (t_din)[1]) * fract)
  
  dtrain = t_din[v_sample[(i_limitc + 1):(length (v_sample))],]
  dtest = t_din[v_sample[1:i_limitc],]
  
  return (list(dtrain, dtest))
  
}

samplingDataNgroup = function (t_din, i_nb_group){
  
  # sample data
  v_sample = sample (dim (t_din)[1])
  
  # ind limit
  i_limitc = as.integer (dim(t_din)[1] / i_nb_group)
  
  output = list ()
  for (g in 1:i_nb_group){
    #print (g)
    # start selct 1
    if (g == 1 ){
      t_group = t_din[v_sample[1:i_limitc],]
    }
    # last end to number of coulumn
    else if (g == i_nb_group){
      #print ("inf")
      #print (i_limitc0 * (g-1) + 1)
      #print (i_limitc1 * (g-1) + 1)
      #print ("sup")
      #print (length (v_sampledc0))
      #print (length (v_sampledc1))
      #print ("**IC**")
      #print ((i_limitc0 * (g-1) + 1):(length (v_sampledc0)))
      #print ((i_limitc1 * (g-1) + 1):(length (v_sampledc1)))
      
      
      t_group = t_din[v_sample[((i_limitc * (g-1) + 1):(length (v_sample)))],]
    }
    else{
      t_group = t_din[(v_sample[(i_limitc * (g-1) + 1):(i_limitc * g)]),]
    }
    # append list
    output[[g]] = t_group
  }
  
  return (output)
}


#################################
#   Control dataset integrity   #
#################################

controlDatasets = function(ldataset, prin){
  
  nbsplit = length(ldataset)
  
  pdf(paste(prin,".pdf", sep = ""), width = 10, height = 10)
  colorrainbow = rainbow(nbsplit)
  i = 1
  colorpoint = NULL
  dPCA = NULL
  dglobal = NULL
  for (d in ldataset){
    dglobal = rbind(dglobal,d)
    h = ggplot(d, aes(x=Aff)) + 
      geom_histogram(aes(y=..density..),     # Histogram with density instead of count on y-axis
                     binwidth=.5,
                     colour="black", fill="white") +
      labs(x = "pAff", y = "Frequencies") + 
      ggtitle(paste("Fold ", i, "-Dim=", dim(d)[1], sep = ""))+
      theme(axis.text.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.text.x = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.x =  element_text(size = 25, hjust = 0.5, vjust =0.1))+
      geom_density(alpha=.2, fill="#FAFFA5")
    print(h)
    
    #points for PCA
    colorpoint = as.factor(d[,dim(d)[2]])#append(colorpoint, rep(colorrainbow[i], (dim(d)[1]))) +> give multicolor
    dPCA = rbind(dPCA, d[,-dim(d)[2]]) # remove affinity coulum
    i = i + 1
  }
  
  h = ggplot(dglobal, aes(x=Aff)) + 
    geom_histogram(aes(y=..density..),     # Histogram with density instead of count on y-axis
                   binwidth=.5,
                   colour="black", fill="white") +
    labs(x = "pAff", y = "Frequencies") + 
    ggtitle(paste("Fold ", i, "-Dim=", dim(dglobal)[1], sep = "")) +
    theme(axis.text.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.text.x = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.x =  element_text(size = 25, hjust = 0.5, vjust =0.1))+
    geom_density(alpha=.2, fill="#FAFFA5")
  print (h)
  
  # plot PCA
  dplot = generatePCAcoords(dPCA)
  var_cap = dplot[[2]]
  data_plot = dplot[[1]]
  #par(mar=c(8,8,8,8))
  plot(data_plot[,1],data_plot[,2], pch=20, col = colorpoint, xlab = paste("CP1: ", signif (var_cap[1], 4), "%", sep = ""), ylab = paste("CP2: ", signif (var_cap[2], 4), "%", sep = ""), cex.lab = 1.5, cex.main = 2, cex.axis = 1.5, cex = 2)
  abline(h=0,v=0)
  dev.off()
  
}

#######################
# order data rownames #
#######################

orderByType = function (d.col, d){
  l_temp = NULL
  
  for (desc in colnames (d.col)){
    if (is.integer0 (which(rownames (d) == desc))== FALSE){
      l_temp = append (l_temp, desc)
    }
  }
  if (length (l_temp) != length (rownames (d))){
    print ("TO DO -> return list desc")
    return (l_temp)
  }
  d = d[l_temp,]
  #names (d) = l_temp
  return (d)
}

orderByTypeOneCol = function (d.col, d){
  l_temp = NULL
  
  for (desc in colnames (d.col)){
    if (is.integer0 (which(rownames (d) == desc))== FALSE){
      l_temp = append (l_temp, desc)
    }
  }
  if (length (l_temp) != length (rownames (d))){
    print ("TO DO -> return list desc")
    return (l_temp)
  }
  d = d[l_temp,]
  names (d) = l_temp
  return (d)
}



####################
# general function #
####################

calculLimListData = function (data1, data2){

	x1 = min (data1)
	x2 = min (data2)
	X1 = max(data1)
	X2 = max (data2)
	l = c (min (x1,x2), max (X1,X2))
	return (l)
}

# check color vector
CheckColorVector = function (l_descriptor, col.des){
	for (d in l_descriptor){
		if (is.integer0(which(d == names(col.des))) == TRUE){
			out = rep (1,length (l_descriptor))
			names (out) = l_descriptor
			return (out)
		}
	}
	return (col.des)
}

#remove SD = 0
delSDnull = function(desc){
  
  sd_desc = apply (desc[,1:(dim(desc)[2])], 2, sd)
  sd_0 = which (sd_desc == 0)
  
  if ( !is.integer0(sd_0)){
    desc=subset(desc,select=-sd_0)
  }
  descout = apply(desc,2,as.double)
  rownames(descout) = rownames(desc)
  return(descout) 
}
