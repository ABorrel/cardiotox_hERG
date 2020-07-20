#!/usr/bin/env Rscript
library(ggplot2)



################
#     MAIN     #
################

args <- commandArgs(TRUE)
p_overlap = args[1]


#p_overlap = "./../../results/ChEMBL_patch_clamp/predict/OverlapModel/overlap_train_test"
  
  
d_overlap = read.csv(p_overlap, sep = "\t")

# make a sum file
nb_overlap = length(which(d_overlap$Include==1))
nb_nooverlap = length(which(d_overlap$Include==0))
per_overlap = round(nb_overlap*100/(nb_overlap + nb_nooverlap),2)

# write sum
c_w = c(nb_overlap, nb_nooverlap, per_overlap)
names(c_w) = c("Overlap", "no_overlap", "per active")

write.csv(c_w, file = paste(p_overlap, ".sum", sep = ""))


d_overlap$AC50.model = as.double(as.character(d_overlap$AC50.model))
d_overlap$AC50.test = as.double(as.character(d_overlap$AC50.test))
d_overlap$Include = as.factor(as.character(d_overlap$Include))

d_overlap = d_overlap[which(d_overlap$AC50.test <= 30),]


ggplot(d_overlap, aes(x=AC50.model, y=AC50.test)) + geom_point()
ggsave(paste(p_overlap, ".png", sep = ""))

