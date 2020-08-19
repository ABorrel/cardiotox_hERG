#!/usr/bin/env Rscript
library(ggplot2)



################
#     MAIN     #
################

args <- commandArgs(TRUE)
p_overlap = args[1]


p_overlap = "../../results/AID_588834/OverlapModel/overlap_train_test"
  
  
d_overlap = read.csv(p_overlap, sep = "\t")
d_overlap$Include = as.double(as.character(d_overlap$Include))

# make a sum file
nb_overlap = length(which(d_overlap$Include==1))
nb_nooverlap = length(which(d_overlap$Include==0))
per_overlap = round(nb_overlap*100/(nb_overlap + nb_nooverlap),2)

# write sum
c_w = c(nb_overlap, nb_nooverlap, per_overlap)
names(c_w) = c("Overlap", "no_overlap", "per active")

write.csv(c_w, file = paste(p_overlap, ".sum", sep = ""))


d_overlap$AC50.model = as.double(as.character(d_overlap$LogAC50.model))
d_overlap$AC50.test = as.double(as.character(d_overlap$LogAC50.test))
d_overlap$Include = as.factor(as.character(d_overlap$Include))


d_overlap = na.omit(d_overlap)
cor_val = round(cor(d_overlap$AC50.model, d_overlap$AC50.test),2)


ggplot(d_overlap, aes(x=LogAC50.model, y=LogAC50.test)) + 
  geom_point()+
  xlab("Log(IC50) model")+
  ylab("Log(IC50) test")+
  theme(axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12))+
  geom_text(x=4.5, y=7.5, label= paste("Cor: ", cor_val, sep = ""), size=6)+
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=1.5)+
ggsave(paste(p_overlap, ".png", sep = ""))

