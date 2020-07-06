#!/usr/bin/env Rscript
library(ggplot2)


#########
# main
#########

args <- commandArgs(TRUE)
p_chemClasses = args[1]

p_chemClasses = "./../../results/dataset/class_hERG_Active_Annotated_listRefChem_July1"

d_chem = read.table(p_chemClasses, sep = "\t", header = TRUE)
#d_chem = na.omit(d_chem)

d_chem <- within(d_chem, 
                   Classes <- factor(Class, 
                                      levels=names(sort(table(Class), 
                                                        decreasing=FALSE))))


# Barplot basique
p<-ggplot(data=d_chem, aes(x=Classes)) +
  geom_bar(stat="count")


# Barplot horizontal
p + coord_flip()

ggsave(paste(p_chemClasses, "_barplot.png", sep = ""),  width = 5, height = 7, dpi = 300, bg="transparent")
