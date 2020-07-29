#!/usr/bin/env Rscript
library(ggplot2)


#########
# main
#########

args <- commandArgs(TRUE)
p_chemClasses = args[1]

#p_chemClasses = "../../results/dataset/class_class_from_interferences"

d_chem = read.table(p_chemClasses, sep = "\t", header = TRUE)
#d_chem = na.omit(d_chem)

l_class = table(d_chem$Class)
write.csv(l_class, paste(p_chemClasses, "_count", sep = ""))

l_class_remove = names(l_class)[which(l_class < 5)]


for(class_chem in l_class_remove){
  d_chem = d_chem[-which(d_chem$Class == class_chem),]
}


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

