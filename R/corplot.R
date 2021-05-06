#!/usr/bin/env Rscript
library(ggplot2)


################
#     MAIN     #
################
args <- commandArgs(TRUE)
p_cor = args[1]
#p_cor = "c://Users/aborr/research/ILS/HERG/results/comparison_CHEMBL_NCAST/comparison"


d_cor = read.csv(p_cor, sep="\t")
d_sum = NULL
d_sum$'Nb chemicals' = dim(d_cor)[1]
d_sum$'Nb active NCAST' = dim(d_cor)[1] - length(which(is.na(d_cor$Aff_NCAST)))
d_sum$'Nb active CHEMBL' = dim(d_cor)[1] - length(which(is.na(d_cor$Aff_CHEMBL)))
write.csv(d_sum, paste(p_cor, ".sum", sep = ""))


d_cor = na.omit(d_cor)
d_path_clamp = d_cor[which(d_cor$Assay_CHEMBL == "Patch clamp"),]

p = ggplot(d_cor, aes(Aff_NCAST, Aff_CHEMBL, color=Assay_CHEMBL))+
  geom_point(size=1.5, shape=19) + 
  theme(axis.text.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.text.x = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.x =  element_text(size = 25, hjust = 0.5, vjust =0.1))+
  labs(x = expression("pIC50 NCAST chemicals"), y =expression("pIC50 CHEMBL chemicals")) + 
  geom_segment(aes(x = 3, y = 3, xend = 9, yend = 9)) + 
  annotate("text", hjust = 0, size = 5, x=3.5, y=8.5, label= paste("r=", round(cor(d_cor$Aff_CHEMBL, d_cor$Aff_NCAST), digits = 2)))+
  annotate("text", hjust = 0, size = 5, x=3.5, y=8.0, label= paste("r=", round(cor(d_path_clamp$Aff_CHEMBL, d_path_clamp$Aff_NCAST), digits = 2), "path clamp"))
print(p)
ggsave(paste(p_cor, "_cor.png",sep=""), width = 9,height = 9, dpi = 300)

