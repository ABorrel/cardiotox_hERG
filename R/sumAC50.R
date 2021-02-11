#!/usr/bin/env Rscript
library(ggplot2)



################
#     MAIN     #
################

args <- commandArgs(TRUE)
p_AC50 = args[1]
pr_out = args[2]


#p_AC50 = "C:/Users/Aborrel/research/ILS/HERG/results/ChEMBL/aff_cleaned.csv"
#pr_out = "C:/Users/Aborrel/research/ILS/HERG/results/ChEMBL/"

dAC50 = read.csv(p_AC50, sep = ",", header = TRUE)
rownames(dAC50) = dAC50[,1]
dAC50 = na.omit(dAC50)

mu = mean(dAC50$Aff)

ggplot(dAC50, aes(x=Aff)) +
  geom_histogram(aes(y=..density..), position="identity", alpha=0.20, color="blue", fill="blue")+
  geom_density(alpha=0.3, fill="blue")+
  theme(text = element_text(size=19))+
  geom_vline(aes(xintercept=mu),
             linetype=c("dashed"))+
  scale_color_manual(values=c("#cde2ff"), labels = c("-log(AC50)")) + 
  labs(title="",x="-log(AC50) (uM)", y = "Density")

ggsave(paste(pr_out, "hist_AC50.png", sep = ""),  width = 8, height = 7, dpi = 300, bg="transparent")
