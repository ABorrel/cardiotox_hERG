#!/usr/bin/env Rscript
library(ggplot2)
library(tools)

################
#     MAIN     #
################
args <- commandArgs(TRUE)
p_perf = args[1]
#p_perf = "c://Users/aborr/research/ILS/HERG/results/NCAST/DNN_model_reg/1/test_pred.csv"
p_out = tools::file_path_sans_ext(p_perf)


dperf = read.csv(p_perf, sep="\t", row.names =  1)
dperf = na.omit(dperf)

cor_perf = cor(dperf$Real, dperf$Pred)

p = ggplot(dperf, aes(Real, Pred, label=rownames(dperf)))+
  geom_point(size=1.5, colour="black", shape=19) + 
  geom_text(x=4.3, y=7.2, label = paste("cor=",round(cor_perf,2), sep = ""), size = 8)+
  labs(x = "pAff", y = "Predicted pAff") +
  theme(axis.text.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.text.x = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.x =  element_text(size = 25, hjust = 0.5, vjust =0.1))+
  xlim (c(4, 7.5))  +
  geom_segment(aes(x = 4, y = 4, xend = 7.5, yend = 7.5), linetype=2, size = 0.5) + 
  ylim (c(4, 7.5))  
#print(p)
ggsave(paste(p_out, "_corplot.png", sep = ""), width = 8,height = 8, dpi = 300)

