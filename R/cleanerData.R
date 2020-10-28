#!/usr/bin/env Rscript
library(Toolbox)


################
#     MAIN     #
################

args <- commandArgs(TRUE)
pdesc = args[1]
pdata = args[2] #to take IC50 or class
prout = args[3]
valcor = args[4]
maxquantile = as.double(args[5])
act = as.integer(args[6])

#pdesc = "../../results/DESC/desc_1D2D.csv"
#pdata = "../../data/AC50_7403.txt"
#prout = "../../results/Cleaned_Data/"
#valcor = 0.9
#maxquantile = 90


##############################
# Process descriptors matrix #
##############################

dglobal = openDataVexcluded(pdesc, valcor, prout, c(1,2))
dglobal = dglobal[[1]]
lCAS = dglobal[,1]
dglobal = dglobal[,-1]
dglobal = dglobal[,-1]

# format in double
dglobal =apply(dglobal, 2, as.numeric)
rownames(dglobal) = lCAS


print("==== Preprocessing ====")
print(paste("Data initial: dim = ", dim(dglobal)[1], dim(dglobal)[2], sep = " "))

##########
# filter #
##########

dglobal = delnohomogeniousdistribution(dglobal, maxquantile)
print(paste("Data after filtering: dim = ", dim(dglobal)[1], dim(dglobal)[2], sep = " "))

#######################
# order with affinity #
#######################
# Opening
daffinity = read.csv(pdata, sep = "\t", header = TRUE)
rownames(daffinity) = daffinity[,1]
colnames(daffinity)[2] = "Aff"
# only extract active
if(act == 1){
    daffinity = na.omit(daffinity)
}

# control same chemicals
linter = intersect(rownames(dglobal), rownames(daffinity))

dclean = dglobal[linter,]
daffinity = daffinity[linter,]


# Write table 
if (act == 1){
    paffout = paste(prout, "AC50_act_cleaned.csv", sep = "")
    pdesout = paste(prout, "desc_act_cleaned.csv", sep = "")
}else{
    paffout = paste(prout, "AC50_cleaned.csv", sep = "")
    pdesout = paste(prout, "desc_cleaned.csv", sep = "")
}
write.csv(daffinity, paffout, col.names = TRUE, row.names = TRUE)
write.csv(dclean, pdesout, col.names = TRUE, row.names = TRUE)
