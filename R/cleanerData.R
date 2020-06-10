#!/usr/bin/env Rscript
source ("dataManager.R")


################
#     MAIN     #
################

args <- commandArgs(TRUE)
pdesc = args[1]
pdata = args[2] #to take IC50 or class
prout = args[3]
valcor = args[4]
maxquantile = as.double(args[5])

pdesc = "../../results/DESC/desc_1D2D.csv"
pdata = "../../data/AC50_7403.txt"
prout = "../../results/Cleaned_Data/"
valcor = 0.9
maxquantile = 90



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

# control same chemicals
linter = intersect(rownames(dglobal), rownames(daffinity))

dclean = dglobal[linter,]
daffinity = daffinity[linter,]

# Write table 
paffout = paste(prout, "AC50_cleaned.csv", sep = "")
write.csv(daffinity, paffout, col.names = TRUE, row.names = TRUE)


pdesout = paste(prout, "desc1D2D_cleaned.csv", sep = "")
write.csv(dglobal, pdesout, col.names = TRUE, row.names = TRUE)
