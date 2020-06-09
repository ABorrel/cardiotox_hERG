import pathFolder
import dataset



# Define folder
################
PR_ROOT = "../../"
PR_DATA = PR_ROOT + "data/"
PR_RESULTS = pathFolder.createFolder(PR_ROOT + "results/")


# MAIN
#######

# 1. process dataset
#####

# file with smi and AC50
p_smi = PR_DATA + "list_chemicals-2020-06-09-09-11-41.csv"
p_AC50 = PR_DATA + "AC50_7403.txt"
pr_dataset = pathFolder.createFolder(PR_RESULTS + "dataset/")

# load dataset
c_dataset = dataset.dataset(p_smi, p_AC50, pr_dataset)
c_dataset.prep_dataset()

# compute desc
pr_desc = pathFolder.createFolder(PR_RESULTS + "DESC/")
c_dataset.computeDesc(pr_desc)


# 2. analysis
#####