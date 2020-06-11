import pathFolder
import dataset
import analysis
import QSAR_modeling


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
p_desc = c_dataset.computeDesc(pr_desc)
#c_dataset.computePNG(pr_desc)


# 2. analysis
#####
COR_VAL = 0.90
MAX_QUANTILE = 90

cAnalysis = analysis.analysis(p_AC50, p_desc, PR_RESULTS, COR_VAL, MAX_QUANTILE)
cAnalysis.prepDesc()

# 2.1. histogram AC50 and summary
#cAnalysis.sumAC50()

# 2.2 PCA
#cAnalysis.PCA_plot()

# 2.3 SOM
#cAnalysis.generate_SOM()
#cAnalysis.analyse_SOM(pr_desc + "PNG/") # have to run !!!!

# 2.4 Hclust
#cAnalysis.HClust_plot()

# 3. QSAR modeling
#######
pr_QSAR = pathFolder.createFolder(PR_RESULTS + "QSAR/")
nb_repetition = 10
n_foldCV = 10
rate_split = 0.15
rate_active = 0.30
cQSAR = QSAR_modeling.QSAR_modeling(cAnalysis.p_desc_cleaned, cAnalysis.p_AC50_cleaned, pr_QSAR, nb_repetition, n_foldCV,rate_active, rate_split)
cQSAR.runQSARClass()

