import analysis
import applyModels
import CHEMBLTable
import dataset
import genericTestSet
import pathFolder
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

# 1.1 load dataset
c_dataset = dataset.dataset(p_smi, p_AC50, pr_dataset)
c_dataset.prep_dataset()

# 1.2 chemical classification
# classification from Interferences
p_classification = PR_DATA + "class_from_interferences.csv"
c_dataset.classChem(p_classification)

#####p_classification = PR_DATA + "hERG_Active_Annotated_listRefChem_July1.csv"
#####c_dataset.classActive(p_classification)


# 1.3 compute desc
pr_desc = pathFolder.createFolder(PR_RESULTS + "DESC/")
p_desc = c_dataset.computeDesc(pr_desc)
#c_dataset.computePNG(pr_desc)

# 1.4 rank active chemical and save in a PDF
#c_dataset.rankActiveChem(pr_desc + "PNG/")


# 2. analysis
#####
COR_VAL = 0.90
MAX_QUANTILE = 90

cAnalysis = analysis.analysis(p_AC50, p_desc, PR_RESULTS, COR_VAL, MAX_QUANTILE)
cAnalysis.prepDesc()

## 2.1. histogram AC50 and summary
#cAnalysis.sumAC50()

## 2.2 PCA
#cAnalysis.PCA_plot()

## 2.3 SOM
size = 15
#cAnalysis.generate_SOM(15)
#cAnalysis.signifDescBySOMCluster()
#cAnalysis.extract_actBySOMCluster(pr_desc + "PNG/") # have to run !!!!
#p_classification = PR_DATA + "class_from_interferences.csv"
#cAnalysis.applySOMtoChemClassification(p_classification, pr_desc + "PNG/")

#p_classification = PR_DATA + "hERG_Active_Annotated_listRefChem_July1.csv"
#cAnalysis.applySOMtoChemClassification(p_classification, pr_desc + "PNG/", filter_active=0)

## 2.4 Hclust
#cAnalysis.HClust_plot()


## 3. QSAR modeling
#######
# 3.1 using a new defintion of test set at each iteration
########################
pr_QSAR = pathFolder.createFolder(PR_RESULTS + "QSAR/")
nb_repetition = 10
n_foldCV = 10
rate_split = 0.15
rate_active = 0.30
cQSAR = QSAR_modeling.QSAR_modeling(cAnalysis.p_desc_cleaned, p_desc, cAnalysis.p_AC50_cleaned, p_AC50, pr_QSAR, nb_repetition, n_foldCV,rate_active, rate_split)
cQSAR.runQSARClassUnderSamplingAllSet()
ee
# pr_RF_models = cQSAR.extractModels(PR_RESULTS, "RF")
#pr_LDA_models = cQSAR.extractModels(PR_RESULTS, "LDA")


# 3.2 using a unique test set 
################################
#pr_QSAR = pathFolder.createFolder(PR_RESULTS + "QSAR_samplingTrain/")
#nb_repetition = 10
#n_foldCV = 10
#rate_split = 0.15
#rate_active = 0.30
#cQSAR = QSAR_modeling.QSAR_modeling(cAnalysis.p_desc_cleaned, p_desc, cAnalysis.p_AC50_cleaned, p_AC50, pr_QSAR, nb_repetition, n_foldCV,rate_active, rate_split)
#cQSAR.runQSARClassUnderSamplingTrain()
#pr_RF_samplingTrain_models = cQSAR.extractModels(PR_RESULTS, "RF_samplingTrain")
#pr_LDA_samplingTrain_models = cQSAR.extractModels(PR_RESULTS, "LDA_samplingTrain")


# 4. Build external test set from ChEMBL
### Target => CHEMBL240
##########################
#P_CHEMBL = PR_DATA + "CHEMBL27-target_chembl240.csv"
#pr_ChEMBL = pathFolder.createFolder(PR_RESULTS + "ChEMBL/")

## 4.1 Prep dataset
#cChEMBL = CHEMBLTable.CHEMBLTable(P_CHEMBL, pr_ChEMBL)
#cChEMBL.parseCHEMBLFile()
#cChEMBL.cleanDataset(l_standard_type=["IC50", "Ki"], l_standard_relation=["'='"])

## 4.2 run descriptor set
#p_desc_ChEMBL = cChEMBL.computeDesc()
#p_aff_ChEMBL = cChEMBL.cleanAff()

## 4.3 run PCA with ChEMBL on PCA tox21
#pr_applyModel = pathFolder.createFolder(PR_RESULTS + "ChEMBL_predict/")
#cApplyModel = applyModels.applyModel(cAnalysis.p_desc_cleaned, cAnalysis.p_AC50_cleaned, p_desc_ChEMBL, p_aff_ChEMBL, pr_RF_models, pr_applyModel)
#cApplyModel.PCACombine()

## 4.4 apply RF models on it
#cApplyModel.predict()
#cApplyModel.mergePrediction()

## 4.5. Apply SOM on it
#p_SOM = PR_RESULTS + "SOM/SOM_model.RData"
#cApplyModel.applySOM(p_SOM)


## 5. use drug external test sest 1
############################

#p_dataset = PR_DATA + "8Candidate_Cipa_drugs.csv"
#pr_TestSet = pathFolder.createFolder(PR_RESULTS + "DrugSet/")
#c_genericSet = genericTestSet.genericTestSet(p_dataset, pr_TestSet)
#c_genericSet.loadDataset()
#p_desc = c_genericSet.computeDesc()
#p_aff = c_genericSet.setAff(allAff=1)

## apply model
#cApplyModel = applyModels.applyModel(cAnalysis.p_desc_cleaned, cAnalysis.p_desc, cAnalysis.p_AC50_cleaned, p_desc, p_aff, pr_RF_models, pr_TestSet)
#cApplyModel.PCACombine()
#cApplyModel.predict()
#cApplyModel.mergePrediction()
#cApplyModel.computeAD()

## SOM
#p_SOM = PR_RESULTS + "SOM/SOM_model.RData"
#cApplyModel.applySOM(p_SOM)

## 6. use drug external test sest 2
#############################

#p_dataset = PR_DATA + "hERG_Active_Annotated_listRefChem_July1.csv"
#pr_TestSet = pathFolder.createFolder(PR_RESULTS + "hERG_Active_Annotated_listRefChem_July1/")
#c_genericSet = genericTestSet.genericTestSet(p_dataset, pr_TestSet)
#c_genericSet.loadDataset()
#p_desc = c_genericSet.computeDesc()
#p_aff = c_genericSet.setAff(allAff=1)


## apply model
#cApplyModel = applyModels.applyModel(cAnalysis.p_desc_cleaned, cAnalysis.p_desc, cAnalysis.p_AC50_cleaned, p_desc, p_aff, pr_RF_models, pr_TestSet)
#cApplyModel.PCACombine()
#cApplyModel.predict()
#cApplyModel.mergePrediction()
#cApplyModel.computeAD()
#cApplyModel.applyAD()


# SOM
#p_SOM = PR_RESULTS + "SOM/SOM_model.RData"
#cApplyModel.applySOM(p_SOM)



## 7. develop test set from Liu 2020 and Zhang 2016
## => https://doi.org/10.1039/C5TX00294J
## => https://doi.org/10.1016/j.snb.2019.127065
######################################
#P_CHEMBL = PR_DATA + "CHEMBL27-target_chembl240.csv"
#pr_ChEMBL_patch_clamp = pathFolder.createFolder(PR_RESULTS + "ChEMBL_patch_clamp/")

## 7.1 Prep dataset
#cChEMBL = CHEMBLTable.CHEMBLTable(P_CHEMBL, pr_ChEMBL_patch_clamp)
#cChEMBL.parseCHEMBLFile()
#cChEMBL.filterOnDescription("patch clamp")
#cChEMBL.cleanDataset(l_standard_type=["IC50"], l_standard_relation=["'='"])


## 7.2 run descriptor set
#p_desc_ChEMBL = cChEMBL.computeDesc()
#p_aff_ChEMBL = cChEMBL.prep_aff(typeAff="class", cutoff_uM=30.0)


## 7.3 run PCA with ChEMBL on PCA tox21
#pr_applyModel = pathFolder.createFolder(pr_ChEMBL_patch_clamp + "predict/")
#cApplyModel = applyModels.applyModel(cAnalysis.p_desc_cleaned, cAnalysis.p_desc, p_AC50, p_desc_ChEMBL, p_aff_ChEMBL, pr_RF_models, pr_applyModel)
#cApplyModel.overlapSet()
#cApplyModel.PCACombine()

## 7.4 apply RF models on it
#cApplyModel.predict()
#cApplyModel.mergePrediction()
#cApplyModel.computeAD()
#cApplyModel.applyAD()


## 7.5. Apply SOM on it
#p_SOM = PR_RESULTS + "SOM/SOM_model.RData"
#cApplyModel.applySOM(p_SOM)


## 8. develop test set from pubchem assay -> 588834
## => https://pubchem.ncbi.nlm.nih.gov/bioassay/588834
########################################################

#p_dataset = PR_DATA + "AID_588834_datatable_all.csv"
#p_mapp = PR_DATA + "AID_588834_ID_mapping.csv"
#p_SMILESComptox = PR_DATA + "CompToxChemicalsDashboard-SMILES.csv"
#pr_TestSet = pathFolder.createFolder(PR_RESULTS + "AID_588834/")
#c_genericSet = genericTestSet.genericTestSet(p_dataset, pr_TestSet)
#c_genericSet.loadDataset(loadDb=p_SMILESComptox, p_mapFile=p_mapp)
#p_desc = c_genericSet.computeDesc()
#p_aff = c_genericSet.setAff(allAff="PUBCHEM_ACTIVITY_OUTCOME")



## apply model
#cApplyModel = applyModels.applyModel(cAnalysis.p_desc_cleaned, cAnalysis.p_desc, cAnalysis.p_AC50, p_desc, p_aff, pr_RF_models, pr_TestSet)
#p_aff = cApplyModel.overlapSetWithID(rm_overlap=1)
#cApplyModel.PCACombine()
#cApplyModel.predict()
#cApplyModel.mergePrediction()
#cApplyModel.computeAD()
#cApplyModel.applyAD()

# SOM
#p_SOM = PR_RESULTS + "SOM/SOM_model.RData"
#cApplyModel.applySOM(p_SOM)


## 9. Shagun test set -> 588834
################################

p_dataset = PR_DATA + "Ext_Test_Set_Final_1014Chem.csv"
p_SMILESComptox = PR_DATA + "CompToxChemicalsDashboard-SMILES.csv"
pr_TestSet = pathFolder.createFolder(PR_RESULTS + "Ext_Test_Set_Final_1014Chem/")
c_genericSet = genericTestSet.genericTestSet(p_dataset, pr_TestSet)
c_genericSet.loadDataset(loadDb=p_SMILESComptox)
p_desc = c_genericSet.computeDesc()
p_aff = c_genericSet.setAff(allAff="PUBCHEM_ACTIVITY_OUTCOME")

## apply model
cApplyModel = applyModels.applyModel(cAnalysis.p_desc_cleaned, cAnalysis.p_desc, cAnalysis.p_AC50, p_desc, p_aff, pr_RF_models, pr_TestSet)
cApplyModel.PCACombine()
cApplyModel.predict()
cApplyModel.mergePrediction()
cApplyModel.computeAD()
cApplyModel.applyAD()

# SOM
p_SOM = PR_RESULTS + "SOM/SOM_model.RData"
cApplyModel.applySOM(p_SOM)
