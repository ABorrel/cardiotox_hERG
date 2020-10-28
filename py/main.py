import analysis
import applyModels
import CHEMBLTable
import dataset
import genericTestSet
import pathFolder
import QSAR_modeling

from os import path
# Define folder
################
PR_ROOT = path.abspath("../../") + "/"
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
l_p_desc = c_dataset.computeDesc(pr_desc)
#c_dataset.computePNG(pr_desc)

# 1.4 rank active chemical and save in a PDF
#c_dataset.rankActiveChem(pr_desc + "PNG/")


# 2. analysis => using only rdkit descriptors
#####

COR_VAL = 0.90
MAX_QUANTILE = 90

cAnalysis = analysis.analysis(p_AC50, c_dataset.p_desc1D2D, c_dataset.p_desc_opera, PR_RESULTS, COR_VAL, MAX_QUANTILE)
#cAnalysis.prepDesc() #only consider 1D2D descriptors

## 2.1. histogram AC50 and summary
#cAnalysis.sumAC50()

## 2.2 PCA
#cAnalysis.PCA_plot()

## 2.3 SOM
#size = 15
#cAnalysis.generate_SOM(15)
#cAnalysis.signifDescBySOMCluster()
#cAnalysis.extract_actBySOMCluster(pr_desc + "PNG/") # have to run !!!!
#p_classification = PR_DATA + "class_from_interferences.csv"
#cAnalysis.applySOMtoChemClassification(p_classification, pr_desc + "PNG/")

#p_classification = PR_DATA + "hERG_Active_Annotated_listRefChem_July1.csv"
#cAnalysis.applySOMtoChemClassification(p_classification, pr_desc + "PNG/", filter_active=0)

## 2.4 Hclust
#cAnalysis.HClust_plot()

## 2.5 desc Significance active vs not active
#cAnalysis.signifDescActInact()


## 3. QSAR modeling - classification
#######
pr_QSAR = pathFolder.createFolder(PR_RESULTS + "QSARclass/")
nb_repetition = 10
n_foldCV = 10
rate_split = 0.15
rate_active = 0.30
# => Prep and combine OPERA and rdkit desc
cAnalysis.pr_out = pr_QSAR # need to redefine the output directory here
cAnalysis.combineDesc()
cAnalysis.prepDesc()


# 3.1 using a new defintion of test set at each iteration
########################
cQSAR = QSAR_modeling.QSAR_modeling(cAnalysis.p_desc_cleaned, c_dataset.p_desc1D2D, cAnalysis.p_AC50_cleaned, p_AC50, pr_QSAR, nb_repetition, n_foldCV,rate_active, rate_split)
#cQSAR.runQSARClassUnderSamplingAllSet(force_run=0)
pr_RF_models = cQSAR.extractModels(PR_RESULTS, "RF")
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


# 4 QSAR modeling regression
#############
pr_QSARreg = pathFolder.createFolder(PR_RESULTS + "QSARreg/")
n_foldCV = 10
rate_split = 0.15
nb_repetition = 10

# 4.1 prep data
#################
# => Prep and combine OPERA and rdkit desc
#COR_VAL = 0.85
#MAX_QUANTILE = 85
#cAnalysis = analysis.analysis(p_AC50, c_dataset.p_desc1D2D, c_dataset.p_desc_opera, PR_RESULTS, COR_VAL, MAX_QUANTILE)
#cAnalysis.prepDesc()
#cAnalysis.combineDesc()

# 4.2 run reg model
####################
#cQSARreg = QSAR_modeling.QSAR_modeling(cAnalysis.p_desc, p_desc, p_AC50, p_AC50, pr_QSARreg, nb_repetition, n_foldCV, 0, rate_split)
#cQSARreg.runQSARReg(cAnalysis.cor_val, cAnalysis.max_quantile)
#cQSARreg.mergeRegResults()


# 5. Build external test set from ChEMBL
### Target => CHEMBL240
##########################
#P_CHEMBL = PR_DATA + "CHEMBL27-target_chembl240.csv"
#pr_ChEMBL = pathFolder.createFolder(PR_RESULTS + "ChEMBL/")

## 5.1 Prep dataset
#cChEMBL = CHEMBLTable.CHEMBLTable(P_CHEMBL, pr_ChEMBL)
#cChEMBL.parseCHEMBLFile()
#cChEMBL.cleanDataset(l_standard_type=["IC50", "Ki"], l_standard_relation=["'='"])

## 5.2 run descriptor set
#p_desc_ChEMBL = cChEMBL.computeDesc()
#p_aff_ChEMBL = cChEMBL.cleanAff()

## 5.3 run PCA with ChEMBL on PCA tox21
#pr_applyModel = pathFolder.createFolder(PR_RESULTS + "ChEMBL_predict/")
#cApplyModel = applyModels.applyModel(cAnalysis.p_desc_cleaned, cAnalysis.p_AC50_cleaned, p_desc_ChEMBL, p_aff_ChEMBL, pr_RF_models, pr_applyModel)
#cApplyModel.PCACombine()

## 5.4 apply RF models on it
#cApplyModel.predict()
#cApplyModel.mergePrediction()

## 5.5. Apply SOM on it
#p_SOM = PR_RESULTS + "SOM/SOM_model.RData"
#cApplyModel.applySOM(p_SOM)


## 6. use drug external test sest 1
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

## 7. use drug external test sest 2
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



## 8. develop test set from Liu 2020 and Zhang 2016
## => https://doi.org/10.1039/C5TX00294J
## => https://doi.org/10.1016/j.snb.2019.127065
######################################
#P_CHEMBL = PR_DATA + "CHEMBL27-target_chembl240.csv"
#pr_ChEMBL_patch_clamp = pathFolder.createFolder(PR_RESULTS + "ChEMBL_patch_clamp/")

## 8.1 Prep dataset
#cChEMBL = CHEMBLTable.CHEMBLTable(P_CHEMBL, pr_ChEMBL_patch_clamp)
#cChEMBL.parseCHEMBLFile()
#cChEMBL.filterOnDescription("patch clamp")
#cChEMBL.cleanDataset(l_standard_type=["IC50"], l_standard_relation=["'='"])


## 8.2 run descriptor set
#p_desc_ChEMBL = cChEMBL.computeDesc()
#p_aff_ChEMBL = cChEMBL.prep_aff(typeAff="class", cutoff_uM=30.0)


## 8.3 run PCA with ChEMBL on PCA tox21
#pr_applyModel = pathFolder.createFolder(pr_ChEMBL_patch_clamp + "predict/")
#cApplyModel = applyModels.applyModel(cAnalysis.p_desc_cleaned, cAnalysis.p_desc, p_AC50, p_desc_ChEMBL, p_aff_ChEMBL, pr_RF_models, pr_applyModel)
#cApplyModel.overlapSet()
#cApplyModel.PCACombine()

## 8.4 apply RF models on it
#cApplyModel.predict()
#cApplyModel.mergePrediction()
#cApplyModel.computeAD()
#cApplyModel.applyAD()


## 8.5. Apply SOM on it
#p_SOM = PR_RESULTS + "SOM/SOM_model.RData"
#cApplyModel.applySOM(p_SOM)


## 9. develop test set from pubchem assay -> 588834
## => https://pubchem.ncbi.nlm.nih.gov/bioassay/588834
########################################################

# path for file for dataset
#p_dataset = PR_DATA + "AID_588834_datatable_all.csv"
#p_mapp = PR_DATA + "AID_588834_ID_mapping.csv"
#p_SMILESComptox = PR_DATA + "CompToxChemicalsDashboard-SMILES.csv"
#pr_TestSet = pathFolder.createFolder(PR_RESULTS + "AID_588834/")

# data from the model
#pr_model_data = pathFolder.createFolder(pr_TestSet + "Data_model/")
#cDataModel = analysis.analysis(p_AC50, c_dataset.p_desc1D2D, c_dataset.p_desc_opera, pr_model_data, COR_VAL, MAX_QUANTILE)
#cDataModel.combineDesc()
#cDataModel.prepDesc()

# data from test set processing
#c_genericSet = genericTestSet.genericTestSet(p_dataset, pr_TestSet)
#c_genericSet.loadDataset(loadDb=p_SMILESComptox, p_mapFile=p_mapp)
#c_genericSet.computeDesc()
## combine OPERA + RDKIT
#c_genericSet.combineDesc()
#c_genericSet.setAff(allAff="PUBCHEM_ACTIVITY_OUTCOME")


## apply classification models ##
#################################
#pr_QSAR_classif = pathFolder.createFolder(pr_TestSet + "QSARClass/")
#cApplyModel = applyModels.applyModel(cDataModel.p_desc_cleaned, cDataModel.p_desc, cDataModel.p_AC50_cleaned, c_genericSet.p_desc, c_genericSet.p_aff, pr_RF_models, pr_QSAR_classif)
#cApplyModel.overlapSetWithID(rm_overlap=1)
#cApplyModel.PCACombine()
#cApplyModel.predict()
#cApplyModel.mergePrediction()
#cApplyModel.computeAD()
#cApplyModel.applyAD()


## apply regression model ##
############################

#p_model_Reg = PR_RESULTS + "QSARreg/7/RFreg/model.RData" # manual selection
#pr_QSAR_reg = pathFolder.createFolder(pr_TestSet + "QSARReg/")
#cApplyModel = applyModels.applyModel(cDataModel.p_desc_cleaned, cDataModel.p_desc, cDataModel.p_AC50_cleaned, c_genericSet.p_desc, c_genericSet.p_aff, p_model_Reg, pr_QSAR_reg)
#cApplyModel.overlapSetWithID(rm_overlap=1, rm_inactive=1)
#cDataModel.extractOnlyActive()
#cApplyModel.computeAD()
#cApplyModel.PCACombine()
#cApplyModel.predictReg()
#cApplyModel.applyAD()



## 10. Shagun test set -> AID588834 no overlap => including 1014 chem
#######################################################################

#p_dataset = PR_DATA + "Ext_Test_Set_Final_1014Chem.csv"
#p_SMILESComptox = PR_DATA + "CompToxChemicalsDashboard-SMILES.csv"
#pr_TestSet = pathFolder.createFolder(PR_RESULTS + "Ext_Test_Set_Final_1014Chem/")
#c_genericSet = genericTestSet.genericTestSet(p_dataset, pr_TestSet)
#c_genericSet.loadDataset(loadDb=p_SMILESComptox)
#p_desc = c_genericSet.computeDesc()
#c_genericSet.combineDesc()
#c_genericSet.setAff(allAff="PUBCHEM_ACTIVITY_OUTCOME")

## 10.1 apply model classification
#####
#cApplyModel = applyModels.applyModel(cAnalysis.p_desc_cleaned, cAnalysis.p_desc, cAnalysis.p_AC50_cleaned, c_genericSet.p_desc, c_genericSet.p_aff, pr_RF_models, pr_TestSet)
#cApplyModel.PCACombine()
#cApplyModel.predict()
#cApplyModel.mergePrediction()
#cApplyModel.computeAD()
#cApplyModel.applyAD()


##10.2 apply model regression
#for i in range(1,11):
#    p_model_Reg = "%sQSARreg/%s/RFreg/model.RData"%(PR_RESULTS, i) # manual selection
#    pr_QSAR_reg = pathFolder.createFolder("%sQSAR-RF-Reg%s/"%(pr_TestSet, i))
#    cApplyModel = applyModels.applyModel(cAnalysis.p_desc_cleaned, cAnalysis.p_desc, cAnalysis.p_AC50_cleaned, c_genericSet.p_desc, c_genericSet.p_aff, p_model_Reg, pr_QSAR_reg)
#    cApplyModel.overlapSetWithID(rm_overlap=1, rm_inactive=1)
#    cApplyModel.computeAD()
#    cApplyModel.predictReg()

#for i in range(1,11):
#    p_model_Reg = "%sQSARreg/%s/NNreg/model.RData"%(PR_RESULTS, i) # manual selection
#    pr_QSAR_reg = pathFolder.createFolder("%sQSAR-NN-Reg%s/"%(pr_TestSet, i))
#    cApplyModel = applyModels.applyModel(cAnalysis.p_desc_cleaned, cAnalysis.p_desc, cAnalysis.p_AC50_cleaned, c_genericSet.p_desc, c_genericSet.p_aff, p_model_Reg, pr_QSAR_reg)
#    cApplyModel.overlapSetWithID(rm_overlap=1, rm_inactive=1)
#    cApplyModel.computeAD()
#    cApplyModel.predictReg("NN")

#for i in range(1,11):
#    p_model_Reg = "%sQSARreg/%s/PCRreg/model.RData"%(PR_RESULTS, i) # manual selection
#    pr_QSAR_reg = pathFolder.createFolder("%sQSAR-PCR-Reg%s/"%(pr_TestSet, i))
#    cApplyModel = applyModels.applyModel(cAnalysis.p_desc_cleaned, cAnalysis.p_desc, cAnalysis.p_AC50_cleaned, c_genericSet.p_desc, c_genericSet.p_aff, p_model_Reg, pr_QSAR_reg)
#    cApplyModel.overlapSetWithID(rm_overlap=1, rm_inactive=1)
#    cApplyModel.computeAD()
#    cApplyModel.predictReg("PCR")

#for i in range(1,11):
#    p_model_Reg = "%sQSARreg/%s/PLSreg/model.RData"%(PR_RESULTS, i) # manual selection
#    pr_QSAR_reg = pathFolder.createFolder("%sQSAR-PLS-Reg%s/"%(pr_TestSet, i))
#    cApplyModel = applyModels.applyModel(cAnalysis.p_desc_cleaned, cAnalysis.p_desc, cAnalysis.p_AC50_cleaned, c_genericSet.p_desc, c_genericSet.p_aff, p_model_Reg, pr_QSAR_reg)
#    cApplyModel.overlapSetWithID(rm_overlap=1, rm_inactive=1)
#    cApplyModel.computeAD()
#    cApplyModel.predictReg("PLS")



## 11. Shagun test set => including 408 chem extracted from the litterature
#######################################

p_dataset = PR_DATA + "Ext_Testset_hERG_review_14Oct.csv"
p_SMILESComptox = PR_DATA + "CompToxChemicalsDashboard-SMILES.csv"
pr_TestSet = pathFolder.createFolder(PR_RESULTS + "Ext_Test_Set_Final_Litt_408Chem/")
c_genericSet = genericTestSet.genericTestSet(p_dataset, pr_TestSet)
c_genericSet.loadDataset()
c_genericSet.computeDesc()
c_genericSet.combineDesc()
c_genericSet.setAff(allAff=1)

## apply model
cApplyModel = applyModels.applyModel(cAnalysis.p_desc_cleaned, cAnalysis.p_desc, cAnalysis.p_AC50_cleaned, c_genericSet.p_desc, c_genericSet.p_aff, pr_RF_models, pr_TestSet)
cApplyModel.PCACombine()
cApplyModel.predict()
cApplyModel.mergePrediction()
cApplyModel.computeAD()
cApplyModel.applyAD()
ss