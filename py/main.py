import NCAST_assays
import CHEMBL_set
import NCAST_CHEMBL_set
import genericTestSet
import pathFolder
import applyModels

from os import path


# Define folder
################
PR_ROOT = path.abspath("../../") + "/"
PR_DATA = PR_ROOT + "data/"
PR_RESULTS = pathFolder.createFolder(PR_ROOT + "results/")


##########
# NCAST DATASET
##################
######################
p_smi_NCAST = PR_DATA + "list_chemicals-2020-06-09-09-11-41.csv"
p_AC50_NCAST = PR_DATA + "AC50_7403.txt"
# chemical class
p_classification_interpred = PR_DATA + "class_from_interferences.csv"
p_classification_shagun = PR_DATA + "hERG_Active_Annotated_listRefChem_July1.csv"
# parameter
COR_VAL = 0.90
MAX_QUANTILE = 90
SOM_size = 15
cutoff_aff = 30 #uM

cNCAST = NCAST_assays.NCAST_assays(p_smi_NCAST, p_AC50_NCAST, p_classification_interpred, cutoff_aff, PR_ROOT, PR_DATA, PR_RESULTS)
cNCAST.main()

############
# BUILD CHEMBL DATASET
########################
############################
p_CHEMBL_raw = PR_DATA + "chembl27/CHEMBL240.csv"

c_CHEMBL_set = CHEMBL_set.CHEMBL_set(p_CHEMBL_raw, "CHEMBL27", PR_RESULTS)
c_CHEMBL_set.main()

####
# Build analysis on merged sets
#############
rate_active = 0

c_NCAST_CHEMBL_set = NCAST_CHEMBL_set.NCAST_CHEMBL_set(c_CHEMBL_set.cCHEMBL.p_merge_sets, c_CHEMBL_set.pr_merge_sets, PR_RESULTS)
c_NCAST_CHEMBL_set.main(c_CHEMBL_set, cNCAST)
sss

######
# EXTERNAL SET OF CHEMICALS
#########


## 6. use drug external test sest 1
############################

#p_dataset = PR_DATA + "8Candidate_Cipa_drugs.csv"
#pr_TestSet = pathFolder.createFolder(PR_RESULTS + "8Candidate_Cipa_drugs/")
#c_genericSet = genericTestSet.genericTestSet(p_dataset, pr_TestSet, PR_RESULTS)
#c_genericSet.main(allAff=1)


#cApplyModel = applyModels.applyModel(cNCAST.cMain.c_analysis.p_desc_cleaned, cNCAST.cMain.c_analysis.p_desc, cNCAST.cMain.c_analysis.p_AC50_cleaned, c_genericSet.p_desc, c_genericSet.p_aff, pr_TestSet)
#cApplyModel.PCACombine()
#cApplyModel.computeAD()

# prediction classif model
#cApplyModel.predict_all_classif()


## 7. use drug external test sest 2
#############################
#p_dataset = PR_DATA + "hERG_Active_Annotated_listRefChem_July1.csv"
#pr_TestSet = pathFolder.createFolder(PR_RESULTS + "listRefChem_July1/")
#c_genericSet = genericTestSet.genericTestSet(p_dataset, pr_TestSet, PR_RESULTS)
#c_genericSet.main(allAff=1)

#cApplyModel = applyModels.applyModel(cNCAST.cMain.c_analysis.p_desc_cleaned, cNCAST.cMain.c_analysis.p_desc, cNCAST.cMain.c_analysis.p_AC50_cleaned, c_genericSet.p_desc, c_genericSet.p_aff, pr_TestSet, PR_RESULTS)
#cApplyModel.PCACombine()
#cApplyModel.computeAD()

# prediction classif model
#cApplyModel.predict_all_classif()



## 9. develop test set from pubchem assay -> 588834
## => https://pubchem.ncbi.nlm.nih.gov/bioassay/588834
########################################################

# path for file for dataset
#p_dataset = PR_DATA + "AID_588834_datatable_all.csv"
#p_mapp = PR_DATA + "AID_588834_ID_mapping.csv"
#p_SMILESComptox = PR_DATA + "CompToxChemicalsDashboard-SMILES.csv"
#pr_TestSet = pathFolder.createFolder(PR_RESULTS + "AID_588834_all/")

# data from test set processing
#c_genericSet = genericTestSet.genericTestSet(p_dataset, pr_TestSet, PR_RESULTS)
#c_genericSet.loadDataset(loadDb=p_SMILESComptox, p_mapFile=p_mapp)
#c_genericSet.main(allAff="PUBCHEM_ACTIVITY_OUTCOME")


#cApplyModel = applyModels.applyModel(cNCAST.cMain.c_analysis.p_desc_cleaned, cNCAST.cMain.c_analysis.p_desc, cNCAST.cMain.c_analysis.p_AC50_cleaned, c_genericSet.p_desc, c_genericSet.p_aff, pr_TestSet, PR_RESULTS)
#cApplyModel.PCACombine()
#cApplyModel.computeAD()
#cApplyModel.overlapSetWithID(rm_overlap=1)

# prediction classif model
#cApplyModel.predict_all_classif()


## 10. Shagun test set -> AID588834 no overlap => including 1014 chem
#######################################################################

p_dataset = PR_DATA + "Ext_Test_Set_Final_1014Chem.csv"
p_SMILESComptox = PR_DATA + "CompToxChemicalsDashboard-SMILES.csv"
pr_TestSet = pathFolder.createFolder(PR_RESULTS + "Ext_Test_Set_Final_1014Chem/")

c_genericSet = genericTestSet.genericTestSet(p_dataset, pr_TestSet, PR_RESULTS)
c_genericSet.loadDataset(loadDb=p_SMILESComptox)
c_genericSet.main(allAff="PUBCHEM_ACTIVITY_OUTCOME")

cApplyModel = applyModels.applyModel(cNCAST.cMain.c_analysis.p_desc_cleaned, cNCAST.cMain.c_analysis.p_desc, cNCAST.cMain.c_analysis.p_AC50_cleaned, c_genericSet.p_desc, c_genericSet.p_aff, pr_TestSet, PR_RESULTS)
cApplyModel.PCACombine()
cApplyModel.computeAD()
cApplyModel.overlapSetWithID(rm_overlap=1)

# prediction classif model
#cApplyModel.predict_all_classif()

ss


## 11. Shagun test set => including 408 chem extracted from the litterature
#######################################

p_dataset = PR_DATA + "Ext_Testset_hERG_review_14Oct.csv"
p_SMILESComptox = PR_DATA + "CompToxChemicalsDashboard-SMILES.csv"
pr_TestSet = pathFolder.createFolder(PR_RESULTS + "Ext_Test_Set_Final_Litt_408Chem/")

#c_genericSet = genericTestSet.genericTestSet(p_dataset, pr_TestSet, PR_RESULTS)
#c_genericSet.loadDataset(loadDb=p_SMILESComptox)
#c_genericSet.main(allAff="PUBCHEM_ACTIVITY_OUTCOME")

#cApplyModel = applyModels.applyModel(cNCAST.cMain.c_analysis.p_desc_cleaned, cNCAST.cMain.c_analysis.p_desc, cNCAST.cMain.c_analysis.p_AC50_cleaned, c_genericSet.p_desc, c_genericSet.p_aff, pr_TestSet, PR_RESULTS)
#cApplyModel.PCACombine()
#cApplyModel.computeAD()
#cApplyModel.overlapSetWithID(rm_overlap=1)

# prediction classif model
#cApplyModel.predict_all_classif()























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

##10.3 comparison with herg-ml
##################
#p_desc = c_genericSet.computeDescForNNComparison("", pr_TestSet) ## need to execute the py report
#p_pred_herg_ml = pr_TestSet + "pred_herg-ml/chemicals_knime_desc_1014Chem_set_pred.csv"
#pr_comparison = pr_TestSet + "pred_herg-ml/"
#c_comparison = hergml.hergml(c_genericSet.d_dataset, p_pred_herg_ml, c_genericSet.p_aff, pr_comparison)
#c_comparison.computePerfAllSet()


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

##11.2 apply model
##################
#cApplyModel = applyModels.applyModel(cAnalysis.p_desc_cleaned, cAnalysis.p_desc, cAnalysis.p_AC50_cleaned, c_genericSet.p_desc, c_genericSet.p_aff, pr_RF_models, pr_TestSet)
#cApplyModel.PCACombine()
#cApplyModel.predict()
#cApplyModel.mergePrediction()
#cApplyModel.computeAD()
#cApplyModel.applyAD()


##11.3 comparison with herg-ml
##################
#p_desc = c_genericSet.computeDescForNNComparison("", pr_TestSet) ## need to execute the py report

p_pred_herg_ml = pr_TestSet + "pred_herg-ml/chemicals_knime_desc_408Chem_set_pred.csv"
pr_comparison = pr_TestSet + "pred_herg-ml/"
c_comparison = hergml.hergml(c_genericSet.d_dataset, p_pred_herg_ml, "1", pr_comparison)
c_comparison.computePerfAllSet()