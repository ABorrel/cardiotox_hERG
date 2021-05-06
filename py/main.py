import NCAST_assays
import CHEMBL_set
import NCAST_CHEMBL_set
import genericTestSet
import pathFolder
import applyModels
import hergml
import toolbox

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


# comparison with herg-ml for NCATS model
#########################
#pr_comparison = pathFolder.createFolder(cNCAST.pr_results + "pred_herg-ml/")
#c_genericSet = genericTestSet.genericTestSet(p_smi_NCAST, pr_comparison, PR_RESULTS)
#c_genericSet.main(allAff=cNCAST.cMain.c_dataset.d_dataset)
#p_desc = c_genericSet.computeDescForNNComparison("", pr_comparison) ## need to execute the py report


#p_pred_herg_ml = pr_comparison + "/chemicals_knime_desc_pred.csv"
#c_comparison = hergml.hergml(c_genericSet.d_dataset, p_pred_herg_ml, c_genericSet.p_aff, pr_comparison)
#c_comparison.computePerfAllSet()

####
# only with test NCAST #
##
#pr_comparison = pathFolder.createFolder(cNCAST.pr_results + "/pred_herg-ml-testset/")
#c_genericSet = genericTestSet.genericTestSet(pr_comparison + "test.csv", pr_comparison, PR_RESULTS)
#d_test = toolbox.loadMatrix(pr_comparison + "test.csv", sep = ",")
#c_genericSet.main(d_test)
#p_desc = c_genericSet.computeDescForNNComparison("", pr_comparison) ## need to execute the py report

#p_pred_herg_ml = pr_comparison + "/chemicals_knime_desc_pred.csv"
#c_comparison = hergml.hergml(c_genericSet.d_dataset, p_pred_herg_ml, c_genericSet.p_aff, pr_comparison)
#c_comparison.computePerfAllSet()




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

######
# EXTERNAL SET OF CHEMICALS
#########


## 6. use drug external test sest 1
############################

#p_dataset = PR_DATA + "8Candidate_Cipa_drugs.csv"
#pr_TestSet = pathFolder.createFolder(PR_RESULTS + "8Candidate_Cipa_drugs/")
#c_genericSet = genericTestSet.genericTestSet(p_dataset, pr_TestSet, PR_RESULTS)
#c_genericSet.main(allAff=1)

# NCAST
#cApplyModel = applyModels.applyModel(cNCAST.cMain.c_analysis.p_desc_cleaned, cNCAST.cMain.c_analysis.p_desc, cNCAST.cMain.c_analysis.p_AC50_cleaned, cNCAST.cMain.c_analysis.p_AC50, c_genericSet.p_desc, c_genericSet.p_aff, pr_TestSet, PR_RESULTS)
#cApplyModel.PCACombine()
#cApplyModel.computeAD()

#NCAST-ChEMBL

# prediction classif model
#cApplyModel.predict_AllClassifModels()



## 7. use drug external test sest 2
#############################
#p_dataset = PR_DATA + "hERG_Active_Annotated_listRefChem_July1.csv"
#pr_TestSet = pathFolder.createFolder(PR_RESULTS + "listRefChem_July1/")
#c_genericSet = genericTestSet.genericTestSet(p_dataset, pr_TestSet, PR_RESULTS)
#c_genericSet.main(allAff=1)

#cApplyModel = applyModels.applyModel(cNCAST.cMain.c_analysis.p_desc_cleaned, cNCAST.cMain.c_analysis.p_desc, cNCAST.cMain.c_analysis.p_AC50_cleaned, cNCAST.cMain.c_analysis.p_AC50, c_genericSet.p_desc, c_genericSet.p_aff, pr_TestSet, PR_RESULTS)
#cApplyModel.PCACombine()
#cApplyModel.computeAD()

# prediction classif model
#cApplyModel.predict_AllClassifModels()


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


#cApplyModel = applyModels.applyModel(cNCAST.cMain.c_analysis.p_desc_cleaned, cNCAST.cMain.c_analysis.p_desc, cNCAST.cMain.c_analysis.p_AC50_cleaned, cNCAST.cMain.c_analysis.p_AC50, c_genericSet.p_desc, c_genericSet.p_aff, pr_TestSet, PR_RESULTS)
#cApplyModel.PCACombine()
#cApplyModel.computeAD()
#cApplyModel.overlapSetWithID(rm_overlap=1)

# prediction classif model
#cApplyModel.predict_AllClassifModels()
#cApplyModel.predict_AllRegModels()


## 10. Shagun test set -> AID588834 no overlap => including 1014 chem
#######################################################################

#p_dataset = PR_DATA + "Ext_Test_Set_Final_1014Chem.csv"
#p_SMILESComptox = PR_DATA + "CompToxChemicalsDashboard-SMILES.csv"
#pr_TestSet = pathFolder.createFolder(PR_RESULTS + "AID588834_filtered/")

#c_genericSet = genericTestSet.genericTestSet(p_dataset, pr_TestSet, PR_RESULTS)
#c_genericSet.loadDataset(loadDb=p_SMILESComptox)
#c_genericSet.main(allAff="PUBCHEM_ACTIVITY_OUTCOME")

#cApplyModel = applyModels.applyModel(cNCAST.cMain.c_analysis.p_desc_cleaned, cNCAST.cMain.c_analysis.p_desc, cNCAST.cMain.c_analysis.p_AC50_cleaned, cNCAST.cMain.c_analysis.p_AC50, c_genericSet.p_desc, c_genericSet.p_aff, pr_TestSet, PR_RESULTS)
#cApplyModel.PCACombine()
#cApplyModel.computeAD()
#cApplyModel.overlapSetWithID(rm_overlap=1)

# prediction classif model
#cApplyModel.predict_AllClassifModels()
# concensus NCATS
#cApplyModel.concensus("NCAST_classif_undersampling", ["DNN", "RF", "LDA"])
#cApplyModel.mergePredictionClass()
# concensus NCATS + ChEMBL
#cApplyModel.concensus("NCAST_CHEMBL_classif_nosampling", ["DNN", "RF"])
#cApplyModel.mergePredictionClass()
#cApplyModel.predict_AllRegModels()

# comparison with herg-ml
#########################
#p_desc = c_genericSet.computeDescForNNComparison("", pr_TestSet) ## need to execute the py report

#pr_comparison = pathFolder.createFolder(pr_TestSet + "pred_herg-ml/")
#p_pred_herg_ml = pr_TestSet + "pred_herg-ml/chemicals_knime_desc_1014Chem_set_pred.csv"
#c_comparison = hergml.hergml(c_genericSet.d_dataset, p_pred_herg_ml, c_genericSet.p_aff, pr_comparison)
#c_comparison.computePerfAllSet()

## 11. Shagun test set => including 408 chem extracted from the litterature
#######################################

#p_dataset = PR_DATA + "Ext_Testset_hERG_review_14Oct.csv"
#p_SMILESComptox = PR_DATA + "CompToxChemicalsDashboard-SMILES.csv"
#pr_TestSet = pathFolder.createFolder(PR_RESULTS + "Ext_Test_Set_Final_Litt_408Chem/")

#c_genericSet = genericTestSet.genericTestSet(p_dataset, pr_TestSet, PR_RESULTS)
#c_genericSet.loadDataset(loadDb=p_SMILESComptox)
#c_genericSet.main(allAff=1)

#cApplyModel = applyModels.applyModel(cNCAST.cMain.c_analysis.p_desc_cleaned, cNCAST.cMain.c_analysis.p_desc, cNCAST.cMain.c_analysis.p_AC50_cleaned, cNCAST.cMain.c_analysis.p_AC50, c_genericSet.p_desc, c_genericSet.p_aff, pr_TestSet, PR_RESULTS)
#cApplyModel.PCACombine()
#cApplyModel.computeAD()

# prediction classif model
#cApplyModel.predict_AllClassifModels()
#cApplyModel.concensus("NCAST_classif_undersampling", ["DNN", "RF"])
#cApplyModel.mergePredictionClass()
#cApplyModel.concensus("NCAST_CHEMBL_classif_nosampling", ["DNN", "RF"])
#cApplyModel.mergePredictionClass()


# comparison with herg-ml
#########################
#p_desc = c_genericSet.computeDescForNNComparison("", pr_TestSet) ## need to execute the py report
#pr_comparison = pathFolder.createFolder(pr_TestSet + "pred_herg-ml/")
#p_pred_herg_ml = pr_TestSet + "pred_herg-ml/chemicals_knime_desc_408Chem_set_pred.csv"
#c_comparison = hergml.hergml(c_genericSet.d_dataset, p_pred_herg_ml, "1", pr_comparison)
#c_comparison.computePerfAllSet()



## 12. external test set from Johnattam llnl
#######################################

p_dataset = PR_DATA + "external_llnl/kcnh2_chembl_testset_base_smiles_union.csv"
p_SMILESComptox = PR_DATA + "CompToxChemicalsDashboard-SMILES.csv"
pr_TestSet = pathFolder.createFolder(PR_RESULTS + "Ext_Test_Set_llnl_chembl_test/")

c_genericSet = genericTestSet.genericTestSet(p_dataset, pr_TestSet, PR_RESULTS)
c_genericSet.loadDataset(loadDb=p_SMILESComptox)
c_genericSet.main(allAff="PUBCHEM_ACTIVITY_OUTCOME")

cApplyModel = applyModels.applyModel(cNCAST.cMain.c_analysis.p_desc_cleaned, cNCAST.cMain.c_analysis.p_desc, cNCAST.cMain.c_analysis.p_AC50_cleaned, cNCAST.cMain.c_analysis.p_AC50, c_genericSet.p_desc, c_genericSet.p_aff, pr_TestSet, PR_RESULTS)
cApplyModel.PCACombine()
cApplyModel.computeAD()

# prediction classif model
cApplyModel.predict_AllClassifModels()
cApplyModel.concensus("NCAST_classif_undersampling", ["DNN", "RF"])
cApplyModel.mergePredictionClass()
cApplyModel.concensus("NCAST_CHEMBL_classif_nosampling", ["DNN", "RF"])
cApplyModel.mergePredictionClass()
main169



