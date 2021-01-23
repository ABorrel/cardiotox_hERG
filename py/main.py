import analysis
import applyModels
import CHEMBLTable
import dataset
import genericTestSet
import pathFolder
import QSAR_modeling
import hergml

from os import path


class main:
    def __init__(self, pr_root, pr_data, name_dataset, pr_results):
        self.pr_root = pr_root
        self.pr_data = pr_data
        self.pr_results = pr_results
        self.name_dataset = name_dataset
        self.pr_desc = pathFolder.createFolder(pr_results + "DESC/")
    
    def setup_NCATS(self, p_NCATS_smi, p_NCAST_AC50, p_class_chemical):
        self.p_smi = p_NCATS_smi
        self.p_AC50 = p_NCAST_AC50
        self.pr_dataset = pathFolder.createFolder(self.pr_results + "dataset_" + self.name_dataset + "/")
        self.p_class_chemical = p_class_chemical

    def load_datasetNCAST(self):
        # 1.1 load dataset
        c_dataset = dataset.dataset(self.p_smi, self.p_AC50, self.pr_dataset)
        c_dataset.prep_dataset()

        self.c_dataset = c_dataset
    
    def classe_chemical_dataset(self):
        self.c_dataset.classChem(self.p_class_chemical)

    def compute_desc(self):
        l_p_desc = self.c_dataset.computeDesc(self.pr_desc)
        self.c_dataset.computePNG(self.pr_desc)
        self.l_p_desc = l_p_desc

    def rankActive(self):
        self.c_dataset.rankActiveChem(self.pr_desc + "PNG/")

    def setup_descAnalysis(self, cor_desc, quantile_desc):
        self.cor_desc = cor_desc
        self.quantile_desc = quantile_desc
        pr_out = pathFolder.createFolder(self.pr_results + "analysis_" + self.name_dataset + "/")
        self.c_analysis = analysis.analysis(self.p_AC50, self.l_p_desc[0], self.l_p_desc[1], pr_out, self.cor_desc, self.quantile_desc)
        self.c_analysis.combineDesc()
        self.c_analysis.prepDesc()
        self.c_analysis.sumAC50()
        self.c_analysis.signifDescActInact()

    def run_structural_PCA(self):
        self.c_analysis.PCA_plot()

    def run_structural_SOM(self, SOM_size):
        self.c_analysis.generate_SOM(SOM_size)
        self.c_analysis.signifDescBySOMCluster()
        self.c_analysis.extract_actBySOMCluster(self.pr_desc + "PNG/") # have to run !!!!
        self.c_analysis.applySOMtoChemClassification(self.p_class_chemical, self.pr_desc + "PNG/")

    def run_structural_Hclust(self):
        self.c_analysis.HClust_plot()


    def prep_QSARClass(self, cor_desc, quantile_desc, nb_repetition, n_foldCV, rate_split, rate_active):
        self.cor_desc = cor_desc
        self.quantile_desc = quantile_desc        
        self.nb_repetition = nb_repetition
        self.n_foldCV = n_foldCV
        self.rate_split = rate_split
        self.rate_active = rate_active

        pr_out = pathFolder.createFolder("%sQSAR_%s__%s-%s-%s-%s-%s-%s/"%(self.pr_results, self.name_dataset, self.cor_desc, self.quantile_desc, self.nb_repetition, self.n_foldCV, self.rate_split, self.rate_active))
        self.c_analysis = analysis.analysis(self.p_AC50, self.l_p_desc[0], self.l_p_desc[1], pr_out, self.cor_desc, self.quantile_desc)
        self.c_analysis.combineDesc()
        self.c_analysis.prepDesc()

    def run_QSARClassif(self, sampling):

        if sampling == "global":
            pr_out = pathFolder.createFolder(self.c_analysis.pr_out + "sampleGlobalSet/")
            self.cQSAR = QSAR_modeling.QSAR_modeling(self.c_analysis.p_desc_cleaned, self.l_p_desc[0], self.c_analysis.p_AC50_cleaned, self.p_AC50, pr_out, self.nb_repetition, self.n_foldCV, self.rate_active, self.rate_split)
            self.cQSAR.runQSARClassUnderSamplingAllSet(force_run=0)
            
            self.cQSAR.pr_RF_models = self.cQSAR.extractModels("%sRF_%s__globalsampling__%s-%s-%s-%s-%s-%s/"%(self.pr_results, self.name_dataset, self.cor_desc, self.quantile_desc, self.nb_repetition, self.n_foldCV, self.rate_split, self.rate_active), "RF")
            self.cQSAR.pr_LDA_models = self.cQSAR.extractModels("%sLDA_%s__globalsampling__%s-%s-%s-%s-%s-%s/"%(self.pr_results, self.name_dataset, self.cor_desc, self.quantile_desc, self.nb_repetition, self.n_foldCV, self.rate_split, self.rate_active), "LDA")

        else:
            pr_out = pathFolder.createFolder(self.c_analysis.pr_out + "sampleTraningSet/")
            self.cQSAR = QSAR_modeling.QSAR_modeling(self.c_analysis.p_desc_cleaned, self.l_p_desc[0], self.c_analysis.p_AC50_cleaned, self.p_AC50, pr_out, self.nb_repetition, self.n_foldCV, self.rate_active, self.rate_split)
            self.cQSAR.runQSARClassUnderSamplingTrain()
            
            self.cQSAR.pr_RF_models = self.cQSAR.extractModels("%sRF_%s__trainnigsampling__%s-%s-%s-%s-%s-%s/"%(self.pr_results, self.name_dataset, self.cor_desc, self.quantile_desc, self.nb_repetition, self.n_foldCV, self.rate_split, self.rate_active), "RF")
            self.cQSAR.pr_LDA_models = self.cQSAR.extractModels("%sLDA_%s__trainnigsampling__%s-%s-%s-%s-%s-%s/"%(self.pr_results, self.name_dataset, self.cor_desc, self.quantile_desc, self.nb_repetition, self.n_foldCV, self.rate_split, self.rate_active), "LDA")

class NCAST_assays:
    def __init__(self, p_smi, p_aff_origin, p_classification,  pr_root, pr_data, pr_results):
        self.p_smi = p_smi
        self.p_aff_origin = p_aff_origin
        self.p_classification = p_classification
        
        self.pr_root = pr_root
        self.pr_data = pr_data
        self.pr_results = pr_results

    def prep_dataset(self):
        self.cMain = main(self.pr_root, self.pr_data, "NCAST", self.pr_results)
        self.cMain.setup_NCATS(self.p_smi, self.p_aff_origin, self.p_classification)
        self.cMain.load_datasetNCAST()
        self.cMain.classe_chemical_dataset()
    
        #compute desc
        self.cMain.compute_desc()
        #rank active
        self.cMain.rankActive()

    def analysis_nopred(self, cor_val, max_quantile, som_size):
        self.cor_val = cor_val
        self.max_quantile = max_quantile
        self.som_size = som_size
        
        self.cMain.setup_descAnalysis(self.cor_val, self.max_quantile)
        self.cMain.run_structural_PCA()
        self.cMain.run_structural_Hclust()
        self.cMain.run_structural_SOM(self.som_size)


    def QSARClassif_builder(self, cor_desc, quantile_desc, nb_repetition, n_foldCV, rate_split, rate_active):

        self.cMain.prep_QSARClass(cor_desc, quantile_desc, nb_repetition, n_foldCV, rate_split, rate_active)
        self.cMain.run_QSARClassif("training")


    def comparison_hergml(self, p_pred_herg_ml, p_test_set):
        
        pr_comparison_hergml = pathFolder.createFolder(self.pr_results + "comparison_hergml_" + self.name_dataset + "/")
        self.cMain.c_dataset.computeDescForNNComparison(self.pr_desc, pr_comparison_hergml)

        self.c_comparison_hergml = hergml.hergml(self.p_smi, p_pred_herg_ml, self.p_aff_origin, pr_comparison_hergml)
        self.c_comparison_hergml.computePerfAllSet()
        self.c_comparison_hergml.computePerfTestSet(p_test_set)

class CHEMBL_set:
    def __init__(self, p_chembl_results, name_dataset, pr_results):
        self.p_chembl_results = p_chembl_results
        self.name_dataset = name_dataset
        self.pr_results = pr_results
    
    def prep_dataset(self, l_standard_type, l_standard_relation):
        pr_dataset = pathFolder.createFolder(self.pr_results + "dataset_" + self.name_dataset + "/")
        self.pr_dataset = pr_dataset
        self.cCHEMBL = CHEMBLTable.CHEMBLTable(self.p_chembl_results, pr_dataset)
        self.cCHEMBL.parseCHEMBLFile()
        self.cCHEMBL.cleanDataset(l_standard_type=["IC50", "Ki"], l_standard_relation=["'='"])
    
    def correlation_dataset(self, p_dataset):
        pr_comparison = pathFolder.createFolder(self.pr_results + "comparison_CHEMBL_NCAST/")
        self.cCHEMBL.correlation_aff(p_dataset, pr_comparison, self.pr_results)


    
#P_CHEMBL = PR_DATA + "CHEMBL27-target_chembl240.csv"
#pr_ChEMBL = pathFolder.createFolder(PR_RESULTS + "ChEMBL/")

## 5.1 Prep dataset
#cChEMBL = CHEMBLTable.CHEMBLTable(P_CHEMBL, pr_ChEMBL)
#cChEMBL.parseCHEMBLFile()
#cChEMBL.cleanDataset(l_standard_type=["IC50", "Ki"], l_standard_relation=["'='"])

## 5.2 run descriptor set
#p_desc_ChEMBL = cChEMBL.computeDesc()
#p_aff_ChEMBL = cChEMBL.cleanAff()



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

#cNCAST = NCAST_assays(p_smi_NCAST, p_AC50_NCAST, p_classification_interpred,  PR_ROOT, PR_DATA, PR_RESULTS)
#cNCAST.prep_dataset()
#cNCAST.analysis_nopred(COR_VAL, MAX_QUANTILE, SOM_size)

#QSAR
nb_repetition = 10
n_foldCV = 10
rate_split = 0.15
rate_active = 0.30
#cNCAST.QSARClassif_builder(COR_VAL, MAX_QUANTILE, nb_repetition, n_foldCV, rate_split, rate_active)

# comparison with herg-ml models
p_pred_herg_ml = PR_ROOT + "comparison_study/pred/list_chemicals-2020-06-09-09-11-41_pred.csv"
p_test_set = PR_RESULTS + "QSARclass/1/test.csv"
#cNCAST.comparison_hergml(p_pred_herg_ml, p_test_set)



############
# BUILD CHEMBL DATASET
########################
############################
p_CHEMBL_raw = PR_DATA + "chembl27/CHEMBL240.csv"
l_standard_type=["IC50", "Ki", "EC50"]
l_standard_relation=["'='"]

c_CHEMBL_set = CHEMBL_set(p_CHEMBL_raw, "CHEMBL27", PR_RESULTS)
c_CHEMBL_set.prep_dataset(l_standard_type, l_standard_relation)

p_NCAST_dataset = PR_RESULTS + "dataset_NCAST/dataset_prep.csv"
c_CHEMBL_set.correlation_dataset(p_NCAST_dataset)
sss










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