import pathFolder
import CHEMBLTable
import analysis
import generic_run

class CHEMBL_set:
    def __init__(self, p_chembl_results, name_dataset, pr_results, pr_root, pr_data):
        self.p_chembl_results = p_chembl_results
        self.name_dataset = name_dataset
        self.pr_results = pr_results
        self.pr_root = pr_root
        self.pr_data = pr_data
        self.pr_out = pathFolder.createFolder(pr_results + name_dataset + "/")
        self.pr_desc = pathFolder.createFolder(pr_results + "DESC/")


    def main(self):
        l_standard_type=["IC50", "Ki", "EC50"]
        l_standard_relation=["'='"]

        self.prep_dataset(l_standard_type, l_standard_relation)
        self.prep_descset()

        ######
        # create model with the chembl set
        ##########
        # general val
        COR_VAL = 0.90
        MAX_QUANTILE = 90
        SOM_size = 15
        cutoff_aff = 30 #uM
        
        # QSAR values
        nb_repetition = 5
        n_foldCV = 10
        rate_split = 0.15

        self.analysis_nopred(COR_VAL, MAX_QUANTILE, SOM_size)

        ###
        # DO QSAR
        rate_active = 0# no undersampling
        self.QSARClassif_builder(COR_VAL, MAX_QUANTILE, nb_repetition, n_foldCV, rate_split, rate_active)



        ######
        # MERGE dataset CHEMBL + NCAST
        ##############
        p_NCAST_aff = self.pr_results + "NCAST/dataset_NCAST/dataset_prep.csv"
        p_NCAST_desc = self.pr_results + "NCAST/dataset_NCAST/desc_1D2D.csv"
        self.merge_dataset(p_NCAST_aff, p_NCAST_desc, "NCAST_CHEMBL")
        self.correlation_dataset(p_NCAST_aff)


    def prep_dataset(self, l_standard_type, l_standard_relation):
        pr_dataset = pathFolder.createFolder(self.pr_out + "dataset_" + self.name_dataset + "/")
        self.pr_dataset = pr_dataset
        self.cCHEMBL = CHEMBLTable.CHEMBLTable(self.p_chembl_results, pr_dataset)
        self.cCHEMBL.parseCHEMBLFile()
        self.cCHEMBL.cleanDataset(l_standard_type=l_standard_type, l_standard_relation=l_standard_relation, assay_type_favorised = "patch clamp")
    
    def prep_descset(self):

        self.cCHEMBL.computeDesc(self.pr_desc)
        self.cCHEMBL.computeOPERADesc(self.pr_desc)
        p_aff_ChEMBL = self.cCHEMBL.prep_aff(typeAff="", cutoff_uM=1)

    def correlation_dataset(self, p_dataset):
        pr_comparison = pathFolder.createFolder(self.pr_results + "comparison_CHEMBL_NCAST/")
        self.cCHEMBL.correlation_aff(p_dataset, pr_comparison, self.pr_results)

    def merge_dataset(self, p_aff_clean_toadd, p_desc1D2D_toadd, name_dataset):
        pr_out = pathFolder.createFolder(self.pr_results + name_dataset + "/")
        self.cCHEMBL.mergeDataset(p_aff_clean_toadd, p_desc1D2D_toadd, pr_out)

        self.pr_merge_sets = pr_out

    def analysis_nopred(self, cor_val, max_quantile, som_size):
        self.cor_val = cor_val
        self.max_quantile = max_quantile
        self.som_size = som_size
        
        self.cMain = generic_run.generic_run(self.pr_root, self.pr_data, "ChEMBL_27", self.pr_out)
        self.cMain.setup_ChEMBL(self.cCHEMBL.p_desc1D2D, self.cCHEMBL.p_descOPERA, self.cCHEMBL.p_aff_clean, self.pr_dataset)

        self.cMain.setup_descAnalysis(self.cor_val, self.max_quantile)
        self.cMain.run_structural_PCA()
        self.cMain.run_structural_Hclust()
        self.cMain.run_structural_SOM(self.som_size)

    def QSARClassif_builder(self, cor_desc, quantile_desc, nb_repetition, n_foldCV, rate_split, rate_active):

        self.cMain.prep_QSARClass(cor_desc, quantile_desc, nb_repetition, n_foldCV, rate_split, rate_active)
        self.cMain.run_QSARClassif("training")

# aff from ncast class

    def QSARReg_builder(self, cor_desc, quantile_desc, nb_repetition, n_foldCV, rate_split, rate_active):

        self.cMain.prep_QSARReg(cor_desc, quantile_desc, nb_repetition, n_foldCV, rate_split, rate_active)
        self.cMain.run_QSARReg()

    def DNN_builder(self, cor_desc, quantile_desc, nb_repetition, n_foldCV, rate_split, rate_active):

        # prep as classiq ML
        self.cMain.prep_QSARClass(cor_desc, quantile_desc, nb_repetition, n_foldCV, rate_split, rate_active)

        # case of no undersampleing
        if rate_split == 0:
            for i in range(0,nb_repetition):
                pr_DNN = pathFolder.createFolder(self.cMain.c_analysis.pr_out + "DNN_model/" + str(i+1) + "/")
                if path.exists(pr_DNN + "combined_perf.csv"):
                    continue
                self.c_DNN = DNN.DNN(pr_DNN, self.cMain.c_analysis.p_desc, self.cMain.c_analysis.p_AC50, self.cMain.c_analysis.p_desc_cleaned, self.cMain.c_analysis.p_AC50_cleaned, nb_repetition, n_foldCV, rate_active, rate_split, "classification")
                self.c_DNN.prepDataset()

                self.c_DNN.GridOptimizeDNN("MCC_train")
                self.c_DNN.evaluateModel()
                self.c_DNN.CrossValidation(n_foldCV)
                self.c_DNN.combineResults()
        
        # already existing undersampling
        else:
            print(self.cMain.c_analysis.pr_out)
            pr_sampling = self.cMain.c_analysis.pr_out + "sampleTraningSet/"
            for i in range(1, nb_repetition + 1):
                pr_rep = pr_sampling + str(i) + "/"
                for i in range(1, nb_repetition + 1):
                    pr_out = pr_rep + str(i) + "/"
                    pr_DNN = pathFolder.createFolder(pr_out + "DNNclass/")
                    self.c_DNN = DNN.DNN(pr_DNN, self.cMain.c_analysis.p_desc, self.cMain.c_analysis.p_AC50, self.cMain.c_analysis.p_desc_cleaned, self.cMain.c_analysis.p_AC50_cleaned, nb_repetition, n_foldCV, rate_active, rate_split, "classification")
                    
                    # only work in case of undersampling on the test set
                    p_train = pr_out + "train.csv"
                    p_test = pr_rep + "test.csv"

                    self.c_DNN.prepDataset(p_train, p_test)
                    self.c_DNN.GridOptimizeDNN("MCC_train")
                    self.c_DNN.evaluateModel()
                    self.c_DNN.CrossValidation(n_foldCV)
                    self.c_DNN.combineResults()

                self.c_DNN.mergeUndersampling(pr_rep)

    def DNN_builder_reg(self, cor_desc, quantile_desc, nb_repetition, n_foldCV, rate_split, rate_active):

        # prep as classiq ML
        self.cMain.prep_QSARClass(cor_desc, quantile_desc, nb_repetition, n_foldCV, rate_split, rate_active)

        # case of no undersampleing
        for i in range(0,nb_repetition):
            pr_DNN = pathFolder.createFolder(self.pr_results + "DNN_model_reg/" + str(i+1) + "/")
            #if path.exists(pr_DNN + "combined_perf.csv"):
            #    continue
                
            # prep data for regression
            p_train = pr_DNN + "trainSet.csv"
            p_test = pr_DNN + "testSet.csv"
            if not path.exists(p_train) or not path.exists(p_test):
                runExternal.prepDataQSARReg(self.cMain.c_analysis.p_desc, self.cMain.c_analysis.p_AC50, pr_DNN, cor_desc, quantile_desc, rate_split,  typeAff="All", logaff=0, nbNA = 10)

            self.c_DNN = DNN.DNN(pr_DNN, "", "", "", "", "", "", "", "", "regression")
            self.c_DNN.prepDataset(p_train, p_test)

            self.c_DNN.GridOptimizeDNN("MSE")
            self.c_DNN.evaluateModel()
            self.c_DNN.CrossValidation(n_foldCV)
            self.c_DNN.combineResults()
