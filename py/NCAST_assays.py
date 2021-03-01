import pathFolder
import generic_run
import DNN
import hergml


class NCAST_assays:
    def __init__(self, p_smi, p_aff_origin, p_classification, cutoff_aff, pr_root, pr_data, pr_results):
        self.p_smi = p_smi
        self.p_aff_origin = p_aff_origin
        self.p_classification = p_classification
        self.cutoff_aff = cutoff_aff

        self.pr_root = pr_root
        self.pr_data = pr_data
        self.pr_results = pathFolder.createFolder(pr_results + "NCAST/")

    def main(self):

        # general val
        COR_VAL = 0.90
        MAX_QUANTILE = 90
        SOM_size = 15
        cutoff_aff = 30 #uM
        
        # QSAR values
        nb_repetition = 5
        n_foldCV = 10
        rate_split = 0.15

        self.prep_dataset()
        self.analysis_nopred(COR_VAL, MAX_QUANTILE, SOM_size)

        #QSAR - classification
        #######################
        rate_active = 0.30
        #self.QSARClassif_builder(COR_VAL, MAX_QUANTILE, nb_repetition, n_foldCV, rate_split, rate_active)

        rate_active = 0# no undersampling
        #self.QSARClassif_builder(COR_VAL, MAX_QUANTILE, nb_repetition, n_foldCV, rate_split, rate_active)

        #QSAR - regression
        #######################
        rate_active = 0
        #self.QSARReg_builder(COR_VAL, MAX_QUANTILE, nb_repetition, n_foldCV, rate_split, rate_active)


        #QSAR - DNN
        #######################
        #rate_active = 0
        #self.DNN_builder(COR_VAL, MAX_QUANTILE, nb_repetition, n_foldCV, rate_split, rate_active)

        
        #comparison - herg-ML
        ########################
        #p_pred_herg_ml = PR_ROOT + "comparison_study/pred/list_chemicals-2020-06-09-09-11-41_pred.csv"
        #p_test_set = PR_RESULTS + "QSARclass/1/test.csv"
        #cNCAST.comparison_hergml(p_pred_herg_ml, p_test_set)


    def prep_dataset(self):
        self.cMain = generic_run.generic_run(self.pr_root, self.pr_data, "NCAST", self.pr_results)
        self.cMain.setup_NCATS(self.p_smi, self.p_aff_origin, self.p_classification)
        self.cMain.load_datasetNCAST(self.cutoff_aff)
        self.cMain.classe_chemical_dataset()
    
        #compute desc
        self.cMain.compute_desc()
        #rank active
        #self.cMain.rankActive()

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

    def QSARReg_builder(self, cor_desc, quantile_desc, nb_repetition, n_foldCV, rate_split, rate_active):

        self.cMain.prep_QSARReg(cor_desc, quantile_desc, nb_repetition, n_foldCV, rate_split, rate_active)
        self.cMain.run_QSARReg()

    def DNN_builder(self, cor_desc, quantile_desc, nb_repetition, n_foldCV, rate_split, rate_active):

        # prep as classiq ML
        self.cMain.prep_QSARClass(cor_desc, quantile_desc, nb_repetition, n_foldCV, rate_split, rate_active)
        pr_DNN = pathFolder.createFolder(self.cMain.c_analysis.pr_out + "DNN_model/")
        for i in range(0,nb_repetition):
            pr_DNN = pathFolder.createFolder(self.cMain.c_analysis.pr_out + "DNN_model/" + str(i+1) + "/")
            if path.exists(pr_DNN + "combined_perf.csv"):
                continue
            self.c_DNN = DNN.DNN(pr_DNN, self.cMain.c_analysis.p_desc, self.cMain.c_analysis.p_AC50, self.cMain.c_analysis.p_desc_cleaned, self.cMain.c_analysis.p_AC50_cleaned, nb_repetition, n_foldCV, rate_active, rate_split)
            self.c_DNN.prepDataset()

            self.c_DNN.GridOptimizeDNN("MCC_train")
            self.c_DNN.evaluateModel()
            self.c_DNN.CrossValidation(n_foldCV)
            self.c_DNN.combineResults()

    def comparison_hergml(self, p_pred_herg_ml, p_test_set):
        
        pr_comparison_hergml = pathFolder.createFolder(self.pr_results + "comparison_hergml_" + self.name_dataset + "/")
        self.cMain.c_dataset.computeDescForNNComparison(self.pr_desc, pr_comparison_hergml)

        self.c_comparison_hergml = hergml.hergml(self.p_smi, p_pred_herg_ml, self.p_aff_origin, pr_comparison_hergml)
        self.c_comparison_hergml.computePerfAllSet()
        self.c_comparison_hergml.computePerfTestSet(p_test_set)
