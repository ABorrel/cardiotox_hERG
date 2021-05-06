import pathFolder
import dataset
import analysis
import QSAR_modeling


class generic_run:
    def __init__(self, pr_root, pr_data, name_dataset, pr_results):
        self.pr_root = pr_root
        self.pr_data = pr_data
        self.pr_results = pr_results
        self.name_dataset = name_dataset
        self.pr_desc = pathFolder.createFolder(pr_root + "results/DESC/")
    
    def setup_NCATS(self, p_NCATS_smi, p_NCAST_AC50, p_class_chemical):
        self.p_smi = p_NCATS_smi
        self.p_AC50 = p_NCAST_AC50
        self.pr_dataset = pathFolder.createFolder(self.pr_results + "dataset_" + self.name_dataset + "/")
        self.p_class_chemical = p_class_chemical


    def load_datasetNCAST(self, cutoff_aff = 30):
        """
        cutoff in uM
        """
        # 1.1 load dataset
        c_dataset = dataset.dataset(self.p_smi, self.p_AC50, self.pr_dataset)
        c_dataset.prep_dataset(cutoff_aff)

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


    def run_structural_PCA_NCATS_CHEMBL(self):
        self.c_analysis.PCA_NCATS_CHEMBL_plot()

    def run_structural_SOM(self, SOM_size):
        self.c_analysis.generate_SOM(SOM_size)
        self.c_analysis.signifDescBySOMCluster()
        self.c_analysis.extract_actBySOMCluster(self.pr_desc + "PNG/") # have to run !!!!
        if "p_class_chemical" in self.__dict__:
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

    def prep_QSARReg(self, cor_desc, quantile_desc, nb_repetition, n_foldCV, rate_split, rate_active):
        self.cor_desc = cor_desc
        self.quantile_desc = quantile_desc        
        self.nb_repetition = nb_repetition
        self.n_foldCV = n_foldCV
        self.rate_split = rate_split
        self.rate_active = rate_active

        pr_out = pathFolder.createFolder("%sQSARReg_%s__%s-%s-%s-%s-%s-%s/"%(self.pr_results, self.name_dataset, self.cor_desc, self.quantile_desc, self.nb_repetition, self.n_foldCV, self.rate_split, self.rate_active))
        self.c_analysis = analysis.analysis(self.p_AC50, self.l_p_desc[0], self.l_p_desc[1], pr_out, self.cor_desc, self.quantile_desc)
        self.c_analysis.combineDesc()
        self.c_analysis.prepDesc()


    def run_QSARReg(self):

        self.cQSARReg = QSAR_modeling.QSAR_modeling(self.c_analysis.p_desc_cleaned, self.l_p_desc[0], self.c_analysis.p_AC50_cleaned, self.p_AC50, self.c_analysis.pr_out, self.nb_repetition, self.n_foldCV, self.rate_active, self.rate_split)
        self.cQSARReg.runQSARReg(self.c_analysis.cor_val, self.c_analysis.max_quantile)
        self.cQSARReg.mergeRegResults()

    def run_QSARClassif(self, sampling):

        if sampling == "global":
            pr_out = pathFolder.createFolder(self.c_analysis.pr_out + "sampleGlobalSet/")
            self.cQSAR = QSAR_modeling.QSAR_modeling(self.c_analysis.p_desc_cleaned, self.l_p_desc[0], self.c_analysis.p_AC50_cleaned, self.p_AC50, pr_out, self.nb_repetition, self.n_foldCV, self.rate_active, self.rate_split)
            self.cQSAR.runQSARClassUnderSamplingAllSet(force_run=0)
            
            self.cQSAR.pr_RF_models = self.cQSAR.extractModels("%sRF_%s__globalsampling__%s-%s-%s-%s-%s-%s/"%(pr_out, self.name_dataset, self.cor_desc, self.quantile_desc, self.nb_repetition, self.n_foldCV, self.rate_split, self.rate_active), "RF")
            self.cQSAR.pr_LDA_models = self.cQSAR.extractModels("%sLDA_%s__globalsampling__%s-%s-%s-%s-%s-%s/"%(pr_out, self.name_dataset, self.cor_desc, self.quantile_desc, self.nb_repetition, self.n_foldCV, self.rate_split, self.rate_active), "LDA")

        elif sampling == "training" and self.rate_active != 0:
            for i in range(1, self.nb_repetition + 1):
                pr_out = pathFolder.createFolder(self.c_analysis.pr_out + "sampleTraningSet/" + str(i) + "/")
                self.cQSAR = QSAR_modeling.QSAR_modeling(self.c_analysis.p_desc_cleaned, self.l_p_desc[0], self.c_analysis.p_AC50_cleaned, self.p_AC50, pr_out, self.nb_repetition, self.n_foldCV, self.rate_active, self.rate_split)
                self.cQSAR.runQSARClassUnderSamplingTrain()
                
                self.cQSAR.pr_RF_models = self.cQSAR.extractModels("%smodels_out/RF_models/"%(pr_out), "RF")
                self.cQSAR.pr_LDA_models = self.cQSAR.extractModels("%smodels_out/LDA_models/"%(pr_out), "LDA")
                self.cQSAR.pr_SVMradial_models = self.cQSAR.extractModels("%smodels_out/SVM-radial_models/"%(pr_out), "SVM-radial")
                self.cQSAR.pr_SVMlinear_models = self.cQSAR.extractModels("%smodels_out/SVM-linear_models/"%(pr_out), "SVM-linear")
                self.cQSAR.pr_DNN_models = self.cQSAR.extractModels("%smodels_out/DNN_models/"%(pr_out), "DNN")

        elif self.rate_active == 0:
            pr_out = pathFolder.createFolder(self.c_analysis.pr_out + "noSampling/")
            self.cQSAR = QSAR_modeling.QSAR_modeling(self.c_analysis.p_desc_cleaned, self.l_p_desc[0], self.c_analysis.p_AC50_cleaned, self.p_AC50, pr_out, self.nb_repetition, self.n_foldCV, 0, self.rate_split)
            self.cQSAR.runQSARClassUnderSamplingAllSet(force_run=0)
            
            self.cQSAR.pr_RF_models = self.cQSAR.extractModels("%sRF_%s__nosampling__%s-%s-%s-%s-%s-%s/"%(pr_out, self.name_dataset, self.cor_desc, self.quantile_desc, self.nb_repetition, self.n_foldCV, self.rate_split, self.rate_active), "RF")
            self.cQSAR.pr_LDA_models = self.cQSAR.extractModels("%sLDA_%s__nosampling__%s-%s-%s-%s-%s-%s/"%(pr_out, self.name_dataset, self.cor_desc, self.quantile_desc, self.nb_repetition, self.n_foldCV, self.rate_split, self.rate_active), "LDA")

