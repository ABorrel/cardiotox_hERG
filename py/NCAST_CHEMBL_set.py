from os import path

import generic_run
import toolbox



class NCAST_CHEMBL_set:
    def __init__(self, p_dataset, pr_results, pr_root):
        self.p_dataset = p_dataset
        self.pr_root = pr_root
        self.pr_results = pr_results

    
    def main(self, c_CHEMBL_set, cNCAST):
        
        # general val
        COR_VAL = 0.90
        MAX_QUANTILE = 90
        SOM_size = 15
        cutoff_aff = 30 #uM
        
        # QSAR values
        nb_repetition = 5
        n_foldCV = 10
        rate_split = 0.15
        rate_active = 0

        self.set_up_for_analysis([c_CHEMBL_set.cCHEMBL.p_desc1D2D, cNCAST.cMain.l_p_desc[0]], [c_CHEMBL_set.cCHEMBL.p_descOPERA, cNCAST.cMain.l_p_desc[1]])
        self.analyseDataset(COR_VAL, MAX_QUANTILE, SOM_size)
        #self.QSARClassif_builder(COR_VAL, MAX_QUANTILE, nb_repetition, n_foldCV, rate_split, rate_active)
        self.QSARReg_builder(COR_VAL, MAX_QUANTILE, nb_repetition, n_foldCV, rate_split, rate_active)
        #self.c_NCAST_CHEMBL_set.QSARClassif_DNN(COR_VAL, MAX_QUANTILE, nb_repetition, n_foldCV, rate_split, rate_active)
    
    
    def set_up_for_analysis(self, l_pdesc1D2D, l_pdescOPERA):

        p_desc1D2D_out = self.pr_results + "desc1D2D.csv"
        p_descOPERA_out = self.pr_results + "descOPERA.csv"
        p_aff = self.pr_results + "aff.csv"

        if path.exists(p_desc1D2D_out) and path.exists(p_descOPERA_out) and path.exists(p_aff) :
            self.p_desc1D2D = p_desc1D2D_out
            self.p_descOPERA = p_descOPERA_out
            self.p_aff = p_aff
            self.cMain = generic_run.generic_run(self.pr_root, "", "NCAST_CHEMBL", self.pr_results)
            self.cMain.p_AC50 = self.p_aff
            self.cMain.l_p_desc = [self.p_desc1D2D, self.p_descOPERA]

        d_chem_dataset = toolbox.loadMatrix(self.p_dataset, sep = "\t")

        # rewrote for formating in analysis
        f_aff = open(p_aff, "w")
        f_aff.write("CASRN\tAff\n")
        for chem in d_chem_dataset.keys():
            if d_chem_dataset[chem]["CASRN"] != "-":
                ID = d_chem_dataset[chem]["CASRN"]
            else:
                ID = d_chem_dataset[chem]["CHEMBLID"]
            
            if d_chem_dataset[chem]["Aff"] == "0":
                aff = "NA"
            else:
                if d_chem_dataset[chem]["paff-add"] != "-":
                    aff = d_chem_dataset[chem]["paff-add"]
                else:
                    aff = d_chem_dataset[chem]["paff-chembl"]
            f_aff.write("%s\t%s\n"%(ID, aff))
        f_aff.close()

        # load OPERA descriptors
        d_opera = {}
        for p_descOPERA in l_pdescOPERA:
            d_temp = toolbox.loadMatrix(p_descOPERA, sep = ",")
            d_opera.update(d_temp)
        
        l_header_opera = list(d_opera[list(d_opera.keys())[0]].keys())
        try:l_header_opera.remove("CASRN")
        except:pass
        try:l_header_opera.remove("CHEMBLID")
        except:pass

        f_descOPERA = open(p_descOPERA_out, "w")
        f_descOPERA.write("CASRN," + ",".join(l_header_opera) + "\n")
        for chem in d_chem_dataset.keys():
            try:
                CASRN = d_chem_dataset[chem]["CASRN"]
                f_descOPERA.write("%s,%s\n"%(CASRN, ",".join([str(d_opera[CASRN][h]) for h in l_header_opera])))
                continue
            except:pass
            try:
                CHEMBL  = d_chem_dataset[chem]["CHEMBLID"]
                f_descOPERA.write("%s,%s\n"%(CHEMBL, ",".join([str(d_opera[CHEMBL][h]) for h in l_header_opera])))
            except: pass
        f_descOPERA.close()


        # load RDKIT descriptor
        d_desc1D2D = {}
        for p_desc1D2D in l_pdesc1D2D:
            d_temp = toolbox.loadMatrix(p_desc1D2D)
            d_desc1D2D.update(d_temp)

        l_header_1D2D = list(d_desc1D2D[list(d_desc1D2D.keys())[0]].keys())
        try:l_header_1D2D.remove("CASRN")
        except:pass
        try:l_header_1D2D.remove("CHEMBLID")
        except:pass
        l_header_1D2D.remove("SMILES")

        f_desc1D2D = open(p_desc1D2D_out, "w")
        f_desc1D2D.write("CASRN\tSMILES\t" + "\t".join(l_header_1D2D) + "\n")
        for chem in d_chem_dataset.keys():
            try:
                CASRN = d_chem_dataset[chem]["CASRN"]
                f_desc1D2D.write("%s\t%s\t%s\n"%(CASRN, d_desc1D2D[CASRN]["SMILES"], "\t".join([str(d_desc1D2D[CASRN][h]) for h in l_header_1D2D])))
                continue
            except:pass
            try:
                CHEMBL  = d_chem_dataset[chem]["CHEMBLID"]
                f_desc1D2D.write("%s\t%s\t%s\n"%(CHEMBL, d_desc1D2D[CHEMBL]["SMILES"], "\t".join([str(d_desc1D2D[CHEMBL][h]) for h in l_header_1D2D])))
                continue
            except:
                pass
        
        f_desc1D2D.close()

        self.p_desc1D2D = p_desc1D2D_out
        self.p_descOPERA = p_descOPERA_out
        self.p_aff = p_aff

        # set up main
        self.cMain = generic_run.generic_run(self.pr_root, "", "NCAST_CHEMBL", self.pr_results)
        self.cMain.p_AC50 = self.p_aff
        self.cMain.l_p_desc = [self.p_desc1D2D, self.p_descOPERA]

    def analyseDataset(self, cor_desc, quantile_desc, som_size):

        self.cMain.setup_descAnalysis(cor_desc, quantile_desc)
        
        self.cMain.run_structural_PCA()
        self.cMain.run_structural_Hclust()
        self.cMain.run_structural_SOM(som_size)

    def QSARClassif_builder(self, cor_desc, quantile_desc, nb_repetition, n_foldCV, rate_split, rate_active):

        self.cMain.prep_QSARClass(cor_desc, quantile_desc, nb_repetition, n_foldCV, rate_split, rate_active)
        self.cMain.run_QSARClassif("training")

    def QSARClassif_DNN(self, cor_desc, quantile_desc, nb_repetition, n_foldCV, rate_split, rate_active):
        
        # prep as classiq ML
        self.cMain.prep_QSARClass(cor_desc, quantile_desc, nb_repetition, n_foldCV, rate_split, rate_active)
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

    def QSARReg_builder(self, cor_desc, quantile_desc, nb_repetition, n_foldCV, rate_split, rate_active):

        self.cMain.prep_QSARReg(cor_desc, quantile_desc, nb_repetition, n_foldCV, rate_split, rate_active)
        self.cMain.run_QSARReg()

