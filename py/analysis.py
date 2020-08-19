from os import path
from statistics import mean, stdev
from shutil import copyfile, rmtree
from re import search

import pathFolder
import runExternal
import toolbox



class analysis:
    def __init__(self, p_AC50, p_desc, p_desc_OPERA, pr_out, cor_val, max_quantile):
        self.p_AC50 = p_AC50
        self.p_desc = p_desc
        self.p_desc_opera = p_desc_OPERA
        self.pr_out = pr_out
        self.cor_val = cor_val
        self.max_quantile = max_quantile



    def prepDesc(self):
        pr_out = pathFolder.createFolder(self.pr_out + "Cleaned_Data/")

        # add shortcut
        p_desc_cleaned = pr_out + "desc1D2D_cleaned.csv"
        p_AC50_cleaned = pr_out + "AC50_cleaned.csv"
        if not path.exists(p_desc_cleaned) and not path.exists(p_AC50_cleaned):
            runExternal.preprocData(self.p_desc, self.p_AC50, pr_out, self.cor_val, self.max_quantile)

        self.p_desc_cleaned = p_desc_cleaned
        self.p_AC50_cleaned = p_AC50_cleaned


    def combineDesc(self, pr_out):

        p_global = pr_out + 'desc_global.csv'

        if not path.exists(p_global):
            runExternal.combineAndPredDesc(self.p_desc, self.p_desc_opera, pr_out)

        self.p_desc = pr_out + 'desc_global.csv'
        self.pr_out = pr_out



    def sumAC50(self):
        pr_out = pathFolder.createFolder(self.pr_out + "Summary_AC50/")

        d_AC50 = toolbox.loadMatrix(self.p_AC50_cleaned, sep=",")

        nb_chem = len(list(d_AC50.keys()))
        l_AC50 = []
        for CASRN in d_AC50.keys():
            if d_AC50[CASRN]["Aff"] != "NA":
                l_AC50.append(float(d_AC50[CASRN]["Aff"]))
        
        p_filout = pr_out + "summary.txt"
        filout = open(p_filout, "w")
        filout.write("Nb chemicals: %s\nNb active: %s\nPercentage active: %s\nAverage log10(AC50): %s +\- %s"%(nb_chem, len(l_AC50), round((len(l_AC50)/nb_chem)*100, 2), round(mean(l_AC50),2), round(stdev(l_AC50),2)))
        filout.close()

        # histogram AC50        
        runExternal.histAC50(self.p_AC50_cleaned, pr_out)


    def PCA_plot(self):

        pr_out = pathFolder.createFolder(self.pr_out + "PCA/")
        runExternal.PCA(self.p_desc_cleaned, self.p_AC50_cleaned, pr_out)

    def signifDescActInact(self):

        pr_out = pathFolder.createFolder(self.pr_out + "SignifDesc/")
        # for rdkit
        runExternal.SignifDesc(self.p_desc, self.p_AC50, pr_out + "Rdkit_")
        # for OPERA
        runExternal.SignifDesc(self.p_desc_opera, self.p_AC50, pr_out + "OPERA_")



    def generate_SOM(self, grid_size):

        pr_out = pathFolder.createFolder(self.pr_out + "SOM/")
        p_model = pr_out + "SOM_model.RData"
        if not path.exists:
            runExternal.SOM(self.p_desc_cleaned, self.p_AC50_cleaned, pr_out, grid_size)


    def extract_actBySOMCluster(self, pr_png):

        p_cluster = pathFolder.createFolder(self.pr_out + "SOM/") + "SOM_Clusters_act"
        if not path.exists(p_cluster):
            return "Error"

        pr_out = pathFolder.createFolder(self.pr_out + "SOM/Active_by_cluster/")

        dcluster = toolbox.loadMatrix(p_cluster, sep = ",")
        for CASRN in dcluster.keys():
            cluster = dcluster[CASRN]["Cluster"]
            pr_cluster = pathFolder.createFolder(self.pr_out + "SOM/Active_by_cluster/" + cluster + "/")
            pathFolder.createFolder(pr_cluster)

            try:copyfile(pr_png + CASRN + ".png", pr_cluster + CASRN + ".png")
            except: pass


    def applySOMtoChemClassification(self, p_classification, pr_png):

        # define exist folder
        name_folder = p_classification.split("/")[-1][:-4]
        pr_out = pathFolder.createFolder(self.pr_out + "SOM/" + name_folder + "/")

        d_classification = toolbox.loadMatrix(p_classification, sep = ",")
        l_classif = []
        for chem in d_classification.keys():
            l_class = d_classification[chem]["Classes"]
            l_class = l_class.split("--")
            for classChem in l_class:
                if not classChem in l_classif and not classChem == "NA":
                    l_classif.append(classChem)

        # need to apply the model for each class and redifine a p_desc file based only on active chemical
        f_desc = open(self.p_desc, "r")
        l_chemdesc = f_desc.readlines()
        f_desc.close()

        d_AC50 = toolbox.loadMatrix(self.p_AC50_cleaned, sep = ",")

        for chemClass in l_classif:
            pr_out_sub = pathFolder.createFolder(pr_out + chemClass.replace(" ", "_").replace("(", "-").replace(")", "") + "/")
            p_fdesc = pr_out_sub + "desc.csv"
            f_desc_sp = open(p_fdesc, "w")
            f_desc_sp.write(l_chemdesc[0])

            c_chem = 0
            l_chem = []
            for chem in d_classification.keys():
                if not chem in list(d_AC50.keys()):
                    continue
                if d_AC50[chem]["Aff"] == "NA":
                    continue
                if search(chemClass, d_classification[chem]["Classes"]):
                    for chemdesc in l_chemdesc:
                        if search("^"+d_classification[chem]["CASRN"] + "\t", chemdesc):
                            if not d_classification[chem]["CASRN"] in l_chem:
                                c_chem = c_chem  + 1
                                #print(chemClass, chemdesc[0:10], d_classification[chem]["CASRN"], d_classification[chem]["Classes"])
                                f_desc_sp.write(chemdesc)
                                l_chem.append(d_classification[chem]["CASRN"])
                            break
            f_desc_sp.close()
            if c_chem > 1:
                runExternal.extractClusterSOM(p_fdesc, self.pr_out + "SOM/SOM_model.RData", pr_out_sub)

            else:
                rmtree(pr_out_sub)

        return 


    def signifDescBySOMCluster(self):

        # check if SOM is computed
        p_cluster = self.pr_out + "SOM/SOM_Clusters"
        if not path.exists(p_cluster):
            print("ERROR: no cluster file existed")
            return 
        
        pr_out = pathFolder.createFolder(self.pr_out + "SOM/DescriptorSignif/")
        runExternal.descSignifByCluster(self.p_desc_cleaned, p_cluster, pr_out)


    def HClust_plot(self):

        pr_out = pathFolder.createFolder(self.pr_out + "HClust/")
        runExternal.HClust(self.p_desc_cleaned, self.p_AC50_cleaned, pr_out)



