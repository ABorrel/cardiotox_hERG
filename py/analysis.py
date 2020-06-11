from os import path
from statistics import mean, stdev
from shutil import copyfile

import pathFolder
import runExternal
import toolbox



class analysis:
    def __init__(self, p_AC50, p_desc, pr_out, cor_val, max_quantile):
        self.p_AC50 = p_AC50
        self.p_desc = p_desc
        self.pr_out = pr_out
        self.cor_val = cor_val
        self.max_quantile = max_quantile



    def prepDesc(self):
        pr_out = pathFolder.createFolder(self.pr_out + "Cleaned_Data/")

        # add shortcut
        p_desc_cleaned = pr_out + "desc1D2D_cleaned.csv"
        p_AC50_cleaned = pr_out + "AC50_cleaned.csv"
        if not path.exists(p_desc_cleaned) and path.exists(p_AC50_cleaned):
            runExternal.preprocData(self.p_desc, self.p_AC50, pr_out, self.cor_val, self.max_quantile)

        self.p_desc_cleaned = p_desc_cleaned
        self.p_AC50_cleaned = p_AC50_cleaned


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


    def generate_SOM(self):

        pr_out = pathFolder.createFolder(self.pr_out + "SOM/")
        runExternal.SOM(self.p_desc_cleaned, self.p_AC50_cleaned, pr_out)


    def analyse_SOM(self, pr_png):

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


    def HClust_plot(self):

        pr_out = pathFolder.createFolder(self.pr_out + "HClust/")
        runExternal.HClust(self.p_desc_cleaned, self.p_AC50_cleaned, pr_out)

