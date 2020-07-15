from os import path, listdir
from numpy import mean, std

import pathFolder
import runExternal
import toolbox


class applyModel:
    def __init__(self, p_desc_model, p_aff_model, p_desc_test, p_aff_test, pr_models, pr_out):
        self.p_desc_model = p_desc_model
        self.p_desc_test = p_desc_test
        self.p_aff_model = p_aff_model
        self.p_aff_test = p_aff_test
        self.pr_models = pr_models
        self.pr_out = pr_out

    
    def PCACombine(self):
        pr_out = pathFolder.createFolder(self.pr_out + "PCA_vs/")
        runExternal.PCAvs(self.p_desc_model, self.p_aff_model, self.p_desc_test, self.p_aff_test, pr_out)


    def predict(self):
        name_model = self.pr_models.split("/")[-2]
        pr_out = pathFolder.createFolder(self.pr_out + "predict_model_" + name_model + "/")

        l_models = listdir(self.pr_models)
        for model in l_models:
            runExternal.predictDataset(self.p_desc_test, self.pr_models + model, name_model + "class", pr_out)

        self.pr_allPredict = pr_out

        
    def mergePrediction(self):

        p_filout = self.pr_allPredict + "Merge_pred.csv"
        if path.exists(p_filout) and path.getsize(p_filout) > 25:
            self.p_pred = p_filout
            return
        
        d_out = {}
        l_fpred = listdir(self.pr_allPredict)
        
        for fpred in l_fpred:
            dpred = toolbox.loadMatrix(self.pr_allPredict + fpred, sep = ",")

            for chemblID in dpred.keys():
                if not chemblID in list(d_out.keys()):
                    d_out[chemblID] = []
                d_out[chemblID].append(float(dpred[chemblID]["Pred"]))
        
        d_aff = toolbox.loadMatrix(self.p_aff_test, sep = ",")

        filout = open(p_filout, "w")
        filout.write("ID\tMpred\tSDpred\tReal\n")
        for chemblID in d_out.keys():
            M_pred = mean(d_out[chemblID])
            SD_pred = std(d_out[chemblID])
            filout.write("%s\t%s\t%s\t%s\n"%(chemblID, round(M_pred, 2), round(SD_pred, 2), d_aff[chemblID]["Aff"]))
        filout.close()
        self.p_pred = p_filout

        # plot prob vs aff
        runExternal.plotAC50VSProb(p_filout)

        return p_filout

    def applySOM(self, p_SOMmodel):

        pr_out = pathFolder.createFolder(self.pr_out + "SOM/")
        runExternal.applySOM(self.p_desc_test, p_SOMmodel, pr_out)

    def computeAD(self):

        pr_out = pathFolder.createFolder(self.pr_out + "AD/")
        
        p_out = pr_out + "AD_zscore.csv"
        if path.exists(p_out):
            self.p_AD = p_out

        else:
            runExternal.AD(self.p_desc_model, self.p_desc_test, pr_out)
            self.p_AD = p_out
    

    def applyAD(self):
        
        pr_out = pathFolder.createFolder(self.pr_out + "predict_AD/")
        runExternal.applyAD(self.p_pred, self.p_AD, pr_out)
