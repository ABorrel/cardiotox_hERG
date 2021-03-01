from os import path, listdir
from numpy import mean, std
from re import search

import pathFolder
import runExternal
import toolbox
import DNN

import CompDesc

class applyModel:
    def __init__(self, p_desc_model, p_desc_global, p_aff_model, p_desc_test, p_aff_test, pr_out, pr_results):
        self.p_desc_model = p_desc_model
        self.p_desc_model_all = p_desc_global
        self.p_desc_test = p_desc_test
        self.p_aff_model = p_aff_model
        self.p_aff_test = p_aff_test
        self.pr_out = pr_out
        self.pr_results = pr_results

    def predict_all_classif(self):

        ## apply model - NCAST RF undersampling
        pr_models = self.pr_results + "NCAST/QSAR_NCAST__0.9-90-5-10-0.15-0.3/sampleTraningSet/1/RF_models/QSARs_models/RF/"
        self.predictClassif(pr_models, "NCAST_classif_RF_undersampling", "RF")
        self.mergePredictionClass()

        ## apply model - NCAST RF nosampling
        pr_models = self.pr_results + "NCAST/QSAR_NCAST__0.9-90-5-10-0.15-0/noSampling/1/RFclass/"
        self.predictClassif(pr_models, "NCAST_classif_RF_nosampling", "RF")
        self.mergePredictionClass()

        ## apply model - NCAST-CHEMBL RF nosampling
        pr_models = self.pr_results + "NCAST_CHEMBL/QSAR_NCAST_CHEMBL__0.9-90-5-10-0.15-0/noSampling/1/RFclass/"
        self.predictClassif(pr_models, "NCAST_CHEMBL_classif_RF_nosampling", "RF")
        self.mergePredictionClass()


        ## apply model DNN - NCAST-CHEMBL DNN nosampling
        pr_models = self.pr_results + "NCAST_CHEMBL/QSAR_NCAST_CHEMBL__0.9-90-5-10-0.15-0/DNN_model/1/"
        self.predictDNN(pr_models, "NCAST_CHEMBL_DNN")
        self.mergePredictionClass()

    def PCACombine(self):
        pr_out = pathFolder.createFolder(self.pr_out + "PCA_vs/")
        runExternal.PCAvs(self.p_desc_model, self.p_aff_model, self.p_desc_test, self.p_aff_test, pr_out)

    def predictClassif(self, pr_models, name_model, ML):
        pr_out = pathFolder.createFolder(self.pr_out + "predict_model_" + name_model + "/")

        l_models = listdir(pr_models)
        print(l_models)
        for model in l_models:
            if search("RData", model) :
                runExternal.predictDataset(self.p_desc_test, pr_models + model, ML + "class", pr_out)
        self.pr_allPredict = pr_out

    def predictDNN(self, pr_models, name_model):
        pr_out = pathFolder.createFolder(self.pr_out + "predict_model_" + name_model + "/")
        
        # prep TEST set with only descriptor from TRAIN
        d_train = toolbox.loadMatrix(pr_models + "train.csv", sep = ",")
        d_test = toolbox.loadMatrix(self.p_desc_test)
        d_test_aff = toolbox.loadMatrix(self.p_aff_test, sep = ",")

        l_desc = list(d_train[list(d_train.keys())[0]].keys())
        p_test = pr_out + "test_cleaned.csv"
        f_test = open(p_test, "w")
        f_test.write(",".join(l_desc) + "\n")
        for chem in d_test.keys():
            l_val = [d_test[chem][desc] for desc in l_desc[1:-2]]
            if "NA" in l_val:
                continue
            else:
                f_test.write("%s,%s,%s\n"%(chem, ",".join(l_val), d_test_aff[chem]["Aff"]))
        f_test.close()

        d_test = toolbox.loadMatrix(p_test, sep = ",")

        l_models = listdir(pr_models)
        for model in l_models:
            if search("h5$", model) :
                c_DNN = DNN.DNN(pr_out, "", "", "", "", "", "", "", "" )
                l_pred = c_DNN.predictDNN(pr_models + model, p_test)
                p_pred = pr_out + "pred_model"
                f_pred = open(p_pred, "w")
                f_pred.write("\"\",\"ID\",\"Pred\"\n")
                i = 0
                for chem in d_test.keys():
                    f_pred.write("\"%s\",\"%s\",\"%s\"\n"%(chem, chem, l_pred[i][0]))
                    i = i + 1 
                f_pred.close()
        self.pr_allPredict = pr_out

    def predictReg(self, ML):
        pr_out = pathFolder.createFolder(self.pr_out + "predict_Reg_model_RF/")
        runExternal.predictRegDataset(self.p_desc_test, self.p_aff_test, self.p_AD, self.pr_models, ML, pr_out)
        
    def mergePredictionClass(self):

        p_filout = self.pr_allPredict + "Merge_pred.csv"
        if path.exists(p_filout) and path.getsize(p_filout) > 25:
            self.p_pred = p_filout
            return
        
        d_out = {}
        l_fpred = listdir(self.pr_allPredict)
        
        for fpred in l_fpred:
            if not search("^pred", fpred):
                continue
            dpred = toolbox.loadMatrix(self.pr_allPredict + fpred, sep = ",")

            for chemblID in dpred.keys():
                if dpred[chemblID]["Pred"] == "NA":
                    continue
                if not chemblID in list(d_out.keys()):
                    d_out[chemblID] = []
                try:d_out[chemblID].append(float(dpred[chemblID]["Pred"]))
                except:pass
        
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

    def applySOM(self, p_SOMmodel):

        pr_out = pathFolder.createFolder(self.pr_out + "SOM/")
        runExternal.applySOM(self.p_desc_test, p_SOMmodel, pr_out)

    def computeAD(self):

        pr_out = pathFolder.createFolder(self.pr_out + "AD/")
        
        p_out = pr_out + "AD_Test_zscore.csv"
        if path.exists(p_out):
            self.p_AD = p_out

        else:
            runExternal.AD(self.p_desc_model, self.p_desc_test, pr_out)
            self.p_AD = p_out  

    def applyAD(self):
        
        pr_out = pathFolder.createFolder(self.pr_out + "predict_AD/")
        runExternal.applyAD(self.p_pred, self.p_AD, pr_out)

    def overlapSetWithInchikey(self):

        pr_out = pathFolder.createFolder(self.pr_out + "OverlapModel/")

        # compare SMILES with a inchkey transformation
        
        # prepre train set #
        ####################
        d_desc_model = toolbox.loadMatrix(self.p_desc_model_all, sep = "\t")

        # prep descriptor for test set
        for chem in d_desc_model.keys():
            if "CASRN" in list(d_desc_model[chem].keys()):
                d_desc_model[chem]["ID"] = d_desc_model[chem]["CASRN"]

            elif "CHEMBLID" in list(d_desc_model[chem].keys()):
                d_desc_model[chem]["ID"] = d_desc_model[chem]["CHEMBLID"]

        # load AC50
        d_AC50_model = toolbox.loadMatrix(self.p_aff_model, sep = "\t")
        # transform -log10 to AC50
        for chem in d_AC50_model.keys():
            if d_AC50_model[chem]["LogAC50"]!= "NA":
                d_AC50_model[chem]["aff"] = 10**(-float(d_AC50_model[chem]["LogAC50"])) * 1000000
            else:
                d_AC50_model[chem]["aff"] = "NA"


        # Add inchkey in dictionarry model
        d_inch_model = {}
        for chem in d_desc_model.keys():
            SMILES = d_desc_model[chem]["SMILES"]
            cChem = CompDesc.CompDesc(SMILES, "", p_salts="./Salts.txt")
            inch = cChem.generateInchiKey()
            if not inch in list(d_inch_model.keys()):
                d_inch_model[inch] = d_desc_model[chem]



        # prep test set #
        ##################
        d_desc_test = toolbox.loadMatrix(self.p_desc_test, sep = "\t")

        # prep descriptor for test set
        for chem in d_desc_test.keys():
            if "CASRN" in list(d_desc_test[chem].keys()):
                d_desc_test[chem]["ID"] = d_desc_test[chem]["CASRN"]

            elif "CHEMBLID" in list(d_desc_test[chem].keys()):
                d_desc_test[chem]["ID"] = d_desc_test[chem]["CHEMBLID"]

        # test set
        d_AC50_test = toolbox.loadMatrix(self.p_aff_test, sep = ",")

        # Test set 
        d_inch_test = {}
        for chem in d_desc_test.keys():
            SMILES = d_desc_test[chem]["SMILES"]
            cChem = CompDesc.CompDesc(SMILES, "", p_salts="./Salts.txt")
            inch = cChem.generateInchiKey()
            if not inch in list(d_inch_test.keys()):
                d_inch_test[inch] = d_desc_test[chem]

        ###################################3


        # path filout
        p_filout = pr_out + "overlap_train_test"
        filout = open(p_filout, "w")
        filout.write("ID\tSMILES\tAC50 model\tAC50 test\tInclude\n")

        for inchTest in d_inch_test.keys():
            if inchTest in list(d_inch_model.keys()):
                filout.write("%s\t%s\t%s\t%s\t1\n"%(d_inch_test[inchTest]["ID"],d_inch_test[inchTest]["SMILES"] ,d_AC50_model[d_inch_model[inchTest]["ID"]]["aff"], d_AC50_test[d_inch_test[inchTest]["ID"]]["AC50-uM"]))
            else:
                filout.write("%s\t%s\tNA\t%s\t0\n"%(d_inch_test[inchTest]["ID"],d_inch_test[inchTest]["SMILES"] , d_AC50_test[d_inch_test[inchTest]["ID"]]["AC50-uM"]))
        filout.close()

    def overlapSetWithID(self, rm_overlap=0, rm_inactive=0):

        pr_out = pathFolder.createFolder(self.pr_out + "OverlapModel/")
        p_filout = pr_out + "overlap_train_test"

        print(self.p_desc_model)
        print(self.p_desc_test)
        print(self.p_aff_model)
        print(self.p_aff_test)

        # desc
        d_desc_model = toolbox.loadMatrix(self.p_desc_model, sep = ",")
        d_desc_test = toolbox.loadMatrix(self.p_desc_test, sep = "\t")

        # AC50
        d_AC50_model = toolbox.loadMatrix(self.p_aff_model, sep = ",")
        d_AC50_test = toolbox.loadMatrix(self.p_aff_test, sep = ",")

        # convert in -log
        for chem in d_AC50_test.keys():
            if d_AC50_test[chem]["LogAC50"] != "NA":
                d_AC50_test[chem]["LogAC50"] = -float(d_AC50_test[chem]["LogAC50"])

        filout = open(p_filout, "w")
        filout.write("ID\tLogAC50 model\tLogAC50 test\tAff\tInclude\n")


        for chem in d_AC50_test.keys():
            if not chem in list(d_desc_test.keys()): # case SMILES is not define
                continue
            if chem in list(d_AC50_model.keys()):
                if d_AC50_model[chem]["Aff"] == "NA":aff = 0
                else:aff = 1
                filout.write("%s\t%s\t%s\t%s\t1\n"%(chem, d_AC50_model[chem]["Aff"], d_AC50_test[chem]["LogAC50"], aff))
            else:
                if d_AC50_test[chem]["LogAC50"] == "NA":aff = 0
                else:aff = 1
                filout.write("%s\tNA\t%s\t%s\t0\n"%(chem, d_AC50_test[chem]["LogAC50"], aff))
        filout.close()
        runExternal.overlapPlot(p_filout)

        # write new pAff in case of rm_overlap = 1
        if rm_overlap == 1 and rm_inactive == 0:
            p_filout = pr_out + "aff_cleaned.csv"
            filout = open(p_filout, "w")
            filout.write("\"ID\",\"LogAC50\",\"Aff\"\n")
            
            for chem in d_AC50_test.keys():
                if not chem in list(d_AC50_model.keys()):
                    filout.write("\"%s\",\"%s\",\"%s\"\n"%(chem, d_AC50_test[chem]["LogAC50"], d_AC50_test[chem]["Aff"]))
            filout.close()
            self.p_aff_test = p_filout

            # rm overlap in the desc set
            p_desc_out = pr_out + "desc_cleaned.csv"

            d_desc = toolbox.loadMatrix(self.p_desc_test, sep ="\t")
            l_h = list(d_desc_test[list(d_desc_test.keys())[0]].keys())
            try:l_h.remove("SMILES")
            except: pass
            try:l_h.remove("CASRN")
            except:pass

            f_desc_out = open(p_desc_out, "w")
            f_desc_out.write("CASRN\t%s\n"%("\t".join(l_h)))
            for chem in d_desc.keys():
                 if not chem in list(d_AC50_model.keys()):
                    f_desc_out.write("%s\t%s\n"%(chem, "\t".join([d_desc[chem][h] for h in l_h])))
            f_desc_out.close()
            self.p_desc_test = p_desc_out
        

        elif rm_overlap == 1 and rm_inactive == 1:
            p_filout = pr_out + "aff_active_cleaned.csv"
            filout = open(p_filout, "w")
            filout.write("\"ID\",\"LogAC50\",\"Aff\"\n")
            
            #print(list(d_AC50_test.keys()))
            #print(list(d_AC50_model.keys()))

            for chem in d_AC50_test.keys():
                if not chem in list(d_AC50_model.keys()):
                    if d_AC50_test[chem]["Aff"] == "0" :
                        continue
                    else:
                        filout.write("\"%s\",\"%s\",\"%s\"\n"%(chem, d_AC50_test[chem]["LogAC50"], d_AC50_test[chem]["Aff"]))
            filout.close()
            self.p_aff_test = p_filout

            # rm overlap in the desc set
            p_desc_out = pr_out + "desc_active_cleaned.csv"

            d_desc = toolbox.loadMatrix(self.p_desc_test, sep ="\t")
            l_h = list(d_desc_test[list(d_desc_test.keys())[0]].keys())
            l_h.remove("CASRN")

            f_desc_out = open(p_desc_out, "w")
            f_desc_out.write("CASRN\t%s\n"%("\t".join(l_h)))
            for chem in d_desc.keys():
                 if not chem in list(d_AC50_model.keys()):
                    if d_AC50_test[chem]["LogAC50"] == "NA":
                        continue
                    else:
                        f_desc_out.write("%s\t%s\n"%(chem, "\t".join([d_desc[chem][h] for h in l_h])))
            f_desc_out.close()
            self.p_desc_test = p_desc_out