import pathFolder
import runExternal
import toolbox

from os import path, listdir, rename, remove
from re import search
from numpy import mean, std
from copy import deepcopy
from shutil import copyfile
from random import shuffle


class QSAR_modeling:
    def __init__(self, p_desc_clean, p_desc_origin, p_AC50, p_AC50_origin, pr_out, nb_repetition, n_foldCV, rate_active, rate_splitTrainTest):
        self.p_desc = p_desc_clean
        self.p_desc_orign = p_desc_origin
        self.p_AC50 = p_AC50
        self.p_AC50_orign = p_AC50_origin
        self.pr_out = pr_out
        self.repetition = nb_repetition
        self.n_foldCV = n_foldCV
        self.rate_active = rate_active
        self.rate_splitTrainTest = rate_splitTrainTest

    def runQSARClassUnderSamplingAllSet(self, force_run = 0):

        # check applicability model
        pr_AD = pathFolder.createFolder(self.pr_out + "AD/")
        l_run = list(range(1, self.repetition + 1))
        shuffle(l_run)

        for i in l_run:###### need to be change
            pr_run = self.pr_out + str(i) + "/"
            #rmtree(pr_run)############################################################################### to remove
            pathFolder.createFolder(pr_run)

            # prepare dataset => split train / test
            self.prepTrainSetforUnderSampling(pr_run)

            # check AD
            pr_AD_run = pathFolder.createFolder(pr_AD + str(i) + "/")
            if len(listdir(pr_AD_run)) < 6:
                runExternal.AD(self.p_train, self.p_test, pr_AD_run)

            # build QSAR
            self.buildQSAR(pr_run, force_run=force_run)

        # merge results
        self.mergeQSARs()

    def runQSARClassUnderSamplingTrain(self):

        # define train and test set here
        self.prepSplitTrainTestSet()

        # check applicability model
        pr_AD = pathFolder.createFolder(self.pr_out + "AD/")
        runExternal.AD(self.p_trainGlobal, self.p_test, pr_AD)

        for i in range(1, self.repetition + 1):
            pr_run = self.pr_out + str(i) + "/"
            #rmtree(pr_run)############################################################################### to remove
            pathFolder.createFolder(pr_run)

            # prepare dataset => split train / test
            self.prepTrainSetforUnderSampling(pr_run, 0)

            # build QSAR
            self.buildQSAR(pr_run)

        # merge results
        self.mergeQSARs()

    def runQSARReg(self, corcoef, maxQuantile):
        l_run = list(range(1, self.repetition + 1))
        shuffle(l_run)

        for i in l_run:
            pr_run = pathFolder.createFolder(self.pr_out + str(i) + "/")
            if not path.exists(pr_run + "trainSet.csv"):
                runExternal.prepDataQSARReg(self.p_desc, self.p_AC50, pr_run, corcoef, maxQuantile, self.rate_splitTrainTest,  typeAff="All", logaff=0, nbNA = 10)
            
            if not path.exists(pr_run + "perfTest.csv"): # control the perf test file is not present
                runExternal.runQSARReg(pr_run + "trainSet.csv", pr_run + "testSet.csv", "0", pr_run, self.n_foldCV)

    def mergeRegResults(self):

        d_result = {}
        pr_result = pathFolder.createFolder(self.pr_out + "mergeResults/")
        l_pr_run = listdir(self.pr_out)

        for pr_run in l_pr_run:
            if pr_run == "mergeResults" or search("csv", pr_run) or pr_run == "Cleaned_Data":
                continue
            else:
                d_result[pr_run] = {}
                p_perfTrain = self.pr_out + pr_run + "/perfTrain.csv"
                p_perfCV = self.pr_out + pr_run + "/perfCV.csv"
                p_perfTest = self.pr_out + pr_run + "/perfTest.csv"

                d_result[pr_run]["train"] =  toolbox.loadMatrix(p_perfTrain, sep = ",")
                d_result[pr_run]["test"] =  toolbox.loadMatrix(p_perfTest, sep = ",")
                d_result[pr_run]["CV"] =  toolbox.loadMatrix(p_perfCV, sep = ",")

        p_filout = pr_result + "results.csv"
        filout = open(p_filout, "w")
        for run in d_result.keys():
            filout.write("Run: %s\n"%(run))
            filout.write("\tCV\t\t\t\t\t\tTrain\t\t\t\t\t\tTest\t\t\t\n")
            filout.write("ML\tR2\tR0\tMAE\tcor\tRMSEP\t\tR2\tR0\tMAE\tcor\tRMSEP\t\tR2\tR0\tMAE\tcor\tRMSEP\n")
            for ML in d_result[run]["CV"].keys():
                filout.write("%s\t%s\t%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s\t%s\n"%(ML, d_result[run]["CV"][ML]["R2"], d_result[run]["CV"][ML]["R02"], d_result[run]["CV"][ML]["MAE"], d_result[run]["CV"][ML]["r"], d_result[run]["CV"][ML]["RMSEP"], d_result[run]["train"][ML]["R2"], d_result[run]["train"][ML]["R02"], d_result[run]["train"][ML]["MAE"], d_result[run]["train"][ML]["r"], d_result[run]["train"][ML]["RMSEP"], d_result[run]["test"][ML]["R2"], d_result[run]["test"][ML]["R02"], d_result[run]["test"][ML]["MAE"], d_result[run]["test"][ML]["r"], d_result[run]["test"][ML]["RMSEP"]))   
            filout.write("\n")

    def prepSplitTrainTestSet(self):

        # define train and test set
        p_train = self.pr_out + "trainGlobal.csv"
        p_test = self.pr_out + "test.csv"

        if path.exists(p_train) and path.exists(p_test):
            self.p_trainGlobal = p_train
            self.p_test = p_test
            return 

        # prep descriptor with classes without ratio of active => define train global
        runExternal.prepDataQSAR(self.p_desc, self.p_AC50, 0, self.pr_out)
        runExternal.SplitTrainTest(self.pr_out + "desc_Class.csv", self.pr_out, self.rate_splitTrainTest)
        # rename train in global train
        rename(self.pr_out + "train.csv", self.pr_out + "trainGlobal.csv")
        self.p_trainGlobal = p_train
        self.p_test = p_test

    def prepTrainSetforUnderSampling(self, pr_run, splitRatio=""):


        # Case where only the train set is created
        if splitRatio == 0: 
            # define train and test set
            p_train = pr_run + "train.csv"

            if path.exists(p_train):
                self.p_train = p_train

            else:
                # prep descriptor with classes without ratio of active
                runExternal.prepDataQSAR(self.p_trainGlobal, self.p_AC50, self.rate_active, pr_run)
                # to apply similar formating and create train.csv
                runExternal.SplitTrainTest(pr_run + "desc_Class.csv", pr_run, 0)
                
                # Rename file in train set file
                self.p_train = p_train

        else:
            # define train and test set
            p_train = pr_run + "train.csv"
            p_test = pr_run + "test.csv"

            if path.exists(p_train) and path.exists(p_test):
                self.p_train = p_train
                self.p_test = p_test
                return 

            # prep descriptor with classes
            if not path.exists(pr_run + "desc_Class.csv"):
                runExternal.prepDataQSAR(self.p_desc, self.p_AC50, self.rate_active, pr_run)

            runExternal.SplitTrainTest(pr_run + "desc_Class.csv", pr_run, self.rate_splitTrainTest)
            self.p_train = p_train
            self.p_test = p_test
    
    def buildQSAR(self, pr_run, force_run = 0):

        # do not run if performance file already exist
        p_perfTrain = pr_run + "perfTrain.csv"
        p_perfCV = pr_run + "perfCV.csv"
        p_perfTest = pr_run + "perfTest.csv"

        if force_run == 0:
            if path.exists(p_perfCV) and path.exists(p_perfTrain) and path.exists(p_perfTest):
                return 
            else:
                runExternal.runRQSAR(self.p_train, self.p_test, self.n_foldCV, pr_run)
        else:
            runExternal.runRQSAR(self.p_train, self.p_test, self.n_foldCV, pr_run)
                 
    def mergeQSARs(self):

        pr_QSAR_average = pathFolder.createFolder(self.pr_out + "Merge_results/")
        pr_QSAR_proba = pathFolder.createFolder(self.pr_out + "Merge_probRF/")# proba of prediction merged
        pr_QSAR_desc_involved = pathFolder.createFolder(self.pr_out + "Merge_involvedDesc/")
        pr_AD =  pathFolder.createFolder(self.pr_out + "Merge_AD/")

        self.mergeResults(pr_QSAR_average)
        self.mergeProbaRF(self.p_AC50_orign, pr_QSAR_proba)
        self.mergeInvolvedDesc("RF", 10, pr_QSAR_desc_involved)
        self.mergeAD(pr_AD)
        return 

    def mergeResults(self, pr_av):

        p_filout = pr_av + "average_perf.csv"
        if path.exists(p_filout):
            return 
        
        l_criteria = ["Acc", "Sp", "Se", "MCC"]
        l_dataset = ["CV", "train", "test"]
        
        # output to write
        d_result = {}
        for dataset in l_dataset:
            d_result[dataset] = {}

        # performance by ML
        d_perf = {}
        for criteria in l_criteria:
            d_perf[criteria] = []

        l_pr_run = listdir(self.pr_out)
        for pr_run in l_pr_run:
            if search("Merge", pr_run) or search("AD", pr_run):
                continue

            p_perfCV = self.pr_out + pr_run + "/perfCV.csv"
            p_perfTrain = self.pr_out + pr_run + "/perfTrain.csv"
            p_perfTest = self.pr_out + pr_run + "/perfTest.csv"

            try:
                M_CV = toolbox.loadMatrix(p_perfCV, sep=",")
                M_train = toolbox.loadMatrix(p_perfTrain, sep=",")
                M_test = toolbox.loadMatrix(p_perfTest, sep=",")
            except:
                continue

            l_ML = list(M_CV.keys())
            
            #build results
            if not l_ML[0] in list(d_result["CV"].keys()):
                for ML in l_ML:
                    d_result["CV"][ML] = deepcopy(d_perf)
                    d_result["train"][ML] = deepcopy(d_perf)
                    d_result["test"][ML] = deepcopy(d_perf)

            for ML in l_ML:
                for criteria in l_criteria:
                    d_result["CV"][ML][criteria].append(float(M_CV[ML][criteria]))
                    d_result["train"][ML][criteria].append(float(M_train[ML][criteria]))
                    d_result["test"][ML][criteria].append(float(M_test[ML][criteria]))


        d_out = deepcopy(d_result)
        for dataset in l_dataset:
            for ML in l_ML:
                for criteria in l_criteria:
                    AV = round(mean(d_result[dataset][ML][criteria]),3)
                    SD = round(std(d_result[dataset][ML][criteria]),3)
                    d_out[dataset][ML][criteria] = [AV, SD]


        # write result
        filout = open(p_filout, "w")
        for dataset in l_dataset:
            filout.write(str(dataset) + "\n")
            filout.write("\t" + "\t".join(["M-" + str(c) + "\t" + "SD-" + str(c) for c in l_criteria]) + "\n")
            for ML in l_ML:
                filout.write(ML)
                for criteria in l_criteria:
                    filout.write("\t" + str(d_out[dataset][ML][criteria][0]) + "\t" + str(d_out[dataset][ML][criteria][1]))
                filout.write("\n")
            filout.write("\n")
        filout.close()

    def mergeAD(self, pr_out):

        pr_AD_all = self.pr_out + "AD/"
        l_pr_run = listdir(pr_AD_all)

        d_Zscore_train = {}
        d_Zscore_test = {}
        for pr_run in l_pr_run:
            try:run = int(pr_run)
            except:continue
            p_Zscore_train = pr_AD_all + pr_run + "/AD_Train_zscore.csv"
            p_Zscore_test = pr_AD_all + pr_run + "/AD_Test_zscore.csv"
            
            if not path.exists(p_Zscore_train) and not path.exists(p_Zscore_test):
                p_Zscore_train = pr_AD_all + "AD_Train_zscore.csv"
                p_Zscore_test = pr_AD_all + "AD_Test_zscore.csv"
                
            d_train = toolbox.loadMatrix(p_Zscore_train, sep = ",")
            d_test = toolbox.loadMatrix(p_Zscore_test, sep = ",")

            for chem_train in d_train.keys():
                if not chem_train in list(d_Zscore_train.keys()):
                    d_Zscore_train[chem_train] = d_train[chem_train]["Zscore"]

            for chem_test in d_test.keys():
                if not chem_test in list(d_Zscore_test.keys()):
                    d_Zscore_test[chem_test] = d_test[chem_test]["Zscore"]
        
        if d_Zscore_test == {}:
            return 

        # write zscore
        p_train = pr_out + "Zscores_train.csv"
        ftrain = open(p_train, "w")
        ftrain.write("CASRN\tZscore\n")
        for chem in d_Zscore_train.keys():
            ftrain.write("%s\t%s\n"%(chem, d_Zscore_train[chem]))
        ftrain.close()

        p_test = pr_out + "Zscores_test.csv"
        ftest = open(p_test, "w")
        ftest.write("CASRN\tZscore\n")
        for chem in d_Zscore_test.keys():
            ftest.write("%s\t%s\n"%(chem, d_Zscore_test[chem]))
        ftest.close()


        # draw histogram for AD and add summary and PCA
        runExternal.mergeADs(p_train, p_test, self.p_desc, pr_AD_all)

    def mergeProbaRF(self, p_AC50, pr_prob):
        # need to change R scripts for 

        if len(listdir(pr_prob)) > 10:
            return 

        d_prob = {}
        l_dataset = ["CV", "train", "test"]

        l_pr_run = listdir(self.pr_out)
        d_AC50 = toolbox.loadMatrix(p_AC50)

        for pr_run in l_pr_run:
            if search("Merge", pr_run) or search("AD", pr_run) or search("Cleaned", pr_run) or search("desc_global.csv", pr_run) or search(".csv", pr_run) or search("merge", pr_run) :
                continue

            d_prob[pr_run] = {}

            # CV
            p_CV = self.pr_out + pr_run + "/RFclass/PerfRFClassCV10.txt"
            d_CV = toolbox.loadMatrix(p_CV, sep = "\t")
            d_prob[pr_run]["CV"] = d_CV

            # train
            p_train = self.pr_out + pr_run + "/RFclass/classTrain.csv"
            d_train = toolbox.loadMatrix(p_train, sep = ",")
            d_prob[pr_run]["train"] = d_train

            # test
            p_test = self.pr_out + pr_run + "/RFclass/classTest.csv"
            d_test = toolbox.loadMatrix(p_test, sep = ",")
            d_prob[pr_run]["test"] = d_test


        d_av = {}
        d_av["CV"] = {}
        d_av["train"] = {}
        d_av["test"] = {}

        for run in l_pr_run:
            if search("Merge", run) or search("AD", run) or search("Cleaned", run) or search("desc", run) or search("csv", run) or search("merge", run):
                continue

            # CV
            for chem in d_prob[run]["CV"].keys():
                if not chem in list(d_av["CV"].keys()):
                    d_av["CV"][chem] = {}
                    try: d_av["CV"][chem]["LogAC50"] = d_AC50[chem]["LogAC50"] 
                    except: d_av["CV"][chem]["LogAC50"] = d_AC50[chem]["Aff"] 
                    d_av["CV"][chem]["predicted"] = []
                d_av["CV"][chem]["predicted"].append(float(d_prob[run]["CV"][chem]["Predict"]))
            
            # train
            for chem in d_prob[run]["train"].keys():
                if not chem in list(d_av["train"].keys()):
                    d_av["train"][chem] = {}
                    try: d_av["train"][chem]["LogAC50"] = d_AC50[chem]["LogAC50"] 
                    except: d_av["train"][chem]["LogAC50"] = d_AC50[chem]["Aff"]
                    d_av["train"][chem]["predicted"] = []
                d_av["train"][chem]["predicted"].append(float(d_prob[run]["train"][chem]["x"]))

            # test
            for chem in d_prob[run]["test"].keys():
                if not chem in list(d_av["test"].keys()):
                    d_av["test"][chem] = {}
                    try:d_av["test"][chem]["LogAC50"] = d_AC50[chem]["LogAC50"] 
                    except:d_av["test"][chem]["LogAC50"] = d_AC50[chem]["Aff"]
                    d_av["test"][chem]["predicted"] = []
                d_av["test"][chem]["predicted"].append(float(d_prob[run]["test"][chem]["x"]))


            
        # train
        p_prob_train = pr_prob + "Prob_train"
        f_prob_train = open(p_prob_train, "w")
        f_prob_train.write("ID\tMpred\tSDpred\tReal\n")
        for IDtrain in list(d_av["train"].keys()):
            f_prob_train.write("%s\t%.3f\t%.3f\t%s\n"%(IDtrain, mean(d_av["train"][IDtrain]["predicted"]), std(d_av["train"][IDtrain]["predicted"]), d_av["train"][IDtrain]["LogAC50"]))
        f_prob_train.close()
        
        runExternal.plotAC50VSProb(p_prob_train)

        # CV
        p_prob_CV = pr_prob + "Prob_CV"
        f_prob_CV = open(p_prob_CV, "w")
        f_prob_CV.write("ID\tMpred\tSDpred\tReal\n")
        for IDCV in list(d_av["CV"].keys()):
            f_prob_CV.write("%s\t%.3f\t%.3f\t%s\n"%(IDCV, mean(d_av["CV"][IDCV]["predicted"]), std(d_av["CV"][IDCV]["predicted"]), d_av["CV"][IDCV]["LogAC50"]))
        f_prob_CV.close()
        runExternal.plotAC50VSProb(p_prob_CV)

        # test
        p_prob_test = pr_prob + "Prob_test"
        f_prob_test = open(p_prob_test, "w")
        f_prob_test.write("ID\tMpred\tSDpred\tReal\n")
        for IDtest in list(d_av["test"].keys()):
            f_prob_test.write("%s\t%.3f\t%.3f\t%s\n"%(IDtest, mean(d_av["test"][IDtest]["predicted"]), std(d_av["test"][IDtest]["predicted"]), d_av["test"][IDtest]["LogAC50"]))
        f_prob_test.close()
        runExternal.plotAC50VSProb(p_prob_test)

    def mergeInvolvedDesc(self, ML, nbdesc, pr_involvedDesc):


        d_importance = {}
        pr_out = pathFolder.createFolder(pr_involvedDesc + ML + "/")

        p_desc_importance = pr_out + "Av_importance"
        if path.exists(p_desc_importance):
            return

        l_pr_run = listdir(self.pr_out)
        for run in l_pr_run:
            if search("Merge", run) or search("AD", run) or search("Cleaned", run) or search("desc", run) or search("csv", run) or search("models", run) or search("mergeDNN", run):
                continue

            if not run in list(d_importance.keys()):
                d_importance[run] = {}

                p_desc_involved = self.pr_out + run + "/" + str(ML) + "class/ImportanceDesc"
                print(p_desc_involved)
                d_desc_involved = toolbox.loadMatrix(p_desc_involved , sep="\t")
                d_importance[run] = d_desc_involved


        # global importance 
        f_desc_importance = open(p_desc_importance, "w")
        f_desc_importance.write("Desc\tRun\tval\n")

        l_desc = list(d_importance["1"].keys())
        
        for desc in l_desc:
            for run in l_pr_run:
                if search("Merge", run) or search("AD", run) or search("Cleaned", run) or search("desc_global.csv", run) or search("desc", run) or search("csv", run) or search("models", run) or search("mergeDNN", run):
                    continue
                try: f_desc_importance.write(desc + "\t" + str(run) + "\t" + str(d_importance[run][desc]["x"]) + "\n")
                except: f_desc_importance.write(desc + "\t" + str(run) + "\t0.0\n")
        f_desc_importance.close()
        
        runExternal.runImportanceDesc(p_desc_importance, nbdesc)

    def extractModels(self, pr_results, ML):

        pr_out = pathFolder.createFolder(pr_results)

        l_pr_run = listdir(self.pr_out)
        l_model = listdir(pr_out)
        if len(l_model) > 2:
            return pr_out

        for pr_run in l_pr_run:
            if search("Merge", pr_run) or search("AD", pr_run) or search("Cleaned", pr_run) or search("desc", pr_run) or search("csv", pr_run) or search("models", pr_run) or search("mergeDNN", pr_run):
                continue

            if search("SVM", ML):
                p_model = "%s%s/SVMclass_%s/model.RData"%(self.pr_out, pr_run, ML.split("-")[-1])
            elif search("DNN", ML):
                l_DNNfile = listdir("%s%s/DNNclass/"%(self.pr_out, pr_run))
                p_model = "out"
                for DNNfile in l_DNNfile:
                    if search(".h5", DNNfile):
                        p_model = "%s%s/DNNclass/%s"%(self.pr_out, pr_run, DNNfile)
                        break
            else:
                p_model = "%s%s/%sclass/model.RData"%(self.pr_out, pr_run, ML)
            if path.exists(p_model) and not path.exists(pr_out + pr_run + ".RData"):
                if p_model[-5:] == "RData":
                    copyfile(p_model, pr_out + pr_run + ".RData")
                else:
                    copyfile(p_model, pr_out + pr_run + ".h5")
        
        return pr_out
        
