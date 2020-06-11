import pathFolder
import runExternal
from os import path


class QSAR_modeling:
    def __init__(self, p_desc_clean, p_AC50, pr_out, nb_repetition, n_foldCV, rate_active, rate_splitTrainTest):
        self.p_desc = p_desc_clean
        self.p_AC50 = p_AC50
        self.pr_out = pr_out
        self.repetition = nb_repetition
        self.n_foldCV = n_foldCV
        self.rate_active = rate_active
        self.rate_splitTrainTest = rate_splitTrainTest


    def runQSARClass(self):

        for i in range(1, self.repetition + 1):
            pr_run = self.pr_out + str(i) + "/"
            #rmtree(pr_run)############################################################################### to remove
            pathFolder.createFolder(pr_run)

            # prepare dataset => split train / test
            self.prepForQSAR(pr_run)

            # build QSAR
            self.buildQSAR(pr_run)

        
        self.mergeQSAR(pr_run)




    def prepForQSAR(self, pr_run):


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
        



    
    def buildQSAR(self, pr_run):

        runExternal.runRQSAR(self.p_train, self.p_test, self.n_foldCV, pr_run)




    def mergeQSAR(self):

        return 