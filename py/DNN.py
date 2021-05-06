import toolbox
import QSAR_modeling
import pathFolder
import runExternal

from numpy import loadtxt, arange
from sklearn import metrics
from sklearn.model_selection import StratifiedKFold, KFold
from sklearn.preprocessing import StandardScaler
sc = StandardScaler()

import keras 
from keras.datasets import mnist 
from keras.models import Sequential 
from keras.layers import Dense, Dropout 
from keras.optimizers import RMSprop 
import tensorflow as tf
import numpy as np 
from os import listdir, remove, path
from re import search
from copy import deepcopy


class DNN:
    def __init__(self, pr_out, p_desc, p_aff,  p_desc_clean, p_aff_clean, nb_repetition, n_foldCV, rate_active, rate_splitTrainTest, typeModel):
        self.pr_out = pr_out
        self.p_desc = p_desc
        self.p_aff = p_aff
        self.p_aff_clean = p_aff_clean
        self.p_desc_clean = p_desc_clean
        self.nb_repetition = nb_repetition
        self.n_foldCV = n_foldCV
        self.rate_active = rate_active
        self.rate_splitTrainTest = rate_splitTrainTest
        self.verbose = 0
        

        # optimization
        self.typeModel = typeModel
        if typeModel == "classification":
            self.kernel_initializer = "random_normal"
            self.loss = "binary_crossentropy"
        else:
            self.kernel_initializer = "normal"
            self.loss = "mean_squared_error"
        self.optimizer = "adam"
        self.l_epochs = [50, 100, 120]
        self.l_batch_size = [32, 64, 128]
        self.l_dense_layer = [3, 4, 5]
        self.l_dense_candidate = [50, 25, 20, 10, 1]
        self.l_activation = ["relu", "selu"]

    def prepDataset(self, p_train = "", p_test = ""):

        self.cQSAR = QSAR_modeling.QSAR_modeling(self.p_desc_clean, self.p_desc, self.p_aff_clean, self.p_aff, self.pr_out, self.nb_repetition, self.n_foldCV, self.rate_active, self.rate_splitTrainTest)
        if path.exists(p_train) and path.exists(p_test):
            self.cQSAR.p_train = p_train
            self.cQSAR.p_test = p_test
        else:
            # split train test
            self.cQSAR.prepTrainSetforUnderSampling(self.pr_out)      

        # descriptor train set 
        with open(self.cQSAR.p_train) as f:
            #determining number of columns from the first line of text
            n_cols = len(f.readline().split(","))
        self.dataset_train = loadtxt(self.cQSAR.p_train, delimiter=",",usecols=arange(1, n_cols-1), skiprows=1)
        self.dataset_train = sc.fit_transform(self.dataset_train)
        self.aff_train = loadtxt(self.cQSAR.p_train, delimiter=",",usecols=arange(n_cols-1, n_cols), skiprows=1)
        self.nb_desc_input = n_cols - 2

        d_train = toolbox.loadMatrix(self.cQSAR.p_train, sep = ",")
        self.l_train = list(d_train.keys())

        # descriptor test set 
        with open(self.cQSAR.p_test) as f:
            #determining number of columns from the first line of text
            n_cols = len(f.readline().split(","))
        self.dataset_test = loadtxt(self.cQSAR.p_test, delimiter=",",usecols=arange(1, n_cols-1), skiprows=1)
        self.dataset_test = sc.fit_transform(self.dataset_test)
        self.aff_test = loadtxt(self.cQSAR.p_test, delimiter=",",usecols=arange(n_cols-1, n_cols), skiprows=1)

        d_test = toolbox.loadMatrix(self.cQSAR.p_test, sep = ",")
        self.l_test = list(d_test.keys())

    #l_batch_size = [32,64,128]
    def GridOptimizeDNN(self, criteria_best_perf):

        # check if model is in the 
        l_filout = listdir(self.pr_out)
        for p_filout in l_filout:
            if search(".h5$", p_filout):
                self.p_model = self.pr_out + p_filout
                return 

        filout = open(self.pr_out + "optimization", "w")
        d_out = {}
        best_perf = 0.0
        for activation in self.l_activation:
            for dense_layer in self.l_dense_layer:
                nb_layer = dense_layer
                self.initModel(activation, self.l_dense_candidate[-nb_layer])
                nb_layer = nb_layer - 1
                while nb_layer > 1:
                    self.addDense(activation, self.l_dense_candidate[-nb_layer])
                    nb_layer = nb_layer - 1
                for epochs in self.l_epochs:
                    for batch_size in self.l_batch_size:
                        l_pref_train = self.fit_compileModel(epochs, batch_size)

                        l_pref_test = self.evaluateOnTest()
                        if self.typeModel == "classification":
                            filout.write("Activation: %s; Nb layers: %s; Epochs: %s; Batch size: %s;\nTRAINNING SET\nAccuracy: %s\nAccuracy-balanced: %s\nMCC: %s\nRecall: %s\nTEST SET\nAccuracy: %s\nAccuracy-balanced: %s\nMCC: %s\nRecall: %s\n[====================]\n"%(activation, dense_layer, epochs, batch_size, l_pref_train[0], l_pref_train[1], l_pref_train[2], l_pref_train[3], l_pref_test[0], l_pref_test[1], l_pref_test[2], l_pref_test[3]))
                            k_out = "%s--%s--%s--%s"%(activation, dense_layer, epochs, batch_size)

                            # save model
                            if criteria_best_perf == "MCC_train":
                                if  float(l_pref_train[2]) > best_perf:
                                    # remove previous model
                                    l_items = listdir(self.pr_out)
                                    for item in l_items:
                                        if item.endswith(".h5"):
                                            remove(path.join(self.pr_out, item))
                                    best_perf = float(l_pref_train[2])
                                    self.model.save(self.pr_out + "best_model_" + k_out + ".h5")
                                    self.p_model = self.pr_out + "best_model_" + k_out + ".h5"
                        
                        else:
                            filout.write("Activation: %s; Nb layers: %s; Epochs: %s; Batch size: %s;\nTRAINNING SET\nMSE: %s\nMAE: %s\nR2: %s\nTEST SET\nMSE: %s\nMAE: %s\nR2: %s\n[====================]\n"%(activation, dense_layer, epochs, batch_size, l_pref_train[3], l_pref_train[0], l_pref_train[1], l_pref_test[3], l_pref_test[0], l_pref_test[1]))
                            k_out = "%s--%s--%s--%s"%(activation, dense_layer, epochs, batch_size)

                            # save model
                            if criteria_best_perf == "MSE":
                                if  float(l_pref_train[3]) > best_perf:
                                    # remove previous model
                                    l_items = listdir(self.pr_out)
                                    for item in l_items:
                                        if item.endswith(".h5"):
                                            remove(path.join(self.pr_out, item))
                                    best_perf = float(l_pref_train[2])
                                    self.model.save(self.pr_out + "best_model_" + k_out + ".h5")
                                    self.p_model = self.pr_out + "best_model_" + k_out + ".h5"
        filout.close()                

    def addDense(self, activation, dense_candidate):
        self.model.add(Dense(int(dense_candidate), activation=activation, kernel_initializer=self.kernel_initializer))
        #self.model.add(Dropout(0.5))

    def initModel(self, activation, dense_candidate):
        self.model = Sequential()
        if dense_candidate > self.nb_desc_input:
            dense_candidate = self.nb_desc_input
        self.model.add(Dense(dense_candidate, input_dim=self.nb_desc_input, activation=activation, kernel_initializer=self.kernel_initializer))
        #self.model.add(Dropout(0.5))
    
    def fit_compileModel(self, epochs, batch_size, l_metric = ["mse"], learning_rate = 0.1):

        # last layer
        if self.typeModel == "classification":
            self.model.add(Dense(1, activation='sigmoid'))
        else:
            self.model.add(Dense(1))

        if self.optimizer == "sgd":
            sgd = tf.keras.optimizers.SGD(learning_rate=learning_rate)
            self.optimizer = sgd
        #elif optimizer == "adam":
        #    adam =  tf.keras.optimizers.Adam(learning_rate=learning_rate)
        #    optimizer = adam
        #    print(adam)

        # compile the keras model
        self.model.compile(loss=self.loss, optimizer=self.optimizer, metrics=l_metric)
        # fit the keras model on the dataset
        self.model.fit(self.dataset_train, self.aff_train, epochs=epochs, batch_size=batch_size), 

        # evaluate the keras model
        self.model.evaluate(self.dataset_train, self.aff_train)
        
        if self.typeModel == "classification":
            # predict on train
            y_pred = self.model.predict(self.dataset_train)
            y_pred = [1. if pred[0] > 0.5 else 0. for pred in y_pred]
            
            acc = metrics.accuracy_score(self.aff_train, y_pred)
            bacc = metrics.balanced_accuracy_score(self.aff_train, y_pred)
            mcc = metrics.matthews_corrcoef(self.aff_train, y_pred)
            recall = metrics.recall_score(self.aff_train, y_pred)
            
            print("======= TRAIN ======")
            print("Acc: ", acc)
            print("b-Acc: ", bacc)
            print("MCC: ", mcc)
            print("Recall: ", recall)
            
            return [acc, bacc, mcc, recall]
        
        else:
            y_pred = self.model.predict(self.dataset_train)

            MAE = metrics.mean_absolute_error(self.aff_train, y_pred)
            R2 = metrics.r2_score(self.aff_train, y_pred)
            EVS = metrics.explained_variance_score(self.aff_train, y_pred)
            MSE = metrics.mean_squared_error(self.aff_train, y_pred)
            MAXERR = metrics.max_error(self.aff_train, y_pred)
            try:MSE_log = metrics.mean_squared_log_error(self.aff_train, y_pred)
            except: MSE_log = 0.0
            MDAE = metrics.median_absolute_error(self.aff_train, y_pred)
            MTD = metrics.mean_tweedie_deviance(self.aff_train, y_pred)
            try:
                MPD = metrics.mean_poisson_deviance(self.aff_train, y_pred)
                MGD = metrics.mean_gamma_deviance(self.aff_train, y_pred)
            except:
                MPD = 0.0
                MGD = 0.0

            print("======= TRAIN ======")
            print("MAE: ", MAE)
            print("R2: ", R2)
            print("Explain Variance score: ", EVS)
            print("MSE: ", MSE)
            print("Max error: ", MAXERR)
            print("MSE log: ", MSE_log)
            print("Median absolute error: ", MDAE)
            print("Mean tweedie deviance: ", MTD)
            print("Mean poisson deviance: ", MPD)
            print("Mean gamma deviance: ", MGD)

            return [MAE, R2, EVS, MSE, MAXERR, MSE_log, MDAE, MTD, MPD, MGD]

    def evaluateOnTest(self):

        y_pred = self.model.predict(self.dataset_test)
        if self.typeModel == "classification":
            y_pred = [1. if pred[0] > 0.5 else 0. for pred in y_pred]
            
            acc = metrics.accuracy_score(self.aff_test, y_pred)
            bacc = metrics.balanced_accuracy_score(self.aff_test, y_pred)
            mcc = metrics.matthews_corrcoef(self.aff_test, y_pred)
            recall = metrics.recall_score(self.aff_test, y_pred)

            print("======= TEST ======")
            print("Acc: ", acc)
            print("b-Acc: ", bacc)
            print("MCC: ", mcc)
            print("Recall: ", recall)
            

            return [acc, bacc, mcc, recall]
        else:

            MAE = metrics.mean_absolute_error(self.aff_test, y_pred)
            R2 = metrics.r2_score(self.aff_test, y_pred)
            EVS = metrics.explained_variance_score(self.aff_test, y_pred)
            MSE = metrics.mean_squared_error(self.aff_test, y_pred)
            MAXERR = metrics.max_error(self.aff_test, y_pred)
            try:MSE_log = metrics.mean_squared_log_error(self.aff_test, y_pred)
            except:MSE_log = 0.0
            MDAE = metrics.median_absolute_error(self.aff_test, y_pred)
            MTD = metrics.mean_tweedie_deviance(self.aff_test, y_pred)
            try:MPD = metrics.mean_poisson_deviance(self.aff_test, y_pred)
            except:MPD = 0.0
            try:MGD = metrics.mean_gamma_deviance(self.aff_test, y_pred)
            except:MGD = 0.0

            print("======= TEST ======")
            print("MAE: ", MAE)
            print("R2: ", R2)
            print("Explain Variance score: ", EVS)
            print("MSE: ", MSE)
            print("Max error: ", MAXERR)
            print("MSE log: ", MSE_log)
            print("Median absolute error: ", MDAE)
            print("Mean tweedie deviance: ", MTD)
            print("Mean poisson deviance: ", MPD)
            print("Mean gamma deviance: ", MGD)

            return [MAE, R2, EVS, MSE, MAXERR, MSE_log, MDAE, MTD, MPD, MGD]
    
    def evaluateModel(self):

        if "p_model" in self.__dict__:
            model = keras.models.load_model(self.p_model)
        else:
            print("No model ready to be loaded")
            return

        if not "d_perf" in self.__dict__:
            self.d_perf = {}

        if self.typeModel == "classification":
            # trainning
            y_train_pred = model.predict(self.dataset_train)
            y_train_pred = [1. if pred[0] > 0.5 else 0. for pred in y_train_pred]
            
            acc_train = metrics.accuracy_score(self.aff_train, y_train_pred)
            bacc_train = metrics.balanced_accuracy_score(self.aff_train, y_train_pred)
            mcc_train = metrics.matthews_corrcoef(self.aff_train, y_train_pred)
            recall_train = metrics.recall_score(self.aff_train, y_train_pred)
            sp_train = (bacc_train * 2) - recall_train


            print("======= TRAIN ======")
            print("Acc: ", acc_train)
            print("b-Acc: ", bacc_train)
            print("MCC: ", mcc_train)
            print("Recall - Se: ", recall_train)
            print("Sp: ", sp_train)
            
            self.d_perf["Train"] = {}
            self.d_perf["Train"]["Acc"] = acc_train
            self.d_perf["Train"]["b-Acc"] = bacc_train
            self.d_perf["Train"]["MCC"] = mcc_train
            self.d_perf["Train"]["Se"] = recall_train
            self.d_perf["Train"]["Sp"] = sp_train


            y_pred = model.predict(self.dataset_test)
            y_pred = [1. if pred[0] > 0.5 else 0. for pred in y_pred]
            
            acc_test = metrics.accuracy_score(self.aff_test, y_pred)
            bacc_test = metrics.balanced_accuracy_score(self.aff_test, y_pred)
            mcc_test = metrics.matthews_corrcoef(self.aff_test, y_pred)
            recall_test = metrics.recall_score(self.aff_test, y_pred)
            sp_test = (bacc_test * 2) - recall_test

            print("======= TEST ======")
            print("Acc: ", acc_test)
            print("b-Acc: ", bacc_test)
            print("MCC: ", mcc_test)
            print("Recall - Se: ", recall_test)
            print("Sp: ", sp_test)


            self.d_perf["Test"] = {}
            self.d_perf["Test"]["Acc"] = acc_test
            self.d_perf["Test"]["b-Acc"] = bacc_test
            self.d_perf["Test"]["MCC"] = mcc_test
            self.d_perf["Test"]["Se"] = recall_test
            self.d_perf["Test"]["Sp"] = sp_test
        
        else:
            # trainning
            y_train_pred = model.predict(self.dataset_train)
            y_train_pred = [pred[0] for pred in y_train_pred]

            MAE_train = metrics.mean_absolute_error(self.aff_train, y_train_pred)
            R2_train = metrics.r2_score(self.aff_train, y_train_pred)
            EVS_train = metrics.explained_variance_score(self.aff_train, y_train_pred)
            MSE_train = metrics.mean_squared_error(self.aff_train, y_train_pred)
            MAXERR_train = metrics.max_error(self.aff_train, y_train_pred)
            try: MSE_log_train = metrics.mean_squared_log_error(self.aff_train, y_train_pred)
            except: MSE_log_train = 0.0
            MDAE_train = metrics.median_absolute_error(self.aff_train, y_train_pred)
            MTD_train = metrics.mean_tweedie_deviance(self.aff_train, y_train_pred)
            try:
                MPD_train = metrics.mean_poisson_deviance(self.aff_train, y_train_pred)
                MGD_train = metrics.mean_gamma_deviance(self.aff_train, y_train_pred)
            except:
                MPD_train = 0.0
                MGD_train = 0.0
            
            print("======= TRAIN ======")
            print("MAE: ", MAE_train)
            print("R2: ", R2_train)
            print("Explain Variance score: ", EVS_train)
            print("MSE: ", MSE_train)
            print("Max error: ", MAXERR_train)
            print("MSE log: ", MSE_log_train)
            print("Median absolute error: ", MDAE_train)
            print("Mean tweedie deviance: ", MTD_train)
            print("Mean poisson deviance: ", MPD_train)
            print("Mean gamma deviance: ", MGD_train)
            
            self.d_perf["Train"] = {}
            self.d_perf["Train"]["MAE"] = MAE_train
            self.d_perf["Train"]["R2"] = R2_train
            self.d_perf["Train"]["Explain Variance score"] = EVS_train
            self.d_perf["Train"]["MSE"] = MSE_train
            self.d_perf["Train"]["Max error"] = MAXERR_train
            self.d_perf["Train"]["MSE log"] = MSE_log_train
            self.d_perf["Train"]["Median absolute error"] = MDAE_train
            self.d_perf["Train"]["Mean tweedie deviance"] = MTD_train
            self.d_perf["Train"]["Mean poisson deviance"] = MPD_train
            self.d_perf["Train"]["Mean gamma deviance"] = MGD_train

            filout = open(self.pr_out + "train_pred.csv", "w")
            filout.write("\tReal\tPred\n")
            i = 0
            imax = len(self.l_train)
            
            while i < imax:
                filout.write("%s\t%s\t%s\n"%(self.l_train[i], self.aff_train[i], y_train_pred[i]))
                i = i + 1
            filout.close()

            runExternal.plotPerfCor(self.pr_out + "train_pred.csv")

            y_pred = model.predict(self.dataset_test)
            y_pred = [pred[0] for pred in y_pred]
            
            MAE_test = metrics.mean_absolute_error(self.aff_test, y_pred)
            R2_test = metrics.r2_score(self.aff_test, y_pred)
            EVS_test = metrics.explained_variance_score(self.aff_test, y_pred)
            MSE_test = metrics.mean_squared_error(self.aff_test, y_pred)
            MAXERR_test = metrics.max_error(self.aff_test, y_pred)
            try:
                MSE_log_test = metrics.mean_squared_log_error(self.aff_test, y_pred)
            except:
                MSE_log_test = 0.0
            MDAE_test = metrics.median_absolute_error(self.aff_test, y_pred)
            MTD_test = metrics.mean_tweedie_deviance(self.aff_test, y_pred)
            try:
                MPD_test = metrics.mean_poisson_deviance(self.aff_test, y_pred)
                MGD_test = metrics.mean_gamma_deviance(self.aff_test, y_pred)
            except:
                MPD_test = 0.0
                MGD_test = 0.0
            
            
            print("======= TEST ======")
            print("MAE: ", MAE_test)
            print("R2: ", R2_test)
            print("Explain Variance score: ", EVS_test)
            print("MSE: ", MSE_test)
            print("Max error: ", MAXERR_test)
            print("MSE log: ", MSE_log_test)
            print("Median absolute error: ", MDAE_test)
            print("Mean tweedie deviance: ", MTD_test)
            print("Mean poisson deviance: ", MPD_test)
            print("Mean gamma deviance: ", MGD_test)
            
            self.d_perf["Test"] = {}
            self.d_perf["Test"]["MAE"] = MAE_test
            self.d_perf["Test"]["R2"] = R2_test
            self.d_perf["Test"]["Explain Variance score"] = EVS_test
            self.d_perf["Test"]["MSE"] = MSE_test
            self.d_perf["Test"]["Max error"] = MAXERR_test
            self.d_perf["Test"]["MSE log"] = MSE_log_test
            self.d_perf["Test"]["Median absolute error"] = MDAE_test
            self.d_perf["Test"]["Mean tweedie deviance"] = MTD_test
            self.d_perf["Test"]["Mean poisson deviance"] = MPD_test
            self.d_perf["Test"]["Mean gamma deviance"] = MGD_test

            filout = open(self.pr_out + "test_pred.csv", "w")
            filout.write("\tReal\tPred\n")
            i = 0
            imax = len(self.l_test)
            while i < imax:
                filout.write("%s\t%s\t%s\n"%(self.l_test[i], self.aff_test[i], y_pred[i]))
                i = i + 1
            filout.close()

            runExternal.plotPerfCor(self.pr_out + "test_pred.csv")

    def combineResults(self):
        if not "d_perf" in self.__dict__:
            print("No data to write")
            return 

        p_filout = self.pr_out + "combined_perf.csv"
        filout = open(p_filout, "w")
        l_perf = ["CV", "Train", "Test"]
        if self.typeModel == "classification":
            filout.write("\tAcc\tSp\tSe\tMCC\tb-Acc\n")
            for perf in l_perf:
                if perf in list(self.d_perf.keys()):
                    filout.write("%s-DNN\t%s\t%s\t%s\t%s\t%s\n"%(perf, self.d_perf[perf]["Acc"],self.d_perf[perf]["Sp"], self.d_perf[perf]["Se"], self.d_perf[perf]["MCC"], self.d_perf[perf]["b-Acc"]))
            filout.close()
        else:
            filout.write("\tMAE\tR2\tExplain Variance score\tMSE\tMax error\tMSE log\tMedian absolute error\tMean tweedie deviance\tMean poisson deviance\tMean gamma deviance\n")
            for perf in l_perf:
                if perf in list(self.d_perf.keys()):
                    filout.write("%s-DNN\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(perf, self.d_perf[perf]["MAE"],self.d_perf[perf]["R2"], self.d_perf[perf]["Explain Variance score"], self.d_perf[perf]["MSE"], self.d_perf[perf]["Max error"], self.d_perf[perf]["MSE log"], self.d_perf[perf]["Median absolute error"], self.d_perf[perf]["Mean tweedie deviance"], self.d_perf[perf]["Mean poisson deviance"], self.d_perf[perf]["Mean gamma deviance"]))
            filout.close()

    def CrossValidation(self, n_folds):
        seed = 7
        if self.typeModel == "classification":
            kfold = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=seed)
        else:
            kfold = KFold(n_splits=n_folds, shuffle=True, random_state=seed)
        # load parameters
        l_files = listdir(self.pr_out)
        for name_file in l_files:
            if search("h5$", name_file):
                l_parameter = name_file.split("_")[-1][:-3].split("--")
                act = l_parameter[0]
                nb_layer = int(l_parameter[1])
                epochs = int(l_parameter[2])
                batch_size = int(l_parameter[3])
                break
        

        if self.typeModel == "classification":
            d_CV = {}
            d_CV["acc"] = []
            d_CV["b-acc"] = []
            d_CV["recall"] = []
            d_CV["mcc"] = []
        else:
            d_CV = {}
            d_CV["MAE"] = []
            d_CV["R2"] = []
            d_CV["Explain Variance score"] = []
            d_CV["MSE"] = []
            d_CV["Max error"] = []
            d_CV["MSE log"] = []
            d_CV["Median absolute error"] = []
            d_CV["Mean tweedie deviance"] = []
            d_CV["Mean poisson deviance"] = []
            d_CV["Mean gamma deviance"] = []

        d_train = deepcopy(self.dataset_train)
        Y_train = deepcopy(self.aff_train)

        for train, test in kfold.split(d_train, Y_train):
            self.dataset_train = d_train[train]
            self.aff_train = Y_train[train]

            self.dataset_test = d_train[test]
            self.aff_test = Y_train[test]

            # create model
            self.initModel(act, self.l_dense_candidate[-nb_layer])
            nb_layer_temp = nb_layer - 1
            while nb_layer_temp > 1:
                self.addDense(act, self.l_dense_candidate[-nb_layer_temp])
                nb_layer_temp = nb_layer_temp - 1

            l_pref_train = self.fit_compileModel(epochs, batch_size)
            l_pref_test = self.evaluateOnTest()
            
            if self.typeModel == "classification":
                # CV performance
                d_CV["acc"].append(l_pref_test[0])
                d_CV["b-acc"].append(l_pref_test[1])
                d_CV["mcc"].append(l_pref_test[2])
                d_CV["recall"].append(l_pref_test[3])
            else:
                d_CV["MAE"].append(l_pref_test[0])
                d_CV["R2"].append(l_pref_test[1])
                d_CV["Explain Variance score"].append(l_pref_test[2])
                d_CV["MSE"].append(l_pref_test[3])
                d_CV["Max error"].append(l_pref_test[4])
                d_CV["MSE log"].append(l_pref_test[5])
                d_CV["Median absolute error"].append(l_pref_test[6])
                d_CV["Mean tweedie deviance"].append(l_pref_test[7])
                d_CV["Mean poisson deviance"].append(l_pref_test[8])
                d_CV["Mean gamma deviance"].append(l_pref_test[9])


        if self.typeModel == "classification":
            filout = open(self.pr_out + "CV_%s.perf"%(n_folds), "w")
            filout.write("Acc\tB-acc\tMCC\tRecall\n")
            filout.write("%.2f+/-%.2f\t%.2f+/-%.2f\t%.2f+/-%.2f\t%.2f+/-%.2f\n"%(np.mean(d_CV["acc"]), np.std(d_CV["acc"]), np.mean(d_CV["b-acc"]), np.std(d_CV["b-acc"]),np.mean(d_CV["mcc"]), np.std(d_CV["mcc"]), np.mean(d_CV["recall"]), np.std(d_CV["recall"])))
            filout.close()
        else:
            filout = open(self.pr_out + "CV_%s.perf"%(n_folds), "w")
            filout.write("MAE\tR2\tExplain Variance score\tMSE\tMax error\tMSE log\tMedian absolute error\tMean tweedie deviance\tMean poisson deviance\tMean gamma deviance\n")
            filout.write("%.2f+/-%.2f\t%.2f+/-%.2f\t%.2f+/-%.2f\t%.2f+/-%.2f\t%.2f+/-%.2f\t%.2f+/-%.2f\t%.2f+/-%.2f\t%.2f+/-%.2f\t%.2f+/-%.2f\t%.2f+/-%.2f\n"%(np.mean(d_CV["MAE"]), np.std(d_CV["MAE"]), np.mean(d_CV["R2"]), np.std(d_CV["R2"]),np.mean(d_CV["Explain Variance score"]), np.std(d_CV["Explain Variance score"]), np.mean(d_CV["MSE"]), np.std(d_CV["MSE"]), np.mean(d_CV["Max error"]), np.std(d_CV["Max error"]), np.mean(d_CV["MSE log"]), np.std(d_CV["MSE log"]), np.mean(d_CV["Median absolute error"]), np.std(d_CV["Median absolute error"]), np.mean(d_CV["Mean tweedie deviance"]), np.std(d_CV["Mean tweedie deviance"]), np.mean(d_CV["Mean poisson deviance"]), np.std(d_CV["Mean poisson deviance"]), np.mean(d_CV["Mean gamma deviance"]), np.std(d_CV["Mean gamma deviance"])))
            filout.close()

        if not "d_perf" in self.__dict__:
            self.d_perf = {}

        self.d_perf['CV'] = {}
        if self.typeModel == "classification":
            self.d_perf["CV"]["Acc"] = np.mean(d_CV["acc"])
            self.d_perf["CV"]["b-Acc"] =  np.mean(d_CV["b-acc"])
            self.d_perf["CV"]["MCC"] = np.mean(d_CV["mcc"])
            self.d_perf["CV"]["Se"] = np.mean(d_CV["recall"])
            self.d_perf["CV"]["Sp"] = (self.d_perf["CV"]["b-Acc"] * 2) - self.d_perf["CV"]["Se"]
        else:
            self.d_perf["CV"]["MAE"] = np.mean(d_CV["MAE"])
            self.d_perf["CV"]["R2"]= np.mean(d_CV["MAE"])
            self.d_perf["CV"]["Explain Variance score"]= np.mean(d_CV["Explain Variance score"])
            self.d_perf["CV"]["MSE"]= np.mean(d_CV["MSE"])
            self.d_perf["CV"]["Max error"]= np.mean(d_CV["Max error"])
            self.d_perf["CV"]["MSE log"]= np.mean(d_CV["MSE log"])
            self.d_perf["CV"]["Median absolute error"]= np.mean(d_CV["Median absolute error"])
            self.d_perf["CV"]["Mean tweedie deviance"]= np.mean(d_CV["Mean tweedie deviance"])
            self.d_perf["CV"]["Mean poisson deviance"]= np.mean(d_CV["Mean poisson deviance"])
            self.d_perf["CV"]["Mean gamma deviance"]= np.mean(d_CV["Mean gamma deviance"])
  
    def predictDNN(self, p_model, p_test, name_out = "test_pred"):
        # descriptor test set 
        with open(p_test) as f:
            #determining number of columns from the first line of text
            n_cols = len(f.readline().split(","))
        dataset_test = loadtxt(p_test, delimiter=",",usecols=arange(1, n_cols-1), skiprows=1)
        dataset_test = sc.fit_transform(dataset_test)
        
        y_test = loadtxt(p_test, delimiter=",",usecols=arange(n_cols-1, n_cols), skiprows=1)
        d_test = toolbox.loadMatrix(p_test)

        model = keras.models.load_model(p_model)
        y_pred = model.predict(dataset_test)

        # compute quality performance
        if self.typeModel == "regression":
            MAE = metrics.mean_absolute_error(y_test, y_pred)
            R2 = metrics.r2_score(y_test, y_pred)
            EVS = metrics.explained_variance_score(y_test, y_pred)
            MSE = metrics.mean_squared_error(y_test, y_pred)
            MAXERR = metrics.max_error(y_test, y_pred)
            try:MSE_log = metrics.mean_squared_log_error(y_test, y_pred)
            except:MSE_log = 0.0
            MDAE = metrics.median_absolute_error(y_test, y_pred)
            MTD = metrics.mean_tweedie_deviance(y_test, y_pred)
            try:MPD = metrics.mean_poisson_deviance(y_test, y_pred)
            except:MPD = 0.0
            try:MGD = metrics.mean_gamma_deviance(y_test, y_pred)
            except:MGD = 0.0

            if self.verbose == 1:
                print("======= TEST SET ======")
                print("MAE:", MAE)
                print("R2:", R2)
                print("Explain Variance score:", EVS)
                print("MSE:", MSE)
                print("Max error:", MAXERR)
                print("MSE log:", MSE_log)
                print("Median absolute error:", MDAE)
                print("Mean tweedie deviance:", MTD)
                print("Mean poisson deviance:", MPD)
                print("Mean gamma deviance:", MGD)
            
            p_filout = self.pr_out + name_out + ".sum"
            filout = open(p_filout, "w")
            filout.write("\tMAE\tR2\tExplain Variance score\tMSE\tMax error\tMSE log\tMedian absolute error\tMean tweedie deviance\tMean poisson deviance\tMean gamma deviance\n")
            filout.write("TEST-DNN\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(MAE, R2, EVS, MSE, MAXERR, MSE_log, MDAE, MTD, MPD, MGD))
            filout.close()

        else:
            print("TEST")
        
        p_filout = self.pr_out + name_out
        filout = open(p_filout, "w")
        filout.write("\"\",\"ID\",\"Pred\"\n")
        i = 0
        l_chem = list(d_test.keys())
        imax = len(l_chem)
        
        while i < imax:
            filout.write("\"\",\"ID\",\"Pred\"\n")
            filout.write("\"%s\",\"%s\",\"%s\"\n"%(l_chem[i], l_chem[i], y_pred[i]))
            i = i + 1
        filout.close()
        
        return y_pred

    def mergeUndersampling(self, pr_rep):

        d_out = {}
        d_out ["train"] = {"Acc":[], "Se":[], "Sp":[], "MCC":[], "b-Acc":[]}
        d_out ["test"] = {"Acc":[], "Se":[], "Sp":[], "MCC":[], "b-Acc":[]}
        d_out ["CV"] = {"Acc":[], "Se":[], "Sp":[], "MCC":[], "b-Acc":[]}
        pr_out = pathFolder.createFolder(pr_rep + "mergeDNN/")
        for i in range(1, self.nb_repetition + 1):
            p_perf = "%s/%s/DNNclass/combined_perf.csv"%(pr_rep, i)
            d_perf = toolbox.loadMatrix(p_perf)
            d_out["train"]["Acc"].append(float(d_perf["Train-DNN"]["Acc"]))
            d_out["train"]["Se"].append(float(d_perf["Train-DNN"]["Se"]))
            d_out["train"]["Sp"].append(float(d_perf["Train-DNN"]["Sp"]))
            d_out["train"]["MCC"].append(float(d_perf["Train-DNN"]["MCC"]))
            d_out["train"]["b-Acc"].append(float(d_perf["Train-DNN"]["b-Acc"]))

            d_out["test"]["Acc"].append(float(d_perf["Test-DNN"]["Acc"]))
            d_out["test"]["Se"].append(float(d_perf["Test-DNN"]["Se"]))
            d_out["test"]["Sp"].append(float(d_perf["Test-DNN"]["Sp"]))
            d_out["test"]["MCC"].append(float(d_perf["Test-DNN"]["MCC"]))
            d_out["test"]["b-Acc"].append(float(d_perf["Test-DNN"]["b-Acc"]))

            d_out["CV"]["Acc"].append(float(d_perf["CV-DNN"]["Acc"]))
            d_out["CV"]["Se"].append(float(d_perf["CV-DNN"]["Se"]))
            d_out["CV"]["Sp"].append(float(d_perf["CV-DNN"]["Sp"]))
            d_out["CV"]["MCC"].append(float(d_perf["CV-DNN"]["MCC"]))
            d_out["CV"]["b-Acc"].append(float(d_perf["CV-DNN"]["b-Acc"]))
        
        l_perf = ["Acc", "Se", "Sp", "MCC", "b-Acc"]
        p_filout = pr_out + "average_perf.csv"
        filout = open(p_filout, "w")
        filout.write("\t" + "\t".join(["M-%s\tSD-%s"%(perf, perf)for perf in l_perf]) + "\n")

        for k in d_out.keys():
            filout.write("%s\t%s\n"%(k, "\t".join(["%s\t%s"%(np.average(d_out[k][perf]), np.std(d_out[k][perf]))for perf in l_perf])))
        filout.close()





#########
### FOR TESTING

    def buildDNNTest(self, l_metric = ["mse"]):
        

        # define the keras model
        model = Sequential()
        model.add(Dense(12, input_dim=self.nb_desc_input, activation='relu'))
        model.add(Dense(8, activation='relu'))
        model.add(Dense(1, activation='sigmoid'))
        # compile the keras model
        model.compile(loss='binary_crossentropy', optimizer='adam', metrics=l_metric)
        # fit the keras model on the dataset
        model.fit(self.dataset_train, self.aff_train, epochs=10, batch_size=10)

        # evaluate the keras model
        _, accuracy = model.evaluate(self.dataset_train, self.aff_train)
        print('ACC: %.2f' % (accuracy))

        y_pred = model.predict(self.dataset_test)
        y_pred = [1. if pred[0] > 0.5 else 0. for pred in y_pred]
        
        acc = metrics.accuracy_score(self.aff_test, y_pred)
        bacc = metrics.balanced_accuracy_score(self.aff_test, y_pred)
        mcc = metrics.matthews_corrcoef(self.aff_test, y_pred)
        recall = metrics.recall_score(self.aff_test, y_pred)

        print("======= TEST ======")
        print("Acc: ", acc)
        print("b-Acc: ", bacc)
        print("MCC: ", mcc)
        print("Recall: ", recall)

    def test(self):

        #dataset = loadtxt('C:/Users/aborr/research/ILS/HERG/data/pima-indians-diabetes.data.csv', delimiter=',')
        dataset = loadtxt('/mnt/c/Users/aborr/research/ILS/HERG/data/pima-indians-diabetes.data.csv', delimiter=',')
        X = dataset[:,0:8]
        y = dataset[:,8]

        print(X)

        # define the keras model
        model = Sequential()
        model.add(Dense(12, input_dim=8, activation='relu'))
        model.add(Dense(8, activation='relu'))
        model.add(Dense(1, activation='sigmoid'))
        # compile the keras model
        model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])
        # fit the keras model on the dataset
        model.fit(X, y, epochs=150, batch_size=10)
        # evaluate the keras model
        _, accuracy = model.evaluate(X, y)
        print('Accuracy: %.2f' % (accuracy*100))



        return 

    def test2(self):

        (x_train, y_train), (x_test, y_test) = mnist.load_data() 

        print(x_train)
        sss

        x_train = x_train.reshape(60000, 784) 
        x_test = x_test.reshape(10000, 784) 
        x_train = x_train.astype('float32') 
        x_test = x_test.astype('float32') 
        x_train /= 255 
        x_test /= 255 

        y_train = keras.utils.to_categorical(y_train, 10) 
        y_test = keras.utils.to_categorical(y_test, 10) 

        model = Sequential() 
        model.add(Dense(512, activation='relu', input_shape = (784,))) 
        model.add(Dropout(0.2)) 
        model.add(Dense(512, activation = 'relu')) 
        model.add(Dropout(0.2)) 
        model.add(Dense(10, activation = 'softmax'))
        model.compile(loss = 'categorical_crossentropy', 
        optimizer = RMSprop(), 
        metrics = ['accuracy']) 

        history = model.fit(x_train, y_train, 
        batch_size = 128, epochs = 20, verbose = 1, validation_data = (x_test, y_test))
