import toolbox
import QSAR_modeling

from numpy import loadtxt, arange
from sklearn import metrics
from sklearn.model_selection import StratifiedKFold
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
    def __init__(self, pr_out, p_desc, p_aff,  p_desc_clean, p_aff_clean, nb_repetition, n_foldCV, rate_active, rate_splitTrainTest):
        self.pr_out = pr_out
        self.p_desc = p_desc
        self.p_aff = p_aff
        self.p_aff_clean = p_aff_clean
        self.p_desc_clean = p_desc_clean
        self.nb_repetition = nb_repetition
        self.n_foldCV = n_foldCV
        self.rate_active = rate_active
        self.rate_splitTrainTest = rate_splitTrainTest

        # optimization
        self.l_epochs = [50, 100, 120]
        self.l_batch_size = [32, 64, 128]
        self.l_dense_layer = [3, 4, 5]
        self.l_dense_candidate = [50, 25, 20, 10, 1]
        self.l_activation = ["relu", "selu"]

    def prepDataset(self):

        # split train test
        self.cQSAR = QSAR_modeling.QSAR_modeling(self.p_desc_clean, self.p_desc, self.p_aff_clean, self.p_aff, self.pr_out, self.nb_repetition, self.n_foldCV, self.rate_active, self.rate_splitTrainTest)
        self.cQSAR.prepTrainSetforUnderSampling(self.pr_out)      

        # define self.cQSAR.p_train // self.cQSAR.p_test
        
        # descriptor train set 
        with open(self.cQSAR.p_train) as f:
            #determining number of columns from the first line of text
            n_cols = len(f.readline().split(","))
        self.dataset_train = loadtxt(self.cQSAR.p_train, delimiter=",",usecols=arange(1, n_cols-1), skiprows=1)
        self.aff_train = loadtxt(self.cQSAR.p_train, delimiter=",",usecols=arange(n_cols-1, n_cols), skiprows=1)
        self.nb_desc_input = n_cols - 2

        # descriptor test set 
        with open(self.cQSAR.p_test) as f:
            #determining number of columns from the first line of text
            n_cols = len(f.readline().split(","))
        self.dataset_test = loadtxt(self.cQSAR.p_test, delimiter=",",usecols=arange(1, n_cols-1), skiprows=1)
        self.aff_test = loadtxt(self.cQSAR.p_test, delimiter=",",usecols=arange(n_cols-1, n_cols), skiprows=1)

    #l_batch_size = [32,64,128]
    def GridOptimizeDNN(self, criteria_best_perf, l_activation = ["relu", "selu"], l_epochs = [50, 100, 120], l_batch_size = [32, 64, 128], l_dense_layer = [3, 4, 5], l_dense_candidate = [50, 25, 20, 10, 1]):

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
                        l_pref_train = self.fit_compileModel(epochs, batch_size, optimizer = "adam")
                        l_pref_test = self.evaluateOnTest()
                        filout.write("Activation: %s; Nb layers: %s; Epochs: %s; Batch size: %s;\nTRAINNING SET\nAccuracy: %s\nAccuracy-balanced: %s\nMCC: %s\nRecall: %s\nTEST SET\nAccuracy: %s\nAccuracy-balanced: %s\nMCC: %s\nRecall: %s\n[====================]\n"%(activation, dense_layer, epochs, batch_size, l_pref_train[0], l_pref_train[1], l_pref_train[2], l_pref_train[3], l_pref_test[0], l_pref_test[1], l_pref_test[2], l_pref_test[3]))
                        k_out = "%s--%s--%s--%s"%(activation, dense_layer, epochs, batch_size)
                        d_out[k_out] = {}
                        d_out[k_out]["Acc"] = l_pref_train[0]
                        d_out[k_out]["b-Acc"] = l_pref_train[1]
                        d_out[k_out]["Recall"] = l_pref_train[3]
                        d_out[k_out]["MCC"] = l_pref_train[2]

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
        filout.close()                

    def addDense(self, activation, dense_candidate):
        self.model.add(Dense(int(dense_candidate), activation=activation, kernel_initializer='random_normal',))
        #self.model.add(Dropout(0.5))

    def initModel(self, activation, dense_candidate):
        self.model = Sequential()
        if dense_candidate > self.nb_desc_input:
            dense_candidate = self.nb_desc_input
        self.model.add(Dense(dense_candidate, input_dim=self.nb_desc_input, activation=activation, kernel_initializer='random_normal'))
        #self.model.add(Dropout(0.5))
    
    def fit_compileModel(self, epochs, batch_size, loss='binary_crossentropy', optimizer = "adam", l_metric = ["mse"], learning_rate = 0.1):

        # last layer
        self.model.add(Dense(1, activation='sigmoid'))

        if optimizer == "sgd":
            sgd = tf.keras.optimizers.SGD(learning_rate=learning_rate)
            optimizer = sgd
        #elif optimizer == "adam":
        #    adam =  tf.keras.optimizers.Adam(learning_rate=learning_rate)
        #    optimizer = adam
        #    print(adam)

        # compile the keras model
        self.model.compile(loss=loss, optimizer=optimizer, metrics=l_metric)
        # fit the keras model on the dataset
        self.model.fit(self.dataset_train, self.aff_train, epochs=epochs, batch_size=batch_size), 

        # evaluate the keras model
        self.model.evaluate(self.dataset_train, self.aff_train)
        

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

    def evaluateOnTest(self):

        y_pred = self.model.predict(self.dataset_test)
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
    
    def evaluateModel(self):

        if "p_model" in self.__dict__:
            model = keras.models.load_model(self.p_model)
        else:
            print("No model ready to be loaded")
            return

        if not "d_perf" in self.__dict__:
            self.d_perf = {}

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

    def combineResults(self):
        if not "d_perf" in self.__dict__:
            print("No data to write")
            return 

        p_filout = self.pr_out + "combined_perf.csv"
        filout = open(p_filout, "w")
        l_perf = ["CV", "Train", "Test"]
        for perf in l_perf:
            if perf in list(self.d_perf.keys()):
                filout.write("%s\n\tAcc\tb-Acc\tSp\tSe\tMCC\nDNN\t%s\t%s\t%s\t%s\t%s\n\n"%(perf, self.d_perf[perf]["Acc"], self.d_perf[perf]["b-Acc"], self.d_perf[perf]["Sp"], self.d_perf[perf]["Se"], self.d_perf[perf]["MCC"]))
        
        filout.close()

    def CrossValidation(self, n_folds):
        seed = 7
        kfold = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=seed)
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
        

        d_CV = {}
        d_CV["acc"] = []
        d_CV["b-acc"] = []
        d_CV["recall"] = []
        d_CV["mcc"] = []

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

            l_pref_train = self.fit_compileModel(epochs, batch_size, optimizer = "adam")
            l_pref_test = self.evaluateOnTest()
            
            # CV performance
            d_CV["acc"].append(l_pref_test[0])
            d_CV["b-acc"].append(l_pref_test[1])
            d_CV["mcc"].append(l_pref_test[2])
            d_CV["recall"].append(l_pref_test[3])


        filout = open(self.pr_out + "CV_%s.perf"%(n_folds), "w")
        filout.write("Acc\tB-acc\tMCC\tRecall\n")
        filout.write("%.2f+/-%.2f\t%.2f+/-%.2f\t%.2f+/-%.2f\t%.2f+/-%.2f\n"%(np.mean(d_CV["acc"]), np.std(d_CV["acc"]), np.mean(d_CV["b-acc"]), np.std(d_CV["b-acc"]),np.mean(d_CV["mcc"]), np.std(d_CV["mcc"]), np.mean(d_CV["recall"]), np.std(d_CV["recall"])))
        filout.close()

        if not "d_perf" in self.__dict__:
            self.d_perf = {}

        self.d_perf['CV'] = {}
        self.d_perf["CV"]["Acc"] = np.mean(d_CV["acc"])
        self.d_perf["CV"]["b-Acc"] =  np.mean(d_CV["b-acc"])
        self.d_perf["CV"]["MCC"] = np.mean(d_CV["mcc"])
        self.d_perf["CV"]["Se"] = np.mean(d_CV["recall"])
        self.d_perf["CV"]["Sp"] = (self.d_perf["CV"]["b-Acc"] * 2) - self.d_perf["CV"]["Se"]
  
    def predictDNN(self, p_model, p_test):
        # descriptor test set 
        with open(p_test) as f:
            #determining number of columns from the first line of text
            n_cols = len(f.readline().split(","))
        dataset_test = loadtxt(p_test, delimiter=",",usecols=arange(1, n_cols-1), skiprows=1)

        model = keras.models.load_model(p_model)
        y_pred = model.predict(dataset_test)

        return y_pred


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
