from os import system, path, remove, chdir, getcwd, listdir



P_RSCRIPTS = "../R/"
P_RQSAR = "./../../../development/QSAR-QSPR/"


######
# Main functions

def runRCMD(cmd, out = 0):

    workdir = getcwd()
    chdir(P_RSCRIPTS)
    wrkdir = getcwd()
    print(cmd)
    if out == 0:
        system(cmd)
        output = 0
    else:
        import subprocess
        output = subprocess.check_output(cmd, shell=True)
    chdir(wrkdir)
    return output


def runRQSARModeling(cmd):

    workdir = getcwd()
    chdir(P_RQSAR)
    print(cmd)
    system(cmd)
    chdir(workdir)



############
# Functions for analysis


def preprocData(p_desc, p_AC50, pr_out, cor_val, max_quantile):

    cmd = "./cleanerData.R %s %s %s %s %s"%(p_desc, p_AC50, pr_out, cor_val, max_quantile)
    runRCMD(cmd)


def histAC50(p_AC50, pr_out):

    cmd = "./sumAC50.R %s %s "%(p_AC50, pr_out)
    runRCMD(cmd)


def PCA(p_desc_cleaned, p_AC50, pr_out):

    cmd = "./PCA_chem.R %s %s %s"%(p_desc_cleaned, p_AC50, pr_out)
    runRCMD(cmd)


def SOM(p_desc_cleaned, p_AC50_cleaned, pr_out, grid_size):

    cmd = "./SOM_chem.R %s %s %s %s"%(p_desc_cleaned, p_AC50_cleaned, pr_out, grid_size)
    runRCMD(cmd)


def HClust(p_desc_cleaned, p_AC50_cleaned, pr_out):

    cmd = "./HClust_chem.R %s %s %s"%(p_desc_cleaned, p_AC50_cleaned, pr_out)
    runRCMD(cmd)

def descSignifByCluster(p_desc, p_cluster, pr_out):

    cmd = "./descSignificantByCluster.R %s %s %s"%(p_desc, p_cluster, pr_out)
    runRCMD(cmd)

def PCAvs(p_desc_model, p_aff_model, p_desc_test, p_aff_test, pr_out):

    cmd = "./PCAvs.R %s %s %s %s %s"%(p_desc_model, p_aff_model, p_desc_test, p_aff_test, pr_out)
    runRCMD(cmd)    


############
# Function for QSAR

def SplitTrainTest(pdescAc50, prout, splitratio):

    pdescAc50 = path.abspath(pdescAc50)
    prout = path.abspath(prout) + "/"

    cmd = "./prepTrainTestSplitClass.R " + pdescAc50 + " " + str(splitratio) + " " + prout
    runRQSARModeling(cmd)


def prepDataQSAR(p_desc, p_AC50, rate_active, pr_run):

    p_desc = path.abspath(p_desc)
    p_AC50 = path.abspath(p_AC50)
    pr_run = path.abspath(pr_run) + "/"

    cmd = "./prepClassDataset.R %s %s %s %s"%(p_desc, p_AC50, rate_active, pr_run)
    runRQSARModeling(cmd)


def runRQSAR(p_train, p_test, n_foldCV, pr_run):

    p_train = path.abspath(p_train)
    p_test = path.abspath(p_test)
    pr_run = path.abspath(pr_run) + "/"

    cmd = "./QSARsClass.R " + p_train + " " + p_test + " 0 " + pr_run + " " + str(n_foldCV) + " > " + pr_run + "perf.txt"
    runRQSARModeling(cmd)


def plotAC50VSProb(p_prob):

    cmd = "./plotAC50vsProb.R %s"%(p_prob)
    runRCMD(cmd)

def runImportanceDesc(p_desc, nb):
    
    p_desc = path.abspath(p_desc)

    cmd = "./importancePlot.R " + str(p_desc) + " " + str(nb)
    runRQSARModeling(cmd)


def predictDataset(p_desc_test, p_model, ML,  pr_out):

    cmd = "./predictTestSet.R %s %s %s %s"%(p_desc_test, p_model, ML, pr_out)
    runRCMD(cmd)