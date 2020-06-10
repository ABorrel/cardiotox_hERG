from os import system, path, remove, chdir, getcwd, listdir



P_RSCRIPTS = "../R/"
P_RQSAR = "/home/borrela2/development/QSARPR/source/"



######
# Main functions

def runRCMD(cmd, out = 0):

    workdir = getcwd()
    chdir(P_RSCRIPTS)
    wrkdir = getcwd()
    print(cmd)
    ddd
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


def SOM(p_desc_cleaned, p_AC50_cleaned, pr_out):

    cmd = "SOM_chem.R %s %s %s"%(p_desc_cleaned, p_AC50_cleaned, pr_out)
    runRCMD(cmd)

    return 
