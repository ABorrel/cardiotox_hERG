import toolbox
import pathFolder
import searchInComptox
import runExternal

from os import path
from random import shuffle
import CompDesc



class genericTestSet:
    def __init__(self, p_dataset, pr_out, sep = ","):
        self.p_dataset = p_dataset
        self.pr_out = pr_out
        self.sep = ","

    def loadDataset(self, loadDb=0, p_mapFile=""):

        p_dataset_preproc = self.pr_out + "dataset_prepoc.csv"
        if path.exists(p_dataset_preproc):
            self.d_dataset = toolbox.loadMatrix(p_dataset_preproc)
            self.p_dataset_preproc = p_dataset_preproc
            return 

        d_dataset = toolbox.loadMatrix(self.p_dataset, sep = self.sep)
        self.d_dataset = d_dataset

        l_h = list(self.d_dataset[list(self.d_dataset.keys())[0]].keys())
        #l_h.append("CARSN")
        #l_h.append("SMILES")

        if p_mapFile != "": # from pubchem
            d_map = toolbox.loadMatrix(p_mapFile, sep = ",")
            for ID_sample in d_map.keys():
                PUBCHEM_SID = d_map[ID_sample]["PUBCHEM_SID"]
                CASRN = d_map[ID_sample]["CAS"]

                for chem in self.d_dataset.keys():
                    if self.d_dataset[chem]["PUBCHEM_SID"] == PUBCHEM_SID:
                        self.d_dataset[chem]["CASRN"] = CASRN
                        break


        if loadDb == 1:
            # load SMILES from comptox

            if not "SMILES" in l_h:
                for chem in d_dataset.keys():
                    if not "CASRN" in list(d_dataset[chem].keys()):
                        continue
                    CASRN = d_dataset[chem]["CASRN"]
                    CASRN = CASRN.split("|")[0]
                    print("CASRN:", CASRN, "LOAD chem")
                    c_search = searchInComptox.loadComptox(CASRN)
                    c_search.searchInDB()
                    if c_search.err != 1:
                        self.d_dataset[chem]["SMILES"] = c_search.SMILES
                    else:
                        self.d_dataset[chem]["SMILES"] = "--"

                    print("Load done:", self.d_dataset[chem]["SMILES"])


        elif path.exists(loadDb):# case of batch search
            d_database = toolbox.loadMatrix(loadDb, sep = ",")
            for chem in self.d_dataset.keys():
                if not "CASRN" in list(d_dataset[chem].keys()) or self.d_dataset[chem]["CASRN"] == "":
                    continue
                CASRN = self.d_dataset[chem]["CASRN"]
                if CASRN in list(d_database.keys()):
                    self.d_dataset[chem]["SMILES"] = d_database[CASRN]["SMILES"]
                else:
                    self.d_dataset[chem]["SMILES"] = "--"

        else:
            if not "SMILES" in l_h:
                for chem in self.d_dataset.keys():
                    if not "CASRN" in list(d_dataset[chem].keys()):
                        continue
                    CASRN = self.d_dataset[chem]["CASRN"]
                    CASRN = CASRN.split("|")[0]
                    self.d_dataset[chem]["SMILES"] = "--"


        f_dataset_cleanned = open(p_dataset_preproc, "w")
        f_dataset_cleanned.write("CASRN\tSMILES\t" + "\t".join(l_h) + "\n")
        for chem in self.d_dataset.keys():
            if "CASRN" in list(self.d_dataset[chem].keys()) and self.d_dataset[chem]["CASRN"] != "":
                f_dataset_cleanned.write("%s\t%s\t%s\n"%(self.d_dataset[chem]["CASRN"], self.d_dataset[chem]["SMILES"], "\t".join([self.d_dataset[chem][h] for h in l_h])))
        f_dataset_cleanned.close()

        self.d_dataset = toolbox.loadMatrix(p_dataset_preproc)
        self.p_dataset_preproc = p_dataset_preproc



    def computeDesc(self):

        if not "d_dataset" in self.__dict__:
            self.loadDataset()

        pr_out = pathFolder.createFolder(self.pr_out + "DESC/")
        self.pr_desc = pr_out

        p_filout_1D2D = pr_out + "desc_1D2D.csv"
        p_filout_OPERA = pr_out + "desc_OPERA.csv"
        if path.exists(p_filout_1D2D) and path.exists(p_filout_OPERA):
            self.p_desc = p_filout_1D2D
            self.p_desc_opera = p_filout_OPERA

        else:
            # 1D2D ================
            # extract descriptor 2D
            if not path.exists(p_filout_1D2D):
                l_desc = Chemical.getLdesc("1D2D")
                l_smi = []
                # compute descriptor
                for chem in self.d_dataset.keys():
                    CASRN = self.d_dataset[chem]["CASRN"]
                    SMILES = self.d_dataset[chem]["SMILES"]

                    cChem = Chemical.Chemical(SMILES, self.pr_desc, p_salts=path.abspath("./Salts.txt"))
                    cChem.prepChem() # prep
                    # case error cleaning
                    if cChem.err == 1:
                        continue
                    if cChem.smi in l_smi:
                        continue
                    else:
                        l_smi.append(cChem.smi)
                    cChem.computeAll2D() # compute
                    cChem.writeMatrix("2D") # write by chem to save time in case of rerun
                    if cChem.err == 1:
                        continue
                    else:
                        # write direcly descriptor
                        filout.write("%s\t%s\t%s\n"%(CASRN, cChem.smi, "\t".join([str(cChem.all2D[desc]) for desc in l_desc])))
                filout.close()
            self.p_desc = p_filout_1D2D


            ### OPERA =============
            pr_OPERA = pathFolder.createFolder(pr_out + "OPERA/")
            p_listchem = pr_OPERA + "listChem.smi"
            f_listchem = open(p_listchem, "w")
            d_desc_1D2D = toolbox.loadMatrix(p_filout_1D2D)
            for chem in d_desc_1D2D.keys():
                f_listchem.write(d_desc_1D2D[chem]["SMILES"] + "\n")
            f_listchem.close()

            if path.exists(pr_OPERA + "desc_OPERA.csv"):
                l_chem_opera = toolbox.loadMatrixToList(pr_OPERA + "desc_OPERA.csv", sep=",")
                l_chem_rdkit = toolbox.loadMatrixToList(p_filout_1D2D, sep = "\t")
                print(len(l_chem_rdkit))
                print(len(l_chem_opera))
                l_header_opera = list(l_chem_opera[0].keys())
                l_header_opera.remove("MoleculeID")
                i = 0
                imax = len(l_chem_opera)
                while i < imax:
                    i_rdkit = int(l_chem_opera[i]["MoleculeID"].split("_")[-1]) - 1
                    casrn = l_chem_rdkit[i_rdkit]["CASRN"]
                    l_chem_opera[i]["CASRN"] = casrn
                    i = i + 1
                
                # write 
                fopera = open(p_filout_OPERA, "w")
                fopera.write("CASRN,%s\n"%(",".join(l_header_opera)))
                for chem in l_chem_opera:
                    fopera.write("%s,%s\n"%(chem["CASRN"], ",".join([str(chem[h]) for h in l_header_opera])))
                fopera.close()
                
                self.p_desc_opera = p_filout_OPERA
            else:
                print("RUN opera with => ", p_listchem)


    def combineDesc(self):

        p_global = self.pr_out + 'desc_global.csv'
        if not path.exists(p_global):
            runExternal.combineAndPredDesc(self.p_desc, self.p_desc_opera, self.pr_out)
        self.p_desc = p_global


    def setAff(self, allAff=0):
        """
        Add here part to open aff file and formate it
        """

        if not "p_desc" in self.__dict__:
            self.computeDesc()

        p_aff = self.pr_out + "aff_formated.csv"
        #if path.exists(p_aff):
        #    return p_aff

        d_desc = toolbox.loadMatrix(self.p_desc, sep = "\t")


        if allAff == "PUBCHEM_ACTIVITY_OUTCOME": # in this case is pubchem format
            filout = open(p_aff, "w")
            filout.write("\"ID\",\"LogAC50\",\"Aff\"\n")
            for chem in self.d_dataset.keys():
                if not "CASRN" in list(self.d_dataset[chem].keys()):
                    continue
                if self.d_dataset[chem][allAff] == "Active":
                    aff = 1
                else:
                    aff = 0
                AC50 = self.d_dataset[chem]["Fit_LogAC50"]
                if AC50 != "":
                    filout.write("\"%s\",\"%s\",\"%s\"\n"%(self.d_dataset[chem]["CASRN"], self.d_dataset[chem]["Fit_LogAC50"], aff))
                else:
                    filout.write("\"%s\",\"NA\",\"%s\"\n"%(self.d_dataset[chem]["CASRN"], aff))
            filout.close()
        else:
            filout = open(p_aff, "w")
            filout.write("\"\",\"ID\",\"Aff\"\n")
            for chem in d_desc.keys():
                filout.write("\"%s\",\"%s\",%s\n"%(chem, chem, allAff))
            filout.close()

        self.p_aff = p_aff
