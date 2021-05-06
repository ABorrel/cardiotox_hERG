import toolbox
import pathFolder
import searchInComptox
import runExternal

from os import path
from random import shuffle
import CompDesc
import rdkit


class genericTestSet:
    def __init__(self, p_dataset, pr_out, pr_result, sep = ","):
        self.p_dataset = p_dataset
        self.pr_out = pr_out
        self.pr_results = pr_result
        self.sep = ","

    def main(self, allAff):

        self.loadDataset()
        self.computeDesc(self.pr_results + "DESC/")
        self.combineDesc()
        self.setAff(allAff=allAff)

    def loadDataset(self, loadDb="0", p_mapFile=""):

        p_dataset_preproc = self.pr_out + "dataset_prepoc.csv"
        if path.exists(p_dataset_preproc) and path.getsize(p_dataset_preproc) > 100:
            self.d_dataset = toolbox.loadMatrix(p_dataset_preproc)
            self.p_dataset_preproc = p_dataset_preproc
            return 

        d_dataset = toolbox.loadMatrix(self.p_dataset, sep = self.sep)

        # load CHEMBL mapping in case of CHEMBL SET
        if p_mapFile != "": # from pubchem
            d_map = toolbox.loadMatrix(p_mapFile, sep = ",")

        # load mapping CASRN
        if path.exists(loadDb) and loadDb != 1:# case of batch search
            d_database = toolbox.loadMatrix(loadDb, sep = ",")

        l_h = list(d_dataset[list(d_dataset.keys())[0]].keys())
        try:l_h.remove("CASRN")
        except:pass
        try:l_h.remove("SMILES")
        except:pass

        f_dataset_cleanned = open(p_dataset_preproc, "w")
        f_dataset_cleanned.write("CASRN\tSMILES\t" + "\t".join(l_h) + "\n")

        for chem in d_dataset.keys():
            if "CASRN" in list(d_dataset[chem].keys()):
                CASRN = d_dataset[chem]["CASRN"]
            else:
                if "INPUT" in list(d_dataset[chem].keys()):
                    CASRN = d_dataset[chem]["INPUT"]
                elif "PUBCHEM_SID" in list(d_dataset[chem].keys()):
                    PUBCHEM_SID = d_dataset[chem]["PUBCHEM_SID"]
                    flag=0
                    for ID_sample in d_map.keys():
                        PUBCHEM_SID_map = d_map[ID_sample]["PUBCHEM_SID"]
                        if PUBCHEM_SID == PUBCHEM_SID_map:
                            CASRN = d_map[ID_sample]["CAS"]
                            flag = 1
                            break
                    if flag == 0:
                        continue
                else:
                    continue
            
            if "SMILES" in list(d_dataset[chem].keys()):
                SMILES = d_dataset[chem]["SMILES"]
            else:
                if loadDb == 1:
                    c_search = searchInComptox.loadComptox(CASRN)
                    c_search.searchInDB()
                    if c_search.err != 1:
                        SMILES = c_search.SMILES
                    else:
                        SMILES = "--"

                elif path.exists(loadDb):# case of batch search
                    try:SMILES = d_database[CASRN]["SMILES"]
                    except:SMILES = "--"

            print(SMILES)

            f_dataset_cleanned.write("%s\t%s\t%s\n"%(CASRN, SMILES, "\t".join([d_dataset[chem][h] for h in l_h])))
        f_dataset_cleanned.close()

        self.d_dataset = toolbox.loadMatrix(p_dataset_preproc)
        self.p_dataset_preproc = p_dataset_preproc

    def computeDesc(self, pr_desc):

        p_filout_RDKIT = self.pr_out + "desc_1D2D.csv"
        p_filout_OPERA = self.pr_out + "desc_OPERA.csv"
        if path.exists(p_filout_RDKIT) and path.exists(p_filout_OPERA):
            self.p_desc1D2D = p_filout_RDKIT
            self.p_desc_opera = p_filout_OPERA
            return 

        # extract descriptor 2D
        if not path.exists(p_filout_RDKIT):
            cChem = CompDesc.CompDesc("", pr_desc)
            l_desc = cChem.getLdesc("1D2D")

            # open filout
            filout = open(p_filout_RDKIT, "w")
            filout.write("CASRN\tSMILES\t%s\n"%("\t".join(l_desc)))

            # compute descriptor
            l_smi = []
            for CASRN in self.d_dataset.keys():
                SMILES = self.d_dataset[CASRN]["SMILES"]
                cChem = CompDesc.CompDesc(SMILES, pr_desc)
                cChem.prepChem() # prep
                # case error cleaning
                if cChem.err == 1:
                    continue
                cChem.computeAll2D() # compute
                cChem.writeMatrix("2D") # write by chem to save time in case of rerun
                if cChem.err == 1:
                    continue
                else:
                    # write direcly descriptor
                    filout.write("%s\t%s\t%s\n"%(CASRN, cChem.smi, "\t".join(["NA" if not desc in list(cChem.all2D.keys()) else str(cChem.all2D[desc])  for desc in l_desc])))
                    l_smi.append(cChem.smi)
            filout.close()

        self.p_desc1D2D = p_filout_RDKIT

        if not path.exists(p_filout_OPERA):
            self.computeDescOPERA(pr_desc)
        else:
            self.p_desc_opera = p_filout_OPERA

    def computeDescOPERA(self, pr_desc):

        if not "pr_desc" in self.__dict__:
            self.pr_desc = pr_desc

        p_filout = self.pr_out + "desc_OPERA.csv"
        if path.exists(p_filout):
            self.p_descOPERA = p_filout
            return p_filout

        # define LDESCOPERA
        L_OPERA_DESC = CompDesc.CompDesc("", "").getLdesc("OPERA")

        # write list of SMILES for OPERA
        pr_OPERA = pathFolder.createFolder(self.pr_desc + "OPERA/", clean=1)
        p_lSMI = pr_OPERA + "listChem.smi"
        flSMI = open(p_lSMI, "w")
        l_w = []
        for CASRN in self.d_dataset.keys():
            SMILES =  self.d_dataset[CASRN]["SMILES"]
            l_w.append(SMILES)
        
        flSMI.write("\n".join(l_w))
        flSMI.close()

        p_desc_opera = pr_OPERA + "desc_opera_run.csv"
        cCompDesc = CompDesc.CompDesc(p_lSMI, pr_OPERA, update=1)
        cCompDesc.computeOperaFromListSMI(p_desc_opera)

        l_chem = list(self.d_dataset.keys())
        l_ddesc_run = toolbox.loadMatrixToList(p_desc_opera, sep = ",")

        fopera = open(p_filout, "w")
        fopera.write("CASRN,%s\n"%(",".join(L_OPERA_DESC)))

        i = 0
        imax = len(l_chem) 
        while i < imax:
            CASRN = l_chem[i]
            for ddesc_run in l_ddesc_run:
                if ddesc_run["MoleculeID"] == "Molecule_%i"%(i+1):
                    fopera.write("%s,%s\n"%(CASRN, ",".join(ddesc_run[desc] for desc in L_OPERA_DESC)))
                    break
            i = i + 1
        fopera.close()

        self.p_desc_opera = p_filout

    def combineDesc(self):

        p_global = self.pr_out + 'desc_global.csv'
        if not path.exists(p_global):
            runExternal.combineAndPredDesc(self.p_desc1D2D, self.p_desc_opera, self.pr_out)
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
                if self.d_dataset[chem][allAff] == "Active" or self.d_dataset[chem][allAff] == "1":
                    aff = 1
                else:
                    aff = 0
                AC50 = self.d_dataset[chem]["Fit_LogAC50"]
                if AC50 != "":
                    filout.write("\"%s\",\"%s\",\"%s\"\n"%(self.d_dataset[chem]["CASRN"], self.d_dataset[chem]["Fit_LogAC50"], aff))
                else:
                    filout.write("\"%s\",\"NA\",\"%s\"\n"%(self.d_dataset[chem]["CASRN"], aff))
            filout.close()
        elif type(allAff) == int:
            filout = open(p_aff, "w")
            filout.write("\"\",\"ID\",\"Aff\"\n")
            for chem in d_desc.keys():
                filout.write("\"%s\",\"%s\",%s\n"%(chem, chem, allAff))
            filout.close()

        else:
            print(allAff)
            d_dataset = allAff
            filout = open(p_aff, "w")
            filout.write("\"ID\",\"LogAC50\",\"Aff\"\n")
            for chem in d_dataset.keys():
                if not "CASRN" in list(d_dataset[chem].keys()):
                    continue
                if d_dataset[chem]["Aff"] == "1":
                    aff = 1
                else:
                    aff = 0
                AC50 = d_dataset[chem]["-log10(AC50)"]
                if AC50 != "":
                    filout.write("\"%s\",\"%s\",\"%s\"\n"%(d_dataset[chem]["CASRN"], d_dataset[chem]["-log10(AC50)"], aff))
                else:
                    filout.write("\"%s\",\"NA\",\"%s\"\n"%(d_dataset[chem]["CASRN"], aff))
            filout.close()
        self.p_aff = p_aff

    def computeDescForNNComparison(self, pr_desc, pr_out):

        l_desc = ["SlogP","SMR","LabuteASA","TPSA","AMW","ExactMW","NumLipinskiHBA","NumLipinskiHBD","NumRotatableBonds","NumHBD","NumHBA","NumAmideBonds","NumHeteroAtoms","NumHeavyAtoms","NumAtoms","NumStereocenters","NumUnspecifiedStereocenters","NumRings","NumAromaticRings","NumSaturatedRings","NumAliphaticRings","NumAromaticHeterocycles","NumSaturatedHeterocycles","NumAliphaticHeterocycles","NumAromaticCarbocycles","NumSaturatedCarbocycles","NumAliphaticCarbocycles","FractionCSP3","Chi0v","Chi1v","Chi2v","Chi3v","Chi4v","Chi1n","Chi2n","Chi3n","Chi4n","HallKierAlpha","kappa1","kappa2","kappa3","slogp_VSA1","slogp_VSA2","slogp_VSA3","slogp_VSA4","slogp_VSA5","slogp_VSA6","slogp_VSA7","slogp_VSA8","slogp_VSA9","slogp_VSA10","slogp_VSA11","slogp_VSA12","smr_VSA1","smr_VSA2","smr_VSA3","smr_VSA4","smr_VSA5","smr_VSA6","smr_VSA7","smr_VSA8","smr_VSA9","smr_VSA10","peoe_VSA1","peoe_VSA2","peoe_VSA3","peoe_VSA4","peoe_VSA5","peoe_VSA6","peoe_VSA7","peoe_VSA8","peoe_VSA9","peoe_VSA10","peoe_VSA11","peoe_VSA12","peoe_VSA13","peoe_VSA14","MQN1","MQN2","MQN3","MQN4","MQN5","MQN6","MQN7","MQN8","MQN9","MQN10","MQN11","MQN12","MQN13","MQN14","MQN15","MQN16","MQN17","MQN18","MQN19","MQN20","MQN21","MQN22","MQN23","MQN24","MQN25","MQN26","MQN27","MQN28","MQN29","MQN30","MQN31","MQN32","MQN33","MQN34","MQN35","MQN36","MQN37","MQN38","MQN39","MQN40","MQN41","MQN42"] 
        # open filout
        p_filout = pr_out + "chemicals_knime_desc.csv"
        filout = open(p_filout, "w")
        filout.write("smiles,%s,%s\n"%(",".join(l_desc), ",".join(["Bit %i"%(i+1) for i in range(0, 1024)])))

        # compute descriptor
        for CASRN in list(self.d_dataset.keys()):
            SMILES = self.d_dataset[CASRN]["SMILES"]

            # define class from CompDESC
            cChem = CompDesc.CompDesc(SMILES, pr_desc)
            cChem.prepChem()
                
            # compute desc
            cChem.update = 1 # to recompute MQNs
            cChem.computeAll2D()
            if cChem.err == 1:
                continue
            
            cChem.convertDesc2DtoKnimeDesc()

            l_FpMorgan = rdkit.Chem.AllChem.GetMorganFingerprintAsBitVect(cChem.mol, radius=2, nBits=1024)
            filout.write("%s,%s,%s\n"%(SMILES, ",".join([str(cChem.all2D[desc]) for desc in l_desc]), ",".join([str(l_FpMorgan[i]) for i in range(0,1024)])))
        filout.close()

        return p_filout

