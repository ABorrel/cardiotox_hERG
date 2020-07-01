import toolbox
import pathFolder
from os import path

# import descriptor computation scripts => precise folder where descriptor are included
import sys
sys.path.insert(0, "./../../../development/molecular-descriptors/")
import Chemical



class genericTestSet:
    def __init__(self, p_dataset, pr_out, sep = ","):
        self.p_dataset = p_dataset
        self.pr_out = pr_out
        self.sep = ","

    def loadDataset(self):

        d_dataset = toolbox.loadMatrix(self.p_dataset, sep = self.sep)
        self.d_dataset = d_dataset


    def computeDesc(self):

        if not "d_dataset" in self.__dict__:
            self.loadDataset()
        
        pr_out = pathFolder.createFolder(self.pr_out + "DESC/")
        self.pr_desc = pr_out

        p_filout = pr_out + "desc_1D2D.csv"
        if path.exists(p_filout):
            self.p_desc = p_filout
            return p_filout

        # extract descriptor 2D
        l_desc = Chemical.getLdesc("1D2D")

        # open filout
        filout = open(p_filout, "w")
        filout.write("CASRN\tSMILES\t%s\n"%("\t".join(l_desc)))

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
        self.p_desc = p_filout
        return p_filout
        
        
    def setAff(self, allAff=0):
        """
        Add here part to open aff file and formate it
        """
        
        if not "p_desc" in self.__dict__:
            self.computeDesc()

        p_aff = self.pr_out + "aff_formated.csv"
        if path.exists(p_aff):
            return p_aff
        
        d_desc = toolbox.loadMatrix(self.p_desc, sep = "\t")

        filout = open(p_aff, "w")
        filout.write("\"\",\"ID\",\"Aff\"\n")
        for chem in d_desc.keys():
            filout.write("\"%s\",\"%s\",%s\n"%(chem, chem, allAff))
        filout.close()

        return p_aff
