from os import path
import toolbox

# import descriptor computation scripts => precise folder where descriptor are included
import sys
sys.path.insert(0, "C:/Users/Aborrel/research/NIEHS/development/descriptor/")
import Chemical



class dataset:
    def __init__(self, p_smi, p_AC50, pr_out):
        self.p_smi = p_smi
        self.p_AC50 = p_AC50
        self.pr_out = pr_out
    

    def prep_dataset(self):

        # define output
        p_filout = self.pr_out + "dataset_prep.csv"
        if path.exists(p_filout):
            self.d_dataset = toolbox.loadMatrix(p_filout)
            return self.d_dataset

        d_chem = toolbox.loadMatrix(self.p_smi, sep = ",")
        d_AC50 = toolbox.loadMatrix(self.p_AC50)

        # rebuild table of chemicals
        d_out = {}
        for CASRN in d_AC50.keys():
            d_temp = {}
            d_temp["CASRN"] = CASRN
            d_temp["log10(AC50)"] = d_AC50[CASRN]["LogAC50"]

            # load on the tox21 chem library
            for DTXID in d_chem.keys():
                if d_chem[DTXID]["CASRN"] == CASRN:
                    d_temp["SMILES"] =  d_chem[DTXID]["SMILES"]
                    d_temp["PREFERRED_NAME"] =  d_chem[DTXID]["PREFERRED_NAME"]
                    d_out[CASRN] = d_temp
                    break
        
        # write the dataset
        l_header = ["CASRN", "SMILES", "PREFERRED_NAME", "log10(AC50)"]
        filout = open(p_filout, "w")
        filout.write("\t".join(l_header) + "\n")
        for CASRN in d_out.keys():
            filout.write("%s\n"%("\t".join([d_out[CASRN][h] for h in l_header])))
        filout.close()
        self.d_dataset = d_out

    
    def computeDesc(self, pr_desc):

        p_filout = pr_desc + "desc_1D2D.csv"
        if path.exists(p_filout):
            return p_filout

        # extract descriptor 2D
        l_desc = Chemical.getLdesc("1D2D")

        # open filout
        filout = open(p_filout, "w")
        filout.write("CASRN\tSMILES\t%s\n"%("\t".join(l_desc)))

        # compute descriptor
        for CASRN in self.d_dataset.keys():
            SMILES = self.d_dataset[CASRN]["SMILES"]
            cChem = Chemical.Chemical(SMILES, pr_desc, p_salts=path.abspath("./Salts.txt"))
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
                filout.write("%s\t%s\t%s\n"%(CASRN, cChem.smi, "\t".join([str(cChem.all2D[desc]) for desc in l_desc])))

        filout.close()

        return p_filout
