from copy import deepcopy
from numpy import mean, std
from os import path
from re import search

import toolbox
import runExternal

# import descriptor computation scripts => precise folder where descriptor are included
import sys
sys.path.insert(0, "./../../../development/molecular-descriptors/")
import Chemical


class CHEMBLTable:
    def __init__(self, p_filin, pr_out, verbose = 1):
        self.p_filin = p_filin
        self.pr_out = pr_out
        self.verbose = 1

    def parseCHEMBLFile(self):

        l_out = []
        filin = open(self.p_filin, "r")
        l_lines = filin.readlines()
        filin.close()

        l_head = l_lines[0].strip().replace("\"", "")
        l_head = l_head.split(";")

        # to key same order for print
        self.l_header = l_head

        i = 1
        imax = len(l_lines)
        while i < imax:
            line_temp = toolbox.formatLineDataset(l_lines[i])
            l_elem = line_temp.split(";")


            d_out = {}
            j = 0
            while j < len(l_head):
                d_out[l_head[j]] = l_elem[j]
                j += 1
            l_out.append(d_out)
            i += 1

        if self.verbose == 1:
            print("Load ChEMBL table - %s lines loaded"%(len(l_out)))

        self.l_ChEMBLin = l_out
        self.l_work = deepcopy(l_out)        

    def cleanDataset(self, l_standard_type, l_standard_relation):

        # load existing data if exist
        p_filout = self.pr_out + "ChEMBL_cleaned.csv"
        if path.exists(p_filout):
            self.l_work = toolbox.loadMatrixToList(p_filout, sep="\t")
        else:
            if not "l_work" in dir(self):
                self.parseCHEMBLFile()

            self.get_standard_type(l_standard_type)
            self.get_standard_relation(l_standard_relation)
            self.filterDuplicates()

            self.writeTable()

    def filterOnDescription(self, search_str):
        
        if not "l_work" in dir(self):
            self.parseCHEMBLFile()

        i = 0
        i_max = len(self.l_work) 

        while i < i_max:
            desc = self.l_work[i]["Assay Description"]
            desc = desc.lower()

            if not search(search_str, desc):
                del self.l_work[i]
                i_max = i_max - 1
            else:
                i = i + 1
            


    def get_standard_relation(self, l_relation):
        """Keep only value with a strict standard relation"""

        i = 0
        imax = len(self.l_work)
        if self.verbose == 1:
            print("Filter by relation => [%s]"%(";".join(l_relation)))
            print("Nb value => %s"%(imax))
        while i < imax:
            row = self.l_work[i]
            if not row["Standard Relation"] in l_relation:
                del self.l_work[i]
                imax = imax - 1
                continue
            else:
                i += 1

        if self.verbose == 1:
            print("After filtering by relation => %s"%(len(self.l_work)))

    def get_standard_type(self, l_type):

        i = 0
        imax = len(self.l_work)
        if self.verbose == 1:
            print("Filter by type => [%s]"%(";".join(l_type)))
            print("Nb value => %s"%(imax))
        while i < imax:
            row = self.l_work[i]
            if not row["Standard Type"] in l_type:
                del self.l_work[i]
                imax = imax - 1
                continue
            else:
                i += 1

        if self.verbose == 1:
            print("After filtering by type => %s"%(len(self.l_work)))

    def filterDuplicates(self):

        d_CHEMBLID = {}
        for row in self.l_work:
            if not row["Molecule ChEMBL ID"] in list(d_CHEMBLID.keys()):
                d_CHEMBLID[row["Molecule ChEMBL ID"]] = []
            d_CHEMBLID[row["Molecule ChEMBL ID"]].append(deepcopy(row))

        if self.verbose == 1:
            print("Nb duplicate => %s"%(len(self.l_work) - len(list(d_CHEMBLID.keys()))))


        # case no copy
        l_CHEMBLid =list(d_CHEMBLID.keys())
        imax = len(l_CHEMBLid)
        if imax == len(self.l_work):
            return

        i = 0
        
        while i < imax:
            # case not problem
            #print(d_CHEMBLID.keys()[i])
            if len(list(d_CHEMBLID[l_CHEMBLid[i]])) == 1:
                i += 1
                continue
            else:
                # Control the published value
                # favorise ki to other affinity
                l_standard_type = [d_CHEMBLID[l_CHEMBLid[i]][k]["Standard Type"] for k in range(0, len(d_CHEMBLID[l_CHEMBLid[i]]))]
                if len(list(set(l_standard_type))) != 1:
                    i_type = 0
                    while i_type < len(l_standard_type):
                        if l_standard_type[i_type] != "Ki":
                            del l_standard_type[i_type]
                            del d_CHEMBLID[l_CHEMBLid[i]][i_type]
                        else:
                            i_type = i_type + 1 

                # remove also unit is not in molar
                l_standard_unit = [d_CHEMBLID[l_CHEMBLid[i]][k]["Standard Units"] for k in range(0, len(d_CHEMBLID[l_CHEMBLid[i]]))]
                if len(list(set(l_standard_unit))) != 1:
                    i_unit = 0
                    while i_unit < len(l_standard_unit):
                        if l_standard_unit[i_unit] == "ug.mL-1":
                            del l_standard_unit[i_unit]
                            del d_CHEMBLID[l_CHEMBLid[i]][i_unit]
                        else:
                            i_unit = i_unit + 1 




                l_PUBLISHED_VALUE = [float(d_CHEMBLID[l_CHEMBLid[i]][k]["Standard Value"]) for k in range(0, len(d_CHEMBLID[l_CHEMBLid[i]]))]
                l_PUBLISHED_UNITS = [d_CHEMBLID[l_CHEMBLid[i]][k]["Standard Units"] for k in range(0, len(d_CHEMBLID[l_CHEMBLid[i]]))]

                if not l_PUBLISHED_UNITS.count(l_PUBLISHED_UNITS[0]) == len(l_PUBLISHED_UNITS):
                    l_PUBLISHED_VALUE = toolbox.convertUnit(l_PUBLISHED_VALUE, l_PUBLISHED_UNITS)

                MPUBLISHED_VALUE = mean(l_PUBLISHED_VALUE)
                SDPUBLISHED_VALUE = std(l_PUBLISHED_VALUE)

                magnitudeval = len(str(int(min(l_PUBLISHED_VALUE))))
                magnitudeSD = len(str(int(SDPUBLISHED_VALUE)))

                if magnitudeval != magnitudeSD:
                    del d_CHEMBLID[l_CHEMBLid[i]]
                    del l_CHEMBLid[i]
                    imax = imax - 1
                    continue

                else:

                    l_STANDARD_VALUE = [float(d_CHEMBLID[l_CHEMBLid[i]][k]["Standard Value"]) for k in
                                        range(0, len(d_CHEMBLID[l_CHEMBLid[i]]))]
                    

                    try: 
                        l_PCHEMBL_VALUE = [float(d_CHEMBLID[l_CHEMBLid[i]][k]["pChEMBL Value"]) for k in
                                       range(0, len(d_CHEMBLID[l_CHEMBLid[i]]))]
                    except:
                        l_PCHEMBL_VALUE = []

                    MSTANDARD_VALUE = mean(l_STANDARD_VALUE)


                    d_CHEMBLID[l_CHEMBLid[i]] = [toolbox.mergeDict(d_CHEMBLID[l_CHEMBLid[i]])]
                    d_CHEMBLID[l_CHEMBLid[i]][0]["Standard Value"] = str(MSTANDARD_VALUE)
                    if l_PCHEMBL_VALUE == []:
                        d_CHEMBLID[l_CHEMBLid[i]][0]["pChEMBL Value"] = ""
                    else:
                        d_CHEMBLID[l_CHEMBLid[i]][0]["pChEMBL Value"] = str(mean(l_PCHEMBL_VALUE))
                    i += 1

        # reformate the table
        self.l_work = []
        for k in d_CHEMBLID.keys():
            self.l_work.append(d_CHEMBLID[k][0])

        if self.verbose:
            print("Nb after duplicate cleaning => %s"%(len(self.l_work)))
       
    def writeTable(self):

        filout = open(self.pr_out + "ChEMBL_cleaned.csv", "w")
        filout.write("\t".join(self.l_header) + "\n")
        for row in self.l_work:
            lw = [row[h] for h in self.l_header]
            filout.write("\t".join(lw) + "\n")
        filout.close()

    def computeDesc(self):

        p_filout = self.pr_out + "desc_1D2D.csv"
        if path.exists(p_filout):
            return p_filout

        # extract descriptor 2D
        l_desc = Chemical.getLdesc("1D2D")

        # open filout
        filout = open(p_filout, "w")
        filout.write("CHEMBLID\tSMILES\t%s\n"%("\t".join(l_desc)))

        l_smi = []
        # compute descriptor
        for d_chem in self.l_work:
            ChEMBLid = d_chem["Molecule ChEMBL ID"]
            SMILES = d_chem["Smiles"]

            cChem = Chemical.Chemical(SMILES, self.pr_out, p_salts=path.abspath("./Salts.txt"))
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
                filout.write("%s\t%s\t%s\n"%(ChEMBLid, cChem.smi, "\t".join([str(cChem.all2D[desc]) for desc in l_desc])))

        filout.close()

        return p_filout

    def prep_aff(self, typeAff="pAff", cutoff_uM=""):

        p_aff_cleaned = "%saff_%s_%s.csv"%(self.pr_out, typeAff, cutoff_uM)
        if path.exists(p_aff_cleaned):
            return p_aff_cleaned

        p_desc = self.computeDesc()
        d_desc = toolbox.loadMatrix(p_desc)

        l_chem = list(d_desc.keys())

        if typeAff == "pAff":
            f_aff_cleaned = open(p_aff_cleaned, "w")
            f_aff_cleaned.write("\"\",\"ChEMBL ID\",\"Aff\"\n")

            for d_chem in self.l_work:
                ChEMBL_ID = d_chem["Molecule ChEMBL ID"]
                if ChEMBL_ID in l_chem:
                    if d_chem["pChEMBL Value"] == "":
                        f_aff_cleaned.write("\"%s\",\"%s\",NA\n"%(ChEMBL_ID, ChEMBL_ID))
                    else:
                        f_aff_cleaned.write("\"%s\",\"%s\",%s\n"%(ChEMBL_ID, ChEMBL_ID, d_chem["pChEMBL Value"]))
            f_aff_cleaned.close()
            runExternal.histAC50(p_aff_cleaned, self.pr_out)
            return p_aff_cleaned
        
        else:
            f_aff_cleaned = open(p_aff_cleaned, "w")
            f_aff_cleaned.write("\"ChEMBLID\",\"AC50-uM\",\"Aff\"\n")
            n_act = 0
            n_inact = 0

            for d_chem in self.l_work:
                ChEMBL_ID = d_chem["Molecule ChEMBL ID"]
                if ChEMBL_ID in l_chem:
                    val_aff =  d_chem["Standard Value"]
                    unit = d_chem["Standard Units"]
                    
                    l_val_converted = toolbox.convertUnit([val_aff], [unit])
                    val_converted = l_val_converted[0]

                    if val_converted <= cutoff_uM:
                        f_aff_cleaned.write("\"%s\",\"%s\",0\n"%(ChEMBL_ID, val_converted))
                        n_inact = n_inact + 1
                    else:
                        f_aff_cleaned.write("\"%s\",\"%s\",1\n"%(ChEMBL_ID, val_converted))
                        n_act = n_act + 1
            f_aff_cleaned.close()

            fsum = open(self.pr_out + "class_sum.txt", "w")
            fsum.write("NB chemical: %s\nNB active: %s\nNB inactive: %s\n"%(len(self.l_work), n_act, n_inact)) 
            fsum.close()

            return p_aff_cleaned

