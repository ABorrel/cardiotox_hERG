from copy import deepcopy
from numpy import mean, std
from os import path

import toolbox



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

        if not "l_work" in dir(self):
            self.parseCHEMBLFile()

        self.get_standard_type(l_standard_type)
        self.get_standard_relation(l_standard_relation)
        self.filterDuplicates()

        self.writeTable()


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

    def checkIdenticSMI(self):

        if not "tableorgafull" in self.__dict__:
            return

        k1 = self.tableorgafull.keys()[0]

        lsmi = []
        i = 0
        imax = len(self.tableorgafull[k1])
        while i < imax:
            smi = self.tableorgafull[k1][i]["CANONICAL_SMILES"]
            smiclean = liganddescriptors.standardizeSMILES(smi)
            self.tableorgafull[k1][i]["SMILES_PREP"] = smiclean
            print("***", i, smiclean, imax, "***")
            if smiclean == 1:
                i += 1
                continue

            if not smiclean in lsmi:
                lsmi.append(smiclean)
            else:
                for k in self.tableorgafull.keys():
                    del self.tableorgafull[k][i]
                imax = imax - 1
                continue

            i += 1



        lorga = self.tableorgafull.keys()
        dw = {}
        for orga in self.tableorgafull.keys():
            dw[orga] = {}
            for chem in self.tableorgafull[orga]:
                chemblID = chem["CMPD_CHEMBLID"]
                print(chemblID)
                print(chem["STANDARD_VALUE"], "MIC")
                print(chem["PUBLISHED_UNITS"], "unit")
                print(chem["MOLWEIGHT"], "weight")
                
                MIC_M =  (float(chem["STANDARD_VALUE"]) / 1000)/ float( chem["MOLWEIGHT"])
                chem["MIC_M"] = str(MIC_M)
                dw[orga][chemblID] = chem
                
        
        filout = open(p_filout, "w")
        filout.write("CMPD_CHEMBLID\tSMILES\t" + "\t".join(lorga) + "\n")
        
        lchemID = dw[lorga[0]].keys()
        for chemID in lchemID:
            filout.write("%s\t%s\t%s\n"%(chemID, dw[lorga[0]][chemID]["SMILES_PREP"] ,"\t".join([dw[orga][chemID]["MIC_M"] for orga in lorga])))
        filout.close()