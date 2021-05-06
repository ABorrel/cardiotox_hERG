from copy import deepcopy
from numpy import mean, std
from os import path
from re import search

import toolbox
import runExternal
import CompDesc
import pathFolder

# selected physico chemical descriptor from OPERA
L_OPERA_DESC = ['LogP_pred', 'MP_pred', 'BP_pred', 'LogVP_pred', 'LogWS_pred', 'LogHL_pred', 'RT_pred', 'LogKOA_pred', 'ionization', 'LogD55_pred', 'LogD74_pred', 'LogOH_pred', 'LogBCF_pred', 'BioDeg_LogHalfLife_pred', 'ReadyBiodeg_pred', 'LogKM_pred', 'LogKoc_pred', 'FUB_pred', 'Clint_pred']



class CHEMBLTable:
    def __init__(self, p_filin, pr_out, verbose = 1):
        self.p_filin = p_filin
        self.pr_out = pr_out
        self.verbose = 1
        self.slog = ""

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

    def cleanDataset(self, l_standard_type, l_standard_relation, assay_type_favorised = ""):

        # load existing data if exist
        p_filout = self.pr_out + "ChEMBL_cleaned.csv"
        if path.exists(p_filout):
            self.l_work = toolbox.loadMatrixToList(p_filout, sep="\t")
            self.p_dataset_cleaned = p_filout
        else:
            if not "l_work" in self.__dict__:
                self.parseCHEMBLFile()

            self.get_standard_type(l_standard_type)
            self.get_standard_relation(l_standard_relation)
            if assay_type_favorised != "":
                self.filter_doubleByAssayType(assay_type_favorised)
            self.filterDuplicates_activity()
            self.filterDuplicate_structure()
            self.writeTable()
            self.p_dataset_cleaned = p_filout

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
        self.slog = self.slog + "Filter by relation => [%s]\nNb value => %s\n"%(";".join(l_relation), imax)
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
        self.slog = self.slog + "After filtering by relation => %s\n"%(len(self.l_work))

    def get_standard_type(self, l_type):

        i = 0
        imax = len(self.l_work)
        if self.verbose == 1:
            print("Filter by type => [%s]"%(";".join(l_type)))
            print("Nb value => %s"%(imax))
        self.slog = self.slog + "Filter by type => [%s]\nNb value => %s\n"%(";".join(l_type), imax)

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
        self.slog = self.slog + "After filtering by type => %s\n"%(len(self.l_work))

    def filter_doubleByAssayType(self, assay_type):

        d_CHEMBLID = {}
        for row in self.l_work:
            if not row["Molecule ChEMBL ID"] in list(d_CHEMBLID.keys()):
                d_CHEMBLID[row["Molecule ChEMBL ID"]] = []
            d_CHEMBLID[row["Molecule ChEMBL ID"]].append(deepcopy(row))

        if self.verbose == 1:
            print("Nb chemicals with multiple activity => %s"%(len(self.l_work) - len(list(d_CHEMBLID.keys()))))
        self.slog = self.slog + "Nb chemicals with multiple activity => %s\n"%(len(self.l_work) - len(list(d_CHEMBLID.keys())))

        nb_del = 0
        for chem_id in d_CHEMBLID.keys():
            if len(d_CHEMBLID[chem_id]) == 1:
                continue

            flag_assay = 0
            for row in d_CHEMBLID[chem_id]:
                if search(assay_type.lower(), row["Assay Description"].lower()):
                    flag_assay = 1
            
            if flag_assay == 1:
                i = 0
                imax = len(d_CHEMBLID[chem_id])
                while i < imax:
                    if not search(assay_type.lower(), d_CHEMBLID[chem_id][i]["Assay Description"].lower()):
                        del d_CHEMBLID[chem_id][i]
                        imax = imax - 1
                        nb_del = nb_del + 1
                    else:
                        i = i + 1
        
        # rebuild work list
        self.l_work = []
        for k in d_CHEMBLID.keys():
            self.l_work = self.l_work + d_CHEMBLID[k]

        if self.verbose:
            print("Nb after favorise %s assays => %s"%(assay_type, len(self.l_work)))
        self.slog = self.slog +  "Nb deleted chemicals: %s\nNb after assay activity cleaning => %s\n"%(nb_del, len(self.l_work))

    def filterDuplicate_structure(self):
        
        if self.verbose == 1:
            print("Filter structural duplicate => %s chemicals\n"%(len(self.l_work)))
        self.slog = self.slog + "Filter structural duplicate => %s chemicals\n"%(len(self.l_work))

        l_smiles_clean = []
        i = 0
        imax = len(self.l_work)
        while i < imax:
            smiles = self.l_work[i]["Smiles"]
            cChem = CompDesc.CompDesc(smiles, self.pr_out)
            cChem.prepChem()
            if cChem.err ==1:
                del  self.l_work[i]
                imax = imax - 1
                continue
            else:
                if cChem.smi in l_smiles_clean:
                    del  self.l_work[i]
                    imax = imax - 1
                    continue
                else:
                    l_smiles_clean.append(cChem.smi)
                    i = i + 1 


        if self.verbose == 1:
            print("After structural duplicate => %s chemicals\n"%(len(self.l_work)))
        self.slog = self.slog + "After structural duplicate => %s chemicals\n"%(len(self.l_work))

    def filterDuplicates_activity(self):

        nb_del = 0

        d_CHEMBLID = {}
        for row in self.l_work:
            if not row["Molecule ChEMBL ID"] in list(d_CHEMBLID.keys()):
                d_CHEMBLID[row["Molecule ChEMBL ID"]] = []
            d_CHEMBLID[row["Molecule ChEMBL ID"]].append(deepcopy(row))

        if self.verbose == 1:
            print("Nb chemicals with multiple activity => %s"%(len(self.l_work) - len(list(d_CHEMBLID.keys()))))
        self.slog = self.slog + "Nb chemicals with multiple activity => %s\n"%(len(self.l_work) - len(list(d_CHEMBLID.keys())))

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
                    nb_del = nb_del + 1
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
            print("Nb after duplicate activity cleaning => %s"%(len(self.l_work)))
        self.slog = self.slog +  "Nb deleted chemicals: %s\nNb after duplicate activity cleaning => %s\n"%(nb_del, len(self.l_work))

    def writeTable(self):

        filout = open(self.pr_out + "ChEMBL_cleaned.csv", "w")
        filout.write("\t".join(self.l_header) + "\n")
        for row in self.l_work:
            lw = [row[h] for h in self.l_header]
            filout.write("\t".join(lw) + "\n")
        filout.close()

        p_log = self.pr_out + "log.txt"
        flog = open(p_log, "w")
        flog.write("%s"%(self.slog))
        flog.close()

    def computeDesc(self, pr_desc):

        if not "pr_desc" in self.__dict__:
            self.pr_desc = pr_desc

        p_filout_RDKIT = self.pr_out + "desc_1D2D.csv"
        if path.exists(p_filout_RDKIT) :
            self.p_desc1D2D = p_filout_RDKIT
            return self.p_desc1D2D

        # extract descriptor 2D
        if not path.exists(p_filout_RDKIT):
            cChem = CompDesc.CompDesc("", self.pr_desc)
            l_desc = cChem.getLdesc("1D2D")

            # open filout
            filout = open(p_filout_RDKIT, "w")
            filout.write("CHEMBLID\tSMILES\t%s\n"%("\t".join(l_desc)))

            # compute descriptor
            l_smi = []
            for row in self.l_work:
                SMILES = row["Smiles"]
                CHEMBL_ID = row["Molecule ChEMBL ID"]
                cChem = CompDesc.CompDesc(SMILES, self.pr_desc)
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
                    filout.write("%s\t%s\t%s\n"%(CHEMBL_ID, cChem.smi, "\t".join(["NA" if not desc in list(cChem.all2D.keys()) else str(cChem.all2D[desc])  for desc in l_desc])))
                    l_smi.append(cChem.smi)
            filout.close()

        self.p_desc1D2D = p_filout_RDKIT
        return p_filout_RDKIT

    def computeOPERADesc(self, pr_desc):

        if not "pr_desc" in self.__dict__:
            self.pr_desc = pr_desc

        p_filout = self.pr_out + "desc_OPERA.csv"
        if path.exists(p_filout):
            self.p_descOPERA = p_filout
            return p_filout

        # write list of SMILES for OPERA
        pr_OPERA = pathFolder.createFolder(self.pr_desc + "OPERA/", clean=1)
        p_lSMI = pr_OPERA + "listChem.smi"
        flSMI = open(p_lSMI, "w")
        l_w = []
        for row in self.l_work:
            SMILES = row["Smiles"]
            CHEMBL_ID = row["Molecule ChEMBL ID"]
            l_w.append(SMILES)
        
        flSMI.write("\n".join(l_w))
        flSMI.close()

        p_desc_opera = pr_OPERA + "desc_opera_run.csv"
        cCompDesc = CompDesc.CompDesc(p_lSMI, pr_OPERA)
        cCompDesc.computeOperaFromListSMI(p_desc_opera)

        l_chem = self.l_work
        l_ddesc_run = toolbox.loadMatrixToList(p_desc_opera, sep = ",")

        fopera = open(p_filout, "w")
        fopera.write("CHEMBLID,%s\n"%(",".join(L_OPERA_DESC)))

        i = 0
        imax = len(l_chem) 
        while i < imax:
            CHEMBLID = l_chem[i]["Molecule ChEMBL ID"]
            for ddesc_run in l_ddesc_run:
                if ddesc_run["MoleculeID"] == "Molecule_%i"%(i+1):
                    fopera.write("%s,%s\n"%(CHEMBLID, ",".join(ddesc_run[desc] for desc in L_OPERA_DESC)))
                    break
            i = i + 1
        fopera.close()

        self.p_descOPERA = p_filout

        return p_filout

    def prep_aff(self, typeAff="pAff", cutoff_uM=30000):

        if type(cutoff_uM) == list:
            p_aff_cleaned = "%saff_%s_%s.csv"%(self.pr_out, typeAff, "-".join([str(i) for i in cutoff_uM]))
        else:
            p_aff_cleaned = "%saff_%s_%s.csv"%(self.pr_out, typeAff, cutoff_uM)
        
        if path.exists(p_aff_cleaned):
            self.p_aff_clean = p_aff_cleaned
            return 

        if typeAff == "pAff":
            f_aff_cleaned = open(p_aff_cleaned, "w")
            f_aff_cleaned.write("\"\",\"ChEMBL ID\",\"Aff\"\n")

            for d_chem in self.l_work:
                ChEMBL_ID = d_chem["Molecule ChEMBL ID"]
                if d_chem["pChEMBL Value"] == "":
                    f_aff_cleaned.write("\"%s\",\"%s\",NA\n"%(ChEMBL_ID, ChEMBL_ID))
                else:
                    f_aff_cleaned.write("\"%s\",\"%s\",%s\n"%(ChEMBL_ID, ChEMBL_ID, d_chem["pChEMBL Value"]))
            f_aff_cleaned.close()
            runExternal.histAC50(p_aff_cleaned, self.pr_out)
            return p_aff_cleaned
        
        else:
            f_aff_cleaned = open(p_aff_cleaned, "w")
            f_aff_cleaned.write("\"ChEMBLID\",\"Aff-uM\",\"Aff\",\"pAff\"\n")
            n_act = 0
            n_inact = 0
            nb_remove = 0

            for d_chem in self.l_work:
                ChEMBL_ID = d_chem["Molecule ChEMBL ID"]
                val_aff =  d_chem["Standard Value"]
                unit = d_chem["Standard Units"]
                pAff = d_chem["pChEMBL Value"]
                    
                l_val_converted = toolbox.convertUnit([val_aff], [unit])
                val_converted = l_val_converted[0]

                
                if type(cutoff_uM) == list:
                    if val_converted <= cutoff_uM[0]:
                        f_aff_cleaned.write("\"%s\",\"%s\",0,%s\n"%(ChEMBL_ID, val_converted, pAff))
                        n_act = n_act + 1
                    elif val_converted >= cutoff_uM[1]:
                        f_aff_cleaned.write("\"%s\",\"%s\",1,%s\n"%(ChEMBL_ID, val_converted, pAff))
                        n_inact = n_inact + 1
                    else:
                        nb_remove = nb_remove + 1
                else:
                    if val_converted <= cutoff_uM:
                        f_aff_cleaned.write("\"%s\",\"%s\",1,%s\n"%(ChEMBL_ID, val_converted, pAff))
                        n_inact = n_inact + 1
                    else:
                        f_aff_cleaned.write("\"%s\",\"%s\",0,%s\n"%(ChEMBL_ID, val_converted, pAff))
                        n_act = n_act + 1
            f_aff_cleaned.close()

            fsum = open(self.pr_out + "class_sum.txt", "w")
            fsum.write("NB chemical: %s\nNB removed chemicals: %s\nNB active: %s\nNB inactive: %s\n"%(len(self.l_work), nb_remove, n_act, n_inact)) 
            fsum.close()

            self.p_aff_clean = p_aff_cleaned

    def correlation_aff(self, p_dataset_to_compare, pr_comparison, pr_root_results):

        l_smiles_chembl = []
        d_chem_chembl = toolbox.loadMatrix(self.p_dataset_cleaned)
        for chem in d_chem_chembl.keys():
            smi = d_chem_chembl[chem]["Smiles"]
            c_ChemDec = CompDesc.CompDesc(smi, pr_root_results + "DESC/")
            c_ChemDec.prepChem()
            if c_ChemDec.err == 1:
                d_chem_chembl[chem]["SMILES_CLEAN"] = "NA"
            else:
                d_chem_chembl[chem]["SMILES_CLEAN"] = c_ChemDec.smi
                l_smiles_chembl.append(c_ChemDec.smi)

        l_smiles_to_compare = []
        d_chem_dataset_to_compare = toolbox.loadMatrix(p_dataset_to_compare)
        for chem_to_compare in d_chem_dataset_to_compare.keys():
            smi = d_chem_dataset_to_compare[chem_to_compare]["SMILES"]
            c_ChemDec = CompDesc.CompDesc(smi, pr_root_results + "DESC/")
            c_ChemDec.prepChem()
            if c_ChemDec.err == 1:
                d_chem_dataset_to_compare[chem_to_compare]["SMILES_CLEAN"] = "NA"
            else:
                d_chem_dataset_to_compare[chem_to_compare]["SMILES_CLEAN"] = c_ChemDec.smi
                l_smiles_to_compare.append(c_ChemDec.smi)

        l_smi_inter = list(set.intersection(set(l_smiles_chembl), set(l_smiles_to_compare)))

        p_fcor = pr_comparison + "comparison"
        if not path.exists(p_fcor):
            fcor = open(p_fcor, "w")
            fcor.write("CASRN\tCHEMBL\tAff_CHEMBL\tAff_NCAST\tAssay_CHEMBL\n")
            for smi_inter in l_smi_inter:
                for chem_CHEMBL in d_chem_chembl.keys():
                    if d_chem_chembl[chem_CHEMBL]["SMILES_CLEAN"] == smi_inter:
                        CHEMBLID = d_chem_chembl[chem_CHEMBL]["Molecule ChEMBL ID"]
                        aff_chembl = d_chem_chembl[chem_CHEMBL]["pChEMBL Value"]
                        assays =  d_chem_chembl[chem_CHEMBL]["Assay Description"]
                        break
                
                for chem_to_compare in d_chem_dataset_to_compare.keys():
                    if d_chem_dataset_to_compare[chem_to_compare]["SMILES_CLEAN"] == smi_inter:
                        CASRN = d_chem_dataset_to_compare[chem_to_compare]["CASRN"]
                        aff_NCAST = d_chem_dataset_to_compare[chem_to_compare]["-log10(AC50)"]
                        break
                
                if search("patch clamp", assays):
                    fcor.write("%s\t%s\t%s\t%s\tPatch clamp\n"%(CASRN, CHEMBLID, aff_chembl, aff_NCAST))
                else:
                    fcor.write("%s\t%s\t%s\t%s\tOther\n"%(CASRN, CHEMBLID, aff_chembl, aff_NCAST))
            fcor.close()
            
            runExternal.corplotAff(p_fcor)


    
    def mergeDataset(self, p_aff_to_merge, p_desc_to_merge, pr_out, rm_inact=1):
        
        p_filout = pr_out + "aff_merged.csv"
        if path.exists(p_filout):
            self.p_merge_sets = p_filout
            return
        
        d_tomerge_aff = toolbox.loadMatrix(p_aff_to_merge)
        d_tomerge_desc = toolbox.loadMatrix(p_desc_to_merge)

        d_chembl_desc = toolbox.loadMatrix(self.p_desc1D2D)
        d_chembl_aff = toolbox.loadMatrix(self.p_aff_clean, sep = ",")

        d_out = {}
        d_duplicate = {}
        for chem_tomerge in d_tomerge_aff.keys():
            try:SMILES = d_tomerge_desc[chem_tomerge]["SMILES"]
            except:pass
            d_out[SMILES] = {}
            d_out[SMILES]["CASRN"] = chem_tomerge
            d_out[SMILES]["CHEMBL"] = "-"
            d_out[SMILES]["pAff-CHEMBL"] = "-"
            d_out[SMILES]["aff-CHEMBL"] = "-"
            d_out[SMILES]["pAff-CASRN"] = str(d_tomerge_aff[chem_tomerge]["-log10(AC50)"])
            d_out[SMILES]["aff-CASRN"] = str(d_tomerge_aff[chem_tomerge]["Aff"])

        
        for chem_chembl in d_chembl_aff.keys():
            SMILES = d_chembl_desc[chem_chembl]["SMILES"]
            if rm_inact == 1 and d_chembl_aff[chem_chembl]["Aff"] == "0":
                continue
            if not SMILES in list(d_out.keys()):
                d_out[SMILES] = {}
                d_out[SMILES]["CASRN"] = "-"
                d_out[SMILES]["aff-CASRN"] = "-"
                d_out[SMILES]["pAff-CASRN"] = "-"
            d_out[SMILES]["CHEMBL"] = chem_chembl
            d_out[SMILES]["pAff-CHEMBL"] = str(d_chembl_aff[chem_chembl]["pAff"])
            d_out[SMILES]["aff-CHEMBL"] = str(d_chembl_aff[chem_chembl]["Aff"])
        
        filout = open(p_filout, "w")
        filout.write("ID\tSMILES\tCASRN\tCHEMBLID\tpaff-add\tpaff-chembl\taff-add\taff-chembl\tAff\n")

        id = 1
        for SMILES in d_out.keys():
            if d_out[SMILES]["aff-CASRN"] != "-":
                Aff = d_out[SMILES]["aff-CASRN"]
            else:
                Aff = d_out[SMILES]["aff-CHEMBL"]
            d_out[SMILES]["Aff"] = Aff
            filout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(id, SMILES, d_out[SMILES]["CASRN"], d_out[SMILES]["CHEMBL"], d_out[SMILES]["pAff-CASRN"], d_out[SMILES]["pAff-CHEMBL"], d_out[SMILES]["aff-CASRN"], d_out[SMILES]["aff-CHEMBL"], Aff))
            id = id + 1
        filout.close()
        self.p_merge_sets = p_filout

        # write summary
        nb_chem = len(list(d_out.keys()))
        nb_active = 0
        nb_chembl = 0
        nb_add = 0
        nb_inactive = 0
        nb_overlap = 0


        for smiles in d_out.keys():
            if d_out[smiles]["Aff"] == "0":
                nb_inactive = nb_inactive + 1
            if d_out[smiles]["Aff"] == "1":
                nb_active = nb_active + 1
            if d_out[smiles]["CASRN"] != "-" and d_out[smiles]["CHEMBL"] != "-":
                nb_overlap = nb_overlap + 1
            if d_out[smiles]["CASRN"] != "-":
                nb_add = nb_add + 1
            if d_out[smiles]["CHEMBL"] != "-":
                nb_chembl = nb_chembl + 1   

        p_filout = pr_out + "merge_dataset.sum"
        filout = open(p_filout, "w")
        filout.write("Nb chemicals: %s\nNb chemicals in NCAST: %s\nNb chemicals in CHEMBL: %s\nNb overlap: %s\nNb active chemical: %s\nNb inactive chemical: %s\n"%(nb_chem, nb_add, nb_chembl, nb_overlap, nb_active, nb_inactive))   
        filout.close()     

        self.p_merge_sets = pr_out + "aff_merged.csv"
            