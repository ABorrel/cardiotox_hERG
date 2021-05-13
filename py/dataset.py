from os import path, rename
from random import shuffle
import rdkit
from rdkit import Chem
from re import search
from math import pow

import toolbox
import pathFolder
import runExternal
import CompDesc

# import descriptor computation scripts => precise folder where descriptor are included

from PIL import Image, ImageFont, ImageDraw
font = ImageFont.truetype("./OpenSans-Regular.ttf", size=24)

class dataset:
    def __init__(self, p_smi, p_AC50, pr_out):
        self.p_smi = p_smi
        self.p_AC50 = p_AC50
        self.pr_out = pr_out
    
    def prep_dataset(self, cutoff_aff):

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
            d_temp["-log10(AC50)"] = d_AC50[CASRN]["LogAC50"]
            if d_AC50[CASRN]["LogAC50"] == "NA":
                d_temp["Aff"] = 0
            else:
                val_uM = pow(10, -float(d_AC50[CASRN]["LogAC50"])) * pow(10,6) 
                if val_uM < cutoff_aff:
                    d_temp["Aff"] = 1
                else:
                    d_temp["Aff"] = 0

            # load on the tox21 chem library
            for DTXID in d_chem.keys():
                if d_chem[DTXID]["CASRN"] == CASRN:
                    d_temp["SMILES"] =  d_chem[DTXID]["SMILES"]
                    d_temp["PREFERRED_NAME"] =  d_chem[DTXID]["PREFERRED_NAME"]
                    d_out[CASRN] = d_temp
                    break
        
        # write the dataset
        l_header = ["CASRN", "SMILES", "PREFERRED_NAME", "-log10(AC50)", "Aff"]
        filout = open(p_filout, "w")
        filout.write("\t".join(l_header) + "\n")
        for CASRN in d_out.keys():
            filout.write("%s\n"%("\t".join([str(d_out[CASRN][h]) for h in l_header])))
        filout.close()
        self.d_dataset = d_out

    def classChem(self, p_classification):
        """
        Return histogram active by chemical class
        """

        if not "d_dataset" in self.__dict__:
            self.prep_dataset()
        
        d_classification = toolbox.loadMatrix(p_classification, sep = ",")

        # for active chem
        d_out_class_active = {}
        d_out_class_all = {}
        for CASRN in self.d_dataset:
            try:l_class_chem = d_classification[CASRN]["Classes"]
            except: l_class_chem = "No defined"

            if l_class_chem == "NA":
                l_class_chem = "No defined"

            if self.d_dataset[CASRN]["-log10(AC50)"] != "NA":
                d_out_class_active[CASRN] = l_class_chem

            d_out_class_all[CASRN] = l_class_chem
        
        
        # for active chemicals
        p_filout = self.pr_out + "class_active_" + str(p_classification.split("/")[-1][0:-4])
        filout = open(p_filout, "w")
        filout.write("CASRN\tClass\n")
        for chem in d_out_class_active.keys():
            l_class_chem = d_out_class_active[chem].split("--")
            for class_chem in l_class_chem:
                if class_chem == "NA":
                    class_chem = "No defined"
                filout.write("%s\t%s\n"%(chem, class_chem))
        filout.close()

        # make hist
        if not path.exists(p_filout + "_barplot.png"):
            runExternal.barplotClass(p_filout)

        # for all
        p_filout = self.pr_out + "class_" + str(p_classification.split("/")[-1][0:-4])
        filout = open(p_filout, "w")
        filout.write("CASRN\tClass\n")
        for chem in d_out_class_all.keys():
            l_class_chem = d_out_class_all[chem].split("--")
            for class_chem in l_class_chem:
                if class_chem == "NA":
                    class_chem = "No defined"
                filout.write("%s\t%s\n"%(chem, class_chem))
        filout.close()

        # make hist
        if not path.exists(p_filout + "_barplot.png"):
            runExternal.barplotClass(p_filout)
        self.d_class = d_out_class_all

    def computeDesc(self, pr_desc):

        p_filout_RDKIT = self.pr_out + "desc_1D2D.csv"
        p_filout_OPERA = self.pr_out + "desc_OPERA.csv"
        if path.exists(p_filout_RDKIT) and path.exists(p_filout_OPERA):
            self.p_desc1D2D = p_filout_RDKIT
            self.p_desc_opera = p_filout_OPERA
            return [self.p_desc1D2D,  self.p_desc_opera]

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
            self.p_desc_opera = self.computeDescOPERA(pr_desc)
        else:
            self.p_desc_opera = p_filout_OPERA

        return [self.p_desc1D2D,  self.p_desc_opera]

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

        return p_filout

    def computePNG(self, pr_desc):

        pr_png = pathFolder.createFolder(pr_desc + "PNG/")

        # compute descriptor
        l_CASRN = list(self.d_dataset.keys())
        shuffle(l_CASRN)
        for CASRN in l_CASRN:
            p_png = pr_png + CASRN + ".png"
            if path.exists(p_png):
                continue
            else:
                SMILES = self.d_dataset[CASRN]["SMILES"]
                cChem = CompDesc.CompDesc(SMILES, pr_desc)#, p_salts=path.abspath("./Salts.txt"), OS="Windows")
                cChem.prepChem() # prep
                p_png_inch = cChem.computePNG()
                if cChem.err == 0:
                    rename(p_png_inch, p_png)

    def rankActiveChem(self, pr_PNG):

        pr_out = pathFolder.createFolder(self.pr_out + "ranking/")
        
        p_filout = pr_out + "chem_rank.pdf"
        if path.exists(p_filout):
            return 


        l_aff = []
        d_temp = {}
        for CASRN in self.d_dataset.keys():
            aff = self.d_dataset[CASRN]["-log10(AC50)"]
            name = self.d_dataset[CASRN]["PREFERRED_NAME"]
            if aff == "NA":
                continue
            else:
                l_aff.append(float(aff))
                d_temp[CASRN] = {}
                d_temp[CASRN]["aff"] = aff
                d_temp[CASRN]["name"] = name
                d_temp[CASRN]["class"] = self.d_class[CASRN]

        l_aff.sort(reverse=True)
        
        lw = []
        l_ppng = []
        nb_page = 0
        i_page = 1
        rank = 1
        nchem_page = 0
        l_pngout = []

        for aff in l_aff:
            i = 0
            l_casrn = list(d_temp.keys())
            imax = len(l_casrn)
            while i < imax:
                if nchem_page == 6:
                    p_image = pr_out + "page_"+ str(i_page) + ".png"
                    l_pngout.append(p_image)
                    imgnew = Image.new("RGB", (1535, 1285), (250, 250, 250))

                    i_image=0
                    i_img_max = len(l_ppng)
                    while i_image < i_img_max:

                        # put png
                        img1 = Image.open(l_ppng[i_image])
                        if i_image < 3:
                            imgnew.paste(img1, (0 + i_image * 510,0))
                                
                            draw = ImageDraw.Draw(imgnew)
                            draw.text((5 + 510 * i_image, 500), lw[0 + (i_image*5)], (0, 0, 0), font=font)
                            draw.text((5 + 510 * i_image, 525), lw[1 + (i_image*5)], (0, 0, 0), font=font)
                            draw.text((5 + 510 * i_image, 550), lw[2 + (i_image*5)], (0, 0, 0), font=font)
                            draw.text((5 + 510 * i_image, 575), lw[3 + (i_image*5)], (0, 0, 0), font=font)
                            draw.text((5 + 510 * i_image, 600), lw[4 + (i_image*5)], (0, 0, 0), font=font)
                        else:
                            imgnew.paste(img1, (0 + (i_image-3) * 510, 650))
                            draw = ImageDraw.Draw(imgnew)
                            draw.text((5 + 510 * (i_image-3), 1150), lw[15 + ((i_image-3)*5)], (0, 0, 0), font=font)
                            draw.text((5 + 510 * (i_image-3), 1175), lw[16 + ((i_image-3)*5)], (0, 0, 0), font=font)
                            draw.text((5 + 510 * (i_image-3), 1200), lw[17 + ((i_image-3)*5)], (0, 0, 0), font=font)
                            draw.text((5 + 510 * (i_image-3), 1225), lw[18 + ((i_image-3)*5)], (0, 0, 0), font=font)
                            draw.text((5 + 510 * (i_image-3), 1250), lw[19 + ((i_image-3)*5)], (0, 0, 0), font=font)
                        i_image = i_image + 1
                    
                    imgnew.save(p_image)
                    # add
                    nchem_page = 0
                    lw = []
                    l_ppng = []
                    i_page = i_page + 1

                if float(d_temp[l_casrn[i]]["aff"]) == aff:
                    if not path.exists(pr_PNG + l_casrn[i] + ".png"):
                        i = i + 1
                        continue
                    l_ppng.append(pr_PNG + l_casrn[i] + ".png")
                    lw.append("%s"%(l_casrn[i]))
                    lw.append("Name: %s"%(d_temp[l_casrn[i]]["name"]))
                    lw.append("pIC50: %s"%(round(aff, 2)))
                    lw.append("Rank: %s"%(rank))
                    lw.append("Classes: %s"%(d_temp[l_casrn[i]]["class"]))
                    nchem_page = nchem_page + 1
                    
                    # remove chem for speed up the process of search
                    del d_temp[l_casrn[i]]
                    del l_casrn[i]
                    imax = imax - 1
                    i = i - 1
                    rank = rank + 1
                i = i + 1

        lpdf = []
        for ppng in l_pngout:
            ppdf = runExternal.pngtopdf(ppng)
            lpdf.append(ppdf)

        # merge pdf sheet
        runExternal.mergepdfs(lpdf, p_filout)

    def computeDescForNNComparison(self, pr_desc, pr_out):


        l_desc = ["SlogP","SMR","LabuteASA","TPSA","AMW","ExactMW","NumLipinskiHBA","NumLipinskiHBD","NumRotatableBonds","NumHBD","NumHBA","NumAmideBonds","NumHeteroAtoms","NumHeavyAtoms","NumAtoms","NumStereocenters","NumUnspecifiedStereocenters","NumRings","NumAromaticRings","NumSaturatedRings","NumAliphaticRings","NumAromaticHeterocycles","NumSaturatedHeterocycles","NumAliphaticHeterocycles","NumAromaticCarbocycles","NumSaturatedCarbocycles","NumAliphaticCarbocycles","FractionCSP3","Chi0v","Chi1v","Chi2v","Chi3v","Chi4v","Chi1n","Chi2n","Chi3n","Chi4n","HallKierAlpha","kappa1","kappa2","kappa3","slogp_VSA1","slogp_VSA2","slogp_VSA3","slogp_VSA4","slogp_VSA5","slogp_VSA6","slogp_VSA7","slogp_VSA8","slogp_VSA9","slogp_VSA10","slogp_VSA11","slogp_VSA12","smr_VSA1","smr_VSA2","smr_VSA3","smr_VSA4","smr_VSA5","smr_VSA6","smr_VSA7","smr_VSA8","smr_VSA9","smr_VSA10","peoe_VSA1","peoe_VSA2","peoe_VSA3","peoe_VSA4","peoe_VSA5","peoe_VSA6","peoe_VSA7","peoe_VSA8","peoe_VSA9","peoe_VSA10","peoe_VSA11","peoe_VSA12","peoe_VSA13","peoe_VSA14","MQN1","MQN2","MQN3","MQN4","MQN5","MQN6","MQN7","MQN8","MQN9","MQN10","MQN11","MQN12","MQN13","MQN14","MQN15","MQN16","MQN17","MQN18","MQN19","MQN20","MQN21","MQN22","MQN23","MQN24","MQN25","MQN26","MQN27","MQN28","MQN29","MQN30","MQN31","MQN32","MQN33","MQN34","MQN35","MQN36","MQN37","MQN38","MQN39","MQN40","MQN41","MQN42"] 
        # open filout
        p_filout = pr_out + "chemicals_desc.csv"
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

            l_out = []
            d_desc_out = {}
            for desc in l_desc:
                if desc in list(cChem.all2D.keys()):
                    d_desc_out[desc] = cChem.all2D[desc]
                else:
                    if search("peoe_", desc):
                        d_desc_out[desc] = cChem.all2D[desc.upper()]
                    elif search("slogp_", desc):
                        d_desc_out[desc] = cChem.all2D["SlogP_" + desc.split("_")[-1]]
                    elif search("smr_", desc):
                        d_desc_out[desc] = cChem.all2D[desc.upper()]
                    elif search("kappa", desc):
                        d_desc_out[desc] = cChem.all2D["s%s"%(desc)]
                    elif search("Chi.v", desc):
                        d_desc_out[desc] = cChem.all2D["Chi%s%s"%(desc[4], desc[3])]  
                    elif search("Chi.n", desc):
                        d_desc_out[desc] = cChem.all2D["Chiv%s"%(desc[3])] 
                    elif desc == 'ExactMW':
                        d_desc_out[desc] = cChem.all2D["ExactMolWt"] 
                    elif desc == 'NumAtoms':
                        d_desc_out[desc] = cChem.all2D["NumAllatoms"]
                    elif desc == 'NumHeteroAtoms':
                        d_desc_out[desc] = cChem.all2D["NumHeteroatoms"]
                    elif desc == 'NumHeavyAtoms':
                        d_desc_out[desc] = cChem.all2D["HeavyAtomCount"]
                    elif desc == 'NumRings':
                        d_desc_out[desc] = cChem.all2D["RingCount"]
                    elif desc == 'NumHBD':
                        d_desc_out[desc] = cChem.all2D["NumHDonors"]
                    elif desc == 'NumHBA':
                        d_desc_out[desc] = cChem.all2D["NumHAcceptors"]
                    elif desc == 'SlogP':
                        d_desc_out[desc] = cChem.all2D["MolLogP"]
                    elif desc == 'SMR':
                        d_desc_out[desc] = cChem.all2D["MolMR"]
                    elif desc == 'AMW':
                        d_desc_out[desc] = cChem.all2D["MolWt"]
                    else:
                        l_out.append(desc)
                
            l_FpMorgan = rdkit.Chem.AllChem.GetMorganFingerprintAsBitVect(cChem.mol, radius=2, nBits=1024)
            filout.write("%s,%s,%s\n"%(SMILES, ",".join([str(d_desc_out[desc]) for desc in l_desc]), ",".join([str(l_FpMorgan[i]) for i in range(0,1024)])))
        filout.close()

