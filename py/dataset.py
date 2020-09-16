from os import path, rename
from random import shuffle

import toolbox
import pathFolder
import runExternal

# import descriptor computation scripts => precise folder where descriptor are included
import sys
sys.path.insert(0, "./../../../development/molecular-descriptors/")
import Chemical

from PIL import Image, ImageFont, ImageDraw
font = ImageFont.truetype("./OpenSans-Regular.ttf", size=24)

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

            if self.d_dataset[CASRN]["log10(AC50)"] != "NA":
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
        runExternal.barplotClass(p_filout)

        self.d_class = d_out_class_all


    def computeDesc(self, pr_desc):

        p_filout_RDKIT = pr_desc + "desc_1D2D.csv"
        p_filout_OPERA = pr_desc + "desc_OPERA.csv"
        if path.exists(p_filout_RDKIT) and path.exists(p_filout_OPERA):
            self.p_desc1D2D = p_filout_RDKIT
            self.p_desc_opera = p_filout_OPERA
        else:

            # create opera descriptor
            pr_OPERA = pathFolder.createFolder(pr_desc + "OPERA/")
            p_listchem = pr_OPERA + "listChem.smi"
            

            # extract descriptor 2D
            if not path.exists(p_filout_RDKIT):
                flistchem = open(p_listchem, "w")
                l_desc = Chemical.getLdesc("1D2D")

                # open filout
                filout = open(p_filout_RDKIT, "w")
                filout.write("CASRN\tSMILES\t%s\n"%("\t".join(l_desc)))

                # compute descriptor
                l_smi = []
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
                        l_smi.append(cChem.smi)
                filout.close()
                flistchem.write("\n".join(l_smi))
                flistchem.close()

            #############
            # RUN OPERA #
            # need to add name in OPERA
            if path.exists(pr_OPERA + "desc_OPERA.csv"):
                l_chem_opera = toolbox.loadMatrixToList(pr_OPERA + "desc_OPERA.csv", sep=",")
                l_chem_rdkit = toolbox.loadMatrixToList(p_filout_RDKIT, sep = "\t")
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

            else:
                print("RUN opera with => ", p_listchem)
                ooo

            self.p_desc1D2D = p_filout_RDKIT
            self.p_desc_opera = p_filout_OPERA


    def computeDescOPERA(self, pr_desc):

        p_filout = self.pr_desc + "desc_OPERA.csv"
        if path.exists(p_filout):
            return p_filout

        # write list of SMILES for OPERA
        pr_OPERA = pathFolder.createFolder(pr_desc + "OPERA/")
        p_lSMI = pr_OPERA + "listChem.smi"
        flSMI = open(p_lSMI, "w")
        l_w = []
        for CASRN in self.d_dataset.keys():
            SMILES = self.d_dataset[CASRN]["SMILES"]
            if SMILES != "--":
                l_w.append(self.d_dataset[CASRN]["SMILES"])
        
        flSMI.write("\n".join(l_w))
        flSMI.close()

        print("RUN OPERA COMMAND LINE")

        print("ERROR - OPERA desc no computed")
        return 0


        return 


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
                cChem = Chemical.Chemical(SMILES, pr_desc)#, p_salts=path.abspath("./Salts.txt"), OS="Windows")
                cChem.prepChem() # prep
                p_png_inch = cChem.computePNG()
                if cChem.err == 0:
                    rename(p_png_inch, p_png)


    def rankActiveChem(self, pr_PNG):

        pr_out = pathFolder.createFolder(self.pr_out + "ranking/")
        
        l_aff = []
        d_temp = {}
        for CASRN in self.d_dataset.keys():
            aff = self.d_dataset[CASRN]["log10(AC50)"]
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
        
        print(l_aff)

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
        runExternal.mergepdfs(lpdf, pr_out + "chem_rank.pdf")
