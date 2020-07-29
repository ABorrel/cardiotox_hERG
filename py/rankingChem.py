import toolbox
from math import log10
from PIL import Image, ImageFont, ImageDraw
import runExternalSoft

font = ImageFont.truetype("./OpenSans-Regular.ttf", size=24)



def rank_chem(p_dataset, p_class, pr_PNG, pr_out):

        
    # load files
    d_MIC = toolbox.loadMatrix(p_dataset)

    d_cluster = toolbox.loadMatrix(p_cluster, sep=",")
    
    
    # rank pMIC
    d_rank = {}

    for chem in d_MIC.keys():
        for orga in d_MIC[chem].keys():
            if orga == "SMILES" or orga =='CMPD_CHEMBLID' :
                continue
            if not orga in d_rank.keys():
                d_rank[orga] = []
            pMIC = -log10(float(d_MIC[chem][orga]))
            d_MIC[chem][orga] = pMIC
            d_rank[orga].append(pMIC)
    
    # order pMIC
    for orga in d_rank.keys():
        d_rank[orga].sort(reverse=True)

    # reorganise d_cluster
    d_cluster_used = {}
    for chem in d_cluster.keys():
        cluster = int(d_cluster[chem]["cluster"])
        if not cluster in d_cluster_used.keys():
            d_cluster_used[cluster] = []
        d_cluster_used[cluster].append(chem)

    # build the png
    lcluster = d_cluster_used.keys()
    lcluster.sort()

    i_page = 1

    l_pngout = []
    for cluster in lcluster:
        nchem = len(d_cluster_used[cluster])
        l_pchempng = []
        lw = []
        i_chem = 0
        nchem_page = 0
        while i_chem <= nchem:
            if nchem_page == 6 or i_chem == nchem:
                p_image = pr_out + "page_"+ str(i_page) + ".png"
                l_pngout.append(p_image)
                imgnew = Image.new("RGBA", (1535, 1285), (250, 250, 250))

                i_image=0
                i_img_max = len(l_pchempng)
                while i_image < i_img_max:

                    # put png
                    img1 = Image.open(l_pchempng[i_image])
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
                l_pchempng = []
                i_page = i_page + 1
                if i_chem == nchem:
                    break
                else:
                    continue

            l_pchempng.append(pr_PNG + d_cluster_used[cluster][i_chem] + ".png")
            lw.append("%s (cluster: %s)"%(d_cluster_used[cluster][i_chem], cluster))
            lw.append("pMIC (E. coli): %.2f (%s)"%(d_MIC[d_cluster_used[cluster][i_chem]]["Escherichia coli"], d_rank["Escherichia coli"].index(d_MIC[d_cluster_used[cluster][i_chem]]["Escherichia coli"])+ 1))
            lw.append("pMIC (P. aeruginosa): %.2f (%s)"%(d_MIC[d_cluster_used[cluster][i_chem]]["Pseudomonas aeruginosa"], d_rank["Pseudomonas aeruginosa"].index(d_MIC[d_cluster_used[cluster][i_chem]]["Pseudomonas aeruginosa"]) + 1))
            lw.append("pMIC (S. aureus): %.2f (%s)"%(d_MIC[d_cluster_used[cluster][i_chem]]["Staphylococcus aureus"], d_rank["Staphylococcus aureus"].index(d_MIC[d_cluster_used[cluster][i_chem]]["Staphylococcus aureus"])+ 1))
            lw.append("pMIC (S. pneumoniae): %.2f (%s)"%(d_MIC[d_cluster_used[cluster][i_chem]]["Streptococcus pneumoniae"], d_rank["Streptococcus pneumoniae"].index(d_MIC[d_cluster_used[cluster][i_chem]]["Streptococcus pneumoniae"])+ 1))            
            nchem_page = nchem_page + 1
            i_chem = i_chem + 1


    # transform png to pdf
    lpdf = []
    for ppng in l_pngout:
        ppdf = runExternalSoft.pngtopdf(ppng)
        lpdf.append(ppdf)

    # merge pdf sheet
    runExternalSoft.mergepdfs(lpdf, pr_out + "chem_pMIC.pdf")