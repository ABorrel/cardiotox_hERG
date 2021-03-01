import pathFolder
import CHEMBLTable


class CHEMBL_set:
    def __init__(self, p_chembl_results, name_dataset, pr_results):
        self.p_chembl_results = p_chembl_results
        self.name_dataset = name_dataset
        self.pr_results = pr_results
        self.pr_desc = pathFolder.createFolder(pr_results + "DESC/")


    def main(self):
        l_standard_type=["IC50", "Ki", "EC50"]
        l_standard_relation=["'='"]

        self.prep_dataset(l_standard_type, l_standard_relation)
        self.prep_descset()


        ######
        # MERGE dataset CHEMBL + NCAST
        ##############

        p_NCAST_aff = self.pr_results + "dataset_NCAST/dataset_prep.csv"
        p_NCAST_desc = self.pr_results + "dataset_NCAST/desc_1D2D.csv"
        self.merge_dataset(p_NCAST_aff, p_NCAST_desc, "NCAST_CHEMBL")
        #self.correlation_dataset(p_NCAST_aff)


    def prep_dataset(self, l_standard_type, l_standard_relation):
        pr_dataset = pathFolder.createFolder(self.pr_results + "dataset_" + self.name_dataset + "/")
        self.pr_dataset = pr_dataset
        self.cCHEMBL = CHEMBLTable.CHEMBLTable(self.p_chembl_results, pr_dataset)
        self.cCHEMBL.parseCHEMBLFile()
        self.cCHEMBL.cleanDataset(l_standard_type=["IC50", "Ki", "EC50"], l_standard_relation=["'='"], assay_type_favorised = "patch clamp")
    
    def prep_descset(self):

        self.cCHEMBL.computeDesc(self.pr_desc)
        self.cCHEMBL.computeOPERADesc(self.pr_desc)
        p_aff_ChEMBL = self.cCHEMBL.prep_aff(typeAff="", cutoff_uM=[1,10])

        return 

    def correlation_dataset(self, p_dataset):
        pr_comparison = pathFolder.createFolder(self.pr_results + "comparison_CHEMBL_NCAST/")
        self.cCHEMBL.correlation_aff(p_dataset, pr_comparison, self.pr_results)

    def merge_dataset(self, p_aff_clean_toadd, p_desc1D2D_toadd, name_dataset):
        pr_out = pathFolder.createFolder(self.pr_results + name_dataset + "/")
        self.cCHEMBL.mergeDataset(p_aff_clean_toadd, p_desc1D2D_toadd, pr_out)

        self.pr_merge_sets = pr_out
