import toolbox
import CompDesc
import runExternal

class hergml:

    def __init__(self, pd_dataset, p_pred_hergml, p_aff, pr_out):
        
        self.p_pred_hergml = p_pred_hergml
        self.p_aff = p_aff
        self.pd_dataset = pd_dataset
        self.pr_out = pr_out


    def computePerfAllSet(self):

        d_pred = toolbox.loadMatrix(self.p_pred_hergml, sep = ",")
        
        if not isinstance(self.pd_dataset, dict) :
            d_dataset = toolbox.loadMatrix(self.p_dataset, sep = "\t")
            if len(list(d_dataset[list(d_dataset.keys())[0]].keys())) == 1:
                d_dataset = toolbox.loadMatrix(self.p_dataset, sep = ",")
        else:
            d_dataset = self.pd_dataset
        
        if self.p_aff != "1":
            d_aff = toolbox.loadMatrix(self.p_aff, sep = "\t")
            if len(list(d_aff[list(d_aff.keys())[0]].keys())) == 1:
                d_aff = toolbox.loadMatrix(self.p_aff, sep = ",")


        d_out = {}
        for chem in d_dataset.keys():
            CASRN = d_dataset[chem]["CASRN"]
            smi = d_dataset[chem]["SMILES"]
            
            d_out[CASRN] = {}
            try:
                d_out[CASRN]["consensus_pred"] = d_pred[smi]["consensus_pred"]
            except:
                continue

            if self.p_aff == "1":
                d_out[CASRN]["aff"] = "1"
            else:
                try:
                    if d_aff[CASRN]["LogAC50"] == "NA" :
                        d_out[CASRN]["aff"] = "0"
                    else:
                        d_out[CASRN]["aff"] = "1"
                except:
                    d_out[CASRN]["aff"] = "NA"
                    continue
        
        p_filout = self.pr_out + "pred_all"
        filout = open(p_filout, "w")
        filout.write("CASRN\taff\tpred\n")
        for CASRN in d_out.keys():
            if d_out[CASRN] == {}:
                continue
            else:
                filout.write("%s\t%s\t%s\n"%(CASRN, d_out[CASRN]["aff"], d_out[CASRN]["consensus_pred"]))
        
        filout.close()
        runExternal.computePerf(p_filout)


    def computePerfTestSet(self, p_test_set):

        d_pred = toolbox.loadMatrix(self.p_pred_hergml, sep = ",")
        d_dataset = toolbox.loadMatrix(self.p_dataset, sep = ",")
        d_aff = toolbox.loadMatrix(self.p_aff)


        # for the test set
        d_test = toolbox.loadMatrix(p_test_set, sep = ",")
        l_CASRN_test = list(d_test)
        
        d_out = {}
        for chem in d_dataset.keys():
            CASRN = d_dataset[chem]["CASRN"]
            smi = d_dataset[chem]["SMILES"]
            
            if not CASRN in l_CASRN_test:
                continue

            d_out[CASRN] = {}
            try:
                d_out[CASRN]["consensus_pred"] = d_pred[smi]["consensus_pred"]
            except:
                continue

            try:
                if d_aff[CASRN]["LogAC50"] == "NA" :
                    d_out[CASRN]["aff"] = "0"
                else:
                    d_out[CASRN]["aff"] = "1"
            except:
                d_out[CASRN]["aff"] = "NA"
                continue
        
        p_filout = self.pr_out + "pred_test"
        filout = open(p_filout, "w")
        filout.write("CASRN\taff\tpred\n")
        for CASRN in d_out.keys():
            if d_out[CASRN] == {}:
                continue
            else:
                filout.write("%s\t%s\t%s\n"%(CASRN, d_out[CASRN]["aff"], d_out[CASRN]["consensus_pred"]))
        
        filout.close()
        runExternal.computePerf(p_filout)

