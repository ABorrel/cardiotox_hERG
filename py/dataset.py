import toolbox




class dataset:
    def __init__(self, p_smi, p_AC50, pr_out):
        self.p_smi = p_smi
        self.p_AC50 = p_AC50
        self.pr_out = pr_out
    

    def build_dataset(self):

        d_chem = toolbox.loadMatrixToDict(self.p_smi)
        d_AC50 = toolbox.loadMatrixToDict(self.p_AC50)

        d_out = {}

        self.d_dataset = self.d_out

        return d_out