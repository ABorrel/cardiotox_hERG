import urllib.request



class loadComptox:
    def __init__(self, search_strg):

        self.input = search_strg
        self.err = 0

    def searchInDB(self):

        try:handle = urllib.request.urlopen('https://comptox.epa.gov/dashboard/dsstoxdb/results?search=' + str(self.input))
        except: 
            self.err = 1
            return
        resp = handle.read()
        resp = str(resp.decode("utf8"))
        handle.close()
    
        try:
            SMILES = resp.split("SMILES: ")[1][0:1000]
            SMILES = SMILES.split("\\nInChI")[0]

            DTXSID = resp.split("dsstox_substance_id&quot;:&quot;")[1][0:1000]
            DTXSID = DTXSID.split("&quot")[0]

            nameChem = resp.split("IUPAC Name: ")[1][0:1000]
            nameChem = nameChem.split("\\n")[0]

            CASRN = resp.split("casrn&quot;:&quot;")[1][0:1000]
            CASRN = CASRN.split("&quot")[0]

            if SMILES == None or nameChem == None or DTXSID == None or CASRN == None:
                self.err = 1
            else:
                self.nameChem = nameChem
                self.SMILES = SMILES
                self.DTXSID = DTXSID
                self.CASRN = CASRN
        except:
            self.err=1
