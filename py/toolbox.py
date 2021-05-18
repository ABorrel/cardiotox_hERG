def loadMatrixToList(pmatrixIn, sep = "\t"):

    filin = open(pmatrixIn, "r", encoding="utf8", errors='ignore')
    llinesMat = filin.readlines()
    filin.close()

    l_out = []
    line0 = formatLine(llinesMat[0])
    line1 = formatLine(llinesMat[1])
    lheaders = line0.split(sep)
    lval1 = line1.split(sep)

    # case where R written
    if len(lheaders) == (len(lval1) -1):
        lheaders = ["ID"] + lheaders


    i = 0
    while i < len(lheaders):
        if lheaders[i] == "":
            lheaders[i] = "ID"
        i += 1

    i = 1
    imax = len(llinesMat)
    while i < imax:
        lineMat = formatLine(llinesMat[i])
        lvalues = lineMat.split(sep)
        j = 0
        if len(lvalues) != len(lheaders):
            print("ERROR - line: ", i)
            print(lvalues)
            print(lheaders)
        jmax = len(lheaders)
        dtemp = {}
        while j < jmax:
            try:dtemp[lheaders[j]] = lvalues[j]
            except:pass
            j += 1
        l_out.append(dtemp)
        i += 1
    
    return l_out



def loadMatrix(pmatrixIn, sep = "\t"):

    filin = open(pmatrixIn, "r", encoding="utf8", errors='ignore')
    llinesMat = filin.readlines()
    filin.close()

    dout = {}
    line0 = formatLine(llinesMat[0])
    line1 = formatLine(llinesMat[1])
    lheaders = line0.split(sep)
    if len(lheaders) == 1:
        lheaders = line0.split(",")
        sep = ","
    lval1 = line1.split(sep)

    # case where R written
    if len(lheaders) == (len(lval1) -1):
        lheaders = ["ID"] + lheaders


    i = 0
    while i < len(lheaders):
        if lheaders[i] == "":
            lheaders[i] = "ID"
        i += 1

    i = 1
    imax = len(llinesMat)
    while i < imax:
        lineMat = formatLine(llinesMat[i])
        lvalues = lineMat.split(sep)
        kin = lvalues[0]
        dout[kin] = {}
        j = 0
        if len(lvalues) != len(lheaders):
            print("ERROR - line: ", i)
            print(lvalues)
            print(lheaders)
        jmax = len(lheaders)
        while j < jmax:
            try:dout[kin][lheaders[j]] = lvalues[j]
            except:pass
            j += 1
        i += 1

    return dout


def formatLine(linein):

    linein = linein.replace("\n", "")
    linenew = ""

    imax = len(linein)
    i = 0
    flagchar = 0
    while i < imax:
        if linein[i] == '"' and flagchar == 0:
            flagchar = 1
        elif linein[i] == '"' and flagchar == 1:
            flagchar = 0

        if flagchar == 1 and linein[i] == ",":
            linenew = linenew + " "
        else:
            linenew = linenew + linein[i]

        i += 1

    linenew = linenew.replace('\"', "")
    return linenew



def writeMatrix(ddesc, pdescAct, sep = "\t"):


    filout = open(pdescAct, "w")
    lheader = list(ddesc[list(ddesc.keys())[0]].keys())

    # put header in first

    if "CAS" in lheader:
        del lheader[lheader.index("CAS")]
        lheader = ["CAS"] + lheader
    elif "CASID" in lheader:
        del lheader[lheader.index("CASID")]
        lheader = ["CASID"] + lheader
    else:
        lheader = ["CASID"] + lheader
        for casID in list(ddesc.keys()):
            ddesc[casID]["CASID"] = casID


    filout.write(sep.join(lheader) + "\n")
    for casID in list(ddesc.keys()):
        lval = [str(ddesc[casID][i]) for i in lheader]
        filout.write(sep.join(lval) + "\n")
    filout.close()


def formatLineDataset(linein):

    linein = linein.replace("\n", "")
    linenew = ""

    imax = len(linein)
    i = 0
    flagchar = 0
    while i < imax:
        if linein[i] == '"' and flagchar == 0:
            flagchar = 1
        elif linein[i] == '"' and flagchar == 1:
            flagchar = 0

        if flagchar == 1 and linein[i] == ",":
            linenew = linenew + " "
        else:
            linenew = linenew + linein[i]
        i += 1

    linenew = linenew.replace('\"', "")
    return linenew




def convertUnit(l_values, l_units):
    """Convert list of affinity in uM"""

    #print l_values
    #print l_units
    lout = []
    i = 0
    imax = len(l_values)
    while i < imax:
        if l_units[i] == "uM" or l_units[i] == "10'-6M" or l_units[i] == "um" or l_units[i] == "":
            lout.append(l_values[i])
            i += 1
        elif l_units[i] == "nM":
            val = float(l_values[i])
            val = val/1000
            lout.append(val)
            i += 1
        else:
            print(l_units[i], "ssss")
            ffff

    #print lout

    return lout



def mergeDict(l_dict):

    dout = {}

    lk = list(l_dict[0].keys())

    for k in lk:
        lval = [l_dict[i][k] for i in range(0,len(l_dict))]
        #identic
        if lval.count(lval[0]) == len(lval):
            dout[k] = lval[0]
        else:
            # remove duplicate
            lval = list(dict.fromkeys(lval))
            dout[k] = "----".join(lval)

    return dout