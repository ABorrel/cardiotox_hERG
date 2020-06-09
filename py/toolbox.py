


def loadMatrix(pmatrixIn, sep = "\t"):

    filin = open(pmatrixIn, "r")
    llinesMat = filin.readlines()
    filin.close()

    dout = {}
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
        kin = lvalues[0]
        dout[kin] = {}
        j = 0
        if len(lvalues) != len(lheaders):
            print(lvalues)
            print(lheaders)
        jmax = len(lheaders)
        while j < jmax:
            try:dout[kin][lheaders[j]] = lvalues[j]
            except:pass
            j += 1
        i += 1

    return dout


def loadMatrixToDict(pmatrixIn, sep ="\t"):

    filin = open(pmatrixIn, "r")
    llinesMat = filin.readlines()
    filin.close()

    dout = {}
    line0 = formatLine(llinesMat[0])
    line1 = formatLine(llinesMat[1])
    lheaders = line0.split(sep)
    lval1 = line1.split(sep)

    # case where R written
    if len(lheaders) == (len(lval1)-1):
        lheaders.append("val")

    i = 1
    imax = len(llinesMat)
    while i < imax:
        lineMat = formatLine(llinesMat[i])
        lvalues = lineMat.split(sep)
        kin = lvalues[0]
        dout[kin] = {}
        j = 0
        if len(lvalues) != len(lheaders):
            print(lineMat)
            print(llinesMat[i])
            print(lvalues)
            print("Error => nb element", i)
            print(len(lvalues))
            print(len(lheaders))
            ddd

        jmax = len(lheaders)
        while j < jmax:
            dout[kin][lheaders[j]] = lvalues[j]
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