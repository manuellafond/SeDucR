from os import system
from getopt import getopt
import sys
import time
import re

#directories
segdupDir = "~/segdup/"
kowhaiDir = "~/software/kowhai/"
multrecDir = "~/software/MultRec/Multrec/"
segrecDir = "~/segrec-bnb/build/"

#kowhai options
nH = 20
nP = 20
rB = 1.0
pC = 0.5
pJ = 0.5

#segdup/multrec/segrec options
d = 10
l = 1
iterations = 10000

#replicates
replicates = 100

opts, args = getopt(sys.argv[1:], "d:l:r:i:",["segdup-dir=","kowhai-dir=","multrec-dir=","segrec-dir=","nH=","nP=","rB=","pC=","pJ=","replicates=","iterations="])

for opt, arg in opts:
    if opt == "-d":
        d = int(arg)
    elif opt == "-l":
        l = int(arg)
    elif opt == "--segdup-dir":
        segdupDir = arg
        if segdupDir[-1] != "\/":
            segdupDir = segdupDir + "\/"
    elif opt == "--kowhai-dir":
        kowhaiDir = arg
        if kowhaiDir[-1] != "\/":
            kowhaiDir = kowhaiDir + "\/"
    elif opt == "--multrec-dir":
        multrecDir = arg
        if multrecDir[-1] != "\/":
            multrecDir = multrecDir + "\/"
    elif opt == "--segrec-dir":
        segrecDir = arg
        if segrecDir[-1] != "\/":
            segrecDir = segrecDir + "\/"
    elif opt == "--nH":
        nH = int(arg)
    elif opt == "--nP":
        nP = int(arg)
    elif opt == "--rB":
        rB = float(arg)
    elif opt == "--pC":
        pC = float(arg)
    elif opt == "--pJ":
        pJ = float(arg)
    elif opt in ("-r", "--replicates"):
        replicates = int(arg)
    elif opt in ("-i", "--iterations"):
        iterations = int(arg)

system("rm summary.csv")

mrResults = []

sdTime = []

for r in range(replicates):
    #run kowhai
    system(kowhaiDir + "kowhai --sim -nH " + str(nH) + " -nP " + str(nP) + " -nR 1 -rB " + str(rB) + " -pC " + str(pC) + " -pJ " + str(pJ) + " --for-segdup --for-multrec --verbose > /dev/null")

    #run segdup
    #print("Running segdup...")
    curTime = time.time()
    system("cat ./for-segdup-from-kowhai.txt | xargs " + segdupDir + "segdup -n " + str(iterations) + " -Tinit 10 -Tfinal 0.0 -d " + str(d) + " -l " + str(l) + " > /dev/null")
    sdTime.append(time.time() - curTime)

    #run multrec
    multrecFile = open("for-multrec-from-kowhai.txt")
    multrecInput = multrecFile.readline()
    multrecFile.close()

    #print("Running multrec...")
    curTime = time.time()
    system(multrecDir + "Multrec -d " + str(d) + " -l " + str(l) + " " + multrecInput[:-3] + "\" -o multrec-output.txt")
    mrTime = time.time() - curTime

    #parse multrec output
    multrecOutput = open("multrec-output.txt")
    for line in multrecOutput:
        if line[1:-2] == "COST":
            mrCost = int(next(multrecOutput).rstrip())
        elif line[1:-2] == "DUPHEIGHT":
            mrDups = int(next(multrecOutput).rstrip())
        elif line[1:-2] == "NBLOSSES":
            mrLosses = int(next(multrecOutput).rstrip())
            break
    mrResults.append((mrDups,mrLosses,mrCost,mrTime))
    multrecOutput.close()

    #parse input for segrec
    srSpecies = open("segrec-sp.txt", "w")
    srGenes = open("segrec-genes.txt", "w")
    sptree = multrecInput[4:multrecInput.find(';')+1]
    srSpecies.write(re.sub("v\d\+", "'&'", sptree))
    srSpecies.close()
    srGenes.close()
    exit(0)

    #run segrec

f = open("summary.csv")
output = open("results.csv", "w")

output.write("nCospec,nIndividualDups,nAllDupEvents,nJointDups,nXtinc,nHostSwitch,nLineageSort,codivs,dups,losses,cost,sdTime,mrDups,mrLosses,mrCost,mrTime\n")

rep = 0
headerRead = False
for line in f:
    if headerRead:
        kowhaiOutput = line.rstrip()
    else:
        kowhaiOutput = next(f).rstrip()
    segdupHeader = next(f).rstrip()

    #check if failure
    if segdupHeader[:6] != "codivs":
        rep = rep + 1
        headerRead = True
        continue

    segdupOutput = next(f).rstrip()
    output.write(kowhaiOutput + segdupOutput + "," + str(sdTime[rep]) + "," + ",".join([str(i) for i in mrResults[rep]]) + "\n")
    rep = rep + 1
    headerRead = False

f.close()
output.close()

#clean
system("rm summary.csv multrec-output.txt")
