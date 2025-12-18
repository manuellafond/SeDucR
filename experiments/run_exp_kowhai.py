from os import system
from getopt import getopt
import sys
import re
import time
import itertools

#directories
segdupDir = ""
kowhaiDir = "~/git/kowhai/"
fastmultrecDir = "/home/manuel/git/FastMultRec/FastMultRec/build/"
insider_dir = "/home/manuel/git/ybrec/build/"
insider_relax_dir = "/home/manuel/git/ybrec2/build/"


#kowhai options
nH_vals = [20, 40]
nP_vals = [20, 100]
rB_vals = [1.0]
pC_vals = [0.5]
pJ_vals = [0.5]

#segdup/multrec options
d_vals = [10]
l = 1
iterations = 20000

#replicates
replicates = 100





#fetches reconciliation costs from a file in insider format (also works for fastmultrec)
#returns (total cost, nbdups, nblosses)
def get_insider_costs(filename):
        
    nbdups = -1
    nblosses = -1
    cost = -1
    
    
    with open(filename) as f:
        prevline = ""
        for line in f:
            line = line.strip()
            
            if prevline == "<COST>":
                cost = float(line)
            elif prevline == "<DUPHEIGHT>":
                nbdups = float(line)
            elif prevline == "<NBLOSSES>":
                nblosses = float(line)
            
            prevline = line
        
    return (cost, nbdups, nblosses)




def get_segdup_costs(filename):
    pattern = r"(\d+)d\+(\d+)l\t(\d+)"

    with open(filename) as f:
        lines = f.readlines()
    

    target_line = lines[-3].rstrip("\n")  # line just before the last

    m = re.search(pattern, target_line)
    if m:
        d, l, c = m.groups()
        return (d, l, c)
    else:
        print("ERROR: segdup out not formed as expected.")        


out = open("stats_kowhai.csv", 'w')

out.write("method,nh,np,rb,pc,pj,dup_cost,solution_cost,solution_nbdups,solution_nblosses,time\n")

for (nH, nP, rB, pC, pJ, d, r) in itertools.product(nH_vals, nP_vals, rB_vals, pC_vals, pJ_vals, d_vals, range(replicates)):


    #do a simulation
    command = kowhaiDir + "kowhai --sim -nH " + str(nH) + " -nP " + str(nP) + " -nR 1 -rB " + str(rB) + " -pC " + str(pC) + " -pJ " + str(pJ) + " --for-segdup --for-multrec --verbose"
    print(command)
    system(command + " > /dev/null")


        
    #-------------------------------------------------------------------------------------
    #segdup
    #-------------------------------------------------------------------------------------
    start = time.perf_counter()
    command = "cat ./for-segdup-from-kowhai.txt | xargs " + segdupDir + "segdup -n " + str(iterations) + " -Tinit 10 -Tfinal 0.0 -d " + str(d) + " -l " + str(l) + " > segdup-output.txt"
    print(command)
    system(command)
    (nbdups, nblosses, cost) = get_segdup_costs("segdup-output.txt")
    elapsed = time.perf_counter() - start 

    out.write(f"segdup,{nH},{nP},{rB},{pC},{pJ},{d},{cost},{nbdups},{nblosses},{elapsed:.3f}\n")


    #-------------------------------------------------------------------------------------
    #fastmultrec
    #-------------------------------------------------------------------------------------
    multrecFile = open("for-multrec-from-kowhai.txt")
    multrecInput = multrecFile.readline()
    multrecFile.close()

    start = time.perf_counter()
    command = fastmultrecDir + "segrec -d " + str(d) + " -l " + str(l) + " " + multrecInput[:-3] + "\" -o multrec-output.txt"
    print(command)
    system(command)
    elapsed = time.perf_counter() - start 

    (mrCost, mrDups, mrLosses) = get_insider_costs("multrec-output.txt")
    out.write(f"fastmultrec,{nH},{nP},{rB},{pC},{pJ},{d},{mrCost},{mrDups},{mrLosses},{elapsed:.3f}\n")


    #-------------------------------------------------------------------------------------
    #insider relax branch
    #-------------------------------------------------------------------------------------
    start = time.perf_counter()
    command = insider_relax_dir + "ybrec -d " + str(d) + " -l " + str(l) + " " + multrecInput[:-3] + "\" -o insider-output.txt"
    print(command)
    system(command)
    elapsed = time.perf_counter() - start 

    (insCost, insDups, insLosses) = get_insider_costs("insider-output.txt")
    out.write(f"insider_relax,{nH},{nP},{rB},{pC},{pJ},{d},{insCost},{insDups},{insLosses},{elapsed:.3f}\n")



    #-------------------------------------------------------------------------------------
    #lca map
    #-------------------------------------------------------------------------------------
    start = time.perf_counter()
    command = fastmultrecDir + "segrec -d " + str(d) + " -l " + str(l) + " " + multrecInput[:-3] + "\" -o lca-output.txt  -al lca"
    print(command)
    system(command)
    elapsed = time.perf_counter() - start 

    (lcaCost, lcaDups, lcaLosses) = get_insider_costs("lca-output.txt")
    out.write(f"lca,{nH},{nP},{rB},{pC},{pJ},{d},{lcaCost},{lcaDups},{lcaLosses},{elapsed:.3f}\n")






