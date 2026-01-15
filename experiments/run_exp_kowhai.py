import os
from os import system
from getopt import getopt
import sys
import re
import time
import itertools
import shutil

#directories
segdupDir = ""
kowhaiDir = "~/git/kowhai/"
fastmultrecDir = "/home/manuel/git/FastMultRec/FastMultRec/build/"
insider_dir = "/home/manuel/git/ybrec/build/"
insider_relax_dir = "/home/manuel/git/ybrec2/build/"

workdir = "./work_kowhai/"
os.makedirs(workdir, exist_ok=True)


#kowhai possible options
nH_vals = [50, 20]	#leaves  [100,20,50]
nP_vals = [20, 50]	#nb gtrees [100,20,50]
rB_vals = [3.0, 1.0, 3.0]
pC_vals = [1]
pJ_vals = [1, 0.5, 0.8,0.2]


#segdup/multrec options
d_vals = [25,5,10,15,20]
l = 1
iterations = 20000	#nb reps that segdup will make

#replicates
replicates = 10


methods = ["segdup", "insider_relax", "fastmultrec", "lca"]




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
        
        
def create_segdup_files_from_kowhai(kowhai_infile, species_outfile, gt_outfile):
    #disclaimer, used chatgpt for this function
    sptree_newick = ""
    gtree_strings = []

    with open(kowhai_infile, "r") as f:
        tokens = f.read().split("-G")


    sptree_newick = tokens[0].split(" ")[1].replace('"', '')

    for i in range(1, len(tokens)):
        if tokens[i].strip() == "":
            continue
        gz = tokens[i].strip().replace('"', "").split(" ", maxsplit=1)
       
        gtree_strings.append(gz[0])
        gtree_strings.append(gz[1])

    with open(species_outfile, "w") as fs:
        fs.write(sptree_newick)

    with open(gt_outfile, "w") as fg:
        fg.write("\n".join(gtree_strings).replace('"', ""))
        









def create_fastmultrec_files_from_kowhai(kowhai_infile, species_outfile, gt_outfile):
    #disclaimer, used chatgpt for this function
    sptree_newick = ""

    with open(kowhai_infile, "r") as f:
        tokens = f.read().split("-g")
    print(tokens)

    sptree_newick = tokens[0].strip().split(" ")[1].replace('"', '')

    gtree_strings = [p.strip() + ";" for p in tokens[1].replace('"', '').strip().split(";")]

    with open(species_outfile, "w") as fs:
        fs.write(sptree_newick)

    with open(gt_outfile, "w") as fg:
        fg.write("\n".join(gtree_strings))
        


out = open("stats_kowhai.csv", 'w')

out.write("method,nh,np,rb,pc,pj,dup_cost,solution_cost,solution_nbdups,solution_nblosses,time\n")

for (nH, nP, rB, pC, pJ, d, r) in itertools.product(nH_vals, nP_vals, rB_vals, pC_vals, pJ_vals, d_vals, range(replicates)):


    #do a simulation
    command = kowhaiDir + "kowhai --sim -nH " + str(nH) + " -nP " + str(nP) + " -nR 1 -rB " + str(rB) + " -pC " + str(pC) + " -pJ " + str(pJ) + " --for-segdup --for-multrec"
    command += " > /dev/null"
    print(command)
    system(command)

    print("Simulation done")
    #filenames will have all the parameter info
    suffix = f"nH{str(nH)}_nP{str(nP)}_nR{1}_rB{str(rB)}_pC{str(pC)}_pJ{str(pJ)}_rep{r}"
    
    kowhai_file_segdup = workdir + "kowhai_for_segdup_" + suffix + ".txt"
    kowhai_file_multrec = workdir + "kowhai_for_multrec_" + suffix + ".txt"
    shutil.move("for-segdup-from-kowhai.txt", kowhai_file_segdup)
    shutil.move("for-multrec-from-kowhai.txt", kowhai_file_multrec)
    
        
    #-------------------------------------------------------------------------------------
    #segdup
    #-------------------------------------------------------------------------------------
    if "segdup" in methods:
        sptree_filename = workdir + "kowhai_sptreesegdup_" + suffix + ".newick"
        gt_segdup_tree_filename = workdir + "kowhai_gtreesegdup_" + suffix + ".txt"
        create_segdup_files_from_kowhai(kowhai_file_segdup, sptree_filename, gt_segdup_tree_filename)
        
        start = time.perf_counter()
        #command = "cat ./for-segdup-from-kowhai.txt | xargs " + segdupDir + "segdup -n " + str(iterations) + " -Tinit 3 -Tfinal 0.0 -d " + str(d) + " -l " + str(l) + " > segdup-output.txt"
        command = f"{segdupDir}segdup -Sfile {sptree_filename} -Gfile {gt_segdup_tree_filename} -d {d} -l {str(l)} -n {iterations} -Tinit 3 -Tfinal 0.0 > {workdir}segdupout_{suffix}.txt"
        print(command)
        system(command)
        (nbdups, nblosses, cost) = get_segdup_costs(f"{workdir}segdupout_{suffix}.txt")
        elapsed = time.perf_counter() - start 

        out.write(f"segdup,{nH},{nP},{rB},{pC},{pJ},{d},{cost},{nbdups},{nblosses},{elapsed:.3f}\n")


    #-------------------------------------------------------------------------------------
    #fastmultrec
    #-------------------------------------------------------------------------------------
    if "fastmultrec" in methods:
        sptree_filename = workdir + "kowhai_sptreemultrec_" + suffix + ".newick"
        gt_multrec_tree_filename = workdir + "kowhai_gtreemultrec_" + suffix + ".txt"
        create_fastmultrec_files_from_kowhai(kowhai_file_multrec, sptree_filename, gt_multrec_tree_filename)

        start = time.perf_counter()
        #command = fastmultrecDir + "segrec -d " + str(d) + " -l " + str(l) + " " + multrecInput[:-3] + "\" -o {workdir}multrec-output.txt"
        command = f"{fastmultrecDir}segrec -d {str(d)} -l {str(l)} -sf {sptree_filename} -gf {gt_multrec_tree_filename} -o {workdir}multrec-output_{suffix}.txt"
        print(command)
        system(command)
        elapsed = time.perf_counter() - start 

        (mrCost, mrDups, mrLosses) = get_insider_costs(f"{workdir}multrec-output_{suffix}.txt")
        out.write(f"fastmultrec,{nH},{nP},{rB},{pC},{pJ},{d},{mrCost},{mrDups},{mrLosses},{elapsed:.3f}\n")

    #-------------------------------------------------------------------------------------
    #insider relax branch
    #-------------------------------------------------------------------------------------
    if "insider_relax" in methods:
        start = time.perf_counter()
        #command = insider_relax_dir + "ybrec -d " + str(d) + " -l " + str(l) + " " + multrecInput[:-3] + "\" -o insider-output.txt"
        command = f"{insider_relax_dir}ybrec -d {str(d)} -l {str(l)} -sf {sptree_filename} -gf {gt_multrec_tree_filename} -o {workdir}insider-output_{suffix}.txt"
        print(command)
        system(command)
        elapsed = time.perf_counter() - start 
    
        (insCost, insDups, insLosses) = get_insider_costs(f"{workdir}insider-output_{suffix}.txt")
        out.write(f"insider_relax,{nH},{nP},{rB},{pC},{pJ},{d},{insCost},{insDups},{insLosses},{elapsed:.3f}\n")



    #-------------------------------------------------------------------------------------
    #lca map
    #-------------------------------------------------------------------------------------
    if "lca" in methods:
        start = time.perf_counter()
        #command = fastmultrecDir + "segrec -d " + str(d) + " -l " + str(l) + " " + multrecInput[:-3] + "\" -o lca-output.txt  -al lca"
        command = f"{fastmultrecDir}segrec -d {str(d)} -l {str(l)} -sf {sptree_filename} -gf {gt_multrec_tree_filename} -al lca -o {workdir}lca-output_{suffix}.txt"
        print(command)
        system(command)
        elapsed = time.perf_counter() - start 

        (lcaCost, lcaDups, lcaLosses) = get_insider_costs(f"{workdir}lca-output_{suffix}.txt")
        out.write(f"lca,{nH},{nP},{rB},{pC},{pJ},{d},{lcaCost},{lcaDups},{lcaLosses},{elapsed:.3f}\n")






