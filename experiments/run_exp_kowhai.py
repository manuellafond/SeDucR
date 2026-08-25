import os
from os import system
from getopt import getopt
import sys
import re
import time
import itertools
import shutil


'''
This script runs some kowhai simulations and then runs each chosen method on them. 
Variables below need to be set.

Note that SeDucR was previously named insider.  Hence you will find the name 
insider scattered around this script.

The script runs the chosen methods (see the variable "methods = [...]" below) on each simphy directory.
It outputs a stats file, including running times and reconciliation cost. 
The script also saves one csv per simulation directory in the work directory.  
If the script fails or is stopped, when it is run again it will not recalculate simulations 
for which the csv file exists.  However, stats.csv will be overwritten.  
The full stats file can be recovered with (this is chatgpt's solution):
head -n 1 "$(ls work_kowhai/*.csv | head -n 1)" > stats_wgd25.csv
tail -n +2 -q work_kowhai/*.csv >> stats_wgd25.csv

'''


#directories
segdupDir = ""
kowhaiDir = "~/git/kowhai/"
fastmultrecDir = "/home/manuel/git/FastMultRec/FastMultRec/build/"
insider_dir = "/home/manuel/git/SeDucR/build/"

insider_relax_dir = "/home/manuel/git/ybrec2/build/"    #old test, please ignore

workdir = "./work_kowhai/"
os.makedirs(workdir, exist_ok=True)



skip_existing_kowhai = True
skip_existing_csv = True


#kowhai possible options, DEFAULT MUST BE FIRST
nH_vals = [50, 20, 30, 40, 60, 70, 80, 90, 100]    #leaves  [100,20,50]
nP_vals = [20, 30, 40, 50, 60, 70, 80, 90, 100]    #nb gtrees [100,20,50]
rB_vals = [2.0, 1.0, 3.0, 4.0, 5.0]
pC_vals = [0.5]
pJ_vals = [0.5, 0.2, 0.8, 1.0]


'''
#for testing
nH_vals = [20]    
nP_vals = [20]    
rB_vals = [2.0]
pC_vals = [0.5]
pJ_vals = [0.5]
'''



#segdup/multrec options
d_vals = [25,5,10,15,20]
l = 1
iterations = 20000    #nb reps that segdup will make

#replicates
replicates = 10


methods = ["segdup", "insider", "fastmultrec", "lca"]







'''
This function receives five arrays and returns the list of 5-tuples, one element per array,
containing the combinations of values where all values are the first element of the arrays, 
except for possibly one array.  This is useful when the first value of the array is the default,
so we keep all fixed to default except one parameter that varies.
'''
def get_param_combinations(a, b, c, d, e):
    #global nH_vals, nP_vals, rB_vals, pC_vals, pJ_vals
    
    base = (a[0], b[0], c[0], d[0], e[0])
    results = [base]

    arrays = (a, b, c, d, e)
    for i, arr in enumerate(arrays):
        for val in arr[1:]:
            combo = list(base)
            combo[i] = val
            results.append(tuple(combo))

    return results



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


header_line = "method,nh,np,rb,pc,pj,repl,dup_cost,solution_cost,solution_nbdups,solution_nblosses,time\n"
out.write(header_line)


param_combos = get_param_combinations(nH_vals, nP_vals, rB_vals, pC_vals, pJ_vals)

#commented line below was testing EVERY param combo
#for (nH, nP, rB, pC, pJ, d, r) in itertools.product(nH_vals, nP_vals, rB_vals, pC_vals, pJ_vals, d_vals, range(replicates)):
for ((nH, nP, rB, pC, pJ), r) in itertools.product(param_combos, range(replicates)):

    #filenames will have all the parameter info
    suffix = f"nH{str(nH)}_nP{str(nP)}_nR{1}_rB{str(rB)}_pC{str(pC)}_pJ{str(pJ)}_rep{r}"

    #each run produces a csv - if it already exists we can skip it
    run_csv_filename = workdir + f"stats_{suffix}.csv"
    if skip_existing_csv and os.path.exists(run_csv_filename):
        print(run_csv_filename + " skipped")
        continue

    run_csv_file = open(run_csv_filename, 'w')
    run_csv_file.write(header_line)

    #do a simulation
    kowhai_file_segdup = workdir + "kowhai_for_segdup_" + suffix + ".txt"
    kowhai_file_multrec = workdir + "kowhai_for_multrec_" + suffix + ".txt"
    
    #unless it should be skipped
    if not (skip_existing_kowhai and os.path.exists(kowhai_file_segdup)):
        command = kowhaiDir + "kowhai --sim -nH " + str(nH) + " -nP " + str(nP) + " -nR 1 -rB " + str(rB) + " -pC " + str(pC) + " -pJ " + str(pJ) + " --for-segdup --for-multrec"
        command += " > /dev/null"
        print(f"Rep {r} : {command}")
    
        system(command)

        print("Simulation done")
    
    

        shutil.move("for-segdup-from-kowhai.txt", kowhai_file_segdup)
        shutil.move("for-multrec-from-kowhai.txt", kowhai_file_multrec)
    
    
    
    #now solve the kohwai instance for every dup cost
    for d in d_vals:
    
        suffix_d = f"{suffix}_d{d}"
        #-------------------------------------------------------------------------------------
        #segdup
        #-------------------------------------------------------------------------------------
        if "segdup" in methods:
            sptree_filename = workdir + "kowhai_sptreesegdup_" + suffix_d + ".newick"
            gt_segdup_tree_filename = workdir + "kowhai_gtreesegdup_" + suffix_d + ".txt"
            create_segdup_files_from_kowhai(kowhai_file_segdup, sptree_filename, gt_segdup_tree_filename)
        
            start = time.perf_counter()
            #command = "cat ./for-segdup-from-kowhai.txt | xargs " + segdupDir + "segdup -n " + str(iterations) + " -Tinit 3 -Tfinal 0.0 -d " + str(d) + " -l " + str(l) + " > segdup-output.txt"
            command = f"{segdupDir}segdup -Sfile {sptree_filename} -Gfile {gt_segdup_tree_filename} -d {d} -l {str(l)} -n {iterations} -Tinit 3 -Tfinal 0.0 > {workdir}segdupout_{suffix_d}.txt"
            print(command)
            system(command)
            (nbdups, nblosses, cost) = get_segdup_costs(f"{workdir}segdupout_{suffix_d}.txt")
            elapsed = time.perf_counter() - start 

            out.write(f"segdup,{nH},{nP},{rB},{pC},{pJ},{r},{d},{cost},{nbdups},{nblosses},{elapsed:.6f}\n")
            run_csv_file.write(f"segdup,{nH},{nP},{rB},{pC},{pJ},{r},{d},{cost},{nbdups},{nblosses},{elapsed:.6f}\n")


        #-------------------------------------------------------------------------------------
        #fastmultrec
        #-------------------------------------------------------------------------------------
        if "fastmultrec" in methods:
            sptree_filename = workdir + "kowhai_sptreemultrec_" + suffix_d + ".newick"
            gt_multrec_tree_filename = workdir + "kowhai_gtreemultrec_" + suffix_d + ".txt"
            create_fastmultrec_files_from_kowhai(kowhai_file_multrec, sptree_filename, gt_multrec_tree_filename)

            start = time.perf_counter()
            #command = fastmultrecDir + "segrec -d " + str(d) + " -l " + str(l) + " " + multrecInput[:-3] + "\" -o {workdir}multrec-output.txt"
            command = f"{fastmultrecDir}segrec -d {str(d)} -l {str(l)} -sf {sptree_filename} -gf {gt_multrec_tree_filename} -o {workdir}multrec-output_{suffix_d}.txt"
            print(command)
            system(command)
            elapsed = time.perf_counter() - start 

            (mrCost, mrDups, mrLosses) = get_insider_costs(f"{workdir}multrec-output_{suffix_d}.txt")
            out.write(f"fastmultrec,{nH},{nP},{rB},{pC},{pJ},{r},{d},{mrCost},{mrDups},{mrLosses},{elapsed:.6f}\n")
            run_csv_file.write(f"fastmultrec,{nH},{nP},{rB},{pC},{pJ},{r},{d},{mrCost},{mrDups},{mrLosses},{elapsed:.6f}\n")

        #-------------------------------------------------------------------------------------
        #insider relax branch
        #-------------------------------------------------------------------------------------
        if "insider_relax" in methods:
            start = time.perf_counter()
            #command = insider_relax_dir + "ybrec -d " + str(d) + " -l " + str(l) + " " + multrecInput[:-3] + "\" -o insider-output.txt"
            command = f"{insider_relax_dir}ybrec -d {str(d)} -l {str(l)} -sf {sptree_filename} -gf {gt_multrec_tree_filename} -o {workdir}insider-relax-output_{suffix_d}.txt"
            print(command)
            system(command)
            elapsed = time.perf_counter() - start 
        
            (insCost, insDups, insLosses) = get_insider_costs(f"{workdir}insider-relax-output_{suffix_d}.txt")
            out.write(f"insider_relax,{nH},{nP},{rB},{pC},{pJ},{r},{d},{insCost},{insDups},{insLosses},{elapsed:.6f}\n")
            run_csv_file.write(f"insider_relax,{nH},{nP},{rB},{pC},{pJ},{r},{d},{insCost},{insDups},{insLosses},{elapsed:.6f}\n")
        
        
        #-------------------------------------------------------------------------------------
        #insider branch (not relax)
        #-------------------------------------------------------------------------------------
        if "insider" in methods:
            start = time.perf_counter()
            #command = insider_relax_dir + "ybrec -d " + str(d) + " -l " + str(l) + " " + multrecInput[:-3] + "\" -o insider-output.txt"
            command = f"{insider_dir}seducr -d {str(d)} -l {str(l)} -sf {sptree_filename} -gf {gt_multrec_tree_filename} -o {workdir}insider-output_{suffix_d}.txt"
            print(command)
            system(command)
            elapsed = time.perf_counter() - start 
        
            (insCost, insDups, insLosses) = get_insider_costs(f"{workdir}insider-output_{suffix_d}.txt")
            out.write(f"insider,{nH},{nP},{rB},{pC},{pJ},{r},{d},{insCost},{insDups},{insLosses},{elapsed:.6f}\n")
            run_csv_file.write(f"insider,{nH},{nP},{rB},{pC},{pJ},{r},{d},{insCost},{insDups},{insLosses},{elapsed:.6f}\n")



        #-------------------------------------------------------------------------------------
        #lca map
        #-------------------------------------------------------------------------------------
        if "lca" in methods:
            #command = fastmultrecDir + "segrec -d " + str(d) + " -l " + str(l) + " " + multrecInput[:-3] + "\" -o lca-output.txt  -al lca"
            command = f"{fastmultrecDir}segrec -d {str(d)} -l {str(l)} -sf {sptree_filename} -gf {gt_multrec_tree_filename} -al lca -o {workdir}lca-output_{suffix_d}.txt"
            print(command)
            start = time.perf_counter()
            system(command)
            elapsed = time.perf_counter() - start 

            (lcaCost, lcaDups, lcaLosses) = get_insider_costs(f"{workdir}lca-output_{suffix_d}.txt")
            out.write(f"lca,{nH},{nP},{rB},{pC},{pJ},{r},{d},{lcaCost},{lcaDups},{lcaLosses},{elapsed:.6f}\n")
            run_csv_file.write(f"lca,{nH},{nP},{rB},{pC},{pJ},{r},{d},{lcaCost},{lcaDups},{lcaLosses},{elapsed:.6f}\n")


    run_csv_file.close()


out.close()


