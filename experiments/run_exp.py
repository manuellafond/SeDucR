import sys
import os
import itertools
import time


#this script assumes the existence of directories of the form 
#"./sims/sim_{[1,2]}WGD_D{[7,10]}/sim_{[1..25]}
#and it outputs stats.csv with everything

wgd_nb = [1,2]
dup_rates = [7,10]
dup_costs = [5, 10]     #values of argument -d to test


#set to directories containing binaries of the programs to test
fastmultrecdir = "FastMultRec/FastMultRec/build"
ybrecdir = "ybrec/build"



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




output_data = []

out = open("stats.csv", 'w')

out.write("method,planted_wgds,duprate,simid,dup_cost,solution_cost,solution_nbdups,solution_nblosses,time\n")

for (simid, wgd, duprate, dup_cost) in itertools.product(range(1, 25 + 1), wgd_nb, dup_rates, dup_costs):
    
    
    directory = f"./sims/sim_{wgd}WGD_D{duprate}/sim_{simid}"
    genetree_file = directory + "/all_genetrees_edited.txt"
    speciestree_file = directory + "/s_tree.newick"
    
    #ybrec/insiser
    cmd = f"{ybrecdir}/ybrec -d {dup_cost} -l 1 -o out.txt -gf '{genetree_file}' -sf '{speciestree_file}'"
    print(cmd)
    
    start = time.perf_counter()
    os.system(cmd)
    end = time.perf_counter()
    
    (cost, nbdups, nblosses) = get_insider_costs("out.txt")
    
    elapsed = time.perf_counter() - start 
    out.write(f"insider,{wgd},{duprate},{simid},{dup_cost},{cost},{nbdups},{nblosses},{elapsed:.3f}\n")
    
    
    #fastmultrec
    cmd = f"{fastmultrecdir}/segrec -d {dup_cost} -l 1 -o out.txt -gf '{genetree_file}' -sf '{speciestree_file}'"
    print(cmd)
    
    start = time.perf_counter()
    os.system(cmd)
    end = time.perf_counter()
    
    (cost, nbdups, nblosses) = get_insider_costs("out.txt")
    
    elapsed = time.perf_counter() - start 
    out.write(f"fastmultrec_greedy,{wgd},{duprate},{simid},{dup_cost},{cost},{nbdups},{nblosses},{elapsed:.3f}\n")
    
    
    
    
    
    #lca mapping
    cmd = f"{fastmultrecdir}/segrec -d {dup_cost} -l 1 -o out.txt -gf '{genetree_file}' -sf '{speciestree_file}' -al lca"
    print(cmd)
    
    start = time.perf_counter()
    os.system(cmd)
    end = time.perf_counter()
    
    (cost, nbdups, nblosses) = get_insider_costs("out.txt")
    
    elapsed = time.perf_counter() - start 
    out.write(f"lcamap,{wgd},{duprate},{simid},{dup_cost},{cost},{nbdups},{nblosses},{elapsed:.3f}\n")

    
out.close()
