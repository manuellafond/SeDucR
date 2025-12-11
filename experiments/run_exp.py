import sys
import os
import re
import itertools
import time
from ete3 import Tree


#this script assumes the existence of directories of the form 
#"{datadir}/sim_{[1,2]}WGD_D{[7,10]}/sim_{[1..25]}
#and it outputs stats.csv with everything

datadir = "/home/manuel/git/wgddata"


wgd_nb = [1,2]
dup_rates = [10, 7]
dup_costs = [5, 10]     #values of argument -d to test
runs = range(1, 10+1)

#wgd_nb = [1]
#dup_rates = [7]
#dup_costs = [10]
#runs = [6]


#set to directories containing binaries of the programs to test
fastmultrecdir = "/home/manuel/git/FastMultRec/FastMultRec/build"
ybrecdir = "/home/manuel/git/ybrec/build"
segdupdir = ""


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





def get_segdup_stree_str(s_tree_filename):
    tree = Tree(s_tree_filename, quoted_node_names=True, format=1);
    for leaf in tree.iter_leaves():
        leaf.name = "x" + leaf.name
    return tree.write(format=9)
        



def get_cleaned_gtree_str(newick, gtree_index):
    tree = Tree(newick, quoted_node_names=True, format=1)
    
    #todo, copy-pasted from get_segdup_gtree_str
    ncopies_per_gene = dict()
    

    first = True
    for leaf in tree.iter_leaves():

        parts = leaf.name.split("_")
        
        sp = parts[0]
        copynum = 0
        if sp in ncopies_per_gene:
            copynum = ncopies_per_gene[sp] + 1
        ncopies_per_gene[sp] = copynum

        leaf.name = f"x{parts[0]}_{copynum}_{parts[2]}_{gtree_index}"
        
    return tree.write(format=9)


#returns leaf association string to be used with segsup
def get_segdup_gtree_str(newick, gtree_index):
    tree = Tree(newick, quoted_node_names=True, format=1)
    
    #so, annoyingly, Reza allowed two genes to have the same id, eg 2_0_0 twice
    #I think it doesn't matter for fastmultrec since it uses internal ids, but it matters for segdup
    #so leaves are renamed according to their copy number
    ncopies_per_gene = dict()
    
    mapout = ""    #leaf assoc map string   
    first = True
    for leaf in tree.iter_leaves():

        parts = leaf.name.split("_")
        
        sp = parts[0]
        copynum = 0
        if sp in ncopies_per_gene:
            copynum = ncopies_per_gene[sp] + 1
        ncopies_per_gene[sp] = copynum

        leaf.name = f"x{parts[0]}_{copynum}_{parts[2]}_{gtree_index}"
        
        if first:
            first = False
        else:
            mapout += " "
            
        mapout += f"{leaf.name}:x{parts[0]}"
        

    #strout = f'"{tree.write(format=9)}" "{mapout}"'
    return (tree.write(format=9), mapout)











output_data = []

out = open("stats.csv", 'w')

out.write("method,planted_wgds,duprate,simid,dup_cost,solution_cost,solution_nbdups,solution_nblosses,time\n")

for (simid, wgd, duprate, dup_cost) in itertools.product(runs, wgd_nb, dup_rates, dup_costs):
    
    
    directory = f"{datadir}/sim_{wgd}WGD_D{duprate}/sim_{simid}"
    genetree_file = directory + "/all_genetrees_edited.txt"
    speciestree_file = directory + "/s_tree.newick"
    
    
    #-------------------------------------------------------------------------------------
    #segdup
    #-------------------------------------------------------------------------------------
    #note: simphy species trees add quotes around their leaf names, they must be removed for segdup
    stree_str = '-S "'
    stree_str += get_segdup_stree_str(speciestree_file)
    stree_str += '"'

    genetree_outfile = open("gfile.txt", "w")
    gtree_str = ""
    with open(genetree_file, "r") as file:
        cnt = 0
        lines = file.readlines()
        for line in lines:
            newick = line.strip()
            if newick != "":
                if cnt < 9999:   #dummy condition, was used to test limits
                    (gstr,mapstr) = get_segdup_gtree_str(newick, cnt)
                    genetree_outfile.write(gstr + "\n")
                    genetree_outfile.write(mapstr + "\n")
                    gtree_str += f' -G {gstr}'
                cnt += 1
    
    genetree_outfile.close()
    
    #segdup_input_filename = "segdup_input.txt"
    #with open(segdup_input_filename, "w") as f:
    #     f.write(f"{stree_str} {gtree_str} -n 100")    
    #command = f"cat {segdup_input_filename} | xargs {segdupdir}segdup"
    #command = f"{segdupdir}segdup {stree_str} {gtree_str} -n 100"
    #print(command)
    
    command = f"{segdupdir}segdup {stree_str} -Gfile gfile.txt -d {dup_cost} -n 50000 -Tinit 2 -Tfinal 0 > segdupout.txt"
    print(command)
    
    start = time.perf_counter()

    os.system(command)
    end = time.perf_counter()

    elapsed = time.perf_counter() - start 


    
    #get segdup's cost line, which has the form nCospec,nSegDup,nLoss,cost
    (nbdups, nblosses, cost) = get_segdup_costs("segdupout.txt")
    
    
    out.write(f"segdup,{wgd},{duprate},{simid},{dup_cost},{cost},{nbdups},{nblosses},{elapsed:.3f}\n")
    

    
    
    #-------------------------------------------------------------------------------------
    #get insider data (aka ybrec)
    #-------------------------------------------------------------------------------------
    genetree_insider_filename = "gfile_insider.txt"
    genetree_insider_file = open(genetree_insider_filename, "w")
    
    snewick = get_segdup_stree_str(speciestree_file)
    sptree_insider_filename = "sfile_insider.txt"
    with open(sptree_insider_filename, "w") as f:
        f.write(snewick)
    
    
    gtree_str = ""
    with open(genetree_file, "r") as file:
        cnt = 0
        lines = file.readlines()
        for line in lines:
            newick = line.strip()
            if newick != "":
                    gstr = get_cleaned_gtree_str(newick, cnt)
                    genetree_insider_file.write(gstr + "\n")
                    cnt += 1
    
    genetree_insider_file.close()
    
    
    #cmd = f"{ybrecdir}/ybrec -d {dup_cost} -l 1 -o out.txt -gf '{genetree_file}' -sf '{speciestree_file}'"
    cmd = f"{ybrecdir}/ybrec -d {dup_cost} -l 1 -o out.txt -gf '{genetree_insider_filename}' -sf '{sptree_insider_filename}'"
    print(cmd)
    
    start = time.perf_counter()
    os.system(cmd)
    end = time.perf_counter()
    
    (cost, nbdups, nblosses) = get_insider_costs("out.txt")
    ybcost = cost
    elapsed = time.perf_counter() - start 
    out.write(f"insider,{wgd},{duprate},{simid},{dup_cost},{cost},{nbdups},{nblosses},{elapsed:.3f}\n")
    
    
    #-------------------------------------------------------------------------------------
    #fastmultrec
    #-------------------------------------------------------------------------------------
    #cmd = f"{fastmultrecdir}/segrec -d {dup_cost} -l 1 -o out.txt -gf '{genetree_file}' -sf '{speciestree_file}'"
    cmd = f"{fastmultrecdir}/segrec -d {dup_cost} -l 1 -o out.txt -gf '{genetree_insider_filename}' -sf '{sptree_insider_filename}'"
    print(cmd)
    
    start = time.perf_counter()
    os.system(cmd)
    end = time.perf_counter()
    
    (cost, nbdups, nblosses) = get_insider_costs("out.txt")
    mrec = cost
    elapsed = time.perf_counter() - start 
    out.write(f"fastmultrec_greedy,{wgd},{duprate},{simid},{dup_cost},{cost},{nbdups},{nblosses},{elapsed:.3f}\n")
    
    
    
    
    #-------------------------------------------------------------------------------------
    #lca mapping
    #-------------------------------------------------------------------------------------
    cmd = f"{fastmultrecdir}/segrec -d {dup_cost} -l 1 -o out.txt -gf '{genetree_insider_filename}' -sf '{sptree_insider_filename}' -al lca"
    print(cmd)
    
    start = time.perf_counter()
    os.system(cmd)
    end = time.perf_counter()
    
    (cost, nbdups, nblosses) = get_insider_costs("out.txt")
    
    #if mrec < cost:
    #    print(simid)
    #    sys.exit()
    
    
    elapsed = time.perf_counter() - start 
    out.write(f"lcamap,{wgd},{duprate},{simid},{dup_cost},{cost},{nbdups},{nblosses},{elapsed:.3f}\n")




    















    
out.close()
