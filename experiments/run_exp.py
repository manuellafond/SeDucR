import sys
import os
import re
import itertools
import time
from ete3 import Tree
from pathlib import Path

#this script assumes the existence of directories of the form 
#"{datadir}/sim_{[1,2]}WGD_D{[7,10]}/sim_{[1..25]}
#and it outputs stats.csv with everything

datadir = "/home/manuel/git/wgddata_full"
output_filename = "stats_wgd_d25.csv"

workdir = "./workwgd"
os.makedirs(workdir, exist_ok=True)

skip_existing_csv = True


segdup_reps = 20000
segdup_initial_temp = 3

wgd_nb = [2, 1]
dup_rates = [7,8,9,10,11]
#dup_costs = [5, 10, 15, 20]     #values of argument -d to test
dup_costs = [25]
runs = range(1, 100+1)

#wgd_nb = [1]
#dup_rates = [7]
#dup_costs = [10]
#runs = [6]


#set to directories containing binaries of the programs to test
fastmultrec_dir = "/home/manuel/git/FastMultRec/FastMultRec/build"
insider_dir = "/home/manuel/git/inSiDeR/build"
insider_relax_dir = "/home/manuel/git/ybrec2/build"
segdup_dir = ""

#available methods here
methods = ["insider", "lca", "fastmultrec", "segdup"] 


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
    print(f"Opening {s_tree_filename}")
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
    ncopies_per_gene =dict()
    
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









out = open(output_filename, 'w')


header_line = "method,planted_wgds,duprate,simid,dup_cost,solution_cost,solution_nbdups,solution_nblosses,time\n"
out.write(header_line)

for (simid, wgd, duprate, dup_cost) in itertools.product(runs, wgd_nb, dup_rates, dup_costs):
    
    
    directory = f"{datadir}/sim_{wgd}WGD_D{duprate}/sim_{simid}"
    
    directory = f"{datadir}/WGD{wgd}Simphy_S1_G100_Duprate-{duprate}/sim_{simid}"
    genetree_file = directory + "/all_genetrees_edited.txt"
    speciestree_file = directory + "/s_tree.newick"
    
    if not os.path.exists(speciestree_file):
        print(f"SPECIES TREE FILE DOES NOT EXIST, SKIPPING {speciestree_file}")
        continue
        
    #for some reason, that input crashes segdup, so we just skip them
    if wgd == 2 and simid == 44 and duprate == 9:
        continue
    if wgd == 2 and simid == 47 and duprate == 11:
        continue
    if wgd == 2 and simid == 55 and duprate == 9:
        continue
    if wgd == 1 and simid == 79 and duprate == 11:
        continue
    if wgd == 2 and simid == 88 and duprate == 7:
        continue
    
    
    
    suffix = f"simid{simid}_wgd{wgd}_duprate{duprate}_d{dup_cost}"
    
    
    #each run produces a csv - if it already exists we can skip it
    run_csv_filename = f"{workdir}/stats_{suffix}.csv"
    if skip_existing_csv and os.path.exists(run_csv_filename):
        print(run_csv_filename + " skipped")
        continue

    run_csv_file = open(run_csv_filename, 'w')
    run_csv_file.write(header_line)
    
    
    
    
    #-------------------------------------------------------------------------------------
    #segdup
    #-------------------------------------------------------------------------------------
    #note: simphy species trees add quotes around their leaf names, they must be removed for segdup
    if "segdup" in methods:
        stree_str = '-S "'
        stree_str += get_segdup_stree_str(speciestree_file)
        stree_str += '"'
        
        segdup_outfilename = f"{workdir}/segdupout_{suffix}.txt"

        genetree_segdup_outfilename = f"{workdir}/gfile_{suffix}.txt"
        genetree_outfile = open(genetree_segdup_outfilename, "w")
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
        #command = f"cat {segdup_input_filename} | xargs {segdup_dir}segdup"
        #command = f"{segdup_dir}segdup {stree_str} {gtree_str} -n 100"
        #print(command)
        
        command = f"{segdup_dir}segdup {stree_str} -Gfile {genetree_segdup_outfilename} -d {dup_cost} -n {segdup_reps} -Tinit {segdup_initial_temp} -Tfinal 0 > {segdup_outfilename}"
        print(command)
        
        #if not os.path.exists(segdup_outfilename):
        start = time.perf_counter()
        os.system(command)
        elapsed = time.perf_counter() - start 


        
        #get segdup's cost line, which has the form nCospec,nSegDup,nLoss,cost
        (nbdups, nblosses, cost) = get_segdup_costs(segdup_outfilename)
        
        
        out.write(f"segdup,{wgd},{duprate},{simid},{dup_cost},{cost},{nbdups},{nblosses},{elapsed:.6f}\n")
        run_csv_file.write(f"segdup,{wgd},{duprate},{simid},{dup_cost},{cost},{nbdups},{nblosses},{elapsed:.6f}\n")

    
    
    #-------------------------------------------------------------------------------------
    #get insider data (aka ybrec), relax or not
    #-------------------------------------------------------------------------------------
    if "insider" or "insider_relax" in methods:
        genetree_insider_filename = f"{workdir}/gfile_insider_{suffix}.txt"
        genetree_insider_file = open(genetree_insider_filename, "w")
        
        snewick = get_segdup_stree_str(speciestree_file)
        sptree_insider_filename = f"{workdir}/sfile_insider_{suffix}.txt"
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
        
        if "insider" in methods:
            outfilename = f"{workdir}/insider_out_{suffix}.txt"
            #cmd = f"{insider_dir}/ybrec -d {dup_cost} -l 1 -o out.txt -gf '{genetree_file}' -sf '{speciestree_file}'"
            cmd = f"{insider_dir}/ybrec -d {dup_cost} -l 1 -o {outfilename} -gf '{genetree_insider_filename}' -sf '{sptree_insider_filename}'"
            print(cmd)
        
            start = time.perf_counter()
            retcode = os.system(cmd)
            end = time.perf_counter()
            elapsed = time.perf_counter() - start 
        
            print(f"inSiDeR return with code {retcode}")
            if retcode < 0:
                print("but it failed...")
        
            #if out file does not exist, we assume that inSiDeR was killed
            if not Path(outfilename).exists():
                (cost, nbdups, nblosses) = (-1, -1, -1)
                elapsed = 9999999999
                ybcost = -1
            else:
                (cost, nbdups, nblosses) = get_insider_costs(outfilename)
                ybcost = cost

            out.write(f"insider,{wgd},{duprate},{simid},{dup_cost},{cost},{nbdups},{nblosses},{elapsed:.6f}\n")
            run_csv_file.write(f"insider,{wgd},{duprate},{simid},{dup_cost},{cost},{nbdups},{nblosses},{elapsed:.6f}\n")
        
        if "insider_relax" in methods:
            #cmd = f"{insider_dir}/ybrec -d {dup_cost} -l 1 -o out.txt -gf '{genetree_file}' -sf '{speciestree_file}'"
            outfilename_relax = f"{workdir}/insiderrelax_out_{suffix}.txt"
            cmd = f"{insider_relax_dir}/ybrec -d {dup_cost} -l 1 -o {outfilename_relax} -gf '{genetree_insider_filename}' -sf '{sptree_insider_filename}'"
            print(cmd)
        
            start = time.perf_counter()
            os.system(cmd)
            end = time.perf_counter()
            elapsed = time.perf_counter() - start 
        
            (cost_relax, nbdups_relax, nblosses_relax) = get_insider_costs(outfilename_relax)
            ybcost_relax = cost_relax

            out.write(f"insider_relax,{wgd},{duprate},{simid},{dup_cost},{cost_relax},{nbdups_relax},{nblosses_relax},{elapsed:.6f}\n")
            run_csv_file.write(f"insider_relax,{wgd},{duprate},{simid},{dup_cost},{cost_relax},{nbdups_relax},{nblosses_relax},{elapsed:.6f}\n")
            
            if "insider" in methods and ybcost != ybcost_relax:
                print(f"Difference insider trie vs relax")
                print(genetree_insider_filename)
                print(f"trie: cost {ybcost}   relax: cost {ybcost_relax}")
                sys.exit()
        
    
    #-------------------------------------------------------------------------------------
    #fastmultrec
    #-------------------------------------------------------------------------------------
    if "fastmultrec" in methods:
        outfilename_fmr = f"{workdir}/fmr_out_{suffix}.txt"
        #cmd = f"{fastmultrec_dir}/segrec -d {dup_cost} -l 1 -o out.txt -gf '{genetree_file}' -sf '{speciestree_file}'"
        cmd = f"{fastmultrec_dir}/segrec -d {dup_cost} -l 1 -o {outfilename_fmr} -gf '{genetree_insider_filename}' -sf '{sptree_insider_filename}'"
        print(cmd)
        
        start = time.perf_counter()
        os.system(cmd)
        end = time.perf_counter()
        
        (cost, nbdups, nblosses) = get_insider_costs(outfilename_fmr)
        mrec = cost
        elapsed = time.perf_counter() - start 
        out.write(f"fastmultrec_greedy,{wgd},{duprate},{simid},{dup_cost},{cost},{nbdups},{nblosses},{elapsed:.6f}\n")
        run_csv_file.write(f"fastmultrec_greedy,{wgd},{duprate},{simid},{dup_cost},{cost},{nbdups},{nblosses},{elapsed:.6f}\n")
    
    
    
    
    #-------------------------------------------------------------------------------------
    #lca mapping
    #-------------------------------------------------------------------------------------
    if "lca" in methods:
        outfilename_lca = f"{workdir}/lca_out_{suffix}.txt"
        cmd = f"{fastmultrec_dir}/segrec -d {dup_cost} -l 1 -o {outfilename_lca} -gf '{genetree_insider_filename}' -sf '{sptree_insider_filename}' -al lca"
        print(cmd)
        
        start = time.perf_counter()
        os.system(cmd)
        end = time.perf_counter()
        
        (cost, nbdups, nblosses) = get_insider_costs(outfilename_lca)
        
        #if mrec < cost:
        #    print(simid)
        #    sys.exit()
        
        
        elapsed = time.perf_counter() - start 
        out.write(f"lcamap,{wgd},{duprate},{simid},{dup_cost},{cost},{nbdups},{nblosses},{elapsed:.6f}\n")
        run_csv_file.write(f"lcamap,{wgd},{duprate},{simid},{dup_cost},{cost},{nbdups},{nblosses},{elapsed:.6f}\n")
        
    run_csv_file.close()




    















    
out.close()
