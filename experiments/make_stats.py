import matplotlib.pyplot as plt
import numpy as np
import argparse

data = dict()
colmap = dict() #key = colname, value = column index, built from 1st row


'''
This script takes your stats file and outputs the average costs and times, 
where averages are grouped by [param1,param2,...,paramk]
where the list of params that form a group is specified by colnames_in_key.
The latter is a dict that contains all column names that are part of a key group, 
specified by key-value pairs of the form (display_name_of_param,colname_of_param)
'''

parser = argparse.ArgumentParser()
parser.add_argument(
    "--kowhai",
    action="store_true",
    help="Set kowhai mode to True"
)

parser.add_argument(
    "--mode",
    type=str,
    default="avg",
    help="Must be avg or error.  avg mode outputs avg per parameter, error outputs number of subopt instances for each method"
)


parser.add_argument(
    "--duprate",
    type=str,
    default=None,
    help="In WGD mode, give this to only get error count for specified duprate."
)



args = parser.parse_args()

kowhai_mode = args.kowhai



stats_file = "stats_wgd.csv"
if kowhai_mode:
    stats_file = "stats_kowhai.csv"


colnames_in_key = dict()

if kowhai_mode:
    colnames_in_key = { "nh" : "nh", "np" : "np" , "rb" : "rb" , "pc" : "pc", "pj" : "pj" , "d" : "dup_cost" }
    if args.mode == "error":
        colnames_in_key["repl"] = "repl"
else:
    colnames_in_key["WGD"] = "planted_wgds"
    colnames_in_key["duprate"] = "duprate"
    colnames_in_key["d"] = "dup_cost"
    
    if args.mode == "error":
        colnames_in_key["simid"] = "simid"


def get_key(vals, colmap):
   global colnames_in_key
   
   key = ""
   for name in colnames_in_key:
       key += f"{name}={vals[colmap[colnames_in_key[name]]]};"
   return key
   


longest_key = 20
seen_head = False
with open(stats_file, "r") as f:
    for line in f:
        
        line = line.strip()
        vals = line.split(",")
        
        if not seen_head:
            for i in range(len(vals)):
                colmap[vals[i]] = i
            seen_head = True 
            continue
        
        line = line.strip()
        vals = line.split(",")
        
        if line == "":
            continue
        
        #col 0 = algo, 1 = wgd, 2 = duprate
        
        algo = vals[colmap["method"]]
        #key = f'WGD={vals[colmap["planted_wgds"]]}, duprate={vals[colmap["duprate"]]}, d={vals[colmap["dup_cost"]]}'
        key = get_key(vals, colmap)
        
        if args.duprate != None and vals[colmap["duprate"]] != args.duprate:
            continue
        
        
        longest_key = max(longest_key, len(key))
        
        if key not in data:
            data[key] = dict()
        
        if algo not in data[key]:
            data[key][algo] = dict()
            data[key][algo]["time"] = []
            data[key][algo]["solution_cost"] = []
            
        data[key][algo]["time"].append( float(vals[colmap["time"]]) )  
        data[key][algo]["solution_cost"].append( float(vals[colmap["solution_cost"]]) )  
        
        
if args.mode == "avg":
    print("Average times (seconds)")
    for key in data:
        print("\nDataset " + key)
        
        for algo in data[key]:
            mean_time = np.mean(data[key][algo]["time"])
            print(f"{algo:20} : {mean_time:0.5f}")
            
    print("\n-----------------------------------------")
    print("Average solution costs")
    for key in data:
        print("\nDataset " + key)
        
        for algo in data[key]:
            mean_cost = np.mean(data[key][algo]["solution_cost"])
            print(f"{algo:20} : {mean_cost}")
elif args.mode == "error":
    
    err_per_method = dict()
    instances_per_method = dict()
    
    for key in data:
        for algo in data[key]:
            
            if algo not in err_per_method:
                err_per_method[algo] = 0
                instances_per_method[algo] = 0
            
            instances_per_method[algo] += 1
            
            opt = data[key]["insider"]["solution_cost"][0]
            algocost = data[key][algo]["solution_cost"][0]
            
            if len(data[key][algo]["solution_cost"]) > 1:
                print(f"Warning: entry {key}.{algo} has multiple solution costs, no supposed to happen in mode error")
            
            if algocost > opt:
                err_per_method[algo] += 1
            elif algocost < opt:
                print(f"WARNING: better than insider: {key}.{algo}")
        
    print("Nb suboptimal solutions per method")
    print(err_per_method)
    
    print("Nb instances per method")
    print(instances_per_method)
        
        
    
    
    
    
    


