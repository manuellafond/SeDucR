import matplotlib.pyplot as plt
import numpy as np

data = dict()
colmap = dict() #key = colname, value = column index, built from 1st row


'''
This script takes your stats file and outputs the average costs and times, 
where averages are grouped by [param1,param2,...,paramk]
where the list of params that form a group is specified by colnames_in_key.
The latter is a dict that contains all column names that are part of a key group, 
specified by key-value pairs of the form (display_name_of_param,colname_of_param)
'''



kowhai_mode = True


stats_file = "stats.csv"
if kowhai_mode:
    stats_file = "stats_kowhai.csv"


colnames_in_key = dict()

if kowhai_mode:
    colnames_in_key = { "nh" : "nh", "np" : "np" , "rb" : "rb" , "pc" : "pc", "pj" : "pj" , "d" : "dup_cost" }
else:
    colnames_in_key["WGD"] = "planted_wgds"
    colnames_in_key["duprate"] = "duprate"
    colnames_in_key["d"] = "dup_cost"


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
        
        longest_key = max(longest_key, len(key))
        
        if key not in data:
            data[key] = dict()
        
        if algo not in data[key]:
            data[key][algo] = dict()
            data[key][algo]["time"] = []
            data[key][algo]["solution_cost"] = []
            
        data[key][algo]["time"].append( float(vals[colmap["time"]]) )  
        data[key][algo]["solution_cost"].append( float(vals[colmap["solution_cost"]]) )  
        
        

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


