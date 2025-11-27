import matplotlib.pyplot as plt
import numpy as np

data = dict()
colmap = dict() #key = colname, value = column index, built from 1st row

seen_head = False
with open("stats.csv", "r") as f:
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
        key = f'WGD={vals[colmap["planted_wgds"]]}, duprate={vals[colmap["duprate"]]}, d={vals[colmap["dup_cost"]]}'
        
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


