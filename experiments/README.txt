
*** Running experiments on the WGD dataset ***

The WGD dataset needs to be obtained from the authors of fastmultrec.
The script run_exp.py can then be used to analyze this dataset.  
In that script, the variable 
datadir
needs to specifiy the location of the WGD dataset.
For example, it is currently set to ML's directory
datadir = "/home/manuel/git/wgddata_full"
and it should be changed with your directory. 

The directories with all the necessary software binaries also needs to be specified.  In the script, these lines must be edited:
fastmultrec_dir = "/home/manuel/git/FastMultRec/FastMultRec/build"
insider_dir = "/home/manuel/git/inSiDeR/build"
segdup_dir = ""
(here, the segdup binary was in /usr/bin/, so no path needed to be specified)


The script goes through every simphy simulation in that directory, and it runs each software on those 
(segdup, fastmultrec, inSiDeR, lca-map).
It saves one csv file per simulation directory, in the work directory (also specified in run_exp.py).  
Each method has its own line in this csv.

The script does not re-analyze a simphy directory if the csv file for it exists in the work directory.
Therefore, even if the script does not terminate, it can be re-run and it will resume where it left off.

Once all experiments are finished, the work directory should contain one csv per simphy directory.
The small csv's can be combined into one large csv that contains everything with (this is chatgpt's solution):
head -n 1 "$(ls workwgdd25/*.csv | head -n 1)" > stats_wgd_all.csv
tail -n +2 -q workwgdd25/*.csv >> stats_wgd_all.csv
where here workwgdd25 is the work directory specified in run_exp.py


*** Running experiments on the kowhai dataset ***

The script run_exp_kowhai.py runs the experiments with the kowhai simulator.  
No dataset needs to be downloaded - the script call kowhai.  
As with the WGD dataset, a work directory needs to be specified in the script, as well as the path to the binary 
of the softwares (including the path to kowhai).

If the script terminates without errors, all the data will be in stats_kowhai.csv.  If the script needs to be run multiple times, 
the csv's in the work directory can be combined in the same manner as the WGD dataset.



*** Computing average running times (WGD dataset) ***

The script make_stats.py outputs average running times for the WGD dataset, for a given duplication rate.
For example, the command
python make_stats.py --duprate=15 
will output a summary of running times for datasets with dup/loss rate 1e-15, for each duplication cost.
--duprate can take values in [8,9,10,11,15].  
These were used to make the running time plots.




*** Computing the number of suboptimal solutions per method (WGD dataset) ***

To obtain a csv containing the number of suboptimal solutions per method:
(stats_wgd_all.csv must exist)

python make_stats.py --mode=error --duprate=8 --outformat=csv > subopt.csv
python make_stats.py --mode=error --duprate=9 --outformat=csvnoheader >> subopt.csv
python make_stats.py --mode=error --duprate=10 --outformat=csvnoheader >> subopt.csv
python make_stats.py --mode=error --duprate=11 --outformat=csvnoheader >> subopt.csv
python make_stats.py --mode=error --duprate=15 --outformat=csvnoheader >> subopt.csv
