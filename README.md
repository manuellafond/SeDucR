# ybrec

A tentative implementation of Yao-Ban's algorithm for reconciliation with segmental duplications.

##To compile:
mkdir build
cd build 
cmake ..
make

Then run ./ybrec

##Input:
```
-sf: file containing species tree newick
-gf: file containing list of gene tree newick, one gene tree per line
-d: dup cost 
-l: loss cost 
```

The format of gene tree leaves should be speciesname_X_Y, where speciesname refers to a leaf name of the species tree, and X, Y are whatever.  

In fact, the software just takes the gene leaf labels, splits them according to the "_" character, and takes the 0-th string.
This can be configured with the "spsep" and "spindex" arguments, no time to explain now, see main.cpp code.  

##Output:
Just the min reconciliation that the algorithm found.  More precisely, it outputs the set of c(x, V) values, where V = all gene tree roots.

The rest is debugging message, which includes the number iterations needed to build the A_x sets from A_{x_l}, A_{x_r}, and the number of iterations to process them. 

##Example:
```
./ybrec -gf "../data/s25/sim_9/all_genetrees_edited.txt" -sf "../data/s25/sim_9/s_tree.newick" -d 5 -l 1
```