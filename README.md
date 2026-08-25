# SeDucR


An exact algorithm to compute most parsimonious reconciliations with segmental duplications and losses.  
The implementation is from the paper :

Yao-Ban Chan and Manuel Lafond, "A fast and exact algorithm for gene-species reconciliation with segmental duplications and losses]{A fast and exact algorithm for gene-species reconciliation with segmental duplications and losses".


## To compile:
```
mkdir build
cd build 
cmake ..
make
```

Then run ./seducr

## Input:
```
-d              Duplication cost [default=5]
-l              Loss cost [default=1]
-o              Output file, contains optimal cost plus all reconciliation info [default=stdout]
-g              String that contains gene trees in newick format, separated by semi-colon
-gf             File that contains gene trees in newick format, one per line
-s              Species tree newick
-sf             File with species tree newick
-spsep  Separator to use in gene leaf labels [default=_]
-spindex        Index in the gene tree label with the species name, after separating with spsep [default=0]

```

The format of gene tree leaves should be speciesname_X_Y, where speciesname refers to a leaf name of the species tree, and X, Y are whatever.  

In fact, the software just takes the gene leaf labels, splits them according to the "_" character, and takes the 0-th string.
This can be configured with the "spsep" and "spindex" arguments, no time to explain now, see main.cpp code.  

## Output:
Specify the output file with -o, otherwise the output will be in stdout along with various messages.

The output contains the optimal cost found, along with the given dup cost, loss cost, species tree, and gene trees.  The internal nodes of the species tree are given a label trees, and the internal nodes of the gene trees as well.  The format of the internal gene tree node labels has the form "[species tree label]_[Dup or Spec]".  From this, the reconciliation can be reconstructed.  The output also contains, for each species, the labels of the gene tree nodes that are in a duplication in that species.

## Example:
```
./seducr -d 10 -l 1 -gf "../data/s25/sim_9/all_genetrees_edited.txt" -sf "../data/s25/sim_9/s_tree.newick"
```
