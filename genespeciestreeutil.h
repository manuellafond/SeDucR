#ifndef GENESPECIESTREEUTIL_H
#define GENESPECIESTREEUTIL_H

#include "treeutils/node.h"
#include "treeutils/util.h"

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <string>
#include <iostream>


//GNode = gene tree node, SNode = species tree node (those all use a Node object, but it helps to distinguish them)
#define GNode Node
#define SNode Node

//GSMap = gene-species map from V(G) to V(S)
#define GSMap unordered_map<GNode*, SNode*>




class GeneSpeciesTreeUtil
{

public:
    
    /**
    Return the gene leaf -> species leaf map.  It uses the leaf labels to make the association.  For each gene leaf g, 
    we read the g->label and extract its species name from it.  To do this, we split g->label split it into substrings
    according to separator, and substrings[species_index] is the species name.  Then the map finds the species with that label
    and maps it to g.
    This is not optimal: the species label is searched every time, so it's O(|V(G)||V(S)|).  Linear version is not hard to implement...
    **/
    static GSMap get_leaf_species_map(Node* gene_tree, Node* species_tree, std::string separator = "_", int species_index = 0);

    /**
    Same as get_leaf_species_map, but with multiple gene trees.  All leaves are in the same returned map.
    **/
    static GSMap get_multi_leaf_species_map(vector<GNode*> gene_trees, SNode* species_tree, string species_separator, int species_index);

};

#endif // GENESPECIESTREEUTIL_H
