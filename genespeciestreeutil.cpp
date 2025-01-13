#include "genespeciestreeutil.h"


unordered_map<Node*, Node*> GeneSpeciesTreeUtil::get_leaf_species_map(Node* gene_tree, Node* species_tree, string separator, int species_index)
{
    unordered_map<Node*, Node*> mapping;

    //compute gene leaf to species leaf mapping.  Could be faster by using a map for the species leaves by name...
    for (auto it = gene_tree->begin(); it != gene_tree->end(); ++it)
    {
        Node* g = *it;

        if (!g->is_leaf())
            continue;

        string lbl = g->label;
        vector<string> sz = Util::Split(g->label, separator);

        if (sz.size() == 0 || (species_index > 0 && sz.size() < species_index - 1))   //there was a bug when speciesIndex = 0
        {
            cout << "Gene label " << lbl << " malformed" << endl << flush;
            throw "Gene label " + lbl + " malformed.";
        }

        // if data is from simphy
        if (!sz[species_index].empty() && g->is_leaf()) {
            // Add single quotes at the beginning and end of sz[speciesIndex]
            sz[species_index].insert(sz[species_index].begin(), std::string::value_type('\''));
            sz[species_index] += std::string::value_type('\'');
        }

        Node* s = species_tree->get_leaf_by_label(sz[species_index]);

        if (s)
        {
            mapping[g] = s;
        }
        else
        {
            string msg = "Could not find species for gene " + lbl;
            cout << msg << endl;
            throw msg;
        }
    }
    

    return mapping;
}




GSMap GeneSpeciesTreeUtil::get_multi_leaf_species_map(vector<GNode*> gene_trees, SNode* species_tree, string species_separator, int species_index)
{
    
    GSMap gsmap;

    for (int i = 0; i < gene_trees.size(); i++)
    {
        GSMap tmpmap = GeneSpeciesTreeUtil::get_leaf_species_map(gene_trees[i], species_tree, species_separator, species_index);

        for (GSMap::iterator it = tmpmap.begin(); it != tmpmap.end(); ++it)
        {
            gsmap[it->first] = it->second;
        }
    }

    return gsmap;
}


