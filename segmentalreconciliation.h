#ifndef RECONCILIATION_H
#define RECONCILIATION_H

#include "genespeciestreeutil.h"
#include "treeutils/newicklex.h"


class SegmentalReconciliation{
public:
	//everything public for now, ftw

	
	SNode* species_tree;
    vector<GNode*> gene_trees;
	GSMap gsmap;
	
	map<SNode*, set<set<GNode*>>> dups_per_species;

    int loss_cost;
    int dup_cost;
	
	
	void init_with_gsmap(SNode* sptree, vector<GNode*> gtrees, GSMap& initial_map){
		species_tree = sptree;
		gene_trees = gtrees;
		gsmap = initial_map;	//Note: makes a copy
	}
	
	
	void build_from_segmental_dups(map<SNode*, set< set<GNode*> > >& _dups_per_species){
		
		this->dups_per_species = _dups_per_species;
		
		//map dup nodes to where they should be (other maps remain unchanged).
		for (auto keyval : dups_per_species){
			
			SNode* x = keyval.first;
			
			for (const set<GNode*>& dupnodes : keyval.second){
				
				for (GNode* g : dupnodes){
					gsmap[g] = x;
					
					if (g->is_leaf() && !x->is_leaf()){
						cout<<"WAT???"<<endl;
					}
				}
			}				
			
		}
	}
	
	
	int get_nb_losses(){
		int nblosses = 0;

		for (int i = 0; i < gene_trees.size(); i++) {
			GNode* genetree = gene_trees[i];

			for (auto it = genetree->begin(); it != genetree->end(); ++it){
				GNode* g = *it;
				if (!g->is_leaf()) {
					bool isdup = this->is_duplication(g);

					int d1 = gsmap[g]->get_distance_to(gsmap[g->get_child(0)]);
					int d2 = gsmap[g]->get_distance_to(gsmap[g->get_child(1)]);

					int losses_tmp = d1 + d2;
					if (!isdup) {
						losses_tmp -= 2;
					}

					nblosses += losses_tmp;
				}
			}
		}

		return nblosses;
	}
	
	
	
	bool is_duplication(GNode* g)
	{
		if (g->is_leaf())
			return false;

		GNode* s = gsmap[g];
		GNode* s1 = gsmap[g->get_child(0)];
		GNode* s2 = gsmap[g->get_child(1)];

		if (s1->has_ancestor(s2) || s2->has_ancestor(s1))
			return true;

		if (s != s1->get_lca_with(s2))
			return true;

		return false;
	}
	
	
	int get_nb_segmental_dups(){
		int cpt = 0;
		for (auto keyval : dups_per_species){
			cpt += keyval.second.size();
		}
		return cpt;
	}
	
	
	
	
	void label_gene_internal_nodes(){
		int cpt = 0;
		for (GNode* gtree : gene_trees){
			for (auto it = gtree->begin(); it != gtree->end(); ++it){
				GNode* g = *it;
				
				string label = Util::ToString(cpt);
				
				label += "_" + gsmap[g]->label;
				
				if (is_duplication(g))
					label += "_Dup";
				else 
					label += "_Spec";
				
				g->label = label;
				
				++cpt;
			}
		}
		
	}
	
	
	
	
	//returns the height of the dup subtree in x rooted at g
	int get_dup_height_under(GNode* g, SNode* x)
	{
		if (gsmap[g] != x || !is_duplication(g))
			return 0;

		int d1 = get_dup_height_under(g->get_child(0), x);
		int d2 = get_dup_height_under(g->get_child(1), x);

		return 1 + max(d1, d2);
	}
	
	
	//returns the sum of the heights of the duplication forests at easch species (so, the number of seg dups)
	//naive and slow algorithm
	int get_dup_height_sum() {
		
		map<SNode*, int> dup_heights;
		
		for (auto it = species_tree->begin(); it != species_tree->end(); ++it) {
			dup_heights[*it] = 0;
		}
		


		for (int i = 0; i < gene_trees.size(); i++) {
			GNode* gtree = gene_trees[i];

			for (auto it = gtree->begin(); it != gtree->end(); ++it) {
				GNode* g = *it;
				if (is_duplication(g)){
					SNode* s = gsmap[g];
					int dupheight = get_dup_height_under(g, s);
					dup_heights[s] = max(dup_heights[s], dupheight);
				}
			}
		}


		
		int dup_height_sum = 0;
		for (auto it = species_tree->begin(); it != species_tree->end(); ++it) {
			dup_height_sum += dup_heights[*it];
		}

		return dup_height_sum;
	}
	
	
	
	
	bool is_time_consistent(){
		for (int i = 0; i < gene_trees.size(); i++) {
			GNode* gtree = gene_trees[i];

			for (auto it = gtree->begin(); it != gtree->end(); ++it) {
				GNode* g = *it;
				
				if (g->is_leaf()){
					if (!gsmap[g]->is_leaf()){
						cout<<"Not time consistent: gene leaf " << g->label << " is mapped to non-leaf " << gsmap[g]->label<<endl;
						return false;
					}
				}
				else{
					for (int i = 0; i < g->get_nb_children(); ++i){
						if (! gsmap[g->get_child(i)]->has_ancestor(gsmap[g]) ){
							cout<<"Not time consistent: gene " << g->label << " has child " << i << " not mapped to a descendant"<<endl;
							return false;
						}
			
					}
				}
			}
		}	//closing brackets galore!
		return true;
		
	}
	
	
	
	
	
	
	string get_output_string(bool label_gene_nodes = true){
		
		if (label_gene_nodes)
			label_gene_internal_nodes();
		
		int nblosses = get_nb_losses();
		int nbdups = get_nb_segmental_dups();
		
		string output = "";
		
		output += "<DUPCOST>\n" + Util::ToString(dup_cost) + "\n<DUPCOST>\n";
		output += "<LOSSCOST>\n" + Util::ToString(loss_cost) + "\n<LOSSCOST>\n";
		output += "<COST>\n" + Util::ToString(nbdups * dup_cost + nblosses * loss_cost) + "\n</COST>\n";
        output += "<DUPHEIGHT>\n" + Util::ToString(nbdups) + "\n</DUPHEIGHT>\n";
        output += "<NBLOSSES>\n" + Util::ToString(nblosses) + "\n</NBLOSSES>\n";
        output += "<SPECIESTREE>\n" + NewickLex::ToNewickString(species_tree) + "\n</SPECIESTREE>\n";

        output += "<GENETREES>\n";
        for (int t = 0; t < gene_trees.size(); t++)
        {
            output += NewickLex::ToNewickString(gene_trees[t]) + "\n";
        }
        output += "</GENETREES>\n";


		
        output += "<DUPS_PER_SPECIES>\n";
		
		map<GNode*, int> node_to_genetree_nb;
		for (int i = 0; i < gene_trees.size(); ++i){
			for (auto git = gene_trees[i]->begin(); git != gene_trees[i]->end(); ++git){
					node_to_genetree_nb[*git] = i;
			}
		}
		
		
		for (auto keyval : dups_per_species){
			
			SNode* x = keyval.first;
			
			output += "[" + x->label + "] ";
			
			set< set<GNode*> >& dups = keyval.second;
			for (const set<GNode*>& dup_nodes : dups){
				for (GNode* g : dup_nodes){
					
					output += g->label + " (G" + Util::ToString(node_to_genetree_nb[g]) + ") ";
				}
			}
			
			output += "\n";
		}
		
        output += "</DUPS_PER_SPECIES>\n";
		
		return output;
	}
	
	
};







#endif