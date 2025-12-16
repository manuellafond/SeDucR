#pragma once

#include "treeutils/treeutil.h"
#include <set>
#include <vector>
#include <map>
#include <unordered_map>
#include <list>
#include <string>
#include <queue>

#include "treeutils/node.h"
#include "treeutils/newicklex.h"

//NOTE: GNode, SNode and GSMap are defined in genespeciestreeutil
#include "genespeciestreeutil.h"

//This defines set_gnodes type as a set of GNode* objects
#include "activesettrie.h"

#include "segmentalreconciliation.h"




/*
**********************************************************************************
* Functionalities related to bitmaps, for compressed set representations
**********************************************************************************
*/
#include "treeutils/ewah/ewah.h"
typedef ewah::EWAHBoolArray<uint32_t> bitmap;

//must return true iff a < b (if equal, must return false)
//a comparator function for bitmaps, so that we can store them in sets and maps
struct BCmp {
    bool operator()(const bitmap& a, const bitmap& b) const {

        if (a.numberOfOnes() > 0 && b.numberOfOnes() == 0)
            return false;
        if (a.numberOfOnes() == 0 && b.numberOfOnes() == 0)
            return false;
        if (a.numberOfOnes() == 0 && b.numberOfOnes() > 0)
            return true;

        auto ita = a.begin();
        auto itb = b.begin();


        while (true) {  //never do a while (true)
            if (ita != a.end()) {
                if (itb != b.end()) {
                    if (*ita < *itb)
                        return true;
                    else if (*itb < *ita)
                        return false;
                    else {
                        ++ita;
                        ++itb;
                    }
                }
                else  //itb == b.end()
                    return false;
            }
            else {
                if (itb != b.end()) //a finished before b
                    return true;
                else
                    return false;    //both ended, a == b
            }
        }
    }
};











/*
**********************************************************************************
* YBRec class
**********************************************************************************
Stores all info necessary to apply the YB dynamic programming algorithm, given
a set of gene trees, a species tree, and d/l costs

todo: for now everything is just public, cleanup later
*/
class YBRec {
public:

    SNode* speciestree;
    vector<GNode*> genetrees;

    int loss_cost;
    int dup_cost;

    YBCost cost_upper_bound;

    //max_remap_dist is the distance that a gene tree node can be remmapped to wrt lca.  Default = 5
    //eg if max_remap_dist = 3, then a gene node can only be moved "up" by 3 from its lca mapping.
    int max_remap_dist;

    //the A_x sets from the YB algorithm.  
    //key = species
    //value = the active sets for that species (each set is represented as a Trie)
    map< SNode*, Trie > active_sets;

    GSMap leafmap;
    GSMap lcamap;

    map< SNode*, set_gnodes > rev_leafmap;     //rev_leafmap[x] = set of gene tree leaves mapped to s
    map<int, GNode*> id_genes;                  //id_genes[i] = the gene tree node with id = i
									   
    map< SNode*, set_gnodes > rev_lcamap;     //same as rev_leafmap, but for LCA
    map< SNode*, int > lca_dupheight;

	map< int, GNode* > gnode_by_id;

    //stores info on deleted active sets, since they may be needed during backtracking
    //to save space, they are stored in bitmap format
    map<SNode*, map<set_gnodes, YBCost> > costs_table;
    

    bool DEBUG;

    YBRec() {
        DEBUG = false;
        dup_cost = 5;
        loss_cost = 1;
        max_remap_dist = 100; //changed by YBC
    }




    
	
	
	//Assumes that cost_exists(x, V) is true
	void set_cost_origin(SNode* x, set_gnodes& V, SNode* orig_x, set_gnodes& orig_V){
		YBCost& ybcost = get_cost(x, V);	//careful not to copy the YBCost object
		ybcost.nb_originators = 1;
		ybcost.origin1 = make_pair(orig_x, orig_V);
	}
	
	//Assumes that cost_exists(x, V) is true	
	void set_cost_origins(SNode* x, set_gnodes& V, SNode* orig_x1, set_gnodes& orig_V1, SNode* orig_x2, set_gnodes& orig_V2){
		YBCost& ybcost = get_cost(x, V);	//careful not to copy the YBCost object
		ybcost.nb_originators = 2;
		ybcost.origin1 = make_pair(orig_x1, orig_V1);
		ybcost.origin2 = make_pair(orig_x2, orig_V2);
	}


    //sets c(x, V).  Updates only if cost is smaller than previous entry, or if such an entry does not exist.
    // We must specify the number of dups and losses.  
	// return true if costs table was updated with new value, false otherwise
    bool set_cost(SNode* x, set_gnodes& gnodes, int nb_dups, int nb_losses) {

        bool update_cost = false;
        int newcost = this->dup_cost * nb_dups + this->loss_cost * nb_losses;
        
        if (cost_exists(x, gnodes)) {
            YBCost& best_cost = get_cost(x, gnodes);
            if (newcost < best_cost.cost)
                update_cost = true;
        }
        else
            update_cost = true;

        if (update_cost) {
			
            YBCost new_ybcost;
            new_ybcost.nb_dups = nb_dups;
            new_ybcost.nb_losses = nb_losses;
            new_ybcost.cost = newcost;
            
            costs_table[x][gnodes] = new_ybcost;
        }
		
		return update_cost;
        
    }




	//returns true iff cost has been set for species/active_set combo
	bool cost_exists(SNode* x, set_gnodes& gnodes){
		return costs_table.count(x) && costs_table[x].count(gnodes);
	}


    //returns a reference to the cost object associated with species/active_set combo.
	//ASSUMES THAT cost_exists(x, gnodes) is true, does not do that verification!
    YBCost& get_cost(SNode* x, const set_gnodes& gnodes) {
         return costs_table[x][gnodes];
    }


    //preprocessing function to fill in rev_leafmap
    void compute_rev_leafmap() {

        for (GNode* g_root : genetrees) {

            auto it = g_root->begin();
            while (it != g_root->end()) {
                Node* g = *it;
                if (g->is_leaf()) {
                    rev_leafmap[leafmap[g]].insert(g);
                }
                ++it;
            }

        }
    }


    //preprocessing to compute lcamap
    //Now also computes rev_lcamap and lca_dupheight
    void compute_lca_map() {
        for (GNode* g_root : genetrees) {
		  for (auto it = g_root->begin(); it != g_root->end(); ++it) {
                Node* g = *it;

                if (g->is_leaf()) {
                    lcamap[g] = leafmap[g];
                }
                else {
                    lcamap[g] = lcamap[g->get_child(0)]->get_lca_with(lcamap[g->get_child(1)]);
                }

			 rev_lcamap[lcamap[g]].insert(g);
            }
        }

	   //Count duplication heights
      for (GNode* g_root : genetrees) {
        for (auto it = g_root->begin(); it != g_root->end(); ++it) {
            Node* g = *it;

			int height = 0;

			//looks like it counts speciations, but it doesn't: the lowest node is always a speciation, so -1 to dup height
			while (!g->is_root() && lcamap[g->get_parent()] == lcamap[g]) {
			  g = g->get_parent();
			  height++;
			}

			if (!lca_dupheight[lcamap[g]])
				lca_dupheight[lcamap[g]] = height;
			else
				lca_dupheight[lcamap[g]] = max(height, lca_dupheight[lcamap[g]]);
        }
      }
    }



    //returns V_l cross V_r
    set_gnodes get_cross_set(const set_gnodes& V_l, const set_gnodes& V_r) {
        set_gnodes ret;

        //pass 1: for every gl in V_l, check if gl's sibling is in V_r.  If so, they get "fused" 
        //and the common parent is inserted.  If not, gl is kept.
        for (GNode* gl : V_l) {
            GNode* gr = gl->get_sibling();

            if (gr && V_r.count(gr)) {
                ret.insert(gl->get_parent());
            }
            else {
                ret.insert(gl);
            }
        }

        //pass 2: add every gr in V_r whose parent was not added.
        for (GNode* gr : V_r) {
            if (!gr->get_parent() || !ret.count(gr->get_parent())) {
                ret.insert(gr);
            }
        }

        return ret;
    }




    //returns the set of nodes whose children are both in V
    set_gnodes get_common_parents(const set_gnodes& V) {
        set_gnodes ret;

        for (GNode* g1 : V) {
            GNode* g2 = g1->get_sibling();

            if (g2 && V.count(g2)) {
                ret.insert(g1->get_parent());
            }
        }

        return ret;
    }




    //get number of edges from lower_sp to higher_sp
    //TODO: make this function faster - also, lower_sp MUST be a descendant of higher_sp
    int get_species_distance(SNode* lower_sp, SNode* higher_sp) {
        int d = 0;

        while (lower_sp != higher_sp) {
            lower_sp = lower_sp->get_parent();
            d++;
        }
        return d;
    }


    //compute cost of LCA map as upper bound
    void compute_lca_upper_bound() {
	    int dups = 0;
	    int losses = 0;


	    //count duplications
         vector<SNode*> species_nodes = speciestree->get_postordered_nodes();
	    for (SNode* x : species_nodes) {
		    dups += lca_dupheight[x];
	    }

	    //count losses
	    for (GNode* g_root : genetrees) 
		    for (auto it = g_root->begin(); it != g_root->end(); ++it) {
			    Node* g = *it;

			    if (!g->is_root()) {
				    //probably not the most efficient
				    int l = get_species_distance(lcamap[g],lcamap[g->get_parent()]);

				    if (l > 0 && lcamap[g->get_sibling()] != lcamap[g->get_parent()])
					    //parent of g is speciation
					    losses += l-1;
				    else
					    losses += l;
			    }
		    }

	    //store
	    cost_upper_bound.nb_dups = dups;
	    cost_upper_bound.nb_losses = losses;
	    cost_upper_bound.cost = dup_cost * dups + loss_cost * losses;
    }




    //preprocessing to assign a unique id to each gnode and snode.  This fills id_genes map.
    void init_indices() {
        auto snodes = speciestree->get_postordered_nodes();
        for (int i = 0; i < snodes.size(); ++i) {
            snodes[i]->id = i;
        }


        int cpt = 0;
        for (GNode* groot : genetrees) {
            auto gnodes = groot->get_preordered_nodes();
            for (GNode* g : gnodes) {
                g->id = cpt;
                id_genes[cpt] = g;
                cpt++;
            }
        }
    }


	/**
	 Returns true if this is descended from V - cannot be IN V
	 **/
    bool is_descendant_to(GNode* g, set_gnodes V) {
	    Node* cur = g;

	    while (!cur->is_root()) {
		    cur = cur->get_parent();

		    if (V.find(cur) != V.end())
			    return true;
	    }

	    return false;
    } 

	//return true if all elements of V are descendant from elements of Vprime
	bool is_below(set_gnodes& V, set_gnodes& Vprime) {
		
		bool is_below = true;

		for (auto n : V) {
			GNode* cur = n;

			while (Vprime.find(cur) == Vprime.end()) {
				if (cur->is_root())
					return false;
				
				cur = cur->get_parent();
			}
		}

		return true;
	}

    //Return true if (x,V) can be removed
    bool bound(set_gnodes& V, SNode* x, bool allow_events_at_x) {
      // return false;

      vector<SNode*> species_nodes = speciestree->get_postordered_nodes();
      set_gnodes V_set = V;	//todo trie:  update this
	  
	  if (!cost_exists(x, V))
		  throw "Trying to find a lower bound for a non-existing active set";
	  
      YBCost lower_bound = get_cost(x, V);	//note: this copies the current cost object, so modifying lower_bound is safe

      //Calculate lower bound for x, V
      //Set feasible reconciliation (LCA completion)
      map< GNode*, SNode* > newrec;
      for (GNode* g_root : genetrees)
      for (auto it = g_root->begin(); it != g_root->end(); ++it) {
        Node *g = *it;

        //start by setting elements of V mapped to x (if events at x allowed) parent of x (if events at x not allowed), then LCA mapping...
        if (V_set.find(g) != V_set.end()) {
          if (allow_events_at_x)
            newrec[g] = x;
          else
            newrec[g] = x->get_parent();
        }
        else if (g->is_leaf())
          newrec[g] = leafmap[g];
        else
          newrec[g] = newrec[g->get_child(0)]->get_lca_with(newrec[g->get_child(1)]);
      }
      //then drop elements of V down to x (if events at x not allowed)
      if (!allow_events_at_x)
        for (GNode* g : V_set)
           newrec[g] = x;


      //Bound duplications
      int dups = 0;
      for (GNode* g_root : genetrees)
      for (auto it = g_root->begin(); it != g_root->end(); ++it) {
        Node *g = *it;

        //g is in V or a leaf not descended from V
        if (V_set.find(g) != V_set.end() || (g->is_leaf() && !is_descendant_to(g, V_set)) ) {
          Node* cur = g;
          int dupheight = 0;

          while (!cur->is_root()) {
			//if cur->get_parent is a duplication in the LCA completion
			if ( newrec[cur->get_parent()] != lcamap[cur->get_parent()] || (newrec[cur->get_parent()] == newrec[cur] || newrec[cur->get_parent()] == newrec[cur->get_sibling()]) )
				dupheight++;

			cur = cur->get_parent();
          }

          dups = max(dups, dupheight);
        }
      }
      lower_bound.nb_dups += dups;
      lower_bound.cost += dups * dup_cost;


      //Bound losses
      //calculate losses, update lower bound
      int losses = 0;
      for (GNode* g_root : genetrees)
      for (auto it = g_root->begin(); it != g_root->end(); ++it) {
        Node *g = *it;

        if (!g->is_root() && !is_descendant_to(g, V_set)) {

          int l = get_species_distance(newrec[g], newrec[g->get_parent()]);

          if (l > 0 && newrec[g->get_sibling()] != newrec[g->get_parent()] && newrec[g->get_parent()] == lcamap[g->get_parent()])
            //parent of g is speciation
            losses += l-1;
          else
            losses += l;
        }
      }
      lower_bound.nb_losses += losses;
      lower_bound.cost += losses * loss_cost;


      //calculate duplications in feasible solution
      map< SNode*, int > feasible_dupheight;

      for (GNode* g_root : genetrees)
      for (auto it = g_root->begin(); it != g_root->end(); ++it) {
        Node *g = *it;

        if (!is_descendant_to(g, V_set)) {
          int height = 0;

          //corner case where a duplication is counted as a speciation mistakenly
          if (newrec[g] == x->get_parent() && newrec[g->get_child(0)] == x && newrec[g->get_child(1)] == x)
            height++;

          while (!g->is_root() && newrec[g->get_parent()] == newrec[g]) {
            g = g->get_parent();
            height++;
          }

          if (!feasible_dupheight[newrec[g]])
            feasible_dupheight[newrec[g]] = height;
          else
            feasible_dupheight[newrec[g]] = max(height, feasible_dupheight[newrec[g]]);
        }
      }


      //use feasible reconciliation to recalculate upper bound
      YBCost upper_bound = get_cost(x, V);
      upper_bound.nb_losses += losses;
      upper_bound.cost += losses * loss_cost;
      int feasible_dups = 0;
      for (SNode* y : species_nodes) {
        if (!y->has_ancestor(x) || (allow_events_at_x && y == x)) {
          feasible_dups += feasible_dupheight[y];
        }
      }
      upper_bound.nb_dups += feasible_dups;
      upper_bound.cost += feasible_dups * dup_cost;
      if (upper_bound.cost < cost_upper_bound.cost) {
        // cout << "updating upper bound from " << cost_upper_bound.cost << " to " << upper_bound.cost << endl;
        // cout << upper_bound.nb_dups << '\t' << upper_bound.nb_losses << endl;
        // if (allow_events_at_x)
        //   cout << "allow x" << endl;
        cost_upper_bound = upper_bound;
      }

      //Remove (return true) if x, V cannot result in optimal reconciliation
      if (lower_bound.cost > cost_upper_bound.cost)
        return true;
      else
        return false;
    }



    // Custom comparator for std::set<int>
    struct set_gnodes_cmp {
        bool operator()(const set_gnodes& s1, const set_gnodes& s2) const {
            if (s1.empty())
                return !s2.empty();
            
            auto it1 = s1.begin();
            auto it2 = s2.begin();
            while (it1 != s1.end() && it2 != s2.end()) {
                if ((*it1)->id < (*it2)->id)
                    return true;
                if ((*it2)->id < (*it1)->id)
                    return false;
                ++it1;
                ++it2;
            }

            return (it1 == s1.end() && it2 != s2.end());
        }
    };


    SegmentalReconciliation reconcile(int bound_option) {

        //max_remap_dist = ceil((float)dup_cost / (float)loss_cost);

        init_indices();
        compute_rev_leafmap();
        compute_lca_map();
	    compute_lca_upper_bound();

	    cout << "LCA upper bound = " << cost_upper_bound.cost << endl;

        vector<SNode*> species_nodes = speciestree->get_postordered_nodes();
        for (SNode* x : species_nodes) {

            if (x->is_leaf()) {
                active_sets[x].insert(rev_leafmap[x]);
                set_cost(x, rev_leafmap[x], 0, 0);
            }
            else {
                SNode* xl = x->get_child(0);
                SNode* xr = x->get_child(1);

				size_t output_counter = 0;
                if (active_sets[xl].size() == 0 || active_sets[xr].size() == 0) {
                    cout << "WARNING: active set has size 0, anything could happen" << endl;
                }
                //build all the possible active sets from those of x's children
                for (auto it = active_sets[xl].begin(); it != active_sets[xl].end(); ++it){
					for (auto it2 = active_sets[xr].begin(); it2 != active_sets[xr].end(); ++it2){
						set_gnodes& V_l = *it;
						set_gnodes& V_r = *it2;
					
                        set_gnodes v_cross = get_cross_set(V_l, V_r);
                        
                        int nb_roots_vl_vr = 0;
                        for (GNode* tempnode : V_l) {
                            if (tempnode->is_root())
                                nb_roots_vl_vr++;
                        }
                        for (GNode* tempnode : V_r) {
                            if (tempnode->is_root())
                                nb_roots_vl_vr++;
                        }



                        //so, v_cross is the set we get if we pair xl and xr guys into a speciation
                        //when they have a common parent, and "raise" all the other ones.
                        //Each such raise causes one loss that must be counted, hence the second line of the cost
                        
                        YBCost& leftcost = get_cost(xl, V_l);
                        YBCost& rightcost = get_cost(xr, V_r);
                        
                        int extra_losses = (2 * v_cross.size() - V_l.size() - V_r.size());
                        extra_losses -= nb_roots_vl_vr;
						
						
                        bool update = set_cost(x, v_cross, leftcost.nb_dups + rightcost.nb_dups,
												leftcost.nb_losses + rightcost.nb_losses + extra_losses);
                        
						
                        if (update) {
                            active_sets[x].insert(v_cross);
                            set_cost_origins(x, v_cross, xl, V_l, xr, V_r);
                        }
                        
            
                        //apply bounding here - YBC 
						//if bound is bad, active set is erased
                        if ( ! (bound_option < 2 || !bound(v_cross, x, true))){
                          active_sets[x].erase(v_cross);
						}

						++output_counter;
                        if (output_counter % 100000 == 0) {
                            cout << "output_counter="<<output_counter <<"  Sp x=" << x->id << "  Ax size is now " << active_sets[x].size();
                            cout << " Axl=" << active_sets[xl].size() << "  Axr=" << active_sets[xr].size() << endl;
                        }
                    }
                }
            }


		    
			//list<set_gnodes> active_sets_queue;
            //for (auto it = active_sets[x].begin(); it != active_sets[x].end(); ++it) {
            //    active_sets_queue.push_back(*it);
            //}
            priority_queue<set_gnodes, vector<set_gnodes>, set_gnodes_cmp> active_sets_queue;
            for (auto it = active_sets[x].begin(); it != active_sets[x].end(); ++it) {
                active_sets_queue.push(*it);
            }


            int nb_iter = 0;

            while (!active_sets_queue.empty()) {
                
				//set_gnodes V = active_sets_queue.front();
				//active_sets_queue.pop_front();
                set_gnodes V = active_sets_queue.top();
                active_sets_queue.pop();
                

				//V was erased, so no point in checking it again
                //TODO: BACKCHECK IF THIS IS CORRECT
				if (!active_sets[x].search(V)){
					continue;
				}

                set_gnodes U = get_common_parents(V);

                if (!U.empty()) {  //if there are actually dups that can be applied

                    set_gnodes Vprime = get_cross_set(V, V);   //not sure that works, I think it does - so it probably doesn't work
                    
					if (!cost_exists(x, V)){
						cout<<"ERROR, cost for "<<x->label<<" does not exist"<<endl;
						continue;
					}
					YBCost& xV_cost = get_cost(x, V);
					int vprime_cost = INT_MAX;

					if (cost_exists(x, Vprime)){
						YBCost& xVprime_cost = get_cost(x, Vprime);
						vprime_cost = xVprime_cost.cost;
					}
						
					//remove previous set if new set is already better (or as good) than previous set
					if (vprime_cost <= xV_cost.cost)
						active_sets[x].erase(V);
					else{
						//add dup unless new set already has a better cost than previous set + dup
						if (xV_cost.cost + dup_cost < vprime_cost) {
							bool updated = set_cost(x, Vprime, xV_cost.nb_dups + 1, xV_cost.nb_losses);
							if (updated){
                                active_sets[x].insert(Vprime);
								set_cost_origin(x, Vprime, x, V);
							}
							
							//apply bounding here - YBC
							//if cost is good, add to the queue, otherwise remove it
							if (bound_option < 2 || !bound(Vprime, x, true)) {
								//active_sets_queue.push_back(Vprime);
                                active_sets_queue.push(Vprime);
							}
							else{
								active_sets[x].erase(Vprime);
                                cout << "Bounding caught it 2" << endl;
							}
						}
					}
							


					bool is_U_forced = false;
					//if someone is trying to go too far, we have to do the dup here
					//also have to do it if it's the root of the species tree
					for (GNode* g : U) {
                        if (x->is_root() || (float)get_species_distance(lcamap[g], x) == max_remap_dist) {
                            is_U_forced = true;

                            /*if ((float)get_species_distance(lcamap[g], x) == max_remap_dist) {
                                cout << "max remap dist is actually useful" << endl;
                            }*/
                        }
					}

					if (U.size() > dup_cost / loss_cost)   //TODO: I'm assuming this ratio is an integer
						is_U_forced = true;


					if (is_U_forced)
						active_sets[x].erase(V);


					nb_iter++;
					if (nb_iter % 10000 == 0) {
						cout << "nb_iter = " << nb_iter << "  Ax.size = " << active_sets[x].size() << "   queue size = " << active_sets_queue.size() << endl;
					}
                }
            }

            //cout << "Cleanup species " << x->label << "  A_x size = " << active_sets[x].size() << endl;
			//Remove active sets that are already worse than something above them
			//Currently quadratic but no better solution for now
			//Also bounding - YBC
            
            /*set<set_gnodes> marked_for_deletion;
            if (!x->is_root()) {
                list<set_gnodes> active_sets_bound_queue;
                for (auto it = active_sets[x].begin(); it != active_sets[x].end(); ++it) {
                    active_sets_bound_queue.push_back(*it);
                }

                for (auto V : active_sets_bound_queue) {
                    bool remove = false;

                    for (auto Vprime : active_sets_bound_queue) {
                        //allow a slightly higher cost for Vprime corresponding to extra losses
                        if (get_cost(x, Vprime).cost < get_cost(x, V).cost + loss_cost * (V.size() - Vprime.size())
                            && is_below(V, Vprime) && V != Vprime)
                        {
                            remove = true;
                            break;
                        }

                    }

                    if (!remove && bound_option >= 1 && bound(V, x, false))   //do not allow further events at x
                        remove = true;

                    if (remove)
                        marked_for_deletion.insert(V);
                }
            }
            for (auto& V : marked_for_deletion) {
                    active_sets[x].erase(V);
            }*/



			if (!x->is_root()) {
                auto it = active_sets[x].begin();
                while (it != active_sets[x].end()) {
                    set_gnodes& V = *it;

                    
                    
                    bool remove = false;
                    if (active_sets[x].has_successor_with_lower_cost(V, costs_table[x]))
                        remove = true;

                    
                    if (!remove && bound_option >= 1 && bound(V, x, false))
                        remove = true;

                    //if (remove)
                    //    marked_for_deletion.insert(V);
                    //++it;

                    if (remove)
                        it = active_sets[x].erase(it);
                    else
                        ++it;
                    
                }

                //for (auto& V : marked_for_deletion) {
                //    active_sets[x].erase(V);
                //}
           
			}
            
            



            cout << "Done with species " << x->id << " " << (x->is_leaf() ? "(leaf " + x->label + ")" : "")
                << " A_x size = " << active_sets[x].size() << endl;

		 //test output for active sets
		   /*for (auto V : active_sets[x]) {
			  for (auto n : V)
					cout << n->label << ",";
			  cout << endl;
			  cout << "cost " << get_cost(x,V).cost << endl;
		   }*/
		}

        set_gnodes all_groots;
        for (GNode* gtree : this->genetrees) {
            all_groots.insert(gtree);
        }
        
        
		GNode* min_species = nullptr;
		int mincost = INT_MAX;
		
        for (auto it = speciestree->begin(); it != speciestree->end(); ++it) {
            SNode* x = *it;

			if (cost_exists(x, all_groots)){
				YBCost& cost = get_cost(x, all_groots);
				
				if (cost.cost < mincost){
					mincost = cost.cost;
					min_species = x;
				}
				
				std::cout<<"Best cost found."<<std::endl;
				std::cout << "Species " << x->label << " (isroot=" << x->is_root() <<")  cost=" << cost.cost
					<< "  nbdups=" << cost.nb_dups
					<< "  nblosses=" << cost.nb_losses << std::endl;

				
			}
        }
		

		SegmentalReconciliation segrec;
		segrec.dup_cost = dup_cost;
        segrec.loss_cost = loss_cost;
		
		
		if (!min_species){
			cout<<"ERROR: no species has all_groots as an active set..."<<endl;
			return segrec;
		}
		
		
		map<SNode*, set< set_gnodes > > dups_per_species;
		build_reconciliation_recursively(min_species, all_groots, dups_per_species);
		
		
		segrec.init_with_gsmap(speciestree, genetrees, lcamap);
		
		map<SNode*, set< set<GNode* > > > converted_dups_per_species;
		for (auto key : dups_per_species){
			for (auto& the_gset : key.second){
				converted_dups_per_species[key.first].insert( set<GNode*>( the_gset.begin(), the_gset.end() ) );
			}
		}
		
		segrec.build_from_segmental_dups(converted_dups_per_species);
		
		//sanity check, remove if slows down things
		YBCost& best_cost = get_cost(min_species, all_groots);
		bool err = false;
		if (best_cost.nb_losses != segrec.get_nb_losses()){
			cout<<"ERROR: DP cost differs from reconciliation cost"<<endl;
			cout<<"ybrec.nb_losses="<<best_cost.nb_losses<<"   segrec.nb_losses="<<segrec.get_nb_losses()<<endl;
			err = true;
		}
		if (best_cost.nb_dups != segrec.get_dup_height_sum()){
			cout<<"ERROR: DP cost differs from reconciliation cost"<<endl;
			cout<<"ybrec.nb_dups="<<best_cost.nb_dups<<"   segrec.nb_dups="<<segrec.get_dup_height_sum()<<endl;
			err = true;
		}
		if (!segrec.is_time_consistent()){
			err = true;
		}
		
		if (!err){
			cout<<"Reconciliation cost agrees with DP, reconciliation is time-consistent, good!"<<endl;
		}
		
		return segrec;
    }






    
	
	
	void build_reconciliation_recursively(SNode* x, set_gnodes& gnodes, map<SNode*, set< set_gnodes > >& dups_per_species){
		
        if (!cost_exists(x, gnodes)) {
            cout << "ERROR: cost at x=" << x->label << " does not exist" << endl;
        }
		YBCost& cost = get_cost(x, gnodes);
		
		//active set created from another in the same species --> all the nodes in active_set but not in originator are dups
		if (cost.nb_originators == 1){	
		
			if (cost.origin1.first != x){
				cout<<"ERROR: active set has 1 origin, but from an active set in another species"<<endl;
			}
		
			set_gnodes desc_gnodes = cost.origin1.second;
		
			set_gnodes dup_nodes;
			for (GNode* g : gnodes){
				if (!desc_gnodes.count(g))
					dup_nodes.insert(g);
			}
		
			dups_per_species[x].insert(dup_nodes);
            
			
			build_reconciliation_recursively(cost.origin1.first, cost.origin1.second, dups_per_species);
		}
		//active set created from two other in lower species
		else if (cost.nb_originators == 2){
			if (cost.origin1.first == x || cost.origin2.first == x || cost.origin1.first == cost.origin2.first){
				cout<<"ERROR: active set has 2 origins, but origins are equal to either themselves or ot x"<<endl;
			}
			build_reconciliation_recursively(cost.origin1.first, cost.origin1.second, dups_per_species);
			build_reconciliation_recursively(cost.origin2.first, cost.origin2.second, dups_per_species);
			
		}
		

	}		



};




