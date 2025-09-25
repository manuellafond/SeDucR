#pragma once

#include "treeutils/treeutil.h"
#include <set>
#include <vector>
#include <map>
#include <unordered_map>
#include <list>
#include <string>

#include "treeutils/node.h"
#include "treeutils/newicklex.h"


//NOTE: GNode, SNode and GSMap are defined in genespeciestreeutil
#include "genespeciestreeutil.h"

//library for compressed set representations
#include "treeutils/ewah/ewah.h"


typedef ewah::EWAHBoolArray<uint32_t> bitmap;







//a comparator function for bitmaps, so that we can store them in sets and maps
struct BCmp {
    bool operator()(const bitmap& a, const bitmap& b) const {
        auto ita = a.begin();
        auto itb = b.begin();


        while (true) {  //never do a while (true)
            if (ita != a.end()) {
                if (itb != a.end()) {
                    if (*ita < *itb)
                        return true;
                    else if (*itb < *ita)
                        return false;
                    else {
                        ++ita;
                        ++itb;
                    }
                }
                else
                    return true;
            }
            else {
                if (itb != b.end())
                    return true;
                else
                    return false;    //both end
            }
        }
    }
};







/**
YBRec class
Stores all info necessary to apply the YB dynamic programming algorithm, given
a set of gene trees, a species tree, and d/l costs

todo: for now everything is just public, cleanup later
**/
class YBRec {
public:

    /*
    Structure to hold cost information in the c table.  So, c(x, V) stores a YBCost 
    object, which remembers the number of losses and segmental dups (and total cost).
    We can also give it debug messages, useful for, well, debugging.
    */
    struct YBCost {
        int nb_losses;
        int nb_dups;
        int cost;

        vector<string> debug;
    };


    Node* speciestree;
    vector<Node*> genetrees;

    int loss_cost;
    int dup_cost;

    YBCost cost_upper_bound;

    //max_remap_dist is the distance that a gene tree node can be remmapped to wrt lca.  Default = 5
    //eg if max_remap_dist = 3, then a gene node can only be moved "up" by 3 from its lca mapping.
    int max_remap_dist;

    //the A_x sets from the YB algorithm.  
    //key = species
    //value = the active sets for that species (each set is represented as a bitmap in an attempt to save memory)
    map< SNode*, set<bitmap, BCmp> > active_sets;

    GSMap leafmap;
    GSMap lcamap;

    map< SNode*, set<GNode*> > rev_leafmap;     //rev_leafmap[x] = set of gene tree leaves mapped to s
    map<int, GNode*> id_genes;                  //id_genes[i] = the gene tree node with id = i
									   
    map< SNode*, set<GNode*> > rev_lcamap;     //same as rev_leafmap, but for LCA
    map< SNode*, int > lca_dupheight;


    //the c(A, x) values from the YB algorithm.  Indexed by species, then by the desired set
    //costs_table[x][my_set] = the cost.  Previous version stored only cost, now we have whole struct
    //map< SNode*, map<bitmap, int, BCmp > > costs_table;    //could be slow to hash keys
    map< SNode*, map<bitmap, YBCost, BCmp > > costs_table;    //could be slow to hash keys


    bool DEBUG;

    YBRec() {
        DEBUG = false;
        dup_cost = 5;
        loss_cost = 1;
        max_remap_dist = 100; //changed by YBC
    }


    //sets the cost c(gnodes, x).  Only updated if entry does not exist of cost is smaller than current cost
    /*void add_cost(SNode* x, const bitmap& gnodes, int cost) {

        int best_cost = INT_MAX;
        if (costs_table.count(x) && costs_table[x].count(gnodes)) {
            best_cost = costs_table[x][gnodes];
        }
        costs_table[x][gnodes] = min(best_cost, cost);
    }*/


    //sets c(x, V).  Updates only if cost is smaller than previous entry, or if such an entry does not exist.
    // We must specify the number of dups and losses.  For debugging, we can pass 
    //debug messages to combine.  If not debugging, just ignore.
    void set_cost(SNode* x, const bitmap& gnodes, int nb_dups, int nb_losses, 
                        vector<string>* dbg1 = nullptr, vector<string>* dbg2 = nullptr, string dbg = "") {

        bool update_cost = false;
        int newcost = this->dup_cost * nb_dups + this->loss_cost * nb_losses;
        
        if (costs_table.count(x) && costs_table[x].count(gnodes)) {
            YBCost& best_cost = costs_table[x][gnodes];
            if (newcost < best_cost.cost) {
                update_cost = true;
            }
        }
        else
            update_cost = true;

        if (update_cost) {
            YBCost new_ybcost;
            new_ybcost.nb_dups = nb_dups;
            new_ybcost.nb_losses = nb_losses;
            new_ybcost.cost = newcost;
            
            if (dbg1)
                new_ybcost.debug.insert(new_ybcost.debug.end(), dbg1->begin(), dbg1->end());
            if (dbg2)
                new_ybcost.debug.insert(new_ybcost.debug.end(), dbg2->begin(), dbg2->end());
            if (dbg != "")
                new_ybcost.debug.push_back(dbg);

            costs_table[x][gnodes] = new_ybcost;

            
        }
        
    }


    //returns c(gnodes, x)
    /*int get_cost(SNode* x, bitmap& gnodes) {

        if (costs_table.count(x) && costs_table[x].count(gnodes)) {
            return costs_table[x][gnodes];
        }

        return INT_MAX;
    }*/


    //returns c(gnodes, x)
    YBCost get_cost(SNode* x, bitmap& gnodes) {

        if (costs_table.count(x) && costs_table[x].count(gnodes)) {
            return costs_table[x][gnodes];
        }

        YBCost dummy_cost;
        dummy_cost.cost = INT_MAX;
        return dummy_cost;
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
    set<GNode*> get_cross_set(const set<GNode*>& V_l, const set<GNode*>& V_r) {
        set<GNode*> ret;

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
    set<GNode*> get_common_parents(const set<GNode*>& V) {
        set<GNode*> ret;

        for (GNode* g1 : V) {
            GNode* g2 = g1->get_sibling();

            if (g2 && V.count(g2)) {
                ret.insert(g1->get_parent());
            }
        }

        return ret;
    }




    //get nuber of edges from lower_sp to higher_sp
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



    //converts a set of GNodes to a bitmap, based on the id of GNodes.
    bitmap nodeset_to_bitmap(const set<GNode*>& gset) {
        bitmap b;

        set<int> ids;
        //NOTE: bitmaps need set to be called in increasing order.  This works ONLY because sets iterate in sorted order.
        for (auto it = gset.begin(); it != gset.end(); ++it) {
            ids.insert((*it)->id);
        }
        for (auto id : ids) {
            b.set(id);
        }

        return b;
    }



    //converts bitmap to set of gnodes, based on id
    set<GNode*> bitmap_to_nodeset(bitmap& b) {
        set<GNode*> gset;

        for (int i : b) {
            gset.insert(id_genes[i]);
        }

        return gset;
    }




    //preprocessing to assign a unique id to each gnode and snode.  This fills id_genes map.
    void init_indices() {
        auto snodes = speciestree->get_postordered_nodes();
        for (int i = 0; i < snodes.size(); ++i) {
            snodes[i]->id = i;
        }


        int cpt = 0;
        for (GNode* groot : genetrees) {
            auto gnodes = groot->get_postordered_nodes();
            for (GNode* g : gnodes) {
                g->id = cpt;
                id_genes[cpt] = g;
                cpt++;
            }
        }
    }

    //Return true if (x,V) can be removed
    bool bound(bitmap V, SNode* x, bool allow_events_at_x) {
      // return false;

      vector<SNode*> species_nodes = speciestree->get_postordered_nodes();
      set<GNode*> V_set = bitmap_to_nodeset(V);
      YBCost lower_bound = get_cost(x, V);

      //Calculate lower bound for x, V
      //Bound duplications
      int dups = 0;
      for (GNode* g_root : genetrees)
      for (auto it = g_root->begin(); it != g_root->end(); ++it) {
        Node *g = *it;

        //g is in V or a leaf not descended from V
        if (V_set.find(g) != V_set.end() || (g->is_leaf() && !g->is_descendant_to(V_set)) ) {
          Node* cur = g;
          int dupheight = 0;

          while (!cur->is_root()) {
            if (lcamap[cur->get_parent()] == lcamap[cur] || lcamap[cur->get_parent()] == lcamap[cur->get_sibling()])
              dupheight++;

            cur = cur->get_parent();
          }

          dups = max(dups, dupheight);
        }
      }
      lower_bound.nb_dups += dups;
      lower_bound.cost += dups * dup_cost;

      //Bound losses
      //Set feasible reconciliation
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

      //calculate losses, update lower bound
      int losses = 0;
      for (GNode* g_root : genetrees)
      for (auto it = g_root->begin(); it != g_root->end(); ++it) {
        Node *g = *it;

        if (!g->is_root() && !g->is_descendant_to(V_set)) {

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

        if (!g->is_descendant_to(V_set)) {
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

    void reconcile() {


        init_indices();
        compute_rev_leafmap();
        compute_lca_map();
	   compute_lca_upper_bound();

	   cout << "LCA upper bound = " << cost_upper_bound.cost << endl;

        vector<SNode*> species_nodes = speciestree->get_postordered_nodes();
        for (SNode* x : species_nodes) {

            if (x->is_leaf()) {

                bitmap bnodes = nodeset_to_bitmap(rev_leafmap[x]);

                set_cost(x, bnodes, 0, 0);
                //add_cost(x, bnodes, 0);

                active_sets[x].insert(bnodes);
            }
            else {
                SNode* xl = x->get_child(0);
                SNode* xr = x->get_child(1);


                //build all the possible active sets from those of x's children
                for (bitmap V_l : active_sets[xl]) {
                    for (bitmap V_r : active_sets[xr]) {
                        set<GNode*> v_cross_set = get_cross_set(bitmap_to_nodeset(V_l), bitmap_to_nodeset(V_r));
                        
                        bitmap v_cross = nodeset_to_bitmap(v_cross_set);
                        
                        int nb_roots_vl_vr = 0;
                        for (int i : V_l) {
                            if (this->id_genes[i]->is_root())
                                nb_roots_vl_vr++;
                        }
                        for (int i : V_r) {
                            if (this->id_genes[i]->is_root())
                                nb_roots_vl_vr++;
                        }

                        //so, v_cross is the set we get if we pair xl and xr guys into a speciation
                        //when they have a common parent, and "raise" all the other ones.
                        //Each such raise causes one loss that must be counted, hence the second line of the cost
                        YBCost leftcost = get_cost(xl, V_l);
                        YBCost rightcost = get_cost(xr, V_r);
                        
                        int extra_losses = (2 * v_cross.numberOfOnes() - V_l.numberOfOnes() - V_r.numberOfOnes());
                        extra_losses -= nb_roots_vl_vr;

                        set_cost(x, v_cross, leftcost.nb_dups + rightcost.nb_dups,
                            leftcost.nb_losses + rightcost.nb_losses + extra_losses, 
                            
                            //passing debugging stuff if activated, horrible stuff
                            (DEBUG ? &leftcost.debug : nullptr), 
                            (DEBUG ? &rightcost.debug : nullptr),
                            (DEBUG ? Util::ToString(extra_losses) + " losses in " + Util::ToString(x->id) + 
                                "(" + Util::ToString(leftcost.nb_losses) + " + " + Util::ToString(rightcost.nb_losses) + ")"
                                : "")
                        );
                        
                        //int cost = get_cost(xl, V_l) + get_cost(xr, V_r) +
                        //    this->loss_cost * (2 * v_cross.numberOfOnes() - V_l.numberOfOnes() - V_r.numberOfOnes());
                        //add_cost(x, v_cross, cost);
            
                        //apply bounding here - YBC
                        //if (!bound(v_cross, x, true))
                          active_sets[x].insert(v_cross);

                        if (active_sets[x].size() % 10000 == 0 && active_sets[x].size() > 0) {
                            cout << "Sp x=" << x->id << "  Ax size is now " << active_sets[x].size();
                            cout << " Axl=" << active_sets[xl].size() << "  Axr=" << active_sets[xr].size() << endl;

                        }
                    }
                }
            }



            list< bitmap > active_sets_queue(active_sets[x].begin(), active_sets[x].end());

            //set< set<GNode*> > active_sets_to_insert;
            //set< set<GNode*> > active_sets_to_delete;
            int nb_iter = 0;

            while (!active_sets_queue.empty()) {
                bitmap V = active_sets_queue.front();
                active_sets_queue.pop_front();

                set<GNode*> U = get_common_parents(bitmap_to_nodeset(V));

                if (!U.empty()) {  //if there are actually dups that can be applied

                    set<GNode*> Vprime_set = get_cross_set(bitmap_to_nodeset(V), bitmap_to_nodeset(V));   //not sure that works, I think it does - so it probably doesn't work
                    bitmap Vprime = nodeset_to_bitmap(Vprime_set);

                    //apply bounding here - YBC
                    YBCost xV_cost = get_cost(x, V);

				//remove previous set if new set is already better
				if (get_cost(x, Vprime).cost < xV_cost.cost + this->dup_cost)
					active_sets[x].erase(V);

                    set_cost(x, Vprime, xV_cost.nb_dups + 1, xV_cost.nb_losses);
                    //if (!bound(Vprime, x, true)) {
                      active_sets[x].insert(Vprime);
                      active_sets_queue.push_back(Vprime);
                    //}

                    bool is_U_forced = false;
                    //if someone is trying to go too far, we have to do the dup here
                    for (GNode* g : U) {
                        if ((float)get_species_distance(lcamap[g], x) == max_remap_dist)
                            is_U_forced = true;
                    }

                    if (U.size() >= dup_cost / loss_cost)   //TODO: I'm assuming this ratio is an integer
                        is_U_forced = true;


                    if (is_U_forced)
                        active_sets[x].erase(V);


                    nb_iter++;
                    if (nb_iter % 5000 == 0) {
                        cout << "nb_iter = " << nb_iter << "  Ax.size = " << active_sets[x].size() << "   queue size = " << active_sets_queue.size() << endl;
                    }
                }
            }

		  
		  //Bounding - YBC
		  /*if (!x->is_root()) {
			  list< bitmap > active_sets_bound_queue(active_sets[x].begin(), active_sets[x].end());

			  for (auto V : active_sets_bound_queue) {
				if (bound(V, x, false))   //do not allow further events at x
				  active_sets[x].erase(V);
			  }
		  }*/



            cout << "Done with species " << x->id << " " << (x->is_leaf() ? "(leaf " + x->label + ")" : "")
                << " A_x size = " << active_sets[x].size() << endl;

		  /*for (auto s : active_sets[x]) {
				set<GNode*> V = bitmap_to_nodeset(s);
				for (auto n : V)
					cout << n->label << ",";
				cout << endl;
				cout << get_cost(x,s).cost << endl;
		  }*/
        }








        std::cout << std::endl << "*** Outputting all c(x, V) entries for V = {all roots of G} ***" << std::endl;

        std::set<GNode*> all_groots_set;
        for (GNode* gtree : this->genetrees) {
            all_groots_set.insert(gtree);
        }
        bitmap all_groots = nodeset_to_bitmap(all_groots_set);
        
        for (auto it = speciestree->begin(); it != speciestree->end(); ++it) {
            SNode* x = *it;

            YBCost cost = get_cost(x, all_groots);
            
            
            if (cost.cost != INT_MAX) {
                std::cout << "Species " << x->label << "  cost=" << cost.cost
                    << "  nbdups=" << cost.nb_dups
                    << "  nblosses=" << cost.nb_losses << std::endl;

                if (DEBUG) {
                    for (string s : cost.debug) {
                        std::cout << s << std::endl;
                    }
                }
            }
            
        }





    }





    



};






//This is an old attempt in which active sets were stored as "set" objects.  This also used an old tree library.
//That took too much memory, but it's there to preserve the history.
/*class YBRec {
public:

    Node* speciestree;
    vector<Node*> genetrees;

    int loss_cost;
    int dup_cost;

    map< SNode*, set<set<GNode*>> > active_sets;

    GSMap leafmap;
    GSMap lcamap;
    map< SNode*, set<GNode*> > rev_leafmap;


    map< pair<SNode*, set<GNode*>>, int > costs_table;    //could be slow to hash keys



    void add_cost(SNode* x, set<GNode*>& gnodes, int cost) {
        auto key = make_pair(x, gnodes);

        int best_cost = 9999999;
        if (costs_table.count(key)) {
            best_cost = costs_table[key];
        }
        costs_table[key] = min(best_cost, cost);
    }


    int get_cost(SNode* x, set<GNode*>& gnodes) {
        auto key = make_pair(x, gnodes);

        if (costs_table.count(key)) {
            return costs_table[key];
        }

        return 999999;
    }


    void compute_rev_leafmap() {

        for (GNode* g_root : genetrees) {

            TreeIterator* it = g_root->GetPostOrderIterator();
            while (Node* g = it->next()) {
                if (g->IsLeaf()) {
                    rev_leafmap[leafmap[g]].insert(g);
                }
            }
            g_root->CloseIterator(it);

        }
    }


    void compute_lca_map() {
        for (GNode* g_root : genetrees) {

            TreeIterator* it = g_root->GetPostOrderIterator();
            while (Node* g = it->next()) {
                if (g->IsLeaf()) {
                    lcamap[g] = leafmap[g];
                }
                else {
                    lcamap[g] = lcamap[g->GetChild(0)]->FindLCAWith(lcamap[g->GetChild(1)]);
                }
            }
            g_root->CloseIterator(it);

        }
    }


    set<GNode*> get_cross_set(set<GNode*>& V_l, set<GNode*>& V_r) {
        set<GNode*> ret;

        //pass 1: for every gl in V_l, check if gl's sibling is in V_l.  If so, they get "fused"
        //and the common parent is inserted.  If not, gl is kept.
        for (GNode* gl : V_l) {
            GNode* gr = gl->GetSibling();

            if (gr && V_r.count(gr)) {
                ret.insert(gl->GetParent());
            }
            else {
                ret.insert(gl);
            }
        }

        //pass 2: add every gr in V_r whose parent was not added.
        for (GNode* gr : V_r) {
            if (!gr->GetParent() || !ret.count(gr->GetParent())) {
                ret.insert(gr);
            }
        }

        return ret;
    }




    //returns the set of nodes whose children are both in V

    set<GNode*> get_common_parents(set<GNode*>& V) {
        set<GNode*> ret;

        for (GNode* g1 : V) {
            GNode* g2 = g1->GetSibling();

            if (g2 && V.count(g2)) {
                ret.insert(g1->GetParent());
            }
        }

        return ret;
    }





    //TODO: make this function faster - also, lower_sp MUST be a descendant of higher_sp
    int get_species_distance(SNode* lower_sp, SNode* higher_sp) {
        int d = 0;

        while (lower_sp != higher_sp) {
            lower_sp = lower_sp->GetParent();
            d++;
        }
        return d;
    }


    void reconcile() {


        int max_remap_dist = 5;

        compute_rev_leafmap();
        compute_lca_map();

        TreeIterator* species_it = speciestree->GetPostOrderIterator();
        while (SNode* x = species_it->next()) {

            if (x->IsLeaf()) {

                add_cost(x, rev_leafmap[x], 0);

                active_sets[x].insert(rev_leafmap[x]);
            }
            else {
                SNode* xl = x->GetChild(0);
                SNode* xr = x->GetChild(1);


                //build all the possible active sets from those of x's children
                for (set<GNode*> V_l : active_sets[xl]) {
                    for (set<GNode*> V_r : active_sets[xr]) {
                        set<GNode*> v_cross = get_cross_set(V_l, V_r);

                        //TODO: if non-merged V_l and V_r is too large, we should abort

                        //so, v_cross is the set we get if we pair xl and xr guys into a speciation
                        //when they have a common parent, and raise all the other ones.
                        //Each such raise causes one loss that must be counted, hence the second line of the cost
                        int cost = get_cost(xl, V_l) + get_cost(xr, V_r) +
                            this->loss_cost * (2 * v_cross.size() - V_l.size() - V_r.size());
                        add_cost(x, v_cross, cost);

                        active_sets[x].insert(v_cross);

                        if (active_sets[x].size() % 10000 == 0) {
                            cout << "Sp x=" << x->GetIndex() << "  Ax size is now " << active_sets[x].size();
                            cout << " Axl=" << active_sets[xl].size() << "  Axr=" << active_sets[xr].size() << endl;

                        }
                    }
                }
            }



            list< set<GNode*> > active_sets_queue(active_sets[x].begin(), active_sets[x].end());

            //set< set<GNode*> > active_sets_to_insert;
            //set< set<GNode*> > active_sets_to_delete;
            int nb_iter = 0;

            while (!active_sets_queue.empty()){
                set<GNode*> V = active_sets_queue.front();
                active_sets_queue.pop_front();

                set<GNode*> U = get_common_parents(V);

                if (!U.empty()) {  //if there are actually dups that can be applied

                    set<GNode*> Vprime = get_cross_set(V, V);   //not sure that works, I think it does - so it probably doesn't work


                    bool isVprimeBad = false;
                    for (GNode* g : Vprime) {
                        if ((float)get_species_distance(lcamap[g], x) > max_remap_dist) {
                            isVprimeBad = true;
                            break;
                        }
                    }

                    if (!isVprimeBad){
                        active_sets[x].insert(Vprime);
                        active_sets_queue.push_back(Vprime);

                        add_cost(x, Vprime, get_cost(x, V) + dup_cost);

                        bool is_U_forced = false;
                        //if someone is trying to go too far, we have to do the dup here
                        for (GNode* g : U) {
                            if ((float)get_species_distance(lcamap[g], x) == max_remap_dist)
                                is_U_forced = true;
                        }

                        if (U.size() >= dup_cost / loss_cost)
                            is_U_forced = true;


                        if (is_U_forced)
                            //active_sets_to_delete.insert(V);
                            active_sets[x].erase(V);
                    }

                    nb_iter++;
                    if (nb_iter % 5000 == 0) {
                        cout << "nb_iter = " << nb_iter << "  Ax.size = " << active_sets[x].size() << "   queue size = " << active_sets_queue.size() << endl;
                    }
                }

            }



            cout << "Done with species " << x->GetIndex() << "  A_x size = " << active_sets[x].size() << endl;
        }
        speciestree->CloseIterator(species_it);
    }



};
*/



