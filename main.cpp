
#include <iostream>

#include <set>
#include <vector>
#include <map>
#include <unordered_map>
#include <list>
#include <string>
#include <cmath>

#include "treeutils/node.h"
#include "treeutils/newicklex.h"
#include "ybrec.h"



using namespace std;


/**
That function generate_worst_case is defined at the end of the .cpp
It generates a full binary species tree of height h = height, so 2^h leaves (and an outgroup just because)
Then for each species leaf x, it generates one gene tree that has dups only in x.  It takes the species tree 
and replaces leaf x with a subtree of height "subheight" leaves, all mapped to x.  
If subheight = 0, you just get the species tree, if subheight = 1, you get one dup with two x children, etc.

The sp tree is written in ../data/sp_[height]_[subheight].newick.  Gene trees are also written in that directory.
Output location is hardcoded for now.

To use this:
./ybrec --mode "worstcase" --height 6 --subheight 1

For curiosity, there's also a "three leaves" mode, where each gene tree has three leaves, two mapped to x, plus an outgroup.
./ybrec --mode "worstcase" --height 6 --subheight 1 --threeleaves true

**/
void generate_worst_case(int height, int subheight, bool three_leaves_mode);


int main(int argc, char** argv)
{
    
    map<string, string> args;

    bool hasHelp = false;
    bool hasTest = false;
	
	string mode = "";

    //BUILD DICTIONARY OF ARGS
    string prevArg = "";
    for (int i = 0; i < argc; i++) {
        if (string(argv[i]) == "-v") {
            prevArg = "";
        }
        else if (string(argv[i]) == "--help") {
            hasHelp = true;
        }
        else {
            if (prevArg != "" && prevArg[0] == '-') {
                args[Util::ReplaceAll(prevArg, "-", "")] = string(argv[i]);
            }

            prevArg = string(argv[i]);
        }
    }


    if (args.find("help") != args.end() || hasHelp)
    {
        cout << "Here is an unhelpful help message" << endl;
    }
    else
    {
        //ML's ad-hoc testing stuff, please ignore.
        /*
        string file_prefix = "C:\\cygwin64\\home\\Manuel\\git\\SegmentalReconciler-Greedy-Version\\supplementaries\\s25";
        args["gf"] = file_prefix + "\\sim_9\\all_genetrees_edited.txt";
        args["sf"] = file_prefix + "\\sim_9\\s_tree.newick";
        args["d"] = "5";
        args["l"] = "1";
        */
        
        


        time_t start, end;
        time(&start);



        vector<GNode*> geneTrees;
        SNode* speciesTree = NULL;

        string species_separator = "_";
        int species_index = 0;
        int dupcost = 5;
        int losscost = 1;
	   int bound_option = 0;
	   
	   


        //parse dup loss cost and max dup height
        if (args.find("d") != args.end()) {
            dupcost = Util::ToInt(args["d"]);
        }
        if (args.find("l") != args.end()) {
            losscost = Util::ToInt(args["l"]);
        }

        string outfile = "";
        if (args.find("o") != args.end()) {
            outfile = args["o"];
        }


        //parse gene trees, either from command line or from file
        if (args.count("g")) {
            vector<string> gstrs = Util::Split(Util::ReplaceAll(args["g"], "\n", ""), ";", false);

            for (int i = 0; i < gstrs.size(); i++) {
                string str = gstrs[i];
                GNode* tree = NewickLex::ParseNewickString(str);

                if (!tree) {
                    cout << "Error: there is a problem with input gene tree " << str << endl;
                    return 0;
                }

                geneTrees.push_back(tree);
            }
        }
        else if (args.find("gf") != args.end()) {
            string gcontent = Util::GetFileContent(args["gf"]);

            vector<string> lines = Util::Split(gcontent, "\n");
            vector<string> gstrs;
            for (int l = 0; l < lines.size(); l++) {
                vector<string> trees_on_line = Util::Split(lines[l], ";", false);
                for (int t = 0; t < trees_on_line.size(); t++) {
				 //hack because one of my systems isn't parsing properly :( - YBC
                    if (trees_on_line[t] != "" && trees_on_line[t].length() > 1)
                        gstrs.push_back(trees_on_line[t]);
                }
            }


            for (int i = 0; i < gstrs.size(); i++) {
                string str = gstrs[i];
                GNode* tree = NewickLex::ParseNewickString(str);

                if (!tree)
                {
                    cout << "Error: there is a problem with input gene tree " << str << endl;
                    return 0;
                }

                geneTrees.push_back(tree);
            }
        }


        //parse species trees, either from command line or from file
        if (args.find("s") != args.end()) {
            speciesTree = NewickLex::ParseNewickString(args["s"]);

            if (!speciesTree) {
                cout << "Error: there is a problem with the species tree." << endl;
                return 0;
            }
        }
        else if (args.find("sf") != args.end()) {
            string scontent = Util::GetFileContent(args["sf"]);
            speciesTree = NewickLex::ParseNewickString(scontent);

            if (!speciesTree)
            {
                cout << "Error: there is a problem with the species tree." << endl;
                return 0;
            }
        }




        //parse species separator and index
        if (args.find("spsep") != args.end()) {
            species_separator = args["spsep"];
        }

        if (args.find("spindex") != args.end()) {
            species_index = Util::ToInt(args["spindex"]);
        }

	   //parse bounding option
	   if (args.find("bound") != args.end()) {
		   bound_option = Util::ToInt(args["bound"]);
	   }
	   
	   
	   
	   if (args.find("mode") != args.end()) {
		   mode = args["mode"];
	   }
	   
	   
	   //special mode checks
	   if (mode == "worstcase"){
			int height = 6;
			if (args.count("height"))
				height = Util::ToInt(args["height"]);
			
			int subheight = 1;
			if (args.count("subheight"))
				subheight = Util::ToInt(args["subheight"]);
			
			bool three_leaves_mode = false;
			if (args.count("threeleaves"))
				three_leaves_mode = true;
			
			generate_worst_case(height, subheight, three_leaves_mode);
			return 0;
		}
	   
	   


		//basic error checks
        if (geneTrees.size() == 0) {
            cout << "No gene tree given.  Program will exit." << endl;
            return 0;
        }
        if (!speciesTree) {
            cout << "No species tree given.  Program will exit." << endl;
            return 0;
        }

        if (dupcost < 0 || losscost <= 0) {
            cout << "dupcost < 0 or losscost <= 0 are prohibited.  Program will exit." << endl;
            return 0;
        }




		

        
        


        GSMap leafmap = GeneSpeciesTreeUtil::get_multi_leaf_species_map(geneTrees, speciesTree, species_separator, species_index);

		if (mode == ""){
			YBRec reconciler;
			reconciler.dup_cost = dupcost;
			reconciler.loss_cost = losscost;
			reconciler.genetrees = geneTrees;
			reconciler.speciestree = speciesTree;
			reconciler.leafmap = leafmap;



			SegmentalReconciliation segrec = reconciler.reconcile(bound_option);
			
			string str = segrec.get_output_string();
			if (outfile == "")
				cout << str;
			else{
				Util::WriteFileContent(outfile, str);
				cout<<"output written to "<<outfile<<endl;
			}
		}
		else if (mode == "lcacost"){
			SegmentalReconciliation segrec;
			segrec.init_with_gsmap(speciesTree, geneTrees, leafmap);
			segrec.apply_lca_mapping();
			int nb_seg_dups = segrec.get_dup_height_sum();
			int nb_losses = segrec.get_nb_losses();
			cout<<"Nb segmental dups = "<<nb_seg_dups<<endl;
			cout<<"Nb losses = "<<nb_losses<<endl;
			cout<<"Cost = "<< (nb_seg_dups * dupcost + nb_losses * losscost)<<endl;
		}
		else{
			cout<<"mode "<<mode<<" does not exists"<<endl;
		}


        for (GNode* g : geneTrees)
            delete g;
        geneTrees.clear();
        delete speciesTree;

        
    }
}



//WATCH OUT: this function creates a "new" and expects the caller to make the delete.  Bad design,
//so to make sure the function name has "to_delete" in it.
Node* get_full_binary_tree_to_delete(int height, string default_leaf_label = ""){
	Node* v = new Node();
	
	if (height > 0){
		Node* t1 = get_full_binary_tree_to_delete(height - 1, default_leaf_label);
		Node* t2 = get_full_binary_tree_to_delete(height - 1, default_leaf_label);
		v->add_subtree(t1);
		v->add_subtree(t2);
	}
	else{
		if (default_leaf_label != "")
			v->label = default_leaf_label;
	}
	
	
	return v;
}


void generate_worst_case(int height, int subheight, bool three_leaves_mode){
	SNode* stree_left = get_full_binary_tree_to_delete(height);
	
	SNode* outgroup = new Node();
	string outgroup_label = Util::ToString(pow(2, height) + 1);
	outgroup->label = outgroup_label;
	
	SNode* stree = new Node();
	stree->add_subtree(stree_left);
	stree->add_subtree(outgroup);
	
	
	
	int cpt = 1;
	int cpt_internalnode = 1;
	for (auto it = stree->begin(); it != stree->end(); ++it){
		SNode* x = *it;
		if (x->is_leaf()){
			x->label = Util::ToString(cpt);
			cpt++;
		}
		else{
			x->label = "SI_" + Util::ToString(cpt_internalnode);
			cpt_internalnode++;
		}
	}
	
	//gene trees 
	//for each sp leaf x, create a gene tree of the form ((x, x), outgroup);
	vector<GNode*> gtrees;
	cpt_internalnode = 1;
	for (auto it = stree->begin(); it != stree->end(); ++it){
		SNode* x = *it;
		if (x->is_leaf()){
			
			
			if (three_leaves_mode){
				string gs = "((" + x->label + ", " + x->label + ")," + outgroup_label + ");";
				cpt_internalnode += 2;
				GNode* g = NewickLex::ParseNewickString(gs);
				gtrees.push_back(g);
			}
			else{	//copy sp tree, add leaves
			
				cpt = 1;
				GNode* gnode = get_full_binary_tree_to_delete(height);
				for (auto itg = gnode->begin(); itg != gnode->end(); ++itg){
					GNode* g = *itg;
					if (g->is_leaf()){
						
						
						if (subheight > 0 && Util::ToString(cpt) == x->label){
							GNode* duptree1 = get_full_binary_tree_to_delete(subheight - 1, x->label);
							GNode* duptree2 = get_full_binary_tree_to_delete(subheight - 1, x->label);
							g->add_subtree(duptree1);
							g->add_subtree(duptree2);
						}
						else{
							g->label = Util::ToString(cpt);
						}
						
						
						cpt++;
					}
					
				}
				gtrees.push_back(gnode);
			}
			
		}
	}
	
	
	
	
	string sstr = NewickLex::ToNewickString(stree);
	Util::WriteFileContent("../data/sp" + Util::ToString(height) + "_" + Util::ToString(subheight) + 
					(three_leaves_mode ? "_threeleaves" : "") + ".newick", sstr);
	delete stree;
	
	cpt_internalnode = 1;
	string gstr = "";
	for (GNode* groot : gtrees){	//i am groot
		for (auto it = groot->begin(); it != groot->end(); ++it){
			GNode* g = *it;
			if (!g->is_leaf()){
				g->label = "GI_" + Util::ToString(cpt_internalnode);
				cpt_internalnode++;
			}
		}
		gstr += NewickLex::ToNewickString(groot) + "\n";
		delete groot;
	}
	Util::WriteFileContent("../data/gt" + Util::ToString(height) + "_" + Util::ToString(subheight) + 
					(three_leaves_mode ? "_threeleaves" : "") + ".newick", gstr);
	
	gtrees.clear();
}