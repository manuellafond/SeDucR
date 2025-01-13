
#include <iostream>

#include <set>
#include <vector>
#include <map>
#include <unordered_map>
#include <list>
#include <string>

#include "treeutils/node.h"
#include "treeutils/newicklex.h"
#include "ybrec.h"



using namespace std;




int main(int argc, char** argv)
{
    
    map<string, string> args;

    bool hasHelp = false;
    bool hasTest = false;

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
                    if (trees_on_line[t] != "")
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

        YBRec reconciler;
        reconciler.dup_cost = dupcost;
        reconciler.loss_cost = losscost;
        reconciler.genetrees = geneTrees;
        reconciler.speciestree = speciesTree;
        reconciler.leafmap = leafmap;




        reconciler.reconcile();

        for (GNode* g : geneTrees)
            delete g;
        geneTrees.clear();
        delete speciesTree;

        
    }
}

