#ifndef NEWICKLEX_H
#define NEWICKLEX_H

#include <string>
#include "node.h"
#include "util.h"

#include <iostream>
#include <set>


class Node;

class NewickLex
{
public:

    /**
      Takes a Newick string and returns the root of a new tree.\n
      Does not try to validate anything, and assumes format correctness.\n
      The accepted Newick format is basic : it doesn't consider special parameters
      such as branch length or anything of the sort - the whole string corresponding to a node
      in the Newick string
      becomes the label of a node in the tree.
      User has to delete returned value. \n
    **/
    static Node* ParseNewickString(string& str);

    /**
      Converts a tree to a Newick string, naming the nodes using Node::GetLabel().
      If addBranchLengthToLabel is true, the branch length will be appended to the outputted label
      for each node (with a "-" between the label and the branch length)
    **/
    static string ToNewickString(Node* root, bool addBranchLengthToLabel = false, bool addInternalNodesLabel = true);


private:
    static int ReadNodeChildren(string& str, int revstartpos, Node* curNode);

    static void WriteNodeChildren(string& str, Node* curNode, bool addBranchLengthToLabel, bool addInternalNodesLabel);

    static void ParseLabel(Node* node, string label);
};



#endif // NEWICKLEX_H