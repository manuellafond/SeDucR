
#include "newicklex.h"




Node* NewickLex::ParseNewickString(string& str)
{
    Node* root = new Node();
    int pos = str.find_last_of(')');
    int lastcolonpos = str.find_last_of(';');


    if (pos != string::npos)
    {
        string label = "";
        if (lastcolonpos > pos)
        {
            label = str.substr(pos + 1, lastcolonpos - pos - 1); 
        }
        else
        {
            label = Util::RTrim(str.substr(pos + 1));
        }
        root->label = label;

        ReadNodeChildren(str, pos - 1, root);
    }
    else
    {
        //this should be the root label
        root->label = str.substr(0, str.length() - 1);
    }

    return root;
}


string NewickLex::ToNewickString(Node* root, bool addBranchLengthToLabel, bool addInternalNodesLabel)
{
    string str;
    WriteNodeChildren(str, root, addBranchLengthToLabel, addInternalNodesLabel);
    str += ";";
    return str;
}



int NewickLex::ReadNodeChildren(string& str, int revstartpos, Node* curNode)
{
    bool ok = true;
    int pos = revstartpos;

    //cout<<"POS="<<pos<<endl;

    while (ok && pos >= 0)
    {
        int closepos = str.find_last_of(')', pos);
        int openpos = str.find_last_of('(', pos);
        int commapos = str.find_last_of(',', pos);

        int maxpos = max(max(closepos, openpos), commapos);

        if (maxpos == string::npos)
        {
            ok = false;
            pos = -1;
        }
        else
        {

            string lbl = "";

            lbl = Util::Trim(str.substr(maxpos + 1, pos - maxpos));

            if (maxpos == closepos)
            {
                //cout<<"Close="<<maxpos<<"      LBL="<<lbl<<endl;
                Node* newNode = curNode->insert_child(0);

                ParseLabel(newNode, lbl);
                //newNode->SetLabel(lbl);
                pos = ReadNodeChildren(str, maxpos - 1, newNode);

                while (str[pos] != ',' && str[pos] != '(')
                    pos--;

                if (str[pos] == ',')
                    pos--;
            }
            else if (maxpos == commapos)
            {
                //cout<<"Comma="<<maxpos<<"      LBL="<<lbl<<endl;

                int ptcomma = str.find_first_of(',', maxpos + 1);
                int ptopen = str.find_first_of('(', maxpos + 1);
                int ptclose = str.find_first_of(')', maxpos + 1);

                if ((ptcomma < ptopen && ptcomma < ptclose) ||
                    (ptopen == string::npos || ptopen > ptclose))
                {
                    Node* newNode = curNode->insert_child(0);
                    ParseLabel(newNode, lbl);
                    //newNode->SetLabel(lbl);
                }
                pos = maxpos - 1;
            }
            else if (maxpos == openpos)
            {
                //cout<<"Open="<<maxpos<<"      LBL="<<lbl<<endl;

                //EDIT ML AUG 2012 : an opening parenthese creates a node only if followed by a ,
                // THIS WAS DONE DURING A PHASE OF TIREDNESS
                // IF SOMETHING IS BUGGY, IT'S PROBABLY AROUND HERE
                //NOTE ML AUG 2013 : it does seem to be holding up, even after extensive use
                int ptcomma = str.find_first_of(',', maxpos + 1);
                int ptopen = str.find_first_of('(', maxpos + 1);
                int ptclose = str.find_first_of(')', maxpos + 1);

                if (ptcomma != string::npos &&
                    (ptcomma < ptclose || ptclose == string::npos) &&
                    (ptcomma < ptopen || ptopen == string::npos))
                {
                    Node* newNode = curNode->insert_child(0);
                    ParseLabel(newNode, lbl);
                    //newNode->SetLabel(lbl);
                }
                pos = maxpos - 1;
                ok = false;
            }
        }

    }

    return pos;
}

void NewickLex::WriteNodeChildren(string& str, Node* curNode, bool addBranchLengthToLabel, bool addInternalNodesLabel)
{
    if (curNode->is_leaf())
    {

        str += curNode->label;
        //str += "_I";
        //str += Util::ToString(curNode->GetIndex());

        if (addBranchLengthToLabel && !curNode->is_root())
            str += ":" + Util::ToString(curNode->branch_length);
    }
    else
    {
        str += "(";
        for (int i = 0; i < curNode->get_nb_children(); i++)
        {
            if (i != 0)
                str += ", ";

            Node* child = curNode->get_child(i);
            WriteNodeChildren(str, child, addBranchLengthToLabel, addInternalNodesLabel);
        }

        str += ")";

        if (addInternalNodesLabel)
            str += curNode->label;
            //str += "_I";
            //str += Util::ToString(curNode->GetIndex());

        if (addBranchLengthToLabel && !curNode->is_root() && !curNode->branch_length == 0.0)
            str += ":" + Util::ToString(curNode->branch_length);
    }


}




void NewickLex::ParseLabel(Node* node, string label)
{
    string sublabel = "";
    int pos = label.find("[");
    if (pos != string::npos)
    {
        sublabel = label.substr(pos);
        label = label.substr(0, pos - 1);
    }

    vector<string> sz = Util::Split(label, ":");

    if (sz.size() == 1)
    {
        //If label is just a double, we consider it's a branch length
        /*if (Util::IsDouble(sz[0]))
        {
            node->SetBranchLength(Util::ToDouble(sz[0]));
        }
        else*/
        {
            node->label = label;
        }
    }
    else
    {
        node->label = sz[0];
        node->branch_length = (Util::ToDouble(sz[1]));
    }

    node->label = (node->label + sublabel);
}


