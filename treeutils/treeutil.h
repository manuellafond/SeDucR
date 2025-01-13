#pragma once

#include "node.h"
#include "util.h"

class TreeUtil {
public:

    static void get_random_binary_tree_rec(Node* v, set<int>& indices) {

        if (indices.size() == 1) {
            v->label = (Util::ToString(*indices.begin()));
            return;
        }

        set<int> left;
        set<int> right;

        bool done = false;

        while (!done) {
            for (auto it = indices.begin(); it != indices.end(); ++it) {
                int x = rand() % 2;
                if (x == 0)
                    left.insert(*it);
                else
                    right.insert(*it);
            }

            //dumb way to ensure no empty child
            if (left.empty() || right.empty()) {
                left.clear();
                right.clear();
            }
            else {
                done = true;
            }
        }


        Node* v1 = new Node();
        get_random_binary_tree_rec(v1, left);
        Node* v2 = new Node();
        get_random_binary_tree_rec(v2, right);

        v->add_subtree(v1);
        v->add_subtree(v2);
    }





    static void get_random_binary_tree(Node* root, int nb_leaves) {
        set<int> labels;
        for (int i = 1; i <= nb_leaves; ++i)
            labels.insert(i);

        return get_random_binary_tree_rec(root, labels);
    }



    static void contract_parent_edge(Node* v) {
        if (v->is_root())
            return;

        Node* p = v->get_parent();
        p->remove_child(v);
        for (int i = 0; i < v->get_nb_children(); ++i) {
            p->add_subtree(v->get_child(i));
        }
        v->remove_all_children(false);
        delete v;
    }




    static void randomize_branch_lengths(Node* v, double min, double max) {
        if (!v->is_root()) {
            double r = (double)rand() / RAND_MAX;
            double b = min + r * (max - min);
            v->branch_length = b;
        }

        for (int i = 0; i < v->get_nb_children(); ++i) {
            randomize_branch_lengths(v->get_child(i), min, max);
        }

    }

    /**
    Creates a degree 2 node between v and its parent, and returns the new node.  If v is the root, does nothing and returns nullptr.
    **/
    static Node* subdivide_parent_edge(Node* v) {
        if (v->is_root())
            return nullptr;

        Node* old_parent = v->get_parent();
        v->get_parent()->remove_child(v);

        Node* w = new Node();
        w->add_subtree(v);
        old_parent->add_subtree(w);

        return w;
    }




    static void reroot_on_node(Node* v) {
        vector<Node*> ancestors;

        Node* cur = v;
        while (cur) {
            ancestors.push_back(cur);
            cur = cur->get_parent();
        }


        for (int i = ancestors.size() - 1; i >= 1; --i){
            Node* w = ancestors[i];

            w->remove_child(ancestors[i - 1]);
            ancestors[i - 1]->add_subtree(w);
        }


    }






};