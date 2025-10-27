#ifndef NODE_H
#define NODE_H

#include "define.h"

#include <vector>
#include <string>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <climits>








/**
  A tree node.
  The root has a NULL parent.
  The leaves have no children.
  Only the root should be created/destroyed by user.  The creation/destruction of descendents is handled by
  the tree (descendents are deleted when the root is deleted).
  **/
class Node
{
private:
    std::vector<Node*> children;
    Node* parent;
public:
    /*
    * Yes I am using public member variables, the greatest sin in programming.  I abide by the principle of "don't use trivial getters and setters".
    https ://github.com/isocpp/CppCoreGuidelines/blob/master/CppCoreGuidelines.md#c131-avoid-trivial-getters-and-setters
    */
    int id; 
    std::string label;
    double branch_length;

    class iterator;

    Node() {
        parent = nullptr;
        
        id = -1;
        label = "";
        branch_length = 0.0;
    }

    ~Node() {
        for (size_t i = 0; i < children.size(); i++){
            delete children[i];
        }
        children.clear();
    }


    //copy constructor
    Node(const Node& src) {
        this->id = src.id;
        this->label = src.label;
        this->branch_length = src.branch_length;

        for (int i = 0; i < src.children.size(); ++i) {
            Node* ch = new Node(*(src.children[i]));
            this->add_subtree(ch);
        }
    }


    /**
      Add a child to the current node, and returns the newly created node.
      **/
    Node* add_child() {
        Node* v = new Node();
        this->add_subtree(v);

        return v;
    }



    Node* insert_child(int index) {
        std::vector<Node*>::iterator it = children.begin();
        Node* v = new Node();
        children.insert(it + index, v);
        v->parent = this;

        return v;
    }



    int get_nb_children() {
        return children.size();
    }


    Node* get_child(int pos) {
        return children[pos];
    }


    Node* get_parent() {
        return parent;
    }



    /**
    Get sibling of the node, which is only well-defined if tree is binary.
    If parent does not have exactly two children, returns nullptr
    **/
    Node* get_sibling() {
        if (!parent || parent->get_nb_children() != 2)
            return nullptr;

        if (this == parent->get_child(0))
            return parent->get_child(1);
        else if (this == parent->get_child(1))
            return parent->get_child(0);

        throw "this->parent doesn't know his child!";
    }




    /**
    Get sibling of the node, which is only well-defined if tree is binary.
    If parent does not have exactly two children, returns nullptr
    **/
    Node* get_right_sibling() {
        if (this->is_root())
            return nullptr;
        
        int pindex = -1;
        //yeah until something better comes up, we have to iterate through all parent's children to find "this"
        for (int i = 0; i < this->parent->get_nb_children(); ++i)
        {
            if (parent->children[i] == this) {
                pindex = i;
                break;
            }
        }

        if (pindex == parent->get_nb_children() - 1)
            return nullptr;

        return parent->children[pindex + 1];
    }




    /**
      Post order iterator
      **/
    iterator begin() {
        return iterator(this, false);
    }

    iterator end() {
        return iterator(this, true);
        
    }








    /**
      Returns true iif node is the root (i.e. if parent is NULL)
      */
    bool is_root() {
        return (parent == nullptr); //yeah just return (!parent) would work
    }

    /**
      Returns true iif node is a leaf (i.e. children.size() == 0)
      */
    bool is_leaf() {
        return children.empty();
    }



    /**
      Add a tree as a new child of the node.
      The added node is the root of the subtree.  It will have its parent
      reaffacted, whether it previously had a parent or not.
      **/
    void add_subtree(Node* v) {
        children.push_back(v);
        v->parent = this;
    }

    /**
      Remove a node that belongs to the children of the node.  This child DOES NOT get deleted, and has its parent set to NULL.
      Since the caller has access to the node, he/she is expected to delete it.
      */
    void remove_child(Node* node) {
        
        std::vector<Node*>::iterator it = children.begin();

        while (it != children.end()){
            if ((*it) == node){
                (*it)->parent = nullptr;
                children.erase(it);
                return;
            }
            else{
                ++it;
            }
        }
    }

    /**
      Clear the children vector of the node.  
      */
    void remove_all_children(bool set_their_parent_to_null) {
        if (set_their_parent_to_null){
            for (size_t i = 0; i < children.size(); i++){
                children[i]->parent = nullptr;
            }
        }
        children.clear();
    }



    /**
      This is the not-so-clever non-constant-but-height-time implementation
      **/
    Node* get_lca_with(Node* v) {
        //NAIVE WAY: list all nodes from this to root
        //then list all from v to root, return first met that was visited previously
        std::unordered_set<Node*> visited;

        Node* cur = this;
        visited.insert(cur);

        while (cur){
            cur = cur->parent;
            visited.insert(cur);
        }

        cur = v;
        while (cur)
        {
            if (visited.find(cur) != visited.end())
                return cur;
            cur = cur->parent;
        }

        //this return should never happen unless the nodes are from different trees
        return nullptr;
    }

    /**
	 Returns true if this is descended from V - cannot be IN V
	 **/
    bool is_descendant_to(std::set<Node*> V) {
	    Node* cur = this;

	    while (!cur->is_root()) {
		    cur = cur->get_parent();

		    if (V.find(cur) != V.end())
			    return true;
	    }

	    return false;
    }


    /**
      Returns true iif ancestor is on the path between current node and the root (inclusively)
      **/
    bool has_ancestor(Node* ancestor) {
        Node* cur = this;
        while (cur){
            if (cur == ancestor)
                return true;
            else
                cur = cur->parent;
        }

        return false;
    }



    /**
     * @brief Returns a list of all the nodes in the tree, ordered by their visiting order in a post-order traversal.
     * The main use is to avoid confusion when the user wants to manipulate the tree while iterating over it.
     * @return
     */
    std::vector<Node*> get_postordered_nodes() {
        std::vector<Node*> ret;
        for (Node* v : *this)
            ret.push_back(v);
        return ret;
    }


    std::vector<Node*> get_preordered_nodes() {
        std::vector<Node*> vec;
        get_preordered_nodes_helper(this, vec);
        return vec;
    }

    void get_preordered_nodes_helper(Node* cur, std::vector<Node*>& to_fill) {
        to_fill.push_back(cur);
        for (Node* v : cur->children)
            get_preordered_nodes_helper(v, to_fill);
    }


    Node* get_leaf_by_label(std::string lbl) {
        
        for (auto it = this->begin(); it != this->end(); ++it) {
            if ((*it)->is_leaf() && (*it)->label == lbl)
                return *it;
        }
        
        return nullptr;
    }





	/**
	If v is an ancestor or descendant of this, returns the distance between this and v  (number of edges).
	If incomparable, returns INT_MAX
	**/
	int get_distance_to(Node* v){
		Node* cur = nullptr;
		Node* end = nullptr;
		if (this->has_ancestor(v)){
			cur = this;
			end = v;
		}
		else if (v->has_ancestor(this)){
			cur = v;
			end = this;
		}
		else{
			return INT_MAX;
		}
		
		int dist = 0;
		while (cur != end){
			cur = cur->get_parent();
			dist++;
		}
		return dist;
	}
    







    class iterator {
    
    private:
        Node* cur;
        Node* root;
        bool isend;
    public:
        iterator(Node* root, bool isend) {
            this->root = root;
            this->isend = isend;

            if (root) {
                cur = root;
                while (cur->get_nb_children() > 0)
                    cur = cur->get_child(0);
            }
        }


        iterator& operator++() {
            if (cur == root)
            {
                cur = nullptr;
                isend = true;
            }
            else
            {
                Node* r = cur->get_right_sibling();
                if (!r)
                {
                    cur = cur->get_parent();
                }
                else
                {
                    cur = r;
                    while (cur->get_nb_children() > 0)
                        cur = cur->get_child(0);
                }
            }

            return *this;
        }



        Node* operator*() {
            return cur;
        }


        bool operator==(const iterator& it) {
            if (isend)
                return it.isend;
            if (it.isend)
                return false;

            return (root == it.root && cur == it.cur);
        }

        bool operator !=(const iterator& it) {
            return !(*this == it);
        }

    };



};





#endif // NODE_H
