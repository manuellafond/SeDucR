#pragma once 

#include <iostream>
#include <unordered_map>
#include <set>
#include <vector>

#include "genespeciestreeutil.h"


struct GNodeCmp {
    bool operator()(const GNode* a, const GNode* b) const {
        return a->id < b->id;
    }
};

typedef std::set<GNode*, GNodeCmp> set_gnodes;




/*
Structure to hold cost information in the c table.  So, c(x, V) stores a YBCost 
object, which remembers the number of losses and segmental dups (and total cost).
We can also give it debug messages, useful for, well, debugging.
*/
struct YBCost {
	int nb_losses;
	int nb_dups;
	int cost;
	
	//for backtracking, remember originator(s) of current cost (there can be 1 or 2, or 0 for leaves)
	size_t nb_originators = 0;
	pair<SNode*, set_gnodes> origin1;
	pair<SNode*, set_gnodes> origin2;
};


class Trie {
private:
    struct TNode {
        bool isEnd = false;
        std::map<GNode*, TNode*, GNodeCmp> children;
		YBCost cost;
    };

    TNode* root;

	size_t nbelem = 0;

public:
    Trie() {
        root = new TNode();
		nbelem = 0;
    }

    ~Trie() {
        clear(root);
    }


	size_t size(){
		return nbelem;
	}
	
	
	
	

    
    void insert(const set_gnodes& key) {
        
		TNode* node = root;
		bool created_something = false;
		
        for (GNode* val : key) {
            if (node->children.find(val) == node->children.end()) {
                node->children[val] = new TNode();
            }
            node = node->children[val];
        }
		
		
		nbelem++;
		
        node->isEnd = true;
    }

    // Search for a full set
    bool search(const set_gnodes& key) const {
        const TNode* node = root;
        for (GNode* val : key) {
            auto it = node->children.find(val);
            if (it == node->children.end()) return false;
            node = it->second;
        }
        return node->isEnd;
    }
	
	
	
    YBCost& operator[](const set_gnodes& key) {
        if (!search(key))
			insert(key);
		
		TNode* node = root;
        for (GNode* val : key) {
            auto it = node->children.find(val);
            if (it == node->children.end()){
				throw "key does not exist, should not happen";
			}				
            node = it->second;
        }
        return node->cost;
    }
	

    // Delete a key from the trie
    bool erase(const set_gnodes& key) {
	
		if (!search(key)){
			return false;
		}
		nbelem--;
		
        return eraseHelper(root, key.begin(), key.end());
    }


    




    // Check if there exists any key that has this prefix
    bool startsWith(const set_gnodes& prefix) const {
        const TNode* node = root;
        for (GNode* val : prefix) {
            auto it = node->children.find(val);
            if (it == node->children.end()) return false;
            node = it->second;
        }
        return true;
    }
	
	
	
	
	
	
	// ===============================
    //          ITERATOR
    // ===============================

    class Iterator {
    private:
        using MapIter = std::map<GNode*, TNode*>::iterator;

        std::vector<std::pair<TNode*, MapIter>> stack;
        std::vector<GNode*> currentPath;
        set_gnodes currentKey;
        TNode* root;

        void goToNext() {
            while (!stack.empty()) {
                auto& [node, it] = stack.back();

                // Backtrack when we've exhausted this node's children
                if (it == node->children.end()) {
                    stack.pop_back();
                    if (!currentPath.empty()) currentPath.pop_back();
                    continue;
                }

                // Visit next child
                GNode* val = it->first;
                TNode* child = it->second;
                ++it; // advance parent's iterator

                currentPath.push_back(val);
                stack.emplace_back(child, child->children.begin());

                if (child->isEnd) {
                    currentKey = set_gnodes(currentPath.begin(), currentPath.end());
                    return;
                }
            }

            currentKey.clear(); // done
        }

    public:
        Iterator(TNode* root_, bool begin) : root(root_) {
            if (begin && root_) {
                stack.emplace_back(root_, root_->children.begin());
                if (root_->isEnd) {
                    currentKey = {};
                } else {
                    goToNext();
                }
            }
        }

        set_gnodes& operator*() { return currentKey; }
        set_gnodes* operator->() { return &currentKey; }

        Iterator& operator++() {
            goToNext();
            return *this;
        }

        bool operator==(const Iterator& other) const {
            return root == other.root && currentKey == other.currentKey;
        }

        bool operator!=(const Iterator& other) const { return !(*this == other); }

        bool isEnd() const { return currentKey.empty(); }
    };

    Iterator begin() { return Iterator(root, true); }
    Iterator end() { return Iterator(root, false); }
	
	Iterator begin() const { return Iterator(root, true); }
	Iterator end()   const { return Iterator(root, false); }
	


    Iterator erase(Iterator it) {
        if (it == end()) return it;
        Iterator next = it;
        ++next;              // advance before deleting
        set_gnodes key = *it;
        erase(key);          // now safe to delete
        return next;
    }


	
    bool has_successor_with_lower_cost(set_gnodes& V) {
        auto it = V.begin();
        set_gnodes C;
        return has_successor_with_lower_cost(V, it, C, root);
    }
	
	
    bool has_successor_with_lower_cost(set_gnodes& V, set_gnodes::iterator v_it, set_gnodes& C, TNode* tcur) {
        //TODO: needs refactoring

        //C is a successor active set if it ends at an end
        if (v_it == V.end()) {
            if (tcur && tcur->isEnd && V != C){
				if ((*this)[C].cost <= (*this)[V].cost + (V.size() - C.size())) {	//TOD: assumes loss_cost = 1
					return true;
				}
			}
                
            return false;
        }
        
        if (C.empty()) {
            GNode* v_cur = *v_it;
            for (auto keyval : root->children) {
                if (v_cur->has_ancestor(keyval.first)) {
                    C.insert(keyval.first);
                    bool res = has_successor_with_lower_cost(V, std::next(v_it), C, keyval.second);
                    C.erase(keyval.first);
					if (res)
						return true;
                }
            }
			return false;
        }
        else {
            GNode* c_last = *std::prev(C.end());
            while (v_it != V.end() && (*v_it)->has_ancestor(c_last)) {
                ++v_it;
            }

            if (v_it == V.end()) {
                //C is a container active set if it ends at an end
                if (tcur && tcur->isEnd && V != C){
					if ((*this)[V].cost <= (*this)[C].cost + (V.size() - C.size())) {	//TOD: assumes loss_cost = 1
						return true;
					}
				}
					
				return false;
            }
            else {
                GNode* v_cur = *v_it;

                //iterate children in order of preorder id.  Once we are greater than v_cur->id we know it cannot be an ancestor
                //of v_cur
                auto it_children = tcur->children.begin();
                while(it_children != tcur->children.end()) {
                    GNode* child = it_children->first;

                    //we got to non-ancestors and haven't returned true yet, so no hope at this point
                    if (child->id > v_cur->id)
                        return false;
                    else {
                        if (v_cur->has_ancestor(child)) {
                            C.insert(child);
                            bool res = has_successor_with_lower_cost(V, std::next(v_it), C, it_children->second);
                            C.erase(child);
                            if (res)
                                return true;
                        }
                    }
                    ++it_children;
                }
                //previous loop that went through every child
                /*for (auto keyval : tcur->children) {
                    if (v_cur->has_ancestor(keyval.first)) {
                        C.insert(keyval.first);
                        bool res = has_successor_with_lower_cost(V, std::next(v_it), C, keyval.second);
                        C.erase(keyval.first);
						if (res)
							return true;
                    }
                }*/
				return false;
            }
        }
        
    }
	
	
	
	
	

private:
    // Recursive deletion of key
    bool eraseHelper(TNode* node, set_gnodes::const_iterator it, set_gnodes::const_iterator end) {
        if (it == end) {
            if (!node->isEnd) return false; // key not found
            node->isEnd = false;
            // if no children, signal that this node can be deleted
            return node->children.empty();
        }

        GNode* val = *it;
        auto childIt = node->children.find(val);
        if (childIt == node->children.end()) return false;

        TNode* child = childIt->second;
        bool shouldDeleteChild = eraseHelper(child, std::next(it), end);

        if (shouldDeleteChild) {
            delete child;
            node->children.erase(childIt);
        }

        // return true if node has no children and is not terminal
        return node->children.empty() && !node->isEnd;
    }

    // Recursively clear the entire trie
    void clear(TNode* node) {
        for (auto& [key, child] : node->children)
            clear(child);
        delete node;
    }
};