/*

MIT License 
Copyright (c) 2023 Junyan Dai, Tobias Rubel, Yunheng Han, Erin Molloy

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

#include <algorithm>
#include <boost/dynamic_bitset/dynamic_bitset.hpp>
#include <chrono>
#include <cstdlib>
#include <cassert>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <list>
#include <iomanip>
#include <climits>
#include <stack>
#include <boost/dynamic_bitset.hpp>
#include <cstdint>
// Allows user to set at compile time
// Use compile flags -DUSE_SHRT or -DUSE_LONG
#ifdef USE_SHRT
    typedef uint16_t index_t;
    #define INDEX_MAX USHRT_MAX
#elif USE_LONG
    typedef uint64_t index_t;
    #define INDEX_MAX ULONG_MAX
#else
    typedef uint32_t index_t;
    #define INDEX_MAX UINT_MAX
#endif  // USE_USHRT or USE_ULONG


#define BYTES_PER_GB 1073741824


int verbose = 0, tid = 0;
std::ofstream logs[8];
size_t MAX_ID = 0;

class Bipartition;
// Adds elements of list 2 to list 1 (does not copy!)
template <typename T> void extend(T &list1, T &list2) {
    typename T::iterator it;
    for (it = list2.begin(); it != list2.end(); ++it) {
        list1.push_back(*it);
    }
}


class Node {
    public:
        Node();
        Node(std::string name);
        ~Node();
        bool is_root();
        bool is_leaf();
        unsigned int num_children();
        Node* get_parent();
        Node* copy_beneath();
        void update_parent(Node * p);
        void add_child(Node *child);
        void remove_child(Node *child);
        void contract();
        void add_children_to_list(std::list<Node*> &nodelist);
        void update_label_list(std::list<std::string> labellist);
        void update_label_list(std::string str);
        std::list<std::string> get_label_list();
        void suppress_unifurcations();
        void compute_c();
        void update_c(Node *child);
        std::list<Node *> get_children() {
            return children;
        }
        unsigned int get_max_id() {
	  return MAX_ID;
	}
        unsigned int get_c();
        std::string newick(bool printindex=false);
        std::string label;  // TODO: we probably want to get rid of this at some point
                            // so we aren't storing so many copies of labels!
                            // this could be done if we read trees onto the same
                            // label set... 
        index_t index;
        index_t size;
        std::list<Node*> children;
        std::list<std::string> label_list;
        index_t c = 0;
        Node *parent;
        size_t ID;
    private:
        unsigned int x, y;
        friend class Tree;
};


namespace Traverse {
struct ToRoot
{
    // constructor that takes in a node
    ToRoot() {
        previous_node = NULL;
        current_node = NULL;
    };

    ToRoot(Node *node) {
        previous_node = NULL;
        current_node = node;
    };

    // incrementing means going to the parent node
    ToRoot &operator++() noexcept
    {
        if (current_node != NULL) {
            previous_node = current_node;
            current_node = current_node->get_parent();
        }
        return *this;
    };

    // post fixing is bad in general but it has it's usages
    ToRoot operator++(int) noexcept
    {
        ToRoot tempIter = *this;   // make a copy of the iterator
        ++*this;                   // increment
        return tempIter;           // return the copy before increment
    };

    // compare nodes
    bool operator!=(const ToRoot &other) const noexcept
    {
        return this->current_node != other.current_node;
    };

    // return the node (dereference operator)
    Node* operator*() const noexcept
    {
        return this->current_node;
    };

    // return a const pointer to the front
    ToRoot begin() const noexcept
    {
        return ToRoot(this->current_node);
    };

    // return a const pointer to the back - the back is always null
    ToRoot end() const noexcept
    {
        return ToRoot();
    };

    private:
        Node *previous_node = NULL;
        Node *current_node = NULL;
};

struct PreOrder
{
    // constructor that takes in a node
    PreOrder() {
        previous_node = NULL;
        current_node = NULL;
    };

    PreOrder(Node *node) {
        previous_node = NULL;
        current_node = node;
        current_node->add_children_to_list(nodelist);
        // TODO: Maybe be able to get rid of storage... ?
    };

    // incrementing means going to the parent node
    PreOrder &operator++() noexcept
    {
        previous_node = current_node;
        if (nodelist.size() == 0) {
            current_node = NULL;
            nodelist.clear();
        } else {
            current_node = nodelist.front();
            nodelist.pop_front();
            current_node->add_children_to_list(nodelist);
        }

        return *this;
    };

    // post fixing is bad in general but it has it's usages
    PreOrder operator++(int) noexcept
    {
        PreOrder tempIter = *this;  // make a copy of the iterator
        ++*this;                    // increment
        return tempIter;            // return the copy before increment
    };

    // compare nodes
    bool operator!=(const PreOrder &other) const noexcept
    {
        return this->current_node != other.current_node;
    };

    // return the node (dereference operator)
    Node* operator*() const noexcept
    {
        return this->current_node;
    };

    // return a const pointer to the front
    PreOrder begin() const noexcept
    {
        return PreOrder(this->current_node);
    };

    // return a const pointer to the back - the back is always null
    PreOrder end() const noexcept
    {
        return PreOrder();
    };

    private:
        Node *previous_node = NULL;
        Node *current_node = NULL;
        std::list<Node*> nodelist;
};


struct PostOrder
{
    // constructor that takes in a node
    PostOrder() {
        previous_node = NULL;
        current_node = NULL;
    };

    PostOrder(Node *node) {
        Node *tmp;

        previous_node = NULL;

        nodelist1.push_back(node);

        while (nodelist1.size() != 0) {
            tmp = nodelist1.back();
            nodelist1.pop_back();
            tmp->add_children_to_list(nodelist1);
            nodelist2.push_back(tmp);
        }

        nodelist1.clear();
        // TO DO: Maybe be able to get rid of stoarge!

        while (nodelist2.size() != 0) {
            tmp = nodelist2.back();
            nodelist2.pop_back();

            current_node = tmp;
            break;
        }
    };

    // incrementing means going to the next node in the postorder traversal
    PostOrder &operator++() noexcept
    {
        Node *node;

        previous_node = current_node;
        if (nodelist2.size() == 0) {
            current_node = NULL;
            nodelist2.clear();
        } else {
            node = nodelist2.back();
            nodelist2.pop_back();
            current_node = node;
        }

        return *this;
    };

    // post fixing is bad in general but it has it's usages
    PostOrder operator++(int) noexcept
    {
        PostOrder tempIter = *this;   // make a copy of the iterator
        ++*this;                      // increment
        return tempIter;              // return the copy before increment
    };

    // compare nodes
    bool operator!=(const PostOrder &other) const noexcept
    {
        return this->current_node != other.current_node;
    };

    // return the node (dereference operator)
    Node* operator*() const noexcept
    {
        return this->current_node;
    };

    // return a const pointer to the front
    PostOrder begin() const noexcept
    {
        return PostOrder(this->current_node);
    };

    // return a const pointer to the back - the back is always null
    PostOrder end() const noexcept
    {
        return PostOrder();
    };

    private:
        Node *previous_node = NULL;
        Node *current_node = NULL;
        std::list<Node*> nodelist1;
        std::list<Node*> nodelist2;
};


struct Leaves
{
    // constructor that takes in a node
    Leaves() {
        previous_node = NULL;
        current_node = NULL;
    };

    Leaves(Node *node) {
        Node *tmp;

        previous_node = NULL;

        nodelist.push_back(node);

        while (nodelist.size() > 0) {
            tmp = nodelist.front();
            nodelist.pop_front();
            tmp->add_children_to_list(nodelist);

            if (tmp->is_leaf()) {
                current_node = tmp;
                break;
            }
        }
    };

    // incrementing means going to the next leaf node
    Leaves &operator++() noexcept
    {
        Node *node;

        previous_node = current_node;

        if (nodelist.size() == 0) {
            current_node = NULL;
            nodelist.clear();
        } else {
            while (nodelist.size() > 0) {
                node = nodelist.front();
                nodelist.pop_front();
                node->add_children_to_list(nodelist);

                if (node->is_leaf()) {
                    current_node = node;
                    break;
                }
            }
        }

        return *this;
    };

    // post fixing is bad in general but it has it's usages
    Leaves operator++(int) noexcept
    {
        Leaves tempIter = *this;   // make a copy of the iterator
        ++*this;                   // increment
        return tempIter;           // return the copy before increment
    };

    // compare nodes
    bool operator!=(const Leaves &other) const noexcept
    {
        return this->current_node != other.current_node;
    };

    // return the node (dereference operator)
    Node* operator*() const noexcept
    {
        return this->current_node;
    };

    // return a const pointer to the front
    Leaves begin() const noexcept
    {
        return Leaves(this->current_node);
    };

    // return a const pointer to the back - the back is always null
    Leaves end() const noexcept
    {
        return Leaves();
    };

    private:
        Node *previous_node = NULL;
        Node *current_node = NULL;
        std::list<Node*> nodelist;
};

}; // namespace traversal

/*
std::vector<Node *> TPO_helper(Node *root, std::vector<Node*> ret) {
    if (root->is_leaf()) {
        ret.push_back(root);
    } else {
        for (auto child: root->children) {
            TPO_helper(child,ret);
        }
        ret.push_back(root);
    }
    return ret;
}

std::vector<Node*> TPO(Node *root) {
    std::vector<Node*> ret;
    return TPO_helper(root,ret);
} 
*/




Node::Node() {
    parent = NULL;
    label = "";
    index = INDEX_MAX;
    size = 0;
    ID = MAX_ID++;
    
}

Node::Node(std::string name) {
    parent = NULL;
    label = name;
    index = INDEX_MAX;
    size = 0;
    ID = MAX_ID++;
    
}


Node::~Node() {
    std::list<Node*>::iterator it;
    for (it = children.begin(); it != children.end(); ++it) {
        this->remove_child(*it);
    }
    if (parent != NULL) {
        parent->remove_child(this);
    }
}

bool Node::is_root() {
    if (parent == NULL) return true;
    return false;
}

bool Node::is_leaf() {
    if (children.size() == 0) return true;
    return false;
}

unsigned int Node::num_children() {
    return children.size();
}

unsigned int Node::get_c() {
    return c;
}

void Node::update_label_list(std::list<std::string> newnodes) {
    label_list.merge(newnodes);
}

void Node::update_label_list(std::string str) {
    label_list.push_back(str);
}

std::list<std::string> Node::get_label_list() {
    return label_list;
}

Node* Node::get_parent() {
    return parent;
}

void Node::update_parent(Node* p) {
    parent = p;

}
void Node::add_child(Node *child) {
    if (child == NULL) return;

    child->parent = this;
    children.push_back(child);

}
void Node::update_c(Node *child) {
    c += child->get_c();
}

void Node::remove_child(Node *child) {
     if (child == NULL) return;

    child->parent = NULL;
    children.remove(child);
}

void Node::compute_c() {
    auto nodeItr = Traverse::PostOrder(this);
    for (; nodeItr != nodeItr.end(); ++nodeItr) {
        if ((*nodeItr)->is_leaf()) {
            (*nodeItr)->c = 1;
        continue;
        }
        std::list<Node*>::iterator it;
        for (it = (*nodeItr)->children.begin(); it != (*nodeItr)->children.end(); ++it) {
            (*nodeItr)->update_c(*it);
        }
    }
}

void Node::contract() {
    if (parent == NULL) return;

    std::list<Node*>::iterator it;
    for (it = children.begin(); it != children.end(); ++it) {
        parent->add_child(*it);
    }
    parent->remove_child(this);

    parent = NULL;
    children.clear();

    delete this;
}


void Node::add_children_to_list(std::list<Node*> &nodelist) {
    extend(nodelist, children);
}

void Node::suppress_unifurcations() {

    auto nodeItr = Traverse::PreOrder(this);
    for (; nodeItr != nodeItr.end(); ++nodeItr) {
        if ((*nodeItr)->num_children() == 1) (*nodeItr)->contract();
    }
}


std::string Node::newick(bool printindex) {
    // TODO: Maybe change to pass by reference?

    std::list<Node*>::iterator it;
    std::string out;

    if (this->is_leaf()) {
        out = "";
        if (!this->label.empty()) {
            out += this->label;
            if (printindex) {
                out += ':';
                out += std::to_string(this->index);
            }         
        }
        return out;
    }

    out = '(';

    for (it = this->children.begin(); it != this->children.end(); ++it) {
        out += (*it)->newick(printindex);
        out += ',';
    }
    out.pop_back();  // Drop trailing comma
    out += ')';

    if (!this->label.empty()) out += this->label;

    return out;
}

Node* Node::copy_beneath() {
    return NULL;
}


class Tree {
    
    public:
        Tree();
        Tree(Node *myroot);
        Tree(const std::string &text);

        Tree(std::vector<std::string> labels, std::unordered_map<Bipartition, unsigned int> &memo, std::unordered_map<Bipartition, std::vector<std::pair<Bipartition, unsigned int>>> &decisions, const Bipartition &state) {
            root = build_tree(labels, memo, decisions, state);
        }

        ~Tree();
        Node* get_root();
        void suppress_unifurcations();
        void compute_c();
        std::string newick(bool printindex=false);
        Tree* get_induced_subtree_copy(std::unordered_set<std::string> taxa);
        Tree* get_induced_subtree_copy_2(std::unordered_set<std::string> taxa);
        void _copy_tree_helper(Node* curr_old, Node* curr_new);
        Tree* copy_tree();

        std::unordered_map<Node *, std::vector<std::string>> get_subtrees() {
            std::unordered_map<Node *, std::vector<std::string>> subtrees;
            get_subtrees(root, subtrees);
            return subtrees;
        }

        std::vector<std::string> get_clade() {
            return get_clade(root);
        }

        unsigned int get_score(std::unordered_map<std::string, index_t> &label2index, const Bipartition &a, const Bipartition &b) {
            return get_score(root, label2index, a, b);
        }
        
        Node *root;
    private:

        std::vector<std::string> get_subtrees(Node *root, std::unordered_map<Node *, std::vector<std::string>> &subtrees) {
            if (root->children.size() == 0) {
                std::vector<std::string> clade({root->label});
                subtrees[root] = clade;
                return clade;
            }
            std::vector<std::string> clade;
            for (Node *child : root->children) {
                std::vector<std::string> labels = get_subtrees(child, subtrees);
                for (std::string label : labels) 
                    clade.push_back(label);
            }
            subtrees[root] = clade;
            return clade;
        }

        std::vector<std::string> get_clade(Node *root) {
            if (root->children.size() == 0) 
                return std::vector<std::string>({root->label});
            std::vector<std::string> clade;
            for (Node *child : root->children) {
                std::vector<std::string> labels = get_clade(child);
                for (std::string label : labels) 
                    clade.push_back(label);
            }
            return clade;
        }

        Node *build_tree(std::vector<std::string> &labels, std::unordered_map<Bipartition, unsigned int> &memo, std::unordered_map<Bipartition, std::vector<std::pair<Bipartition, unsigned int>>> &decisions, const Bipartition &state);
        unsigned int get_score(Node *root, std::unordered_map<std::string, index_t> &label2index, const Bipartition &a, const Bipartition &b);
};


Tree::Tree() {
    root = new Node();
}


Tree::Tree(Node *myroot) {
    root = myroot;
}


Tree::Tree(const std::string &text) {
    Node *node, *child;
    std::string label;
    unsigned int i;

    root = new Node();
    node = root;

    i = 0;
    while (i < text.size()) {
        if (text[i] == ';') {
            // End of Newick string
            if ((i != text.size() - 1) || (node != root)) {
                std::cerr << "Not a valid newick!" << std::endl;
                exit(1);
            }  
        } else if (text[i] == '(') {
            // Go to new child
            child = new Node();
            node->add_child(child);
            node = child;
        } else if (text[i] == ')') {
            // Go to parent
            node = node->get_parent();
        } else if (text[i] == ',') {
            // Go to new sibling
            node = node->get_parent();
            child = new Node();
            node->add_child(child);
            node = child;
        } else if (text[i] == ':') {
            // Parse edge length
            i += 1;
            while ((text[i] != ',') && (text[i] != ')') && (text[i] != ';'))
                i += 1;
            i -= 1;
        } else if ((text[i] == ' ') || (text[i] == '\t')) {
            continue;
        } else {
            // Parse node label
            label.clear();
            while ((text[i] != ':') && (text[i] != ',') && 
                    (text[i] != ')') && (text[i] != ';')) {
                label += text[i];
                i += 1;
            }
            i -= 1;
            node->label = label;
        }
        i += 1;
    }

    root->suppress_unifurcations();
}


Tree::~Tree() { 
    auto nodeItr = Traverse::PostOrder(root);
    for (; nodeItr != nodeItr.end(); ++nodeItr) {
        delete *nodeItr;
    }
}


Node* Tree::get_root() {
    return root;
}


void Tree::suppress_unifurcations() {
    std::vector<Node*> original_children;
    std::list<Node*>::iterator it;
    for (it = root->children.begin(); it != root->children.end(); ++it) {
        original_children.push_back(*it);
    }
    for (Node* n : original_children) {
        n->suppress_unifurcations();
    }
}

void Tree::compute_c() {
    root->compute_c();
}


std::string Tree::newick(bool printindex) {
    std::string out = this->root->newick(printindex);
    out += ";";
    return out;
}

void Tree::_copy_tree_helper(Node* curr_old, Node* curr_new) {
    for (Node* chld: curr_old->children) {
        Node *c_new = new Node();
        c_new->update_parent(curr_new);
        curr_new->add_child(c_new);
        c_new->label = chld->label;
        _copy_tree_helper(chld,c_new);
    }    
}

Tree* Tree::copy_tree() {
    Tree *out = new Tree();
    Node* curr_old = this->root;
    Node* curr_new = out->root;
    _copy_tree_helper(curr_old,curr_new);
    return out;
   
}

Tree* Tree::get_induced_subtree_copy(std::unordered_set<std::string> taxa) {
    Tree * out = copy_tree();
     auto leafItr = Traverse::Leaves(out->root);
     for (; leafItr != leafItr.end(); ++leafItr) {
         auto leaf = *leafItr;
         if (taxa.count(leaf->label) == 0) {
            leaf->parent->remove_child(leaf);
            delete leaf;
            out->suppress_unifurcations();
         }
     }

    out->suppress_unifurcations();

    return out;
}

class Forest {
    public:
        Forest();
        Forest(std::vector<Tree*> &mytrees);
        ~Forest();
        unsigned int num_trees();
        unsigned int num_labels();
        void compute_c();
        std::vector<Tree*> fetch_trees();
        Forest get_induced_subforest_copy(std::unordered_set<std::string> taxa);
        std::vector<Tree*> trees;
        std::vector<std::string> index2label;
        std::unordered_map<std::string, index_t> label2index;
};


Forest::Forest(std::vector<Tree*> &mytrees) {
    std::unordered_map<std::string, index_t>::const_iterator mapIter;
    Node *root, *leaf;
    std::string label;
    index_t index;

    trees = mytrees;

    for (unsigned int i = 0; i < trees.size(); i++) {
        auto leafItr = Traverse::Leaves(trees[i]->get_root());
        for (; leafItr != leafItr.end(); ++leafItr) {
            leaf = *leafItr;
            label = leaf->label;

            mapIter = label2index.find(label);
            if ( mapIter == label2index.end() ) {
                index = index2label.size();
                leaf->index = (index_t) index;
                label2index.insert({label, index});
                index2label.push_back(label);
            }
            else {
                leaf->index = label2index.at(label);
            }
        }
    }
    /*
    std::cout << "Displaying label2index" << std::endl;
    for (auto l: label2index) {
        std::cout << l.first << ": " << l.second <<std::endl; 
    }
    std::cout << "Displaying index2label" << std::endl;
    for (int i=0; i < index2label.size(); i++) {
        
        std::cout << i << ": " << index2label[i] <<std::endl; 
    }
    std::cout << "searching for something fishy in the trees:" <<std::endl;
    //for (unsigned int i = 0; i < trees.size(); i++) {
    for (unsigned int i = 0; i < 1; i++) {
        auto leafItr = Traverse::Leaves(trees[i]->get_root());
         for (; leafItr != leafItr.end(); ++leafItr) {
            leaf = *leafItr;
            label = leaf->label;
            index = leaf->index;
            std::cout << label << ": " << index << std::endl;
            if (label2index[label] != index || index2label[index] != label) {
                std::cout << "found something fucked up!" << std::endl;
            }
         }
    }
    */


}


Forest::~Forest() {

}


unsigned int Forest::num_trees() {
    return trees.size();
}

std::vector<Tree*> Forest::fetch_trees() {
    return trees;
}


unsigned int Forest::num_labels() {
    return index2label.size();
}

void Forest::compute_c() {
    for (unsigned int i = 0; i < num_trees(); i++ ) {
        trees[i]->compute_c();
    }
}

Forest Forest::get_induced_subforest_copy(std::unordered_set<std::string> taxa) {
    std::vector<Tree*> tree_vec;
    for (Tree* t : trees) {
        auto nt = t->get_induced_subtree_copy(taxa);
        tree_vec.push_back(nt);
    }
    auto out = Forest(tree_vec);
    return out;
}

