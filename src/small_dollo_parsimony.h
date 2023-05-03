#include "tree-lib.hpp"
#include "read_characters.h"
#include<queue>
#include<algorithm>
#include<string>
#include<cstdint>
#include<vector>
#include<fstream>
#include<list>
#include<array>


// Assume in the input tree file there is only contain one newick string.
std::string get_newick(std::string input_tree) {
  std::ifstream file(input_tree);
  
  if (file.fail()) {
    std::cout << "Failed to open the input tree file" << std::endl;
    exit(EXIT_FAILURE);
  }
  std::string line;

  std::getline(file, line);
   
  file.close();
  return line;
}



unsigned int one_step(int i, uint8_t** C, unsigned int k, boost::unordered_map<std::string, unsigned int> &label2index, Node* r, size_t max_id) {
  unsigned int res = 0;
  
  int *below = new int[max_id];
  int *above = new int[max_id];
  int *character = new int[max_id];
  
  // compute below for every vertex
  for (auto j = Traverse::PostOrder(r); j != j.end(); j++) {
          
    if ((*j)->is_leaf()) {
      std::string leaf_label = (*j)->label;
      //std::cout << "leaf_label: " << leaf_label << std::endl;
      below[(*j)->ID] = C[i][label2index[leaf_label]];
      //std::cout << "below[" << (*j)->ID <<"]" << ": " << below[(*j)->ID] << std::endl;
    } else {
      std::list<Node*> children = (*j)->get_children();
      int missing_num = 0;
      int one_num = 0;
      for (auto child : children) {
	if (below[child->ID] == 2) {
	  missing_num++;
	} else if (below[child->ID] == 1) {
	  one_num++;
	}
      }
      if (missing_num == children.size()) {
	below[(*j)->ID] = 2;
      } else if (one_num >= 1) {
	below[(*j)->ID] = 1;
      } else {
	below[(*j)->ID] = 0;
      }
    }
  }
  //std::cout << "Computed below" << std::endl;
  // compute above for every vertex
  above[r->ID] = -1;
  std::list<Node*> root_children = r->get_children();
  auto it = root_children.begin();
  
  Node* root_left = root_children.front();
  //std::cout << "computed root left" << std::endl;
  std::advance(it, 1);
  Node* root_right = *it;
  above[root_left->ID] = below[root_right->ID];
  //std::cout << "computed root right" << std::endl;
  above[root_right->ID] = below[root_left->ID];
  //std::cout << "computed root left right above" << std::endl;

  for (auto j = Traverse::PreOrder(r); j != j.end(); j++) {

    if ((*j)->ID == (r)->ID || (*j)->ID == root_left->ID || (*j)->ID == root_left->ID) {
      continue;
    }
    
    Node* parent = (*j)->get_parent();
    std::list<Node*> siblings = parent->get_children();
    
    int tot_size_to_check = siblings.size();
    int miss_num = 0;
    int one_num = 0;
    for (Node* sibling : siblings) {

      if (sibling->ID != (*j)->ID) {
	
	if (below[sibling->ID] == 2) {
	  miss_num++;
	} else if (below[sibling->ID] == 1) {
	  one_num++;
	}
      }
    }
    if (above[parent->ID] == 2) {
      miss_num++;
    } else if (above[parent->ID] == 1) {
      one_num++;
    }
    //std::cout << "Compute above:" << (*j)->ID << std::endl;
    if (miss_num == tot_size_to_check) {
      above[(*j)->ID] = 2;
    } else if (one_num >= 1) {
      above[(*j)->ID] = 1;
    } else {
      above[(*j)->ID] = 0;
    }
    //std::cout << "above[" << (*j)->ID << "]: " << above[(*j)->ID] << std::endl;
  }

  //std::cout << "Computed all above" << std::endl;
  // computing character for every vertex
  for (auto j = Traverse::PostOrder(r); j != j.end(); j++) {
    if ((*j)->is_leaf()) {
      character[(*j)->ID] = C[i][label2index[(*j)->label]];
      //std::cout << "computed leaf character[" << (*j)->ID << "]: " << character[(*j)->ID] << std::endl;
    } else {
      int one_num = 0;
      int miss_num = 0;
      std::list<Node*> children = (*j)->get_children();
      //int size_to_check = children.size();

      for (auto child : children) {
	if (below[child->ID] == 2) {
	  miss_num++;
	} else if (below[child->ID] == 1) {
	  one_num++;
	}
      }
      
      if (above[(*j)->ID] == 2) {
	miss_num++;
      } else if (above[(*j)->ID] == 1) {
	one_num++;
      }
      
      if (miss_num >= 2) {
	character[(*j)->ID] = 2;
      } else if (one_num >= 2) {
	character[(*j)->ID] = 1;
      } else {
	character[(*j)->ID] = 0;
      }
      // std::cout << "character[" << (*j)->ID << "]: " << character[(*j)->ID] << std::endl;
    }
  }
  //std::cout << "computed all character" << std::endl;
  for (auto j = Traverse::PostOrder(r); j != j.end(); j++) {
    if ((*j)->ID == r->ID) continue;
    if (character[(*j)->ID] == 0 && character[((*j)->get_parent())->ID] == 1) res++;
  }
  delete[] character;
  delete[] above;
  delete[] below;
  return res;
}


unsigned int Dollo_parsimony_score(uint8_t** C, unsigned int k,std::string input_tree, boost::unordered_map<std::string, unsigned int> &label2index) {
  unsigned int tot = 0;
  std::string newick_str = get_newick(input_tree);
  Tree *T = new Tree(newick_str);
  
  size_t max_id = (T->get_root())->get_max_id();
  //std::cout << "K: " << k << std::endl;
  for (int i = 0; i < k; i++) {
    tot += one_step(i, C, k, label2index, T->get_root(), max_id);
    //if (i % 100 == 0) std::cout << "Finished scoring " << i + 1<< " characters" << std::endl;
    std::cout << "\rComputed Score " << i + 1 << " Characters";
    std::cout.flush();
  }
  std::cout << std::endl;
  delete T;

  for (int i = 0; i < k; i++) {
    delete[] C[i];
  }
  delete[] C;
  
  return tot;
}
