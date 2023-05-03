#include "read_clades.h"
#include<tuple>
#include<string>
#include<vector>
#include<algorithm>
#include<queue>
#include<unordered_set>
#include <cstdint>

typedef Bipartition clades;
typedef std::unordered_set<Bipartition> clades_set;
typedef uint8_t** character_matrix;
typedef std::tuple<std::string, std::string> pp;
typedef boost::unordered_map<std::string, std::string> character_record;


struct TupleHash {
    template <typename... Args>
    std::size_t operator()(const std::tuple<Args...>& tuple) const {
        std::size_t seed = 0;
        boost::hash_combine(seed, tuple);
        return seed;
    }
};

// Define a custom equality function for tuples
struct TupleEqual {
    template <typename... Args>
    bool operator()(const std::tuple<Args...>& lhs, const std::tuple<Args...>& rhs) const {
        return lhs == rhs;
    }
};


typedef boost::unordered_map<pp, unsigned int, TupleHash, TupleEqual> opt;
typedef boost::unordered_map<std::string, std::vector<std::string>> states_map;
typedef boost::unordered_map<pp, std::string, TupleHash, TupleEqual> state_record;
typedef boost::unordered_map<pp, clades, TupleHash, TupleEqual> child_record;

const unsigned int INF = 0x7f7f7f7f;


std::tuple<unsigned int,std::string> bottom_up(clades_set X, character_record &record, character_matrix C, opt &f, states_map &St, state_record &st_l, state_record &st_r, child_record &g, unsigned int k, clades &S, std::vector<std::string> &labels, boost::unordered_map<std::string, unsigned int> &label2index) {
  auto start = std::chrono::high_resolution_clock::now();
  std::vector<Bipartition> X_vec(X.begin(), X.end());
  sort(X_vec.begin(), X_vec.end(), [](clades &s1, clades &s2) {
    return s1.count() < s2.count();
  });

  std::cout << std::endl;
  std::cout << "size of X: " << X_vec.size() << std::endl;
  size_t subp_num = 1;
  unsigned int score_fun_num = 0;
  for (auto u : X_vec) {
    std::string u_str = u.to_string();
    if (u.count() == 1) {
      std::string taxon = u.to_labels(labels);
      
      St[u_str] = {record[taxon]};

      f[pp(u_str, record[taxon])] = 0;
      if (subp_num % 100 == 0) std::cout << std::setw(12) << subp_num << " subsubproblems computed." << std::endl;
      std::cout << "\rComputing subsubproblem:" << std::setw(12) << subp_num++;
      std::cout.flush();
    } else if (u.count() >= 2) {
      clades_set memo;
      memo.clear();
      for (auto& a : X_vec) {

	if (a.count() >= u.count()) {
	  break;
	}
	
	if ((a.get_bitset().is_subset_of(u.get_bitset()))) {
	  clades a_comp = a.complement(u);
	  
	  auto ptr_to_memo = memo.find(a_comp);

	  if (ptr_to_memo != memo.end()) {
	    continue;
	  } else {
	    memo.insert(a);
	  }
	  
	  auto it1 = X.find(a_comp);
	 
	  if (it1 != X.end()) {
	    
	    std::string t = get_state(u, a, k, C);
	    
	    auto it2 = St.find(u_str);
	    if (it2 != St.end()) {
	      St[u_str].push_back(t);
	    } else {
	      St[u_str] = {t};
	    }
	    
	    unsigned int tmp = INF;
	    clades h;
	    std::string i, j;
	    unsigned int v;
	    std::string a_str = a.to_string();
	    std::string a_comp_str = a_comp.to_string();
	    for (const auto &s1 : St[a_str]) {
	      
	      for (const auto &s2 : St[a_comp_str]) {
		
		unsigned int score_tmp = score(t, s1, s2, k);
		
		v = f[pp(a_str,s1)] + f[pp(a_comp_str,s2)] + score_tmp;
		
		
		if (tmp > v) {
		  
		  tmp = v;
		  h = a;
		  i = s1;
		  j = s2;
		}
	      }
	    }
	    
	    auto it3 = f.find(pp(u_str,t));
        
	    if (it3 == f.end() || f[pp(u_str,t)] > tmp) {

	      if (subp_num % 100 == 0) std::cout << std::setw(12) << subp_num << " subsubproblems computed." << std::endl;
	      std::cout << "\rComputing subsubproblem:" << std::setw(12) << subp_num++;
	      std::cout.flush();

	      f[pp(u_str,t)] = tmp;
	      g[pp(u_str,t)] = h;
	        
	      st_l[pp(u_str,t)] = i;
	      st_r[pp(u_str,t)] = j;
	    }
	  }
	}
      }
    }
  }
  unsigned int res = INF;
  std::string res_state = "Error State";
  
  std::string S_str = S.to_string();
  
  for (const auto &st : St[S_str]) {
    if (res > f[(pp(S_str,st))]) {
      res = f[(pp(S_str,st))];
      res_state = st;
    } 
  }
  std::cout << std::endl;
  std::cout << "The Dollo parsimony score: " << res << std::endl;
  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  std::cout << "execution time for DP part: " << duration.count() << "ms" << std::endl;
  return std::tuple<int, std::string>{res, res_state};
 
}

Bipartition get_bip(std::string str, unsigned int size) {
  boost::dynamic_bitset<> bs(size);
  for (unsigned int i = 0; i < size; i++) {
    if (str[i] == '1') {
      bs.set(size - 1 - i);
    }
  }
  
  return Bipartition(bs);
}

Node* back_solve(opt &f, state_record &st_l, state_record &st_r, child_record &g, pp start_pair, std::vector<std::string> &labels) {
  auto start = std::chrono::high_resolution_clock::now();
  std::queue<std::tuple<pp, Node*>> que;
  Node *r = new Node();
  que.push(std::tuple<pp, Node*>{start_pair, r});
  while (!que.empty()) {
    std::tuple<pp,Node*> t = que.front();
    que.pop();
    pp tp = std::get<0>(t);
    clades full = get_bip(std::get<0>(tp), labels.size());
    

    Node *parent_ptr = std::get<1>(t);
    if (full.count() == 1) {
      Node *new_child_ptr = new Node(full.to_labels(labels)); 
      parent_ptr->add_child(new_child_ptr);
      
    } else if (full.count() >= 2) {
      
      clades left = g[tp];
      
      clades right = left.complement(full);
      
      pp new_pair1 = pp{left.to_string(), st_l[tp]};
      pp new_pair2 = pp{right.to_string(), st_r[tp]};
      Node *left_ptr = new Node();
      Node *right_ptr = new Node();
      parent_ptr->add_child(left_ptr);
      parent_ptr->add_child(right_ptr);
      que.push(std::tuple<pp, Node*>{new_pair1, left_ptr});
      que.push(std::tuple<pp, Node*>{new_pair2, right_ptr});
    }
    
    
  }
  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  std::cout << "execution time to build tree: " << duration.count() << "ms" << std::endl;
  return r;
  
}

Tree constraint_large_dollo_parsimony(clades_set X, character_record record, character_matrix C, int k, std::vector<std::string> &labels, boost::unordered_map<std::string, unsigned int> &label2index) {
  opt f;
  state_record st_l;
  state_record st_r;
  states_map St;
  child_record g;

  // compute the full set S
  
  boost::dynamic_bitset<> tbs(labels.size());

  tbs.flip();
  
  Bipartition S(tbs);

  
  std::tuple<int, std::string> vp = bottom_up(X, record, C, f, St, st_l, st_r, g, k, S, labels, label2index);

  for (int i = 0; i < k; i++) {
    delete[] C[i];
  }
  delete[] C;
  
  std::string start_state = std::get<1>(vp);
  
  Node* r = back_solve(f, st_l, st_r, g, pp{S.to_string(), start_state}, labels);
  return Tree(r);
}
