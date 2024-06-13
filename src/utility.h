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

#include "bipartition.hpp"
#include<string>
#include<vector>
#include<algorithm>
#include<string.h>
#include<map>
#include<set>
#include<sstream>
#include "small_dollo_parsimony.h"
#include<cstdint>

unsigned int score(std::string s1, std::string s2, std::string s3, int k) {
  int R = 0;
  for (int i = 0; i < k; i++) {
    if (s1[i] == '1' && s2[i] == '0')
      R++;
    if (s1[i] == '1' && s3[i] == '0')
      R++;
  }
  return R;

}

// If there is some taxon h
bool get_state_aux(Bipartition A, int i, uint8_t** C) {
  
  for (int j = 0; j < A.size(); j++) {
    if (C[i][j] == 1 && A.contain_index(j))
      return true;
  }

  return false;
  
}

bool all_missing_data(Bipartition A, int i, uint8_t** C) {
  int count = 0;
  for (int j = 0; j < A.size(); j++) {
    if (C[i][j] == 2 && A.contain_index(j)) {
      count++;
    }
  }

  return count == A.count();

}
// Input: A, B, k, C, S
// output: the state when B and A/B is the two children of A in the arbitrary opt tree using k characters.
std::string get_state(Bipartition A, Bipartition B, int k, uint8_t** C) {
  int V[k];
  boost::dynamic_bitset<> tbs(A.size());
  tbs.flip();
  Bipartition S(tbs); 
  memset(V, 0, sizeof(V));
  for (int i = 0; i < k; i++) {
    int flag = 0;

    if (all_missing_data(A, i, C)) {
      V[i] = 2;
      continue;
    }

    Bipartition D = B.complement(A);

    if (get_state_aux(B, i, C)) flag++;

    if (get_state_aux(D, i, C)) flag++;

    Bipartition E = A.complement(S);

    if (get_state_aux(E, i, C)) flag++;

    if (flag >= 2) V[i] = 1;
    
  }
  std::string R;
  R.reserve(k);
  for (int i = 0; i < k; i++)
    R.push_back(V[i] + '0');
  return R;

}
