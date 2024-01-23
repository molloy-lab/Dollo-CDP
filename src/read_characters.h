#include<fstream>
#include<cstdio>
#include<iostream>
#include<regex>
#include<string>
#include<boost/unordered_map.hpp>
#include <cstdint>
//#include<boost/algorithm/string.hpp>
#include<cstring>

uint8_t** read_characters(std::string filename, unsigned int &k, boost::unordered_map<std::string, unsigned int> &label2index, std::vector<std::string> &labels, boost::unordered_map<std::string, std::string> &C) {
  std::ifstream file;

  file.open(filename);
  
  if (file.fail()){
    std::cout << "Failed to open the character matrix file" << std::endl;
    exit(EXIT_FAILURE);
  }
  //std::regex nex("#NEXUS");
  std::regex ntax_reg("ntax=(\\d*)");
  std::regex nchar_reg("nchar=(\\d*)");
  //std::regex taxon_reg("Taxon([A-Z]+)|Out");
  //std::regex taxon_char_reg("(\\d)");
  unsigned int taxons_num = 0;
  //unsigned int chars_num = 0;
  unsigned int total_taxons_num = 0;
  // boost::unordered_map<std::string, std::string> taxon_chars;
  uint8_t** taxon_chars;
  std::string line;
  std::smatch match;
  unsigned int line_num = 1;
  unsigned int end;

  while (!file.eof()) {
    std::getline(file, line);
  	
    if (std::regex_search(line, match, ntax_reg)) {
      total_taxons_num = stoi(match[1]);
      end = total_taxons_num + 7;
      
    }
    
    if (std::regex_search(line, match, nchar_reg)) {
      k = stoi(match[1]);
      
      taxon_chars = new uint8_t*[k];
      //std::cout << "taxons num is " << total_taxons_num << std::endl;
      for (int i = 0; i < k; i++)
	taxon_chars[i] = new uint8_t[total_taxons_num];
    std::cout << "finished creating transposed matrix " << std::endl;
    }

    if (line_num >= 7 && line_num < end) {
      /*
      if (std::regex_search(line, match, taxon_reg)) {
	std::string tax = match[0];
	label2index[tax] = taxons_num;
	labels.push_back(tax);
      
	if (std::regex_search(line, match, taxon_char_reg)) {
	  int index = match.position();
	  C[tax] = line.substr(index);
	  for (int i = index; i < line.length(); i++) {
	    //std::cout << "index is " << index << "i is " << i << std::endl;
	    if (line[i] == '1')
	      taxon_chars[i - index][taxons_num] = 1;
	    if (line[i] == '0')
	      taxon_chars[i - index][taxons_num] = 0;
	  }
	
	}
      */
      unsigned int curr_pos = 0;
      while (!std::isblank(line[curr_pos])) {
	curr_pos++;
      
      }
    
      std::string tax = line.substr(0,curr_pos);
     

      while (std::isblank(line[curr_pos])) {
	curr_pos++;
	
      }
      
      C[tax] = line.substr(curr_pos);
      label2index[tax] = taxons_num;
      labels.push_back(tax);
      
      
      for (int i = curr_pos; i < line.length(); i++) {
	if (line[i] == '1')
	  taxon_chars[i - curr_pos][taxons_num] = 1;
	if (line[i] == '0')
	  taxon_chars[i - curr_pos][taxons_num] = 0;    
	if (line[i] == '?')
	  taxon_chars[i - curr_pos][taxons_num] = 2;
      }
      
      
	taxons_num++;
	
      }
    line_num++;
    
  }
  //std::cout << "The size of label2index: " << label2index.size() << std::endl;
  file.close();
  /*
  for (const auto& pair : label2index) {
    std::cout << pair.first << ": " << pair.second << std::endl;
  }
  */
  std::cout << "finished reading characters matrix " << std::endl;
  return taxon_chars;
}


void write_newick_from_C(std::ostream &os, unsigned int &k, std::vector<std::string> &labels, uint8_t** C) {
  for (unsigned int i = 0; i < k; i++) {
    unsigned int num_of_one_taxons = 0;
    unsigned int num_of_zero_taxons = 0;
    std::string str1 = "(";
    std::string str0 = "(";
    for (unsigned int j = 0; j < labels.size(); j++) {
      if (C[i][j] == 1) {
	if (num_of_one_taxons == 0) {
	  str1 += labels[j];
	} else {
	  str1 += "," + labels[j];
	}
	num_of_one_taxons++;
      } else if (C[i][j] == 0) {
	if (num_of_zero_taxons == 0) {
	  str0 += labels[j];
	} else {
	  str0 += "," + labels[j];
	}
	num_of_zero_taxons++;
      }
    
   }
    if (num_of_one_taxons > 1 && num_of_zero_taxons > 1) {
      os << str0 << "," << str1 << "))" << std::endl;
      
    }
  }
}

