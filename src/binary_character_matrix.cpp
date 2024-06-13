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


/**
 * @file binary_character_matrix.cpp
 * @brief Binary Character Matrix Class Functions
 */


#include "binary_character_matrix.hpp"

using phylotools::BinaryCharacterMatrix;


size_t BinaryCharacterMatrix::countColumns(std::istream &is) {
    std::string buff;
    std::string::iterator iter;
    size_t nchr = 0;
    char c;

    while (std::getline(is, buff)) {
        if (buff[0] == '>') {
            break;
        }
        for (iter = buff.begin(); iter != buff.end(); ++iter) {
            c = *iter;
            switch (c) {
                case '0':
                    ++nchr;
                    break;
                case '1':
                    ++nchr;
                    break;
                default:
                    break;
            }
        }
    }
    return nchr;
}


void BinaryCharacterMatrix::readPhylip(std::istream &is) {
    std::string line, word;
    size_t nseq, nchr, last, i, j;
    char c;

    // Try to read as PHYLIP file
    std::getline(is, line);
    std::stringstream ss(line);
    ss >> word;
    try {
        nseq = std::stoi(word);
    }
    catch (const std::invalid_argument& ia) {
        std::cerr << "Unable to read PHYLIP file!" << std::endl;
        exit(1);
    }

    ss >> word;
    try {
        nchr = std::stoi(word);
    }
    catch (const std::invalid_argument& ia) {
        std::cerr << "Unable to read PHYLIP file!" << std::endl;
        exit(1);
    }

    labels_.resize(nseq);
    splits_.resize(nseq);
    for (i = 0; i < nseq; i++) {
        splits_[i].reserve(nchr);
    }

    // Read file
    i = 0;
    while (std::getline(is, line)) {
        if (i == nseq) {
            std::cerr << "Wrong number of sequences!" << std::endl;
            exit(1);
        }

        // Read label
        for (j = 0; j < line.length(); j++) {
            c = line[j];
            if ((c == ' ') || (j == 10)) {
                last = j;
                break;
            }            
            labels_[i].push_back(c);
        }

        // Read sequence
        for (j = last; j < line.length(); j++) {
            c = line[j];
            if ((c == '0') || (c == '1') || (c == '?')) {
                    splits_[i].push_back(c);
            }
        }

        if (splits_[i].size() != nchr) {
            std::cerr << "Sequence " << labels_[i - 1]
                      << " is the wrong length!" << std::endl;
            exit(1);
        }

        ++i;
    }
}


void BinaryCharacterMatrix::readFasta(std::istream &is) {
    std::string line;
    size_t nchr, last, i;
    int nseq;
    char c;

    // Try to read as FASTA file
    std::getline(is, line);
    if (line[0] != '>') {
        std::cerr << "Unable to read FASTA file!" << std::endl;
        exit(1);
    }

    // Count number of splits
    nseq = -1;
    last = 0;
    nchr = this->countColumns(is);
    is.seekg(0, std::ios::beg);

    // Read file
    while (std::getline(is, line)) {
        if (line[0] == '>') {
            // Read label
            if (nseq > -1) {
                if (last != nchr) {
                    std::cerr << "Sequence " << labels_[nseq - 1]
                              << " is the wrong length!" << std::endl;
                    exit(1);
                }
                
            }
            ++nseq;
            labels_.push_back(line.substr(1));
            splits_.push_back("");
            splits_[nseq].reserve(nchr);
            last = 0;
        } else {
            // Read sequence
            for (i = 0; i < line.length(); i++) {
                if (last == nchr) {
                    std::cerr << "Sequence " << labels_[nseq - 1]
                              << " is the wrong length!" << std::endl;
                    exit(1);
                }
                c = line[i];
                if ((c == '0') || (c == '1') || (c == '?')) {
                    splits_[nseq].push_back(c);
                    ++last;
                }
            }
        }
    }
}


void BinaryCharacterMatrix::readNexus(std::istream &is) {
    std::cerr << "No support for NEXUS format; convert to PHYLIP format!"
              << std::endl;
}


void BinaryCharacterMatrix::readMatrix(std::istream &is) {
    std::string buffer;
    getline(is, buffer);
    is.seekg(0, std::ios::beg);

    if (buffer[0] == '#') {
        this->readNexus(is);
    } else if (buffer[0] == '>') {
        this->readFasta(is);
    } else {
        this->readPhylip(is);
	}
}


void BinaryCharacterMatrix::writeNewick(std::ostream &os) {
    std::string str0, str1;
    size_t nseq, nchr;
    int num0, num1;
    char c;

    nseq = labels_.size();
    nchr = splits_[0].size();

    for (size_t j = 0; j < nchr; j++) {
        str0 = "(";
        str1 = "(";
        num0 = 0;
        num1 = 0;
        for (size_t i = 0; i < nseq; i++) {
            c = splits_[i][j];
            switch (c) {
                case '0' :
                    if (num0 == 0) {
                        str0 += labels_[i];
                    } else {
                        str0 += "," + labels_[i];
                    }
                    ++num0;
                    break;
                case '1' :
                    if (num1 == 0) {
                        str1 += labels_[i];
                    } else {
                        str1 += "," + labels_[i];
                    }
                    ++num1;
                default :
                    break;
            }
        }
        if ((num0 > 1) && (num1 > 1)) {
            os << str0 << "," << str1 << "))" << std::endl;
        }
    }
}
