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
 * @file binary_character_matrix.hpp
 * @brief Binary Character Matrix Class
 */
#ifndef _PHYLOTOOLS_SPLITS_HPP
#define _PHYLOTOOLS_SPLITS_HPP


#include <iostream>
#include <stdint.h>
#include <sstream>
#include <string>
#include <vector>


namespace phylotools {


/**
 *
 */
class BinaryCharacterMatrix {
 public:
    /**
     * Splits Constructor
     */
    BinaryCharacterMatrix() {}

    /**
     * Splits Destructor
     */
    ~BinaryCharacterMatrix() {}

    void readPhylip(std::istream &is);
    void readFasta(std::istream &is);
    void readNexus(std::istream &is);
    void readMatrix(std::istream &is);
    void writeNewick(std::ostream &os);

 protected:
    size_t countColumns(std::istream &is);  /**< */
    std::vector<std::string> splits_;       /**< */
    std::vector<std::string> labels_;       /**< */
};


}  // namespace phylotools


#endif  // _PHYLOTOOLS_SPLITS_HPP
