/* Copyright 2019 Erin K. Molloy
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors
 * may be used to endorse or promote products derived from this software without
 * specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
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
