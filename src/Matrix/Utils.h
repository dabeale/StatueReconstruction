#ifndef UTILS_H
#define UTILS_H

/* ----------------------------------------------------------------------
 * Copyright (C) 2016 Daniel Beale. All rights reserved.
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 * ---------------------------------------------------------------------- */

#include <string>
#include <vector>
#include <sstream>

#include "Matrix.h"

namespace Utils
{ // The namespace is necessary to avoid conflicts with another library, so it is also quite short
    void print_matrix(const Math::Matrix&, const std::string&);
    void print_matrix_octave(const Math::Matrix&, const std::string&);
    void print_matrix_octave(const std::vector< Math::Matrix >&A, const std::string& name);
    void print_matrix_octave(const std::vector<std::vector<double>>& A, const std::string& name);
    void print_matrix_octave(const double* data, const uint32_t M, const uint32_t N, const std::string& name, bool binary=false);
    void print_matrix_octave(const uint32_t* data, const uint32_t M, const uint32_t N, const std::string& name, bool binary=false);
    void print_matrix_octave(const uint32_t* data, const uint32_t M, const uint32_t N, const uint32_t K, const std::string& name);
    void print_matrix_octave(const double* data, const uint32_t M, const uint32_t N, const uint32_t K, const std::string& name, bool binary=false);

    typedef std::pair< std::vector<double> , std::vector<uint32_t> > DataDims;
    DataDims read_matrix_octave(const std::string &filename, std::string& name);

    Math::Matrix read_matrix_bundler( const std::string& filename );

}

#endif
