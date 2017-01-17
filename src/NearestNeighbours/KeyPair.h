#ifndef KEYPAIR_H
#define KEYPAIR_H

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

#include <stdint.h>
#include <limits>

namespace KD
{

/**
 * @brief The KeyPair struct
 * A struct for the elements of a KD tree. It contains the index in to the
 * data matrix and the distance from the test point.
 */
struct KeyPair
{
public:
    /**
     * @brief KeyPair
     * Construct an empty KeyPair
     */
    KeyPair() : distance(std::numeric_limits<double>::max()), index(0) {}

    /**
     * @brief KeyPair
     * Construct a KeyPair with a distance and an index.
     * @param d The distance
     * @param ind The index
     */
    KeyPair( double d , uint32_t ind) : distance(d), index(ind) {}

    double distance;  ///< A distance
    uint32_t index;   ///< An index
};
}
#endif // KEYPAIR_H

