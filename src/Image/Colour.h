#ifndef COLOUR_H
#define COLOUR_H

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

#include <cmath>
#include <vector>
namespace Colour
{
    /**
     * @brief rgb2xyz
     * Convert rgb [0-1] to xyz. This implementation is from Octave.,
     * who take it from Burger, Burge "Digitale Bildverarbeitung", 3rd edition  (2015).
     */
    void rgb2xyz( double r, double g, double b, double& x, double& y, double& z );
    void xyz2rgb( double x, double y, double z, double& r, double& g, double& b );
    void xyz2lab( double x, double y, double z, double& l, double& a, double& b );
    void lab2xyz( double l, double a, double b , double& x, double& y, double& z);

    void rgb2lab(double r, double g, double b, double& l, double& a, double& bp );
    void lab2rgb(double l, double a, double bp, double& r, double& g, double& b );
}

#endif
