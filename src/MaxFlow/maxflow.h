#ifndef MAXFLOW_H
#define MAXFLOW_H

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

#include <vector>
#include <cmath>
#include <array>

#include "Matrix.h"
#include "graph.h"
#include "Message.h"

namespace MaxFlow
{

    struct HyperParameters
    {
        std::vector<double> Alpha; ///< Parameter of the dirichlet prior
        uint32_t N; ///< Draw N samples for parameters estimation
        double edgeWeight;
        double terminalWeight;
        Stream::Message msg;
    };

    typedef std::array<double, 2> Arr2d;

    inline double LogMultinomial( const Arr2d& histogram, const Arr2d& parameter )
    {
        double p1 = histogram[0]*std::log(parameter[0]);
        double p2 = histogram[1]*std::log(parameter[1]);
        return std::lgamma( histogram[0] + histogram[1] + 1 ) - std::lgamma(histogram[0] + 1) - std::lgamma(histogram[1] + 1) + p1 + p2;
    }

    Math::Matrix EvaluateFGBGFromProbabilities( const Math::Matrix& fgprobs, const Math::Matrix& bgprobs, HyperParameters& hypes );
}


#endif // MAXFLOW_H
