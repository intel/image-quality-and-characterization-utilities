/*
* Copyright (c) 2019, Intel Corporation
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*     * Redistributions of source code must retain the above copyright
*       notice, this list of conditions and the following disclaimer.
*     * Redistributions in binary form must reproduce the above copyright
*       notice, this list of conditions and the following disclaimer in the
*       documentation and/or other materials provided with the distribution.
*     * Neither the name of the Intel Corporation nor the
*       names of its contributors may be used to endorse or promote products
*       derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
* DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/


#pragma once

// In alphabetical order
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <stdexcept>
#include <vector>

/// \brief This namespace defines the public interfaces of
/// the \ref libcalc_module module
namespace Teisko
{
    /// \brief Generates *n* points between *a* and *b*.
    ///
    /// The *n* points will be evenly distributed and will include
    /// the starting point *a* and the ending point *b*. If *a < b*
    /// then the points will be increasing in value. If *a > b* then
    /// the points will be decreasing in value. If *a == b* then the
    /// value is repeated *n* times.
    ///
    /// There are two special cases. If *n = 0* then the vector will
    /// be empty. If *n = 1* then the vector will contain just the
    /// first parameter *a*.
    ///
    /// \param a Starting point
    /// \param b Ending point
    /// \param n Number of points
    ///
    /// \returns Vector of *n* linearly spaced points in range [*a*, *b*]
    ///
    /// See also: \ref libcalc_specs_linear_space
    inline std::vector<double> linear_space(double a, double b, unsigned int n)
    {
        std::vector<double> points(n, a);

        if (n <= 1 || std::fabs(a - b) < std::numeric_limits<double>::epsilon())
        {
            return points;
        }

        // Remember to skip the last point as it is specified to be
        // exactly the same as the parameter *b*.
        for (decltype(n) k = 0; k < n - 1; ++k)
        {
            points[k] = a + ((b - a) * k) / (n - 1);
        }

        points[n - 1] = b;

        return points;
    }
}
