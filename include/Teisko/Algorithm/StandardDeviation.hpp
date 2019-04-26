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
    // Utility class to accumulate variance / standard deviation
    //  - this is stable enough when mean is close to zero
    //  - improvement is to accumulate (item - K), where K is close to mean
    //    or any sample in the set e.g. the first sample
    struct std_s
    {
        double sum = 0;
        double sum2 = 0;
        size_t n = 0;
        std_s() : sum(0), sum2(0), n(0) { }
        std_s& operator+=(const double sample)
        {
            sum += sample;
            sum2 += sample * sample;
            n++;
            return *this;
        }
        // returns unbiased std (as in matlab), scaling by 1/(n-1) or by 1/n (when false)
        double operator() (bool unbiased = true)
        {
            auto bias = unbiased ? 1.0 : 0.0;
            return n < 2
                ? 0
                : std::sqrt((sum2 - (sum*sum) / (double)n) / (n - bias));
        }
        // Returns average (so far)
        double mean() {
            return n == 0 ? 0.0 : sum / (double)n;
        }
    };

    /// \brief calculates unbiased (or biased) standard deviation of a vector
    /// \param vec  Vector to calculate
    /// \param unbiased  true for matlab compliance (sqrt(sigma/(n-1))), false to calculate sqrt(sigma/n)
    /// \returns standard deviation
    inline double vector_std(std::vector<double> &vec, bool unbiased = true)
    {
        std_s variance;
        for (auto &x : vec)
            variance += x;
        return variance(unbiased);
    }
}
