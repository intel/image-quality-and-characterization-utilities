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

#include "Teisko/Algorithm/LinearSpace.hpp"
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
    /// \brief Computes the mean of the given *data* excluding highest and lowest samples
    /// The mean is calculated *with* modifying the data
    /// \param  data    Data whose average without outliers is to be estimated
    /// \param  trim    Proportion of good samples -- (1-trim)/2 items to be excluded on both sides
    /// \returns    Average of data without outliers
    template <typename T>
    double trimmed_mean(std::vector<T> &data, double trim)
    {
        auto size = data.size();
        auto first_idx = std::min(size, static_cast<size_t>(std::max(0.0, 0.5 + 0.5 * (1.0 - trim) * size)));
        auto last_idx = size - first_idx;
        if (first_idx >= last_idx)
            return static_cast<T>(0);

        std::sort(data.begin(), data.end());

        decltype (T(1)*T(1)) sum = 0;
        auto divisor = static_cast<decltype(sum)>(last_idx - first_idx);
        for (size_t i = first_idx; i < last_idx; i++)
            sum += data[i];
        return (double)sum / (double)divisor;
    }

    /// \brief trimmed mean calculated from a histogram
    /// \param histogram    count of each sample
    /// \param skip         number of items to skip
    /// \param count        number of items to count
    /// \param centers      optional center (or weight) of each histogram bin
    template <typename T, typename U>
    double trimmed_mean(std::vector<T> &histogram, size_t skip, size_t count, std::vector<U> &centers)
    {
        auto h_size = histogram.size();
        decltype(h_size) i = 0;
        if (h_size == 0)
            return 0.0;

        if (centers.size() > 0 && centers.size() != h_size)
            throw std::runtime_error("Expecting empty, or equal sized weight vector");

        // make a weight vector conforming to 0,1,2,3...N-1 or use the provided one
        auto w = centers.size() == 0
            ? linear_space(0, (double)(h_size - 1), (unsigned int)h_size)
            : std::vector<double>(centers.begin(), centers.end());

        auto h = static_cast<size_t>(histogram[i]);      // number of items to be counted in current bin
        auto n = 0.0;
        auto sum = 0.0;

        // Ignore 'skip' entries
        while (h <= skip)
        {
            skip -= h;
            // Make sure we don't index out of range, because skip >= sum(histogram)
            if (++i >= h_size)
                return 0.0;
            h = static_cast<size_t>(histogram[i]);
        }
        h -= skip;

        while (h <= count)
        {
            count -= h;
            n += h;
            sum += w[i] * h;
            if (++i >= h_size)
            {
                i = 0;
                count = 0;
                break;
            }
            h = static_cast<size_t>(histogram[i]);
        }
        sum += count * w[i];
        n += count;
        return n > 0 ? sum / n : 0.0;
    }
}
