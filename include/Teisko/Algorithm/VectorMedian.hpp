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
    /// \brief Computes the median of the given *values*
    ///
    /// The median is calculated without changing the given input
    /// vector. The median is a value for which half of the data
    /// is bigger than it and half of the data is smaller than it.
    ///
    /// \param values Input values
    ///
    /// \returns The median value of the given *values*
    ///
    /// \throws std::runtime_error is thrown if *values* is empty
    ///
    /// See also: \ref libcalc_specs_vector_median
    template<typename T> T vector_median(std::vector<T> const & values)
    {
        if (values.empty())
        {
            throw std::runtime_error("Teisko::vector_median - the input vector is empty");
        }

        // Copy the input values to a temporary variable as the
        // applied algorithms modify their input data
        std::vector<T> temp(values.begin(), values.end());

        // Partial sort using the "middle index" as pivot point
        auto middle = temp.size() / 2;
        std::nth_element(temp.begin(), temp.begin() + middle, temp.end());

        T median = temp[middle];

        if (middle * 2 == temp.size())
        {
            // If there are even number of values, find the maximum of
            // the lower half of the values. Please note that due to
            // the partial sorting we know that the *median* value is
            // the next element after the lower half and therefore
            // *lower_half* <= *median*.
            //
            // For decimal values this returns the arithmetic average
            // of the *lower_half* and *median*.
            //
            // For integral values this returns the truncated average
            // of the *lower_half* and *median* where the reminder is
            // rounded towards zero.

            T lower_half = *std::max_element(temp.begin(), temp.begin() + middle);
            return (lower_half + (median - lower_half) / 2);
        }
        else
        {
            return median;
        }
    }
}
