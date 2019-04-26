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
#include <cstdint>
#include <vector>
#include <map>
#include <array>
#include <string>
#include <algorithm>  // std::find

namespace Teisko
{
    /// Averaging functor -- this is core math / statistics
    /// Similar functors exists for root mean calculation or min, max, etc.
    template <typename T>
    struct average_f
    {
        T sum;
        size_t count;

        average_f() : sum(0), count(0) {};

        /// Functor form of operation -- e.g in image.foreach(average_f{})
        void operator() (const T add)
        {
            sum += add;
            count++;
        }

        /// Adds a value to statistics
        average_f& operator += (const T add) {
            operator()(add);
            return *this;
        }

        /// Evaluation of the average is integrated to cast operation
        /// Thus:  auto average_f<int> tmp;
        /// tmp += 3; tmp += 5;
        /// tmp == 4
        operator T() const
        {
            return count == 0 ? 0 : sum / (T)count;
        }
    };

    template <typename T>
    struct weighted_average_f
    {
        double sum;
        double weighted_count;

        weighted_average_f() : sum(0), weighted_count(0) {};

        /// Functor form of operation
        void operator() (const T add, const double w)
        {
            sum += add * w;
            weighted_count += w;
        }

        /// Evaluation of the average is integrated to cast operation
        operator T() const
        {
            return weighted_count == 0 ? 0 : (T)(sum / weighted_count);
        }
    };

    template <typename T>
    struct maximum_f
    {
        T maximum{ 0 };
        /// Functor form of operation -- e.g in image.foreach(average_f{})
        void operator() (const T value)
        {
            maximum = std::max(maximum, value);
        }
        operator T() const
        {
            return maximum;
        }
    };
}