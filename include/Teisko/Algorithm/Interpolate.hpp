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

namespace Teisko
{
    /// Interpolate arbitrary length vectors in arbitrary number of dimensions
    /// This belongs more like to a math library
    /// \todo -- decide if the `keys` should be in reverse order or if it should be a list
    ///          for better performance. Any domain specific implementation can and should
    ///          implement the order reversal and/or naming of the keys
    /// \todo -- Add Interpolation method per key -- nearest, linear, bilinear, logarithmic etc.
    ///          Add Extrapolation method per key -- nearest, linear, bilinear, logarithmic etc.
    template <typename T>
    struct interpolate_x_dim
    {
        using type = T;
        std::vector<type> leaf;                    // Final values to interpolate
        std::map<type, interpolate_x_dim> tree;    // X-values to interpolate
    public:
        interpolate_x_dim() { }
        std::vector<type> get(std::vector<type> keys) { return operator() (keys); }
        std::vector<type> operator() (std::vector<type> keys)
        {
            // End of recursion
            if (keys.size() == 0 || tree.size() == 0)
                return leaf;

            // Consume the first key
            // Consuming last key would be more efficient but we want
            // tight control of evaluation order -- first arguments first
            auto key = keys.front();
            keys.erase(keys.begin());

            auto y = tree.begin();
            // If there's only one branch in the tree, just delve in
            if (tree.size() == 1)
                return y->second(keys);

            // There are at least two branches -- binary search the next biggest branch
            auto x = tree.lower_bound(key);
            if (x != tree.end())
            {
                // Item is exactly represented by the <key,value>
                if (x->first == key)
                    return x->second(keys);

                if (x == y)
                    ++y;
                else
                    y = x--;
            }
            else  // the value to be searched is larger than any item in the list -- extrapolate
            {
                y = --x;
                --x;
            }
            if (y == tree.end())
                return leaf;       // should not happen -- but this is to satisfy klocwork

            // recursively get two vectors and interpolate/extrapolate the results
            // Notice that both of these branches must work on own copy of keys
            auto vec1 = x->second(keys);
            auto vec2 = y->second(keys);
            if (vec1.size() != vec2.size() || x->first == y->first)
                return leaf;        // Misconfigured tree (too many dimensions or leaf size conflict)

            auto factor = (key - x->first) / (y->first - x->first);
            for (unsigned int i = 0; i < vec1.size(); i++)
                vec1[i] += (vec2[i] - vec1[i]) * factor;
            return vec1;
        }

        // Recursively adds(or replaces) knee point to tree
        void set(std::vector<type> keys, std::vector<type> &data)
        {
            if (keys.size() == 0)
                leaf = data;
            else
            {
                // Converts first item from the keys vector as knee point
                auto key = keys.front();
                keys.erase(keys.begin());
                tree[key].set(keys, data);
            }
        }
    };

    /// <summary>
    /// Linearly interpolates input data in sampling points as in Matlab.
    /// x_data must be strictly monotonically increasing.
    /// Values outside of input range are set to zero.
    /// </summary>
    template <typename T>
    inline std::vector<T> interp1d(
        const std::vector<T> &x_data,
        const std::vector<T> &y_data,
        const std::vector<T> &sampling_points)
    {
        auto current_iterator = x_data.begin();
        auto interp = [&current_iterator, &x_data, &y_data](T sampling_point)
        {
            // Find if sampling point is out of range of the input data
            if (x_data.front() > sampling_point || x_data.back() < sampling_point)
                return 0.0f;
            if (sampling_point == x_data.back())
                return y_data.back();

            while (*current_iterator <= sampling_point)
                current_iterator++;

            auto lower_iterator = current_iterator - 1;
            auto upper_index = current_iterator - x_data.begin();
            auto lower_index = upper_index - 1;

            // (1-t)*v0 + t*v1
            auto t = (sampling_point - *lower_iterator) / (*current_iterator - *lower_iterator);
            return (1 - t) * y_data[lower_index] + t * y_data[upper_index];
        };

        std::vector<T> output;
        output.reserve(sampling_points.size());
        for (auto point : sampling_points)
            output.push_back(interp(point));

        return output;
    }
}
