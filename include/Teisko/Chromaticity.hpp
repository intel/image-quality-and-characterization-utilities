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
#include "Teisko/Algorithm/Functors.hpp"
#include <cstdint>
#include <vector>
#include <map>
#include <array>
#include <string>
#include <algorithm>  // std::find

namespace Teisko
{
    /// \brief  Container of chromaticity information
    struct chromaticity
    {
        double _r_per_g;
        double _b_per_g;
        double _i_per_g;

        /// \brief  Constructs chromaticity point from three scalars
        chromaticity(double r = 0.0, double b = 0.0, double i = 0.0)
            : _r_per_g(r)
            , _b_per_g(b)
            , _i_per_g(i) { }

        /// \brief  Returns squared euclidian distance of chromaticity points
        double dist2(const chromaticity &other) const
        {
            double x = _r_per_g - other._r_per_g;
            double y = _b_per_g - other._b_per_g;
            double z = _i_per_g - other._i_per_g;
            return x*x + y*y + z*z;
        }

        bool operator ==(const chromaticity &other) { return 0.0 == dist2(other); }
    };

    inline chromaticity& operator +=(chromaticity &a, const chromaticity &b)
    {
        a._r_per_g += b._r_per_g;
        a._b_per_g += b._b_per_g;
        a._i_per_g += b._i_per_g;
        return a;
    }

    inline chromaticity& operator /=(chromaticity &a, int i)
    {
        if (i != 0)
        {
            a._r_per_g /= (double)i;
            a._b_per_g /= (double)i;
            a._i_per_g /= (double)i;
        }
        return a;
    }

    /// \brief chromaticity_factory_f   Calculates chromaticity incrementally
    struct chromaticity_factory_f
    {
        average_f<double> _channels[4]{};  // R, G, B, I

        /// \brief  Returns average functor by color index
        average_f<double>& operator[](color_info_e channel)
        {
            return _channels[static_cast<int>(channel)];
        }

        /// \brief  returns chromaticity from statistics
        operator chromaticity()
        {
            double g_inv = _channels[(int)color_info_e::green];
            if (g_inv != 0.0)
                g_inv = 1.0 / g_inv;
            double r = _channels[(int)color_info_e::red];
            double b = _channels[(int)color_info_e::blue];
            double i = _channels[(int)color_info_e::ir];
            return chromaticity(r * g_inv, b * g_inv, i * g_inv);
        }
    };
}