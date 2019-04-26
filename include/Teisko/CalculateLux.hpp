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


#include <algorithm>
#include <exception>
#include <cstdint>
#include <cmath>
#include <limits>
#include <vector>
#include <tuple>

namespace Teisko
{
    inline int pseudo_inverse(double const* A, double* I)
    {
        // For 2x2 matrix the inverse can be calculated analytically
        //
        // [a  b] -1       1     [ d  -b]
        // [    ]    = --------- [      ]
        // [c  d]       ad - bc  [-c   a]

        // This implementation is sensitive to matrices which determinant
        // is very small number.
        auto determinant = A[0] * A[3] - A[1] * A[2];

        if (determinant == 0.0) return 1;

        I[0] =  A[3] / determinant;
        I[1] = -A[1] / determinant;
        I[2] = -A[2] / determinant;
        I[3] =  A[0] / determinant;

        return 0;
    }

    inline void matrix_multiply(double const* I, double const* b, double* x)
    {
        // I (2x2) * b (2x1) = x (2x1)
        //
        // [ I0, I1 ]   [ b0 ]   [ I0*b0 + I1*b1 ]   [ x0 ]
        // [        ] * [    ] = [               ] = [    ]
        // [ I2, I3 ]   [ b1 ]   [ I2*b0 + I3*b1 ]   [ x1 ]

        x[0] = I[0] * b[0] + I[1] * b[1];
        x[1] = I[2] * b[0] + I[3] * b[1];
    }

    // The code values are set to same as they were in LIBLAB implementation
    const int32_t CALCULATELUX_SUCCESS = 0;
    const int32_t CALCULATELUX_ERROR_INVALID_PARAMETER = -6;
    const int32_t CALCULATELUX_ERROR_INTERNAL_FAILURE = -12;

    inline int32_t calculate_lux(
        const double* exposure, const double* total_gain, const uint32_t input_arrays_size,
        const double base_iso, const double low_limit_pct, const double high_limit_pct,
        const double low_limit_lux, const double high_limit_lux, const double dark_adaptation_factor,
        double* lux_average_array, double* lux_lower_array, double* lux_upper_array)
    {
        if (exposure == nullptr || total_gain == nullptr || lux_average_array == nullptr
            || lux_upper_array == nullptr || lux_lower_array == nullptr)
            return CALCULATELUX_ERROR_INVALID_PARAMETER;

        if (input_arrays_size == 0)
            return CALCULATELUX_ERROR_INVALID_PARAMETER;

        if (base_iso <= 0)
            return CALCULATELUX_ERROR_INVALID_PARAMETER;

        // Create sorted total gain array
        std::vector<double> total_gain_sorted;
        total_gain_sorted.assign(total_gain, total_gain + input_arrays_size);
        std::sort(total_gain_sorted.begin(), total_gain_sorted.end());

        // Constant 5486000 is needed because the exposure time is in milliseconds.
        // Target is the desired normalized average raw level for lux estimation.
        // The sqrt(2) factor removes the headroom from the Base ISO which is
        // defined in the saturation.This equation is based on the Base ISO in ISO
        // standard.
        const int constant = 5486000;
        double target = (low_limit_pct + high_limit_pct) / 2;

        for (uint32_t i = 0; i < input_arrays_size; ++i)
        {
            lux_average_array[i] = constant / (exposure[i] * total_gain_sorted[i] * base_iso) * target * sqrt(2);

            if (lux_average_array[i] < low_limit_lux)
            {
                /// \todo Update calculate_lux specs to test this branch
                /// if (lux_average_array[i] < low_limit_lux) is true
                lux_average_array[i] *= dark_adaptation_factor;
            }
            else if (lux_average_array[i] < high_limit_lux)
            {
                // Matlab reference for dark adaptation lux factor
                /* A = [highLimitLux 1; lowLimitLux 1];
                b = [1; darkAdaptorFactor];
                x = A\b;
                c = [average(k) 1] * x; % dark adaption lux factor
                average(k) = average(k) * c; */

                double A[] = { high_limit_lux, 1, low_limit_lux, 1 };
                double I[] = { 1, 0, 0, 1 };
                auto res = Teisko::pseudo_inverse(A, I);

                if (res != 0) return CALCULATELUX_ERROR_INTERNAL_FAILURE;

                const int height = 2;
                double b[] = { 1, dark_adaptation_factor };
                double x[height];
                Teisko::matrix_multiply(I, b, x);

                double c = lux_average_array[i] * x[0] + x[1];

                lux_average_array[i] *= c;
                lux_lower_array[i] = lux_average_array[i] / target * low_limit_pct;
                lux_upper_array[i] = lux_average_array[i] / target * high_limit_pct;
            }

            // In case multiple exposures, first values are always inf
            if (input_arrays_size > 1)
            {
                lux_average_array[0] = std::numeric_limits<double>::infinity();
                lux_lower_array[0] = std::numeric_limits<double>::infinity();
                lux_upper_array[0] = std::numeric_limits<double>::infinity();
            }

            // Round values
            lux_average_array[i] = round(lux_average_array[i]);
            lux_lower_array[i] = round(lux_average_array[i] / target * low_limit_pct);
            lux_upper_array[i] = round(lux_average_array[i] / target * high_limit_pct);
        }

        return CALCULATELUX_SUCCESS;
    }
}
