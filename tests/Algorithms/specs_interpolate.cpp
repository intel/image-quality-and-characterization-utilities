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


#include "Teisko/Algorithm/Interpolate.hpp"
#include "catch.hpp"
#include <memory>
#include <random>
#include <algorithm>
#include <vector>

using namespace Teisko;

namespace teisto_interpolate_tests
{

    SCENARIO("Interpolation in two dimensions")
    {
        GIVEN("A few values in a uniform 2x2 grid")
        {
            std::vector<float> A = { 4 };               //       y1   4     5
            std::vector<float> B = { 5 };               //       y2   6 (7) 8
            std::vector<float> C = { 6 };               //        |  (7)
            std::vector<float> D = { 8 };               //        +--x1--- x2----x3---> x

            // key/knee points in horizontal axis
            const float x1 = 1.0f;
            const float x2 = 2.0f;
            // key points in vertical axis
            const float y1 = 3.0f;
            const float y2 = 2.0f;

            WHEN("The interpolator is setup with fixed data")
            {
                auto foo = interpolate_x_dim<float>{};
                foo.set({ x1, y1 }, A);
                foo.set({ x2, y1 }, B);
                foo.set({ x1, y2 }, C);
                foo.set({ x2, y2 }, D);
                {
                    THEN("The return value at each grid point equals the configured value")
                    {
                        CHECK(A == foo.get({ x1, y1 }));
                        CHECK(B == foo.get({ x2, y1 }));
                        CHECK(C == foo.get({ x1, y2 }));
                        CHECK(D == foo.get({ x2, y2 }));
                    }
                    AND_THEN("We can interpolate the middle value along the bottom row")
                    {
                        std::vector<float> middle_of_C_AND_D = { (6.0 + 8.0) / 2 };  // reference value
                        CHECK(middle_of_C_AND_D == foo({ (x1 + x2) / 2, y2 }));
                    }
                    AND_THEN("We can extrapolate by 50% along the left column")
                    {
                        std::vector<float> outside_of_A_AND_C = { 7.0 };
                        CHECK(outside_of_A_AND_C == foo({ x1, y2 + 0.5f * (y2 - y1) }));
                    }
                    AND_THEN("The middle point returned by the interpolator should be average of all A,B,C and D")
                    {
                        std::vector<float> average = { (4 + 5 + 6 + 8) / 4.0f };
                        CHECK(average == foo({ (x1 + x2) * 0.5f, (y1 + y2) * 0.5f }));
                    }
                }
            }
        }
    }

    SCENARIO("Interpolation in Matlab style", "[interp1d]")
    {
        GIVEN("A few sample values and sampling points")
        {
            std::vector<float> x = { 1.0f, 2.0f, 3.0f, 4.0f };
            std::vector<float> y = { 1.0f, 2.0f, 3.0f, 4.0f };

            std::vector<float> sampling_points = { 0.0f, 1.5f, 3.0f, 4.0f, 4.5f };

            WHEN("data is interpolated")
            {
                auto sampling_values = interp1d(x, y, sampling_points);

                THEN("data outside input range are zeros")
                {
                    CHECK(sampling_values[0] == 0.0f);
                    CHECK(sampling_values[4] == 0.0f);
                }
                AND_THEN("data at input points are same then input values")
                {
                    CHECK(sampling_values[2] == x[2]);
                    CHECK(sampling_values[3] == x[3]);
                }
                AND_THEN("sampling point at interpolated location matches calculated one")
                {
                    CHECK(sampling_values[1] == (x[1] + x[0]) / 2);
                }
            }
        }
    }
}