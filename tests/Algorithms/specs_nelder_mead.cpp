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

#include "Teisko/Algorithm/NelderMead.hpp"
#include "catch.hpp"


/// \page libcalc_specs_convex_hull Specs: Teisko convex_hull
///
/// \snippet this snippet-specs-convex_hull

using namespace Teisko;

namespace libcalc_nelder_mead_tests
{
    SCENARIO("Nelder Mead simplex converges to Rosenbrock banana local minimum")
    {
        GIVEN("The banana function")
        {
            auto banana = [](double *p) {
                return 100 * (p[1] - p[0] * p[0])*(p[1] - p[0] * p[0]) + (1.0 - p[0]) * (1.0 - p[0]);
            };
            WHEN("this function is minimized starting from [-1.2, 1.0]")
            {
                auto result = nelder_mead_simplex({ -1.2, 1.0 }, banana, { 1500 });
                THEN("The resulting vector is very close to [1.0, 1.0]")
                {
                    CHECK(result.second[0] == Approx(1.0).epsilon(1e-4));
                    CHECK(result.second[1] == Approx(1.0).epsilon(1e-4));
                }
            }
        }
    }

    SCENARIO("Nelder Mead simplex finds average locus")
    {
        GIVEN("A custom function of single variable and a set of points to optimize to")
        {
            auto func = [](double x, double *p) { return (p[0] + p[3] * pow(x - p[1], p[2])); };
            auto xy = std::vector<point_xy>{
                point_xy{ 0.5945, 0.5925 },
                point_xy{ 0.9367, 0.3937 },
                point_xy{ 0.4340, 0.8805 },
                point_xy{ 0.4270, 0.8750 },
                point_xy{ 0.8430, 0.4113 }
            };
            xy_func_optimizer locus(func, xy);

            WHEN("The optimizer class is processed with Nelder Mead optimizer")
            {
                //  for (int k = 0; k < 9999;k++)
                //      nelder_mead_simplex({ 0, 0, -1, 1 }, locus, { 800, 1e-10, 1e-10, 800 });
                auto vec = nelder_mead_simplex({ 0, 0, -1, 1 }, locus, { 800, 1e-10, 1e-10, 800 });
                THEN("The result quality and local minimum point match to Matlab with high accuracy")
                {
                    double reference_quality = 4.070174179481758e-04;
                    double ref_params[4] = { 0.259456533522921, -0.185094026577870, -2.632094144573014, 0.172299346450623 };

                    CHECK(vec.first == Approx(reference_quality).epsilon(1e-6));
                    REQUIRE(vec.second.size() == 4);
                    for (auto i = 0; i < 4; i++)
                    {
                        CHECK(vec.second[i] == Approx(ref_params[i]).epsilon(1e-6));
                    }
                }
            }
        }
    }
}