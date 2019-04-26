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

#include "Teisko/Algorithm/ConvexHull.hpp"
#include "Teisko/Algorithm/PointXY.hpp"
#include "catch.hpp"
#include <random>

/// \page libcalc_specs_convex_hull Specs: Teisko convex_hull
///
/// \snippet this snippet-specs-convex_hull

using namespace Teisko;

namespace libcalc_convex_hull_tests
{
    struct xy
    {
        double x;
        double y;
    };

    bool are_equal(std::vector<xy> container_a, std::vector<xy> container_b)
    {
        auto is_item_equal = [](const xy &a, const xy &b) {
            return a.x == b.x && a.y == b.y;
        };

        return std::equal(
            container_a.begin(), container_a.end(),
            container_b.begin(),
            is_item_equal);
    };

    bool in_range(size_t value, size_t min_value, size_t max_value)
    {
        return value >= min_value && value <= max_value;
    };

    SCENARIO("Convex hull test case")
    {
        GIVEN("A predefined 2D input")
        {
            // Use input validated in Matlab to produce precalculated output
            auto input = std::vector<xy>{
                { 0.751267059305653, 0.351659507062997 },
                { 0.255095115459269, 0.830828627896291 },
                { 0.505957051665142, 0.585264091152724 },
                { 0.699076722656686, 0.549723608291140 },
                { 0.890903252535799, 0.917193663829810 },
                { 0.959291425205444, 0.285839018820374 },
                { 0.547215529963803, 0.757200229110721 },
                { 0.138624442828679, 0.753729094278495 },
                { 0.149294005559057, 0.380445846975357 },
                { 0.257508254123736, 0.567821640725221 },
                { 0.840717255983663, 0.075854289563064 },
                { 0.254282178971531, 0.053950118666607 },
                { 0.814284826068816, 0.530797553008973 },
                { 0.243524968724989, 0.779167230102011 },
                { 0.929263623187228, 0.934010684229183 },
                { 0.349983765984809, 0.129906208473730 },
                { 0.196595250431208, 0.568823660872193 },
                { 0.251083857976031, 0.469390641058206 },
                { 0.616044676146639, 0.011902069501241 },
                { 0.473288848902729, 0.337122644398882 }
            };

            WHEN("we calculate convex hull from input data")
            {
                THEN("calculated convex hull should match predefined convex hull")
                {
                    auto convex_hull = std::vector<xy>{
                        { 0.138624442828679, 0.753729094278495 },
                        { 0.255095115459269, 0.830828627896291 },
                        { 0.929263623187228, 0.934010684229183 },
                        { 0.959291425205444, 0.285839018820374 },
                        { 0.840717255983663, 0.075854289563064 },
                        { 0.616044676146639, 0.011902069501241 },
                        { 0.254282178971531, 0.053950118666607 },
                        { 0.149294005559057, 0.380445846975357 }
                    };

                    auto output = convex_hull::hull(input);
                    CHECK(are_equal(output, convex_hull));

                    AND_THEN("convex hull can be pruned to smaller feature set")
                    {
                        convex_hull::prune_hull(output, 4, 6, 1e-4);
                        CHECK(in_range(output.size(), 4, 6));
                    }
                }
            }
        }
    }

    SCENARIO("Convex hull calculation of a known convex hull produces an exact match")
    {
        GIVEN("A polygon of 8 points with all points residing in the convex hull")
        {
            auto input = std::vector<point_xy>{
                { 0.611489559164733, 0.469958236658933 }, // smallest x, (smallest y in case of tie)
                { 0.611489559164733, 0.485958236658933 },
                { 0.612947630922693, 0.486802992518703 },
                { 0.628947630922693, 0.486802992518703 },
                { 0.634566416040100, 0.486696741854637 },
                { 0.634566416040100, 0.470696741854637 },
                { 0.628758483033932, 0.469045908183633 },
                { 0.612758483033932, 0.469045908183633 }
            };

            WHEN("The input polygon is randomized")
            {
                auto shuffled = input;
                std::random_device rd;
                std::mt19937 g(rd());
                std::shuffle(shuffled.begin(), shuffled.end(), g);

                auto output = convex_hull::hull(shuffled);
                THEN("The output matches the input polygon which was already in proper order")
                {
                    CHECK(output == input);
                }
            }
        }
    }
}