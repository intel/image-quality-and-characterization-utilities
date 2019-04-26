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


#include "Teisko/ValidImageArea.hpp"
#include "catch.hpp"
#include <cstdint>
#include <map>
#include <string>
#include <vector>

namespace valid_image_area_tests
{
    using namespace Teisko;
    SCENARIO(
        "Valid Image Area output grid is valid when no characterization images are given. "
        "[cmc_32]"
        )
    {
        GIVEN("Positive grid dimension")
        {
            const int width = 5;
            const int height = 3;
            WHEN("An ValidImageArea is instantiated with positive dimensions")
            {
                auto via = ValidImageArea(width, height);

                THEN("The class contains a non-zero (i.e. all valid) grid of 5*3 elements")
                {
                    auto ref = std::vector<uint16_t>(width * height, 1);
                    CHECK(via.data == ref);
                }
            }
        }
    }

    SCENARIO("Valid Image Area finds all valid grid from synthetic paraboloid data")
    {
        GIVEN("A largish image of VGA resolution")
        {
            auto center = roi_point(640 / 2, 480 / 2);
            auto image = Teisko::image<uint16_t>(center * 2);
            for (auto xy : image.size())
            {
                auto pt = xy - center;
                // paraboloid with max value of 6000, min value of 6000-(320*320 + 240*240 + 1000) / 32  = 1000
                image(xy) = static_cast<uint16_t>(6000 - ((pt._x * pt._x + pt._y * pt._y) / 32));
            }

            auto bayer = bayer_image_s<uint16_t>(image, bayer_pattern_e::bgir);

            const int width = 40;
            const int height = 40;

            WHEN("The image is recognized for valid image area")
            {
                auto via = ValidImageArea(width, height);
                via.Add(bayer);
                THEN("The class contains a valid image area of 40x40 pixels")
                {
                    auto ref = std::vector<uint16_t>(width * height, 1);
                    CHECK(via.data == ref);
                }
            }

            AND_WHEN("The image is recognized for valid image area with a 2000 unit black level removed")
            {
                auto via = ValidImageArea(width, height);

                image.foreach([](uint16_t &x) { x = x > 2000 ? x - 2000 : 0; });

                via.Add(bayer);
                THEN("The class contains a valid image area of 40x40 pixels with corners detected as invalid and center detected as valid")
                {
                    auto grid = via.get_grid();
                    CHECK(grid.size() == roi_point(width, height));
                    CHECK(grid(0, 0) == 0);
                    CHECK(grid(0, width - 1) == 0);
                    CHECK(grid(height - 1, 0) == 0);
                    CHECK(grid(height - 1, width - 1) == 0);
                    CHECK(grid(height / 2, width / 2) == 1);
                }
            }
        }
    }
}
