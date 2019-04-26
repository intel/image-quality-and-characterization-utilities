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


#include "Teisko/Image/RGB.hpp"
#include "catch.hpp"

using namespace Teisko;

SCENARIO("Libimage_rgb allows creation of self contained RGB images with access to image as whole and access to color channels")
{
    GIVEN("Image dimensions of 3 x 2")
    {
        roi_point size(3, 2);
        WHEN("An RGB image is created with default layout and color order")
        {
            auto img = rgb_image_s<int>(size);
            AND_WHEN("We initialize the three planes with various methods")
            {
                // Initialize channel zero with constant value
                auto red = img[0].fill(7);

                // Initialize channel one with custom generator functor
                int generator = 0;
                auto green = img[1]
                    .generate([&generator](int &x) { x = ++generator; });

                // Initialize channel two from external buffer
                int buffer[6] = { 10, 20, 30, 40, 50, 60 };
                auto blue = img[2].init_from(buffer);

                THEN("The serialized image as a whole matches a planar or paged memory layout")
                {
                    auto vec = img.full_image().to_vector();
                    auto ref = std::vector<int>({
                        7, 7, 7,
                        7, 7, 7,
                        1, 2, 3,
                        4, 5, 6,
                        10, 20, 30,
                        40, 50, 60
                    });
                    CHECK(vec == ref);
                }
            }
        }
    }
}

SCENARIO("Libimage_rgb allows mapping external contiguous data as RGB images")
{
    GIVEN("External buffers to three planes of RGB data")
    {
        roi_point size(3, 2);
        auto tiled_data = std::vector<int>({
            // BLUE            RED              GREEN
            10, 10, 10, /* */ 20, 20, 20, /* */ 30, 30, 30,
            10, 10, 10, /* */ 20, 20, 20, /* */ 30, 30, 30
        });

        WHEN("An RGB image is created with BRG color order with side-by-side layout")
        {
            auto img = rgb_image_s<int>(size, tiled_data.data(), rgb_layout_e::tiled, rgb_order_e::brg);
            THEN("The individual channels of the data matches the expected constant values")
            {
                auto red = img[rgb_color_e::red].to_vector();
                auto green = img[rgb_color_e::green].to_vector();
                auto blue = img[rgb_color_e::blue].to_vector();

                CHECK(blue == std::vector<int>(3 * 2, 10));      // vector of six tens
                CHECK(red == std::vector<int>(3 * 2, 20));       // vector of six twenties
                CHECK(green == std::vector<int>(3 * 2, 30));     // vector of six thirties
            }
        }
    }

    GIVEN("External buffer to interleaved rgb data with excess data at each row")
    {
        roi_point size(3, 2);
        uint8_t r = 10;
        uint8_t g = 20;
        uint8_t b = 30;
        auto inter_leaved_data = std::vector<uint8_t>({
            r, g, b, r, g, b, r, g, b, 0, 0, 0,
            r, g, b, r, g, b, r, g, b, 0, 0, 0,
        });

        WHEN("An RGB image is created with BRG color order with interleaved scheme and padding")
        {
            auto img = rgb_image_s<uint8_t>(size, inter_leaved_data.data(), rgb_layout_e::interleaved_stride_4, rgb_order_e::brg);
            THEN("The individual channels of the data matches the expected constant values")
            {
                auto red = img[rgb_color_e::red].to_vector();
                auto green = img[rgb_color_e::green].to_vector();
                auto blue = img[rgb_color_e::blue].to_vector();

                CHECK(blue == std::vector<uint8_t>(3 * 2, r));      // vector of six tens
                CHECK(red == std::vector<uint8_t>(3 * 2, g));       // vector of six twenties
                CHECK(green == std::vector<uint8_t>(3 * 2, b));     // vector of six thirties
            }
        }
    }

}

SCENARIO("Error handling -- creating an image with unknown layout throws")
{
    CHECK_THROWS(rgb_image_s<int>(roi_point(1, 1), static_cast<rgb_layout_e>(54)));
}

SCENARIO("Error handling -- creating an image with unknown rgb order throws")
{
    CHECK_THROWS(rgb_image_s<int>(roi_point(1, 1), rgb_layout_e::interleaved, static_cast<rgb_order_e>(7)));
}

SCENARIO("Error handling -- accessing 3rd channel of rgb image throws")
{
    auto img = rgb_image_s<int>(roi_point(1, 1));
    CHECK_NOTHROW(img[2]);
    CHECK_THROWS(img[3]);
}
