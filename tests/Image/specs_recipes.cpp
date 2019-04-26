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

#include "Teisko/Image/API.hpp"
#include "catch.hpp"

// This is the cookbook for libimage recipes
using namespace Teisko;

SCENARIO("In the first example we allocate one image containing integers (and destroy it)")
{
    auto const rows = 15;
    auto const columns = 20;
    auto size = roi_point(columns, rows);       // this class allows "x", "y" order of parameters

    CHECK_NOTHROW(image<int>(size));
}

SCENARIO("In the second example we allocate one image containing floats (and destroy it)")
{
    auto const rows = 15;
    auto const columns = 20;

    CHECK_NOTHROW(image<float>(rows, columns)); // this matches the C style: array[rows][columns]
}

SCENARIO(
    "In this example we allocate an image of fixed size, then initialize the image data to all zeros. "
    "Then we show that it is all zeros by converting the image to a vector.")
{
    GIVEN("An image initialized to all zeros")
    {
        auto my_next_image = image<float>(10, 8).fill(0);
        WHEN("We convert the image to a vector [of floats]")
        {
            auto my_image_as_vector = my_next_image.to_vector();

            THEN("The vector matches a newly made vector of 80 zeros")
            {
                auto vector_of_80_zeros = std::vector<float>(10 * 8);
                CHECK(vector_of_80_zeros == my_image_as_vector);
            }
        }
    }
}

SCENARIO("Converting a vector to an image")
{
    GIVEN("A vector of ten values (all ones)")
    {
        auto initial_vector = std::vector<int>(10, 1);
        WHEN("We initialize an image of 5x2 pixels from the vector -- as a view --")
        {
            auto my_image = image<int>(2, 5, initial_vector.data());
            AND_WHEN("We convert that image back to another vector")
            {
                auto my_image_as_vector = my_image.to_vector();
                THEN("These two vectors are identical")
                {
                    CHECK(initial_vector == my_image_as_vector);
                }
            }
        }
    }
}

SCENARIO("Copying a vector to image")
{
    GIVEN("A vector of values")
    {
        auto initial_vector = std::vector<int>(10, 1);      // just ten items (of all ones)
        WHEN("We initialize an image from a vector")
        {
            auto my_image = image<int>(2, 5).init_from(initial_vector.data());
            AND_WHEN("We convert that image back to another vector")
            {
                auto my_image_as_vector = my_image.to_vector();
                THEN("These two vectors are identical")
                {
                    CHECK(initial_vector == my_image_as_vector);
                }
            }
        }
    }
}

/// Internal test case to check that images of various width and heights are iterated w*h times
/// including zero width/height cases and that all and not other pixels in the image are visited
template <int W, int H>
void test_case_for_loop(int random_init_value)
{
    auto count = W * H;
    auto img = image<int>(H, W).init(random_init_value);
    auto expected_result_after_increment = std::vector<int>(count, random_init_value + 1);

    for (auto &p : img)
    {
        REQUIRE(count-- > 0);
        REQUIRE(p == random_init_value);
        p++;
    }
    CHECK(count == 0);

    // And we check that the pixels `p` visited and incremented are actually from the image
    CHECK(img.to_vector() == expected_result_after_increment);
}

SCENARIO("Checking that we can create images of various lengths and that our iterators return ")
{
    test_case_for_loop<0, 0>(1);
    test_case_for_loop<0, 1>(2);
    test_case_for_loop<1, 0>(3);
    test_case_for_loop<1, 1>(4);

    test_case_for_loop<3, 4>(5);
    test_case_for_loop<4, 4>(6);
    test_case_for_loop<4, 3>(7);

}

SCENARIO("Looping over image pixels with ranged for")
{
    GIVEN("An image of size 4x3 initialized to all zeros")
    {
        const int width = 4;
        const int height = 3;
        auto my_image = image<int>(height, width).fill(0);
        WHEN("The pixel data is iterated one by one adding to each pixel")
        {
            for (auto &pixel : my_image)
                ++pixel;

            THEN("The image contains all ones")
            {
                auto all_ones_vector = std::vector<int>(width * height, 1);
                CHECK(my_image.to_vector() == all_ones_vector);
            }
        }
    }
}

SCENARIO("Looping over flipped, mirrored and/or subsampled images with ranged for")
{
    GIVEN("An image of size 4x3 initialized with values from 0 to 11")
    {
        const int width = 4;
        const int height = 3;
        auto data = std::vector<int>({
            0, 1, 2, 3,
            4, 5, 6, 7,
            8, 9, 10, 11
        });

        auto my_image = image<int>(height, width, data.data());

        WHEN("The mirrored (left to right) image is iterated with for loop")
        {
            auto result = std::vector<int>();

            for (auto &pix : my_image.mirror())
                result.emplace_back(pix);

            THEN("The result confirms to order from right to left")
            {
                auto data_mirrorer = std::vector<int>({
                    3,2,1,0,
                    7,6,5,4,
                    11,10,9,8
                });
                CHECK(result == data_mirrorer);
            }
        }

        WHEN("The flipped image is iterated with for loop")
        {
            auto result = std::vector<int>();

            for (auto &pix : my_image.flip())
                result.emplace_back(pix);

            THEN("The result confirms to order from bottom to top")
            {
                auto data_flipped = std::vector<int>({
                    8, 9, 10, 11,
                    4, 5, 6, 7,
                    0, 1, 2, 3
                });
                CHECK(result == data_flipped);
            }
        }

        WHEN("The image is subsampled by a factor of 2x2 and iterated")
        {
            auto result = std::vector<int>();

            for (auto &pix : my_image.subview(2,2))
                result.emplace_back(pix);

            THEN("The result has gathered the points 0,2,8,10")
            {
                auto data_subsampled = std::vector<int>({ 0, 2, 8, 10 });
                CHECK(result == data_subsampled);
            }
        }

        WHEN("The transposed image is subsampled (by 2x2) and iterated")
        {
            auto result = std::vector<int>();

            for (auto &pix : my_image.transpose().subview(2, 2))
                result.emplace_back(pix);

            THEN("The result consists of points 0,8,2,10")
            {
                auto data_subsampled = std::vector<int>({ 0, 8, 2, 10 });
                CHECK(result == data_subsampled);
            }
        }
    }
}

SCENARIO("Looping over an image coordinates")
{
    auto get_reference_striped_image_vector = [](int width, int height)
    {
        auto result = std::vector<int>(width * height);
        for (int i = 0; i < height; i++)
        {
            for (int j = 0; j < width; j++)
            {
                result[j + i * width] = i + j;
            }
        }
        return result;
    };

    GIVEN("An image of size 4x3")
    {
        const int width = 4;
        const int height = 3;
        auto my_image = image<int>(height, width);
        WHEN("We loop over the image coordinates using `auto` keyword with single `&` suggested for performance")
        {
            for (auto &p : my_image.size())
            {
                // roi_point iterator has ALSO members _x and _y, so the API looks and feels the same
                // however the xy_iterator carries the extra `width and height`
                my_image(p) = p._x + p._y;
            }

            THEN("The resulting image equals the reference image with diagonal stripes")
            {
                auto expected_image_data = get_reference_striped_image_vector(width, height);
                CHECK(my_image.to_vector() == expected_image_data);
            }
        }

        AND_WHEN("We loop over the image coordinates converting to roi_point")
        {
            for (roi_point p : my_image.size())
            {
                // here `p` has been explicitly constructed from the xy_iterator
                my_image(p) = p._x + p._y;
            }

            THEN("The resulting image equals the reference image with diagonal stripes")
            {
                auto expected_image_data = get_reference_striped_image_vector(width, height);
                CHECK(my_image.to_vector() == expected_image_data);
            }
        }
    }
}

SCENARIO("Looping over an image and accessing individual pixels by coordinate")
{
    GIVEN("An image of size 4x3")
    {
        const int width = 4;
        const int height = 3;
        auto my_image = image<int>(height, width);
        WHEN("The image is using two nested for loops writing to each pixel location the sum of the coordinates")
        {
            for (auto y : my_image._height)
            {
                for (auto x : my_image._width)
                {
                    my_image(y, x) = y + x;
                }
            }
            THEN("The resulting image as vector matches the expected")
            {
                auto expected_image_data = std::vector<int>({
                    0, 1, 2, 3,
                    1, 2, 3, 4,
                    2, 3, 4, 5
                });
                CHECK(my_image.to_vector() == expected_image_data);
            }
        }
    }
}
