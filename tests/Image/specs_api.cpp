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
#include <memory>
#include <random>
#include <algorithm>
#include <vector>

using namespace Teisko;
#define ENABLE_FILE_READ_WRITE_SCENARIOS

/// Convenience function to convert stack allocated arrays to Teisko images
template <int HEIGHT, int WIDTH, typename T>
image<T> to_image(T(&array2d)[HEIGHT][WIDTH])
{
    return image<T>(HEIGHT, WIDTH, (T*)array2d);
}

template <typename T>
uint32_t count_non_zero_elements(image<T> &src)
{
    uint32_t result = 0;
    src.foreach([&result](T &src){ if (src != 0) result++; });
    return result;
}

SCENARIO("Libimage provides a method to allocate a 2d image with uninitialized data")
{
    GIVEN("Type and image dimensions")
    {
        using type = uint8_t;
        const int width = 77;
        const int height = 13;
        WHEN("An image in constructed with the given parameters only")
        {
            auto img = image<type>(height, width);
            THEN("The image dimensions match the given dimensions")
            {
                CHECK(img._width == width);
                CHECK(img._height == height);
                AND_THEN("The image is allocated as contiguous")
                {
                    REQUIRE(img._begin != nullptr);
                    CHECK(img.is_contiguous() == true);
                    auto address_of_first_pixel = &img(0, 0);
                    auto address_of_last_pixel = &img(height-1, width-1);
                    auto distance_of_first_and_last_pixels = address_of_last_pixel - address_of_first_pixel;
                    CHECK(distance_of_first_and_last_pixels + 1 == height * width);
                }
            }
        }
    }
}

SCENARIO("Libimage provides a method to construct an image view to pre-allocated data")
{
    GIVEN("A vector of data")
    {
        const int width = 4;
        const int height = 3;
        std::vector<int32_t> pixel_values({ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 });
        WHEN("We map that data to 3x4 image")
        {
            auto img = image<int32_t>(height, width, pixel_values.data());
            THEN("By accessing the pixel values by the img provided methods, the vector changes as well")
            {
                img.at(0, 0)++;                     // increasing the first pixel
                img.at(height-1, width-1) >>= 1;                 // dividing the last pixel by two
                CHECK(pixel_values.front() == 2);
                CHECK(pixel_values.back() == 6);
            }
        }
    }
}

SCENARIO("Libimage size can be passed and returned as roi_point")
{
    GIVEN("A roi_point of x=40, y=30")
    {
        auto size = roi_point(40, 30);
        WHEN("A Teisko image is constructed with this image size")
        {
            auto i = image<int>(size);
            THEN("The image width is 40 and the image height is 30")
            {
                CHECK(i._width == 40);
                CHECK(i._height == 30);
                AND_THEN("The image class returns the size as a roi point")
                {
                    CHECK(i.size() == roi_point(40, 30));
                }
            }
        }
    }
}

SCENARIO(
    "Libimage can be indexed using roi_points. The order of indexing with roi_points is (column,row). "
    "In contrast, the order when indexing by individual coordinates is image(row, column). ")
{
    GIVEN("An image of 3x2 pixels")
    {
        float z[6] = {
            0.0f, 8.1f, 8.1f,       // this row has [two of] the maximum
            8.0f, 0.3f, -0.1f };    // and this row has the minimum
        auto i = image<float>(2, 3, z);

        THEN("The pixels can be accessed with zero based indexing using the roi_point")
        {
            CHECK(&i(roi_point{}) == i._begin);
            CHECK(i(roi_point{ 1, 0 }) == 8.1f);
            CHECK(i(roi_point{ 2, 1 }) == -0.1f);

            // NOTICE the difference in order of parameters
            CHECK(&i(roi_point{ 2, 1 }) == &i(1, 2));
        }

        WHEN("The image maximum is queried")
        {
            auto max_point = i.max_element();
            THEN("The coordinate matches to first of the maximums x=1, y=0")
            {
                CHECK(max_point._x == 1);
                CHECK(max_point._y == 0);
                AND_THEN("The image addressed by the returned point equals 8.1f")
                {
                    CHECK(i(max_point) == 8.1f);
                }
            }
        }

        WHEN("The image minimum is queried")
        {
            auto min_point = i.min_element();
            THEN("The coordinate matches to x=2, y=1")
            {
                CHECK(min_point._x == 2);
                CHECK(min_point._y == 1);
                AND_THEN("The image addressed by the returned point equals -0.1f")
                {
                    CHECK(i(min_point) == -0.1f);
                }
            }
        }
    }
}

SCENARIO("Libimage allows regions to be accessed with roi_points")
{
    GIVEN("An image of W40 x H30")
    {
        auto size = roi_point(40, 30);
        auto i = image<int>(size);
        WHEN("We extract a view of 50% (truncated) offset at 25% (also truncated) of the image size")
        {
            auto view = i.region(size / 2, size / 4);
            THEN("The size of the view is 20 x 15")
            {
                CHECK(view.size() == roi_point(20, 15));
                AND_THEN("The region top left corner is aligned to position X=10, Y=7")
                {
                    CHECK(&view(0, 0) == &i(7, 10));
                }
            }
        }
    }
}

#ifdef ENABLE_FILE_READ_WRITE_SCENARIOS
#if (defined(WIN32) || defined(_WIN32))
static const std::string temp_file_name = ".\\temp.raw";
#else
// Linux
static const std::string temp_file_name = "./temp.raw";
#endif
SCENARIO("Libimage provides methods to read and write raw data to file system")
{
    GIVEN("A vector of data mapped as an image")
    {
        const int width = 4;
        const int height = 3;
        std::vector<int32_t> pixel_values({ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 });
        auto img = image<int32_t>(height, width, pixel_values.data());
        WHEN("We write that data to a temporary file")
        {
            img.write(temp_file_name);
            THEN("We can read it to another image of same size")
            {
                auto img2 = image<int32_t>(height, width).read(temp_file_name);
                CHECK(img2.to_vector() == img.to_vector());
            }
        }
    }
}
#endif

SCENARIO("ImagesView class provides methods to access sub regions")
{
    // Ad hoc converter from string literal to vector
    // Catch can compare std::vectors
    auto to_vector = [](const char *literal)
    {
        std::vector<char> a;
        while (*literal)
            a.push_back(*literal++);
        return a;
    };

    GIVEN("A 2d Buffer filled with 24 letters A-X")
    {
        char letters[] =
            "ABCDEFGH"      // In reality the strings are concatenated
            "IJKLMNOP"
            "QRSTUVWX\0";

        auto view = image<char>(3, 8).init_from(letters);

        THEN("The iterator for Columns 1 and 2 'BC' 'JK' 'RS'")
        {
            auto output = view.columns(2, 1).to_vector();
            CHECK(output == to_vector("BC" "JK" "RS"));
        }

        THEN("A subsampler can access every third item such as 'J M P'")
        {
            auto output = view.subview(3, 3, 1, 1);  // moving origin to 'J', then subsampling
            CHECK(output.to_vector() == to_vector("JMP"));
            AND_THEN("The subsampler can be further mirrored for order 'PMJ'")
            {
                CHECK(output.mirror().to_vector() == to_vector("PMJ"));
            }
        }

        THEN("The iterator for Row 2 returns last row QRSTUVWX")
        {
            auto output = view.rows(1, 2).to_vector();
            CHECK(output == to_vector("QRSTUVWX"));
        }

        THEN("The mirrored view serializes all rows right to left")
        {
            auto mirrored_output = view.mirror().to_vector();
            CHECK(mirrored_output == to_vector("HGFEDCBA" "PONMLKJI" "XWVUTSRQ"));
        }

        THEN("The flipped view serializes all rows bottom to top")
        {
            auto flipped_output = view.flip().to_vector();
            CHECK(flipped_output == to_vector("QRSTUVWX" "IJKLMNOP" "ABCDEFGH"));
        }

        THEN("The view to array flipped and mirrored reverses the original vector")
        {
            auto reversed = view.to_vector();
            std::reverse(reversed.begin(), reversed.end());

            CHECK(reversed == view.flip().mirror().to_vector());
        }

        THEN("A Transposed view iterates columns first as in 'AIQ ... 'HPX'")
        {
            CHECK(view.transpose().to_vector() ==
                to_vector("AIQ" "BJR" "CKS" "DLT" "EMU" "FNV" "GOW" "HPX"));
        }

        THEN("2x2 sized subregion starting at offset 1,1 contains the characters 'JKRS'")
        {
            CHECK(to_vector("JKRS") == view.region(2, 2, 1, 1).to_vector());
        }
    }
}

SCENARIO("region() interface allows extraction of rectangular view")
{
    GIVEN("An image of size 5x4")
    {
        auto size = roi_point(5, 4);
        auto image5x5 = image<float>(size);

        WHEN("A region of size 3x2 is extracted from the original image")
        {
            auto newsize = roi_point(3, 2);
            auto region = image5x5.region(newsize);
            THEN("The new region is sized 3x2 with the origins pointing to same element")
            {
                CHECK(region.size() == newsize);
                CHECK(region._begin == image5x5._begin);
            }
        }

        WHEN("A region of size 3x2 is extracted from the original image giving an OFFSET of x=2,y=1")
        {
            auto newsize = roi_point(3, 2);
            auto offset = roi_point(2, 1);
            auto region = image5x5.region(newsize, offset);
            THEN("The new region is sized 3x2 with the origin of the new image pointing to x=2, y=1 in the source image")
            {
                CHECK(region.size() == newsize);
                CHECK(region._begin == &image5x5(offset));
            }
        }

        WHEN("A region with non-positive coordinate(s) is extracted")
        {
            auto relative_size = roi_point(-1, 0);
            auto region = image5x5.region(relative_size);
            THEN("The new region is sized 4 x 4")
            {
                CHECK(region.size() == roi_point(4, 4));
            }
        }

        WHEN("A region with one positive and one non-positive coordinate is extracted")
        {
            auto offset_x = -1;
            auto offset_y = -2;
            auto region_a = image5x5.region(roi_point(offset_x, 3));
            auto region_b = image5x5.region(roi_point(3, offset_y));
            THEN("The axis with positive new dimension matches the given positive value")
            {
                CHECK(region_a.size()._y == 3);
                CHECK(region_b.size()._x == 3);
                AND_THEN("The size of the axis with non-positive (or zero) dimension is relative to the original")
                {
                    CHECK(region_a.size()._x == (size._x + offset_x));
                    CHECK(region_b.size()._y == (size._y + offset_y));
                }
            }
        }
    }
}

SCENARIO("Image buffers that manage internally the data can be created by the image class")
{
    auto unique_addresses = [](image<short> &img)
    {
        std::vector<short *> addresses;
        img.foreach([&addresses](short &a) { addresses.push_back(&a); });

        std::sort(addresses.begin(), addresses.end());
        return std::distance(addresses.begin(), std::unique(addresses.begin(), addresses.end()));
    };
    const int width = 4;
    const int height = 3;

    GIVEN("An image buffer with initially of size Width x Height")
    {
        auto buffer = image<short>(height, width).fill(7);

        THEN("The buffer has given dimensions and is initialized to 7")
        {
            REQUIRE(buffer._width == width);
            REQUIRE(buffer._height == height);
            for (int i = 0; i < height; i++)
            {
                for (int j = 0; j < width; j++)
                {
                    CHECK(buffer(i, j) == 7);
                }
            }

            AND_THEN("There are width*height unique addresses used to store the elements")
            {
                CHECK(unique_addresses(buffer) == width * height);
            }
        }
    }
}

SCENARIO("A view can be copied to using convert_to()")
{
    GIVEN("A 2d array and it's Image View representation")
    {
        const int width = 6;
        const int height = 4;
        char local_array[height][width] = {
            // Column
            // 0    1    2    3    4    5
            { 'A', 'B', 'C', 'D', 'E', 'F' },        // Row 0
            { 'G', 'H', 'I', 'J', 'K', 'L' },        // Row 1
            { 'M', 'N', 'O', 'P', 'Q', 'R' },        // Row 2
            { 'S', 'T', 'U', 'V', 'W', 'X' } };      // Row 3
        auto view = to_image(local_array);
        WHEN("A new buffer is initialized from the view")
        {
            auto buffer = view.convert_to();

            THEN("The contents of the two buffers match while the addresses are distinct")
            {
                REQUIRE(buffer._width == width);
                REQUIRE(buffer._height == height);

                for (int i = 0; i < height; i++)
                {
                    for (int j = 0; j < width; j++)
                    {
                        CHECK(buffer(i, j) == view(i, j));
                        CHECK(&buffer(i, j) != &view(i, j));
                    }
                }
            }
        }
    }
}

SCENARIO("Benchmarking 3x3 filter using the filter method [.perf]")
{
    /// X-a X-b X-c    X == width * y + x
    /// X-d  X  X+d    Sum of the 9 terms == 9*X, since due to symmetry
    /// X+c X+b X+a    terms marked with +- cancel each other
    auto analytical = [](int y, int x, int h, int w)
    {
        (void)(h);
        return x + y*w;
    };

    typedef int customtype;
    const int w = 4000;
    const int h = 3000;
    GIVEN("An 12M pixel image initialized from 0 to 12e6-1")
    {
        customtype i = customtype(0);
        auto buffer = image<customtype>(h, w)
            .generate([&i](customtype &dst){ dst = i++; });
        WHEN("The image is filtered with custom sum filter")
        {
            // The analytical solution to this kernel should equal
            // to the middle element (given that the generator produces a gradient)
            auto filtered = buffer.filter<3, 3>([](support<3, 3, customtype> &arr) {
                return
                    (arr[0][0] + arr[0][1] + arr[0][2] +
                    arr[1][0] + arr[1][1] * 8 + arr[1][2] +
                    arr[2][0] + arr[2][1] + arr[2][2]) / 16;
            });
            THEN("The few calculated values match the analytically calculated values")
            {
                customtype ref = (customtype)(analytical(1, 1, h, w));
                customtype ref_last = (customtype)(analytical(h - 2, w - 2, h, w));

                CHECK(filtered(1, 1) == ref);
                CHECK(filtered(h - 2, w - 2) == ref_last);
            }
        }
    }
}

SCENARIO("A small 3x3 filter output matches matlab reference value")
{
    GIVEN("A random 4x6 image and a 3x3 filter functor")
    {
        short img[4][6] = {
            { 53, 41, 63, 63, 27, 43 },
            { 59, 6, 63, 32, 60, 2 },
            { 8, 18, 10, 52, 52, 56 },
            { 60, 36, 64, 9, 63, 61 } };

        struct filter3x3func
        {
            short operator()(support<3,3,short> &input) const
            {
                return
                    input[0][0] + input[0][1] + input[0][2] +
                    input[1][0] + input[1][1] + input[1][2] +
                    input[2][0] + input[2][1] + input[2][2];
            }
        } custom;

        AND_THEN("The SUM filter in 3x3 kernel with borders replicated matches the reference")
        {
            auto filtered = to_image(img).filter<3,3>(custom, REPLICATE);

            short ref[4][6] = {
                { 418, 442, 435, 461, 360, 290 },
                { 305, 321, 348, 422, 387, 341 },
                { 314, 324, 290, 405, 387, 413 },
                { 346, 356, 298, 386, 426, 534 } };
            CHECK(filtered.to_vector() == to_vector(ref));
        }
    }
}

SCENARIO("Image class provides automatic border generation functions")
{
    GIVEN("A 3x3 matrix and it's associated view")
    {
        short tmp[3][3] = {
            { 1, 2, 3 },
            { 4, 5, 6 },
            { 7, 8, 9 } };

        auto original = to_image(tmp);
        THEN("The left border replicated by 2 columns matches the hard coded 3x5 matrix")
        {
            auto new_buffer = original.make_borders(0, 2, 0, 0, REPLICATE);

            short left_replicated[3][5] = {
                { 1, 1, 1, 2, 3 },
                { 4, 4, 4, 5, 6 },
                { 7, 7, 7, 8, 9 } };
            CHECK(new_buffer.to_vector() == to_vector(left_replicated));
        }
        THEN("The left border odd mirrored by two matches the hard coded 3x5 matrix")
        {
            auto new_buffer = original.make_borders(0, 2, 0, 0, MIRROR_ODD);

            short left_odd_mirrored[3][5] = {
                { 3, 2, 1, 2, 3 },
                { 6, 5, 4, 5, 6 },
                { 9, 8, 7, 8, 9 } };

            CHECK(new_buffer.to_vector() == to_vector(left_odd_mirrored));
        }
        THEN("The left border even mirrored by two matches the hard coded 3x5 matrix")
        {
            auto new_buffer = original.make_borders(0, 2, 0, 0, MIRROR_EVEN);

            short left_odd_mirrored[3][5] = {
                { 2, 1, 1, 2, 3 },
                { 5, 4, 4, 5, 6 },
                { 8, 7, 7, 8, 9 } };

            CHECK(new_buffer.to_vector() == to_vector(left_odd_mirrored));
        }

        THEN("The right border replicated by 3 columns matches the pre-calculated matrix")
        {
            auto new_buffer = original.make_borders(0, 0, 0, 3, REPLICATE);

            short right_replicated[3][6] = {
                { 1, 2, 3, 3, 3, 3 },
                { 4, 5, 6, 6, 6, 6 },
                { 7, 8, 9, 9, 9, 9 } };

            CHECK(new_buffer.to_vector() == to_vector(right_replicated));
        }
        THEN("The right border odd mirrored by TWO matches the hard coded 3x5 matrix")
        {
            auto new_buffer = original.make_borders(0, 0, 0, 2, MIRROR_ODD);

            short right_odd_mirrored[3][5] = {
                { 1, 2, 3, 2, 1 },
                { 4, 5, 6, 5, 4 },
                { 7, 8, 9, 8, 7 } };

            CHECK(new_buffer.to_vector() == to_vector(right_odd_mirrored));
        }
        THEN("The right border even mirrored by two matches the hard coded 3x5 matrix")
        {
            auto new_buffer = original.make_borders(0, 0, 0, 2, MIRROR_EVEN);

            short right_even_mirrored[3][5] = {
                { 1, 2, 3, 3, 2 },
                { 4, 5, 6, 6, 5 },
                { 7, 8, 9, 9, 8 } };

            CHECK(new_buffer.to_vector() == to_vector(right_even_mirrored));
        }
        THEN("The Top Border can be replicated by e.g. four")
        {
            auto new_buffer = original.make_borders(4, 0, 0, 0, REPLICATE & TOP);

            short top_replicated[7][3] = {
                { 1, 2, 3 },
                { 1, 2, 3 },
                { 1, 2, 3 },
                { 1, 2, 3 },
                { 1, 2, 3 },
                { 4, 5, 6 },
                { 7, 8, 9 } };

            CHECK(new_buffer.to_vector() == to_vector(top_replicated));
        }
        THEN("The Top border odd mirrored by two matches reference")
        {
            auto new_buffer = original.make_borders(2, 0, 0, 0, MIRROR_ODD);

            short top_ref[5][3] = {
                { 7, 8, 9 },
                { 4, 5, 6 },
                { 1, 2, 3 },
                { 4, 5, 6 },
                { 7, 8, 9 } };

            CHECK(new_buffer.to_vector() == to_vector(top_ref));
        }
        THEN("The Top border even mirrored by two matches reference")
        {
            auto new_buffer = original.make_borders(2, 0, 0, 0, MIRROR_EVEN);

            short top_ref[5][3] = {
                { 4, 5, 6 },
                { 1, 2, 3 },
                { 1, 2, 3 },
                { 4, 5, 6 },
                { 7, 8, 9 } };

            CHECK(new_buffer.to_vector() == to_vector(top_ref));
        }

        THEN("The bottom border can be replicated by three rows")
        {
            auto new_buffer = original.make_borders(0, 0, 3, 0, REPLICATE & BOTTOM);

            short bottom_replicated[6][3] = {
                { 1, 2, 3 },
                { 4, 5, 6 },
                { 7, 8, 9 },
                { 7, 8, 9 },
                { 7, 8, 9 },
                { 7, 8, 9 } };

            CHECK(new_buffer.to_vector() == to_vector(bottom_replicated));
        }

        THEN("The bottom border odd mirrored by two rows matches reference")
        {
            auto new_buffer = original.make_borders(0, 0, 2, 0, MIRROR_ODD & BOTTOM);

            short bottom_reference[5][3] = {
                { 1, 2, 3 },
                { 4, 5, 6 },
                { 7, 8, 9 },
                { 4, 5, 6 },
                { 1, 2, 3 } };

            CHECK(new_buffer.to_vector() == to_vector(bottom_reference));
        }
        THEN("The bottom border even mirrored by two rows matches reference")
        {
            auto new_buffer = original.make_borders(0, 0, 2, 0, MIRROR_EVEN & BOTTOM);

            short bottom_reference[5][3] = {
                { 1, 2, 3 },
                { 4, 5, 6 },
                { 7, 8, 9 },
                { 7, 8, 9 },
                { 4, 5, 6 } };

            CHECK(new_buffer.to_vector() == to_vector(bottom_reference));
        }
        THEN("Arbitrary amounts of borders can be added using replicate scheme")
        {
            const int left = 6;
            const int right = 7;
            const int bottom = 9;
            const int top = 11;
            auto new_buffer = original.make_borders(top, left, bottom, right, REPLICATE);

            AND_THEN("The new image dimensions and number of non-zero elements matches the expected")
            {
                CHECK(new_buffer._width == 3 + left + right);
                CHECK(new_buffer._height == 3 + top + bottom);
                CHECK(count_non_zero_elements(new_buffer) == (3 + left + right) * (3 + top + bottom));
            }
        }
        THEN("Arbitrary amounts of borders can be added using zero padding")
        {
            const int left = 6;
            const int right = 7;
            const int bottom = 9;
            const int top = 11;
            auto new_buffer = original.make_borders(top, left, bottom, right, ZERO_PAD);

            AND_THEN("The new image dimensions and number of non-zero elements matches the expected")
            {
                CHECK(new_buffer._width == 3 + left + right);
                CHECK(new_buffer._height == 3 + top + bottom);
                CHECK(count_non_zero_elements(new_buffer) == count_non_zero_elements(original));
            }
        }
    }
}

SCENARIO(
    "Benchmarking 3x3 filter using foreach method with"
    " manual setup of stride parameters in capture list [.perf]"
    )
{
    /// X-a X-b X-c    X == width * y + x
    /// X-d  X  X+d    Sum of the 9 terms == 9*X, since due to symmetry
    /// X+c X+b X+a    terms marked with +- cancel each other
    auto analytical = [](int y, int x, int h, int w)
    {
        (void)(h);
        return x + y*w;
    };

    typedef short customtype;
    const int w = 4000;
    const int h = 3000;
    GIVEN("An 12M pixel image initialized from 0 to 12e6-1")
    {
        customtype i = customtype(0);
        auto buffer = image<customtype>(h, w)
            .generate([&i](customtype &dst){ dst = i++; });
        WHEN("The image is filtered with custom sum filter")
        {
            auto filtered = image<customtype>(h, w);
            // We have to pass the relative offsets of src[y +- 1][] and src[][x +- 1]
            const int skip_x = buffer._skip_x;
            const int skip_y = buffer._skip_y;
            // Also we have to restrict both the input and output buffers to read/write
            // from the center of the image
            auto center = buffer.region(-2, -2, 1, 1);
            filtered.region(-2,-2,1,1).
                foreach([skip_x, skip_y](customtype &dst, customtype &src)
            {
                customtype *s = &src;
                dst = (s[-skip_y - skip_x] + s[-skip_y] + s[-skip_y + skip_x]
                    + s[-skip_x] + 8 * s[0] + s[skip_x]
                    + s[skip_y - skip_x] + s[skip_y] + s[skip_y + skip_x]) / 16;
            }, center);
            THEN("The few calculated values match the analytically calculated values")
            {
                customtype ref = (customtype)(analytical(1, 1, h, w));
                customtype ref_last = (customtype)(analytical(h - 2, w - 2, h, w));

                CHECK(filtered(1, 1) == ref);
                CHECK(filtered(h - 2, w - 2) == ref_last);
            }
        }
    }
}

SCENARIO("Libimage can filter source image to target image using the support class. "
    "The function `target.filter<W,H>(src_image, ...)` allows full control of target image layout, "
    "and using different element type between input and output images. ")
{
    GIVEN("A source image filled with data and an empty target image of same size")
    {
        typedef short customtype;
        const int w = 10;
        const int h = 15;
        customtype value = customtype(0);
        auto source = image<customtype>(h, w)
            .generate([&value](customtype &dst){ dst = value++; });

        auto in_place_target = image<double>(h, w);

        WHEN("We apply a filter between the target and source images with automatic border generation")
        {
            in_place_target.template filter<3, 3>(source, [](support<3, 3, customtype> &data) -> double
            {
                return static_cast<double>(data.median());
            });
            auto x = in_place_target.to_vector();
            AND_WHEN("We apply the same kind of filter with the non in-place version")
            {
                auto out_place_target = source.filter<3, 3>([](support<3, 3, customtype> &data)
                {
                    return data.median();
                });
                auto y = out_place_target.to_vector();
                THEN("We notice that the content is the same, except for the type")
                {
                    auto vec_first = in_place_target.to_vector();
                    auto vec_second = out_place_target.to_vector();
                    REQUIRE(vec_first.size() == vec_second.size());
                    for (size_t i = 0; i < vec_first.size(); i++)
                    {
                        CHECK(vec_first[i] == static_cast<double>(vec_second[i]));
                    }
                    AND_THEN("We notice that the in-place layout is contiguous, while the out-place layout is not")
                    {
                        CHECK(in_place_target.is_contiguous() == true);
                        CHECK(out_place_target.is_contiguous() == false);
                    }
                }
            }
        }
    }
}

SCENARIO("Libimage allows foreach only when other images are at least as large as the calling object. ")
{
    GIVEN("One image of 8x10")
    {
        auto size = roi_point(10, 8);
        auto img = image<uint16_t>(size);
        auto img_bigger = image<float>(size + 5);
        THEN("The it's possible to `foreach` the image with 0,1,2 (and more) images of at least same size")
        {
            CHECK_NOTHROW(img.foreach([](uint16_t &i) { ++i; }));
            CHECK_NOTHROW(img.foreach([](uint16_t &i, float j) { i += static_cast<uint16_t>(j); }, img_bigger.flip()));
            CHECK_NOTHROW(img.foreach([](uint16_t &i, float &j, uint16_t &k) { i += (static_cast<uint16_t>(j)-k); }, img_bigger, img.mirror()));
        }
        AND_THEN("It's not possible to `foreach` the bigger image with a smaller image included")
        {
            CHECK_THROWS(img_bigger.foreach([](float &i, uint16_t &j) { i += j; }, img));
            CHECK_THROWS(img_bigger.foreach([](float &i, float j, uint16_t &k) { i += j * k; }, img_bigger.mirror(), img.flip()));
        }
    }
}