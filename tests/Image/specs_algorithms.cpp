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


#include "Teisko/Preprocessing.hpp"
#include "Teisko/Algorithm/Bit.hpp"
#include "catch.hpp"
#include <memory>
#include <random>
#include <algorithm>
#include <vector>

using namespace Teisko;

namespace teisko_median7x7_tests
{
    /// Conversion function from fixed size array to std::vector
    ///  - CATCH framework can more easily compare vectors for equality
    template <int ROWS, int COLS, typename T>
    std::vector<T> to_vector(T(&arr)[ROWS][COLS])
    {
        return std::vector<T>(&arr[0][0], &arr[0][0] + ROWS*COLS);
    }

    /// Convenience function making a Teisko image (view) from a fixed sized array
    template <int HEIGHT, int WIDTH, typename T>
    image<T> to_image(T(&array2d)[HEIGHT][WIDTH])
    {
        return image<T>(HEIGHT, WIDTH, (T*)array2d);
    }

    /// data used in our tests -- items between 10..99 give nicer look
    /// %% r = (10 + fix(rand(8, 8) * 90));
    /// %% sprintf('{%d, %d, %d, %d, %d, %d, %d, %d },\n', r')
    int8_t data[8][8] = {
        { 60, 86, 30, 98, 33, 18, 38, 39 },
        { 89, 68, 44, 75, 63, 63, 57, 19 },
        { 70, 43, 62, 40, 12, 52, 68, 64 },
        { 27, 27, 32, 62, 48, 72, 46, 80 },
        { 43, 48, 36, 19, 38, 72, 83, 48 },
        { 51, 53, 65, 91, 24, 67, 74, 18 },
        { 98, 20, 33, 89, 26, 13, 97, 33 },
        { 24, 63, 84, 83, 48, 16, 57, 23 }
    };

    /// %% medfilt2(r, [7 7]);
    int8_t matlab_median_zero_padded[8][8] = {
        { 0, 0, 0, 27, 27, 0, 0, 0 },
        { 0, 12, 30, 38, 38, 32, 12, 0 },
        { 0, 30, 40, 48, 46, 38, 24, 0 },
        { 27, 33, 43, 52, 48, 40, 38, 18 },
        { 27, 36, 44, 53, 52, 48, 38, 18 },
        { 0, 26, 33, 46, 46, 36, 23, 0 },
        { 0, 19, 24, 33, 33, 24, 13, 0 },
        { 0, 0, 0, 20, 19, 0, 0, 0 }
    };

    int8_t matlab_median_symmetric[8][8] = {
        { 62, 62, 60, 57, 52, 48, 52, 52 },
        { 60, 60, 60, 57, 48, 46, 48, 48 },
        { 53, 51, 52, 52, 48, 46, 48, 48 },
        { 51, 48, 51, 52, 48, 48, 52, 52 },
        { 51, 48, 51, 53, 52, 48, 57, 57 },
        { 51, 48, 48, 51, 48, 48, 52, 48 },
        { 48, 48, 48, 48, 48, 48, 48, 48 },
        { 53, 51, 51, 53, 48, 38, 48, 38 }
    };

    SCENARIO("Preprocessor median filtering supports multiple padding methods")
    {
        GIVEN("Input data as Teisko image")
        {
            auto inp = to_image(data);
            WHEN("The image is median filtered as zero padded")
            {
                auto outp = median7x7(inp, ZERO_PAD);
                THEN("The output matches exactly the matlab reference output")
                {
                    CHECK(outp.to_vector() == to_vector(matlab_median_zero_padded));
                }
            }
            AND_WHEN("The image is filtered with evenly symmetric replication scheme")
            {
                auto outp = median7x7(inp, SYMMETRIC);
                THEN("The output matches exactly the matlab reference output")
                {
                    CHECK(outp.to_vector() == to_vector(matlab_median_symmetric));
                }
            }
        }
    }

    // Generic test case operating on the same input/output data as different types
    template <typename type>
    void run_median7x7_test_case_as()
    {
        auto inp = to_image(data).convert_to<type>();
        auto med = median7x7(inp, SYMMETRIC);

        auto out_ref = to_image(matlab_median_symmetric).convert_to<type>();
        CHECK(med.to_vector() == out_ref.to_vector());
    };

    SCENARIO("Median filter works with a variety of data types")
    {
        run_median7x7_test_case_as<int8_t>();
        run_median7x7_test_case_as<uint8_t>();
        run_median7x7_test_case_as<int16_t>();
        run_median7x7_test_case_as<uint16_t>();
        run_median7x7_test_case_as<int32_t>();
        run_median7x7_test_case_as<uint32_t>();
        run_median7x7_test_case_as<float>();
        // compare_median7x7_as<int64_t>();    /* currently gives static_assert */
        // compare_median7x7_as<uint64_t>();   /* here as well */
    }

    template <int plane_count>
    void array_to_bitplanes(const int(&data)[8][8], uint64_t (&planes)[plane_count])
    {
        for (int i = 0; i < 64; i++)
        {
            for (int j = 0; j < plane_count; j++)
                if (data[i/8][i%8] & (1 << (plane_count - 1 - j)))
                    planes[j] |= (1ULL << i);
        }
    }

    SCENARIO("We show that the median filter kernel is fully configurable by showing "
        "that by varying the threshold between 1..MaskSize the operation is equal to sorting")
    {
        int small_matrix[8][8] = {
            { 1, 1, 7, 5, 2, 5, 6, 6 },
            { 3, 0, 5, 2, 2, 0, 7, 6 },
            { 4, 3, 3, 5, 7, 0, 0, 6 },
            { 6, 4, 5, 2, 5, 7, 2, 4 },
            { 7, 2, 5, 5, 7, 3, 5, 1 },
            { 5, 3, 5, 2, 2, 0, 3, 0 },
            { 1, 7, 2, 0, 7, 6, 7, 4 },
            { 5, 4, 5, 2, 5, 5, 6, 5 }
        };

        // Helper function to select items from the 8x8 matrix by the mask --> sorted
        auto get_elements_as_sorted = [](int(&matrix)[8][8], uint64_t mask)
        {
            std::vector<int> sorted;
            for (int i = 0; i < 64; i++)
            {
                if ((mask >> i) & 1)
                    sorted.emplace_back(matrix[i / 8][i % 8]);
            }
            std::sort(sorted.begin(), sorted.end());
            return sorted;
        };

        GIVEN("A bitplane serialization of a matrix")
        {
            const int number_of_bits = 3;
            uint64_t planes[number_of_bits] = { 0 };

            array_to_bitplanes(small_matrix, planes);

            WHEN("The median kernel is run for every value of threshold between 1..64")
            {
                std::vector<int> items;
                const uint64_t all_values = ~0ULL;

                for (int threshold = 1; threshold <= 64; threshold++)
                    items.emplace_back((int)median8x8_iteration(planes, all_values, threshold));

                THEN("The algorithm has sorted the 64-element input vector")
                {
                    auto ref = get_elements_as_sorted(small_matrix, all_values);
                    CHECK(items == ref);
                }
            }

            AND_WHEN("The initial mask is randomized to contain fewer than 64 set bits "
                "and we run the median kernel for all thresholds between 1..MaskSize")
            {
                uint64_t mask = 0xd81409534fecaULL;
                int mask_size = (int)popcount(mask);
                std::vector<int> items;

                for (int threshold = 1; threshold <= mask_size; threshold++)
                    items.emplace_back((int)median8x8_iteration(planes, mask, threshold));

                THEN("The algorithm has sorted elements with corresponding bit set in the mask")
                {
                    auto ref = get_elements_as_sorted(small_matrix, mask);
                    CHECK(items == ref);
                }
            }
        }
    }

    SCENARIO("Regression test for horizontal filtering - The coefficients are pre-calculated for "
        "a given image width (10), replication scheme (Zero Pad) and filter kernel (1,2,4,2,1)")
    {
        std::vector<double> coeffs({
            4, 2, 1, 0, 0,  // Precalculating these first two rows allow calculating the first two
            2, 4, 2, 1, 0,  // values with the kernel being partially outside the input image data
            1, 2, 4, 2, 1,  // <-- this row marks the kernel
            0, 1, 2, 4, 2,  // Likewise the last two rows allow calculating output data for columns
            0, 0, 1, 2, 4   // with the filter kernel tail partially outside the input image data
        });
        auto one_row = std::vector<int>({ 1, 0, 0, 0, 0, 0, 1, 1, 0, 0 });
        auto dst_row = std::vector<int>(one_row.size());
        int ref_row[1][10] = { { 4, 2, 1, 0, 1, 3, 6, 6, 3, 1 } };

        auto in = image<int>(1, 10, one_row.data());
        auto out = image<int>(1, 10, dst_row.data());
        filter1d_horizontal(coeffs, 5, in, out);
        CHECK(out.to_vector() == to_vector(ref_row));
    }

    SCENARIO("Testing horizontal filtering with coeff generation")
    {
        std::vector<double> coeffs({ 0, 3, 0, 4 });     // some random 4-element kernel 'h'
        std::vector<double> imagedata({ 6, 2, 1, 3, 1, 4, 8, 5 });  // some random 8-element data
        std::vector<double> reference({ 22, 18, 7, 25, 35, 32, 24, 15 });   // imfilter2(h, data)

        auto coffs4x4 = replicate_filter_kernel(coeffs.data(), 4, (int)imagedata.size(), ZERO_PAD);
        auto image_in = image<double>(1, 8, imagedata.data());
        auto image_out = image<double>(1, 8);

        filter1d_horizontal(coffs4x4, 4, image_in, image_out);
        CHECK(reference == image_out.to_vector());
    }

    SCENARIO("Gaussian filtering can be done for multiple types")
    {
        auto tmpimage = image<int>(200, 300).fill(0);
        tmpimage(100, 150) = 10000;
        CHECK_NOTHROW(gaussian35x35(tmpimage, ZERO_PAD));

        auto tmpimage2 = image<uint16_t>(200, 300).fill(113);
        CHECK_NOTHROW(gaussian35x35(tmpimage2, REPLICATE));
        //tmpimage.write("E:\\tempraw.raw");
    }

    SCENARIO("Trimmed mean calculates averages ignoring high and low values")
    {
        GIVEN("A trimmed mean functor")
        {
            trimmed_mean_f<double> func;
            WHEN("A single element is averaged with 0 percent outliers removed")
            {
                double arr[1] = { 7.0 };
                double trim = 0;
                THEN("The average is the only element 7.0")
                {
                    CHECK(func(image<double>(1, 1, arr), trim) == 7.0);
                }
            }
            WHEN("An array of values from -100, 2, 3, 4, 9 is averaged with 40% of items excluded")
            {
                double arr[5] = { -100, 2, 3, 4, 9 };
                double trim = 0.4;
                THEN("The average is the middle element 3.0")
                {
                    CHECK(func(image<double>(1, 5, arr), trim) == 3.0);
                }
            }
            WHEN("An array is averaged with 100% of items excluded")
            {
                double arr[8] = { 1, 2, 3, 4, 5, 6, 7, 8 };
                double trim = 1.0;
                THEN("The average is zero")
                {
                    CHECK(func(image<double>(2, 4, arr), trim) == 0.0);
                }
            }
        }
    }

    SCENARIO("Image Integrals can be calculated")
    {
        GIVEN("A magic 6x6 square in Teisko image")
        {
            // magic(6); %% matlab command
            auto magic = std::vector<int>({
                35, 1, 6, 26, 19, 24,
                3, 32, 7, 21, 23, 25,
                31, 9, 2, 22, 27, 20,
                8, 28, 33, 17, 10, 15,
                30, 5, 34, 12, 14, 16,
                4, 36, 29, 13, 18, 11
            });
            auto img = image<int>(6, 6, magic.data());
            WHEN("The image is integrated")
            {
                summed_area_table(img);
                THEN("The result matches the matlab calculated result")
                {
                    // cumsum(cumsum(a')'); %% matlab command
                    auto ref = std::vector<int>({
                        35, 36, 42, 68, 87, 111,
                        38, 71, 84, 131, 173, 222,
                        69, 111, 126, 195, 264, 333,
                        77, 147, 195, 281, 360, 444,
                        107, 182, 264, 362, 455, 555,
                        111, 222, 333, 444, 555, 666
                    });
                    CHECK(magic == ref);
                }
            }
        }
    }

    SCENARIO("Images can be box-filtered")
    {
        auto magic = std::vector<int>({
            35, 1, 6, 26, 19, 24,
            3, 32, 7, 21, 23, 25,
            31, 9, 2, 22, 27, 20,
            8, 28, 33, 17, 10, 15,
            30, 5, 34, 12, 14, 16,
            4, 36, 29, 13, 18, 11
        });
        GIVEN("A magic 6x6 image")
        {
            auto copy = magic;
            auto img = image<int>(6, 6, copy.data());
            WHEN("The image is box filtered with kernel <= 1")
            {
                box_filter(img, 1, 0);
                THEN("The image is unaltered")
                {
                    CHECK(copy == magic);
                }
            }
            WHEN("The image is box filtered with kernel == 2")
            {
                box_filter(img, 2, ZERO_PAD);

                THEN("The image confirms to matlab filtered result")
                {
                    // % round(filter2(ones(2), magic(6))/4)
                    auto ref_vector = std::vector<int>({
                        18, 12, 15, 22, 23, 12,
                        19, 13, 13, 23, 24, 11,
                        19, 18, 19, 19, 18, 9,
                        18, 25, 24, 13, 14, 8,
                        19, 26, 22, 14, 15, 7,
                        10, 16, 11, 8, 7, 3
                    });
                    CHECK(copy == ref_vector);
                }
            }
            WHEN("The image is box filtered with kernel == 5 with replicated scheme")
            {
                box_filter(img, 5, REPLICATE);
                THEN("The image confirms to matlab filtered result")
                {
                    // round(imfilter(magic(6),ones(5)/25, 'replicate'))
                    auto ref_vector = std::vector<int>({
                        20, 19, 18, 17, 20, 23,
                        18, 18, 18, 18, 19, 21,
                        19, 19, 18, 18, 19, 20,
                        18, 18, 19, 19, 18, 18,
                        19, 19, 19, 19, 17, 15,
                        18, 19, 20, 20, 17, 14
                    });
                    CHECK(copy == ref_vector);
                }
            }
        }
    }

    SCENARIO("Region properties are extracted")
    {
        GIVEN("A bwlabelled image")
        {
            int labels[4 * 8] = {
                0, 1, 1, 0,
                0, 1, 0, 1,
                1, 0, 0, 1,
                0, 0, 0, 1,
                0, 1, 1, 1,
                1, 0, 0, 1,
                0, 1, 1, 1,
                0, 0, 0, 0
            };
            // % in matlab: regionprops(imfill(labels), 'Area', 'Centroid', 'PerimeterOld')

            auto i = image<int>(8, 4, labels);
            WHEN("The region props are run for any '1' labelled pixel in the image")
            {
                roi_point p = { 2, 0 };
                REQUIRE(i(p) == 1);
                auto result = get_regionprops(i, p);
                CHECK(result.perimeter == Approx(23.313708498984763).epsilon(1e-7));
                CHECK(result.area == 17.0);
                CHECK(result.centroid_x == Approx(1.823529411764706).epsilon(1e-7));
                CHECK(result.centroid_y == Approx(3.470588235294118).epsilon(1e-7));
            }
        }
    }

    SCENARIO("Otsu's method finds a threshold from a bimodal histogram")
    {
        GIVEN("A bimodal (i.e. two-peak) histogram")
        {
            auto twopeak = std::vector<int> { 5, 15, 30, 20, 15, 10, 8, 6, 5, 4, 3, 2, 5, 6, 10, 30, 80, 50, 31};
            WHEN("The histogram is passed to otsu_threshold")
            {
                auto idx = calculate_otsu_threshold(twopeak);
                THEN("The found threshold must be somewhere between the two peaks")
                {
                    int last_index = static_cast<int>(twopeak.size() - 1);
                    CHECK(idx > 2);
                    CHECK(idx < last_index - 2);
                }
            }
        }
    }

    SCENARIO("A quantizer class finds dynamically the two tone feature")
    {
        GIVEN("A bitonal image and a quantizer")
        {
            const uint16_t A = 13213;
            const uint16_t B = A + 1;
            auto image_data = std::vector<uint16_t>(
            {
                A, B, A, A,
                B, B, B, B,
                A, B, A, A
            });

            auto q = quantizer{};

            WHEN("The image is quantized")
            {
                auto result = q(image<uint16_t>(3, 4, image_data.data()));
                THEN("The result matches the expected two-tone image with values of zeros and ones")
                {
                    auto expected_vector = std::vector<uint16_t>(
                    {
                        0, 1, 0, 0,
                        1, 1, 1, 1,
                        0, 1, 0, 0
                    });
                    CHECK(result.to_vector() == expected_vector);
                }
            }
        }
    }

    SCENARIO("Image can be downscaled by scaling ratio")
    {
        GIVEN("A 20x15 image containing one square in the middle")
        {
            auto size = roi_point(20, 16);
            auto img = image<float>(size).fill(100);
            img.region(size / 2, size / 4).fill(200);

            WHEN("The image is downscaled by ratio of 0.5")
            {
                auto half_sized = resize_image(img, point_xy(0.5, 0.5));
                THEN("The rescaled image corresponds to matlab rescaled matrix up to last digit")
                {
                    auto const max_diff = 0.0001f;
                    auto matlab_result = std::vector<float>({
                        100.0000f, 100.0549f, +99.4141f, +98.7732f, +98.8281f, +98.8281f, +98.7732f, +99.4141f, 100.0549f, 100.0000f,
                        100.0000f, +99.6887f, 103.3203f, 106.9519f, 106.6406f, 106.6406f, 106.9519f, 103.3203f, +99.6887f, 100.0000f,
                        100.0000f, +95.6238f, 146.6797f, 197.7356f, 193.3594f, 193.3594f, 197.7356f, 146.6797f, +95.6238f, 100.0000f,
                        100.0000f, +95.2576f, 150.5859f, 205.9143f, 201.1719f, 201.1719f, 205.9143f, 150.5859f, +95.2576f, 100.0000f,
                        100.0000f, +95.2576f, 150.5859f, 205.9143f, 201.1719f, 201.1719f, 205.9143f, 150.5859f, +95.2576f, 100.0000f,
                        100.0000f, +95.6238f, 146.6797f, 197.7356f, 193.3594f, 193.3594f, 197.7356f, 146.6797f, +95.6238f, 100.0000f,
                        100.0000f, +99.6887f, 103.3203f, 106.9519f, 106.6406f, 106.6406f, 106.9519f, 103.3203f, +99.6887f, 100.0000f,
                        100.0000f, 100.0549f, +99.4141f, +98.7732f, +98.8281f, +98.8281f, +98.7732f, +99.4141f, 100.0549f, 100.0000f
                    });

                    auto teisko_result = half_sized.to_vector();
                    REQUIRE(teisko_result.size() == matlab_result.size());
                    for (size_t i = 0; i < teisko_result.size(); i++)
                    {
                        CHECK(std::fabs(teisko_result[i] - matlab_result[i]) < max_diff);
                    }
                }
            }

            AND_WHEN("The image is downscaled by ratio of 0.2")
            {
                auto fifth_sized = resize_image(img, point_xy(0.2, 0.2));
                THEN("The rescaled image corresponds to matlab rescaled matrix up to last digit")
                {
                    auto const max_diff = 0.0001f;
                    auto matlab_result = std::vector<float>({
                        101.7060f, 120.3325f, 120.3325f, 101.7060f,
                        108.2038f, 197.7756f, 197.7756f, 108.2038f,
                        103.1611f, 137.6750f, 137.6750f, 103.1611f,
                        +99.6864f, 096.2624f, +96.2624f, +99.6864f
                    });
                    auto teisko_result = fifth_sized.to_vector();
                    REQUIRE(teisko_result.size() == matlab_result.size());
                    for (size_t i = 0; i < teisko_result.size(); i++)
                    {
                        CHECK(std::fabs(teisko_result[i] - matlab_result[i]) < max_diff);
                    }
                }
            }
        }
    }

    SCENARIO("Image can be upscaled by resize_image")
    {
        GIVEN("A 6x6 image containing one square in the middle")
        {
            auto input_data = std::vector<float>({
                101, 102, 103, 104, 105, 106,
                107, 108, 109, 110, 111, 112,
                113, 114, 666, 666, 115, 116,
                117, 118, 666, 666, 119, 120,
                121, 122, 123, 124, 125, 126,
                127, 128, 129, 130, 131, 132
            });
            auto input_img = image<float>(6, 6, input_data.data());
            WHEN("The image is upscaled by ratio of 1.16")
            {
                auto rescaled_img = resize_image(input_img, point_xy(1.16, 1.16));
                THEN("The resulting image matches the matlab generated one")
                {
                    CHECK(rescaled_img._width == 7);
                    CHECK(rescaled_img._height == 7);
                    auto const max_diff = 0.0001f;
                    auto matlab_result = std::vector<float>({
                        100.7908f, 101.5968f, 102.4758f, 103.3379f, 104.2000f, 105.0841f, 105.8622f,
                        105.6269f, 108.7659f, +82.5598f, +67.8736f, +85.8559f, 112.5452f, 110.8285f,
                        111.0492f, +87.0886f, 375.4950f, 541.4132f, 360.5325f, +87.4777f, 114.7402f,
                        115.1596f, +75.7881f, 543.1090f, 811.7306f, 517.7555f, +74.2387f, 117.9814f,
                        118.3980f, +96.1192f, 365.0022f, 519.7095f, 351.1641f, +96.7068f, 122.1713f,
                        122.4121f, 125.6903f, +97.8679f, +82.2542f, 101.2586f, 129.4881f, 127.6225f,
                        127.2196f, 128.0256f, 128.9046f, 129.7667f, 130.6288f, 131.5129f, 132.2910f
                    });
                    auto teisko_result = rescaled_img.to_vector();
                    REQUIRE(teisko_result.size() == matlab_result.size());
                    for (size_t i = 0; i < teisko_result.size(); i++)
                    {
                        CHECK(std::fabs(teisko_result[i] - matlab_result[i]) < max_diff);
                    }
                }
            }
            AND_WHEN("The image is upscaled to predetermined size of 7x9")
            {
                auto dst_size = roi_point(7, 9);
                auto destination_img = image<float>(dst_size);
                resize_image(input_img, destination_img);
                THEN("The resulting image matches the matlab generated one")
                {
                    CHECK(destination_img._width == 7);
                    CHECK(destination_img._height == 9);
                    auto const max_diff = 0.0001f;
                    auto matlab_result = std::vector<float>({
                        100.6220f, 101.4205f, 102.2956f, 103.1528f, 104.0099f, 104.8851f, 105.6836f,
                        103.5942f, 106.6705f, +81.9953f, +67.4102f, +83.8016f, 110.2643f, 108.7808f,
                        107.9924f, 104.4008f, 154.5184f, 185.1372f, 156.0555f, 107.6167f, 112.8132f,
                        112.0849f, +80.8095f, 441.4541f, 659.7517f, 441.8736f, +82.4568f, 115.3873f,
                        115.0942f, +75.0403f, 534.1529f, 811.9609f, 534.2113f, +76.1805f, 117.9058f,
                        117.6127f, +86.6007f, 444.2918f, 660.8003f, 444.7113f, +88.2479f, 120.9151f,
                        120.1868f, 116.6332f, 166.3250f, 196.6858f, 167.8621f, 119.8492f, 125.0076f,
                        124.2192f, 127.2790f, 102.7893f, +88.3164f, 104.5955f, 130.8728f, 129.4058f,
                        127.3164f, 128.1149f, 128.9901f, 129.8472f, 130.7044f, 131.5795f, 132.3780f
                    });
                    auto teisko_result = destination_img.to_vector();
                    REQUIRE(teisko_result.size() == matlab_result.size());
                    for (size_t i = 0; i < teisko_result.size(); i++)
                    {
                        CHECK(std::fabs(teisko_result[i] - matlab_result[i]) < max_diff);
                    }
                }
            }
        }
    }

    SCENARIO("Image can be rotated with nearest neighbor interpolation. "
        "We allow currently small deviation from Matlab result. "
        "The few pixels that do not match matlab results must be a close neighbor of the matlab result."
        "This is tested by a regular pattern with left/right/top/bottom neighbor calculated analytically.")
    {
        GIVEN("A 6x6 input with values from 100..135")
        {
            int generator = 100;
            auto img = image<int>(6, 6).generate([&generator](int &x) { x = generator++; });

            /// The input is a grid of form
            /// 100 101 102 ...
            /// 106 107 108 ...
            /// 112 113 114 ...

            WHEN("The input is rotated by 30 degrees with the out-of-image pixels set to 999")
            {
                // We set the out of image pixels to 999 to align the
                auto res = rotate_image(img, 30.0, 999);
                auto teisko_result = res.to_vector();
                THEN("All the resulting values either belong to the union of input values and 999")
                {
                    for (auto &x : teisko_result)
                    {
                        CHECK(((x == 999) || (x >= 100 && x <= 135)) == true);
                    }
                }
                AND_THEN("The image matches the matlab version in most positions")
                {
                    CHECK(res.size() == roi_point(9, 9));
                    // % k = imrotate(reshape(100:135,[6 6])', 30.0);
                    // % k(k==0) = 999
                    auto reference_vec = std::vector<int>({
                        999, 999, 999, 999, 999, 105, 999, 999, 999,
                        999, 999, 999, 999, 104, 105, 999, 999, 999,
                        999, 999, 102, 103, 110, 110, 117, 999, 999,
                        100, 100, 107, 108, 115, 116, 123, 999, 999,
                        999, 106, 107, 114, 115, 121, 128, 129, 999,
                        999, 999, 112, 119, 120, 127, 128, 135, 135,
                        999, 999, 118, 125, 126, 132, 133, 999, 999,
                        999, 999, 999, 130, 131, 999, 999, 999, 999,
                        999, 999, 999, 130, 999, 999, 999, 999, 999 });
                    REQUIRE(teisko_result.size() == reference_vec.size());
                    int mismatches = 0;
                    for (size_t i = 0; i < teisko_result.size(); i++)
                    {
                        if (teisko_result[i] != reference_vec[i])
                        {
                            mismatches++;
                            auto diff = std::abs(teisko_result[i] - reference_vec[i]);
                            CHECK(((diff == 1) || (diff == 6)) == true);
                        }
                    }
                    CHECK(mismatches < 5);
                }
            }
            AND_WHEN("The image is rotated to existing array")
            {
                auto img2 = image<float>(7, 7).fill(888);   // initial value is not included in original values
                rotate_image(img, img2, -22.5, 0.f);
                THEN("The image only contains the original pixel values (or 0.f)")
                {
                    for (auto &x : img2.to_vector())
                    {
                        CHECK(((x == 0.f) || (x >= 100.f && x <= 135.f)) == true);
                    }
                }
            }
        }
    }

#if 0
    SCENARIO("Median filter performance test")
    {
        struct {
            int a;
            void operator()(uint16_t &dst)
            {
                dst = (uint16_t)++a;
            }
        } rnd = { 0 };
        auto img = image<uint16_t>(3104, 4192).generate(rnd);
        img = median7x7(img, SYMMETRIC);
    }
#endif

#if 0
    SCENARIO("LSC Shading correction performance and comparison test")
    {
        auto img = image<uint16_t>(3104, 4192).generate(rnd).read("E:\\temp\\1044_4192x3104.raw");
        img = median7x7(img, SYMMETRIC);
        gaussian35x35(img, REPLICATE);
        img.write("E:\\processed.raw");

        /* Matlab code to verify -- cut and paste
            a = fopen('E:\temp\1044_4192x3104.raw'); b0=reshape(fread(a, 'short'),[4192 3104])'; fclose(a);
            tic;b1 = medfilt2(b0, [7 7], 'symmetric');
            h=fspecial('gaussian',[35 35], 12);
            b2=imfilter(b1,h,'replicate', 'same'); toc  %% gives ~ 6.7 seconds
            a = fopen('E:\processed.raw'); c0=reshape(fread(a, 'short'),[4192 3104])'; fclose(a);
            isequal(c0, round(b2))  % gives '1'
        */
    }
#endif

#if 0
    // ~0.03s per 1668 items vs 0.3 s on matlab vs 0.6 s for combined regionprops(imfill('holes'))
    SCENARIO("Feature extraction performance test")
    {
        GIVEN("An labeled image")
        {
            auto img = image<uint16_t>(3456, 4608).read("C:\\work\\label.raw");
            auto size = img.size();
            WHEN("The image is classified")
            {
                int labels = 0;
                for (int j = 0; j < size._y; j++)
                {
                    for (int i = 0; i < size._x; i++)
                    {
                        if (img(j, i) == labels + 1)
                        {
                            ++labels;
                            auto p = roi_point(i, j);
                            auto props = get_regionprops(img, p);
                            CHECK(props.area > 0);
                        }
                    }
                }
                THEN("There are an expected number of labels")
                {
                    CHECK(labels == 1668);
                }
            }
        }
    }
#endif
}
