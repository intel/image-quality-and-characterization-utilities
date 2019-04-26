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


#include "Teisko/Image/Conversion.hpp"
#include "catch.hpp"
#include <memory>
#include <random>
#include <algorithm>
#include <functional>
#include <vector>

using namespace Teisko;

template <typename T>
void fill_with_random_data(std::vector<T> &data, T lower = 1, T upper = 9)
{
    std::random_device random_device;
    std::mt19937 random_generator(random_device());
    std::uniform_int_distribution<T> distribution(lower, upper);
    auto gen = std::bind(distribution, random_generator);
    std::generate(data.begin(), data.end(), gen);
}

template <typename T>
image<T> random_image(int height, int width, int bits)
{
    auto result = image<T>(height, width);
    std::random_device random_device;
    std::mt19937 random_generator(random_device());
    std::uniform_int_distribution<T> distribution(0, (1 << bits) - 1);
    auto gen = std::bind(distribution, random_generator);
    std::generate(result._begin, result._begin + height * width, gen);
    return result;
}

// This is reference implementation of the Chroma Upsample -- as copied from SC_ChromaUpsample.cpp
struct CSC_ChromaUpsample
{
    int m_UvWidth;
    int m_UvHeight;
    int m_P_EE[2][2];
    int m_P_OE[2][2];
    int m_P_EO[2][2];
    int m_P_OO[2][2];

    std::vector<int> input_data;
    std::vector<int> output_data;
    std::vector<int *> rows_input;
    std::vector<int *> rows_output;

    int lin_mul_444(int** arr, int P[2][2], int x, int y){
        int tmp = 0;

        for (int k_x = 0; k_x < 2; ++k_x){
            for (int k_y = 0; k_y < 2; ++k_y){
                tmp += arr[height_check(y + k_y)][width_check(x + k_x)] * P[k_x][k_y];
            }
        }
        return tmp;
    }
    int width_check(int x) { return std::min(std::max(x, 0), m_UvWidth - 1); }
    int height_check(int x){ return std::min(std::max(x, 0), m_UvHeight - 1); }
    void chromaUpsample_Reference(int **pUout, int **pUin)
    {
        for (int x = 0; x < m_UvWidth; ++x) {
            for (int y = 0; y < m_UvHeight; ++y) {
                pUout[2 * y][2 * x] = (lin_mul_444(pUin, m_P_EE, x - 1, y - 1) + (1 << 3)) >> 4;
                pUout[2 * y + 1][2 * x] = (lin_mul_444(pUin, m_P_OE, x - 1, y) + (1 << 3)) >> 4;
                pUout[2 * y][2 * x + 1] = (lin_mul_444(pUin, m_P_EO, x, y - 1) + (1 << 3)) >> 4;
                pUout[2 * y + 1][2 * x + 1] = (lin_mul_444(pUin, m_P_OO, x, y) + (1 << 3)) >> 4;
            }
        }
    }

    static void init_array(int(&m_P)[2][2], int a, int b, int c, int d)
    {
        m_P[0][0] = a;
        m_P[0][1] = b;
        m_P[1][0] = c;
        m_P[1][1] = d;
    }

    void make_ate_image(int width, int height, std::vector<int> &data, std::vector<int *> &rows)
    {
        data.resize(width * height);
        for (int i = 0; i < height; i++)
            rows.push_back(data.data() + width * i);
    }

    CSC_ChromaUpsample(int width, int height) : m_UvWidth(width), m_UvHeight(height)
    {
        init_array(m_P_EE, 1, 3, 3, 9);
        init_array(m_P_EO, 3, 9, 1, 3);
        init_array(m_P_OE, 3, 1, 9, 3);
        init_array(m_P_OO, 9, 3, 3, 1);
        make_ate_image(width, height, input_data, rows_input);
        make_ate_image(width * 2, height * 2, output_data, rows_output);
    }

    void upsample()
    {
        chromaUpsample_Reference(rows_output.data(), rows_input.data());
    }
};

// Reference implementation
struct CSC_RGBConversion
{
    int m_nA11;
    int m_nA12;
    int m_nA13;
    int m_nA21;
    int m_nA22;
    int m_nA23;
    int m_nA31;
    int m_nA32;
    int m_nA33;

    int m_nB1;
    int m_nB2;
    int m_nB3;

    int m_Bit_Precision;
    int m_MAX_RGB;
    int m_MIN_RGB;

    CSC_RGBConversion(int *RGB2YUV, int *bias, int pixel_resolution)
        : m_nA11(RGB2YUV[0])
        , m_nA12(RGB2YUV[1])
        , m_nA13(RGB2YUV[2])
        , m_nA21(RGB2YUV[3])
        , m_nA22(RGB2YUV[4])
        , m_nA23(RGB2YUV[5])
        , m_nA31(RGB2YUV[6])
        , m_nA32(RGB2YUV[7])
        , m_nA33(RGB2YUV[8])
        , m_nB1(bias[0])
        , m_nB2(bias[1])
        , m_nB3(bias[2])
        , m_Bit_Precision(pixel_resolution)
        , m_MAX_RGB((1 << pixel_resolution) - 1)
        , m_MIN_RGB(0)
    {}

    void Yuv2Rgb(int inY, int inU, int inV, int* outR, int* outG, int* outB)
    {
        const int ROUND_AFTER_MULT = 128;
        const int SHIFT_AFTER_MULT = 8;
        inY = inY - m_nB1;
        inU = inU - m_nB2;
        inV = inV - m_nB3;

        *outR = (inY * m_nA11 + inU * m_nA12 + inV * m_nA13 + ROUND_AFTER_MULT) >> SHIFT_AFTER_MULT;
        *outG = (inY * m_nA21 + inU * m_nA22 + inV * m_nA23 + ROUND_AFTER_MULT) >> SHIFT_AFTER_MULT;
        *outB = (inY * m_nA31 + inU * m_nA32 + inV * m_nA33 + ROUND_AFTER_MULT) >> SHIFT_AFTER_MULT;

        // Cliping
        *outR = (*outR < m_MIN_RGB)*m_MIN_RGB + (*outR > m_MAX_RGB)*m_MAX_RGB + ((*outR <= m_MAX_RGB) && (*outR >= m_MIN_RGB))* (*outR);
        *outG = (*outG < m_MIN_RGB)*m_MIN_RGB + (*outG > m_MAX_RGB)*m_MAX_RGB + ((*outG <= m_MAX_RGB) && (*outG >= m_MIN_RGB))* (*outG);
        *outB = (*outB < m_MIN_RGB)*m_MIN_RGB + (*outB > m_MAX_RGB)*m_MAX_RGB + ((*outB <= m_MAX_RGB) && (*outB >= m_MIN_RGB))* (*outB);
    }

    void convert_vector(std::vector<int> &inout)
    {
        size_t n = inout.size();
        int *data = inout.data();
        for (size_t i = 0; i + 2 < n; i+=3)
        {
            Yuv2Rgb(data[i], data[i+1], data[i+2], data + i, data + i + 1, data + i + 2);
        }
    }
};

SCENARIO("Libimage can upscale (chroma) channels by factor of two")
{
    GIVEN("A 2x2 Chroma Image")
    {
        int arr[4] = {
            5, 9,
            3, 8
        };
        auto img = image<int>(2, 2, arr);
        WHEN("The image is upscaled")
        {
            auto upscaled = chroma_upscale(img);
            THEN("The upscaled image matches the reference data")
            {
                auto ref = std::vector<int>({
                    5, 6, 8, 9,
                    5, 6, 8, 9,
                    4, 5, 7, 8,
                    3, 4, 7, 8
                });
                CHECK(upscaled._height == 4);
                CHECK(upscaled._width == 4);
                auto foo = upscaled.to_vector();
                CHECK(foo == ref);
            }
        }
    }
}

SCENARIO("Upsampling comparison test - simple")
{
    GIVEN("Original image of 59x31 stored in CSC_ChromaUpsample class")
    {
        const int max_value = 0xffffffffU >> 4;
        const int height = 31;
        const int width = 59;
        auto original = CSC_ChromaUpsample(width, height);
        fill_with_random_data(original.input_data, -max_value, max_value);

        WHEN("The image is upscaled with the reference method -- stored within the class")
        {
            original.upsample();
            AND_WHEN("The image is upscaled with the Teisko library method -- returning a new image")
            {
                auto view = image<int>(height, width, original.input_data.data());
                auto teisko_result = chroma_upscale(view);
                THEN("The upscaled reference image matches with Teisko implementation")
                {
                    auto &expected = original.output_data;
                    auto actual = teisko_result.to_vector();
                    CHECK(actual == expected);
                }
            }
        }
    }
}

SCENARIO("Libimage can convert images equally to reference implementation using 8-bit matrices. "
    "The original implementation uses int32 as the elementary type due to dependency to ATE format. "
    "The Teisko implementation uses uint16_t type as this is used internally in IQS")
{
    int biases_reference[3] = { 0, 1024, 1024 };
    int coefficients_8_bit[9] = { 256, 0, 403, 256, -48, -120, 256, 475, 0 };
    int bits = 11;
    double coeffs_new[9] = {
        coefficients_8_bit[0] / 256.0,
        coefficients_8_bit[1] / 256.0,
        coefficients_8_bit[2] / 256.0,
        coefficients_8_bit[3] / 256.0,
        coefficients_8_bit[4] / 256.0,
        coefficients_8_bit[5] / 256.0,
        coefficients_8_bit[6] / 256.0,
        coefficients_8_bit[7] / 256.0,
        coefficients_8_bit[8] / 256.0
    };

    CSC_RGBConversion original(coefficients_8_bit, biases_reference, bits);

    GIVEN("A small interleaved image of 51x79 pixels")
    {
        const int width = 51;
        const int height = 79;
        auto image = std::vector<uint16_t>(width * height * 3);
        fill_with_random_data(image, (uint16_t)0, (uint16_t)((1 << bits) - 1));  // uniform random data
        WHEN("The image data is converted (in-place) with the original implementation")
        {
            auto image_as_int = std::vector<int>(image.begin(), image.end());
            original.convert_vector(image_as_int);

            AND_WHEN("The original image is converted in-place with Teisko library")
            {
                auto view = Teisko::image<uint16_t>(height, width * 3, image.data());
                yuv_to_rgb_interleaved(view, bits, coeffs_new);
                THEN("The results match exactly")
                {
                    auto own_as_int = std::vector<int>(image.begin(), image.end());
                    CHECK(image_as_int == own_as_int);
                }
            }
        }
    }
}

SCENARIO("All SIMD/Non simd chroma upscaling variants produce equal results")
{
    GIVEN("An array of input resolutions and bitdepth")
    {
        auto widths = std::vector<int>({ 1, 2, 3, 4, 5, 6, 7, 8, 10, 32, 55 });
        auto heights = std::vector<int>({ 1, 2, 3, 4 });
        auto bit_depths = std::vector<int>({ 8, 12, 16 });
        WHEN("All scalers are exercised with the given dimensions, the outputs match")
        {
            for (auto bits : bit_depths)
            {
                for (auto width : widths)
                {
                    for (auto height : heights)
                    {
                        auto raw_data = std::vector<uint16_t>(width*height);
                        fill_with_random_data(raw_data, (uint16_t)0, (uint16_t)((1 << bits) - 1));

                        auto original = CSC_ChromaUpsample(width, height);
                        for (auto i = 0; i < width*height; i++)
                            original.input_data[i] = raw_data[i];
                        original.upsample();

                        auto img = image<uint16_t>(height, width, raw_data.data());
                        auto reference = chroma_upscale(img).to_vector();
                        auto fast = chroma_upscale_asm(img, bits).to_vector();
                        auto fast_32 = chroma_upscale_asm_32(img).to_vector();
                        auto ref_as_int = std::vector<int>(reference.begin(), reference.end());
                        CHECK(reference == fast);
                        CHECK(reference == fast_32);
                        CHECK(ref_as_int == original.output_data);
                    }
                }
            }
        }
    }
}

SCENARIO("RGB inline assembly color conversion produces equal results to reference version")
{
    static const double ycc2rgb_bt709[9] = {
        1.00000000000000000000, 0.00000000000000000000, 1.57472635664535240000,
        0.99999999999999989000, -0.18728134594285878000, -0.46819459633465549000,
        1.00000000000000000000, 1.85563960703714900000, 0.00000000000000000000
    };

    GIVEN("Small contiguous original image filled with random data")
    {
        int bits = 11;
        const int width = 157;
        const int height = 95;
        auto image_data = std::vector<uint16_t>(width * height * 3);
        fill_with_random_data(image_data, (uint16_t)0, (uint16_t)((1 << bits) - 1));  // uniform random data
        auto image_data_copy = image_data;

        WHEN("The image is converted with inline assembler")
        {
            auto red = image<uint16_t>(height, width, image_data.data());
            auto green = image<uint16_t>(height, width, image_data.data() + height * width);
            auto blue = image<uint16_t>(height, width, image_data.data() + 2 * height * width);
            yuv_to_rgb_planar_asm(red, green, blue, bits, ycc2rgb_bt709);
            AND_WHEN("The copy is converted using the reference method")
            {
                auto r = image<uint16_t>(height, width, image_data_copy.data());
                auto g = image<uint16_t>(height, width, image_data_copy.data() + height * width);
                auto b = image<uint16_t>(height, width, image_data_copy.data() + 2 * height * width);
                yuv_to_rgb_planar(r, g, b, bits, ycc2rgb_bt709);
                THEN("All the image planes are equal")
                {
                    auto rv0 = red.to_vector();
                    auto rv1 = r.to_vector();
                    auto gv0 = green.to_vector();
                    auto gv1 = g.to_vector();
                    auto bv0 = blue.to_vector();
                    auto bv1 = b.to_vector();
                    CHECK(rv0 == rv1);
                    CHECK(gv0 == gv1);
                    CHECK(bv0 == bv1);
                }
            }
        }
    }
}

// Current statistics
// -- RGB inline assembler:  21.5 x speed vs SC_RGBConversion           (5.4x on GCC)
//    RGB inline assembler:  15.5 x speed vs Teisko RGB conversion      (6.8x on GCC)
//    RGB Teisko conversion:  1.4 x speed vs SC_RGBConversion           (0.8x on GCC)
// -- Chroma Upscale inline  14.5 x speed vs SC_ChromaUpscale           (19x on GCC)
//    Chroma Upscale inline   2.4 x speed vs Teisko no SIMD             (2.3x on GCC)
//    Chroma Upscale noSIMD   5.9 x speed vs SC_ChromaUpscale           (8.3x on GCC)

// #define PERFORMANCE_COMPARISON
#ifdef PERFORMANCE_COMPARISON

SCENARIO("rgb conversion inline assembly performance test 12MP")
{
    static const double ycc2rgb_bt709[9] = {
        1.00000000000000000000, 0.00000000000000000000, 1.57472635664535240000,
        0.99999999999999989000, -0.18728134594285878000, -0.46819459633465549000,
        1.00000000000000000000, 1.85563960703714900000, 0.00000000000000000000
    };

    GIVEN("Original image of 4000x3000")
    {
        int bits = 11;
        const int width = 4000;
        const int height = 3000;
        auto image_data = std::vector<uint16_t>(width * height * 3);
        fill_with_random_data(image_data, (uint16_t)0, (uint16_t)((1 << bits) - 1));

        auto red = image<uint16_t>(height, width, image_data.data());
        auto green = image<uint16_t>(height, width, image_data.data() + height * width);
        auto blue = image<uint16_t>(height, width, image_data.data() + 2 * height * width);

        THEN("The image is converted using SIMD instructions")  // This is about 0.01 seconds
        {
            yuv_to_rgb_planar_asm(red, green, blue, bits, ycc2rgb_bt709);
        }
    }
}

SCENARIO("Reference implementation rgb converts 12MP image performance test")
{
    const int bits = 10;
    int biases_reference[3] = { 0, 512, 512 };
    int coefficients[9] = { 256, 0, 403, 256, -48, -120, 256, 475, 0 };
    CSC_RGBConversion ref_impl(coefficients, biases_reference, bits);

    GIVEN("Original image of 4000x3000")
    {
        const int width = 4000;
        const int height = 3000;
        auto image_as_int = std::vector<int>(width * height * 3);
        fill_with_random_data(image_as_int, 0, ((1 << bits) - 1));

        THEN("We convert the image in-place - reference RGB")   // This is about 0.2 seconds
        {
            ref_impl.convert_vector(image_as_int);
        }
    }
}

SCENARIO("Teisko library rgb converts 12MP image performance test")
{
    const int bits = 10;
    double coeffs_new[9] = {
        1.000000000000000, 0, 1.574218750000000,
        1.000000000000000, -0.187500000000000, -0.468750000000000,
        1.000000000000000, 1.855468750000000, 0 };

    GIVEN("Original image of 4000x3000")
    {
        const int width = 4000;
        const int height = 3000;
        auto image_as_short = std::vector<uint16_t>(width * height * 3);
        fill_with_random_data(image_as_short, (uint16_t)0, (uint16_t)((1 << bits) - 1));

        THEN("We convert the interleaved image in-place")   // This is about 0.125 seconds (~40% speedup)
        {
            auto test_img = image<uint16_t>(height, width * 3, image_as_short.data());
            yuv_to_rgb_interleaved(test_img, bits, coeffs_new);
        }
    }
}

SCENARIO("Teisko library rgb converts 12MP planar image performance test")
{
    const int bits = 10;
    double coeffs_new[9] = {
        1.000000000000000, 0, 1.574218750000000,
        1.000000000000000, -0.187500000000000, -0.468750000000000,
        1.000000000000000, 1.855468750000000, 0 };

    GIVEN("Original image of 4000x3000")
    {
        const int width = 4000;
        const int height = 3000;
        auto image_as_short = std::vector<uint16_t>(width * height * 3);
        fill_with_random_data(image_as_short, (uint16_t)0, (uint16_t)((1 << bits) - 1));

        THEN("We convert the image in-place - Teisko RGB planar")   // about 0.11 seconds
        {
            auto y = image<uint16_t>(height, width, image_as_short.data());
            auto u = image<uint16_t>(height, width, image_as_short.data() + height * width);
            auto v = image<uint16_t>(height, width, image_as_short.data() + height * width * 2);
            yuv_to_rgb_planar(y,u,v,bits, coeffs_new);
        }
    }
}

SCENARIO("Upsampling performance test original")
{
    GIVEN("Original image of 2000x1500")
    {
        auto ref = CSC_ChromaUpsample(2000, 1500);
        THEN("The image is upscaled - reference")
        {
            ref.upsample();
            CHECK(ref.rows_output.size() == 1500 * 2);
        }
    }
}


SCENARIO("Upsampling performance test - inline asm")
{
    GIVEN("Original image of 2000x1500")
    {
        auto ref = random_image<uint16_t>(1500, 2000, 10);
        THEN("The image is upscaled - Teisko inline asm")
        {
            auto tmp = chroma_upscale_asm(ref,10);
            CHECK(tmp._height == 1500 * 2);
            CHECK(tmp._width == 2000 * 2);
        }
    }
}

SCENARIO("Upsampling performance test - Teisko reference")
{
    GIVEN("Original image of 2000x1500")
    {
        auto ref = random_image<uint16_t>(1500, 2000, 10);
        THEN("The image is upscaled - Teisko reference")
        {
            auto tmp = chroma_upscale(ref);
            CHECK(tmp._height == 1500 * 2);
            CHECK(tmp._width == 2000 * 2);
        }
    }
}

SCENARIO("Upsampling performance test - inline asm 32 bit")
{
    GIVEN("Original image of 2000x1500")
    {
        auto ref = random_image<uint16_t>(1500, 2000, 10);
        THEN("The image is upscaled - Teisko inline asm version 2")
        {
            auto tmp = chroma_upscale_asm_32(ref);
            CHECK(tmp._height == 1500 * 2);
            CHECK(tmp._width == 2000 * 2);
        }
    }
}
#endif
