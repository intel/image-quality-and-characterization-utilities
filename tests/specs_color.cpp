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

#include "Teisko/Color.hpp"
#include "catch.hpp"
#include <vector>

using namespace Teisko;

SCENARIO(
    "It is possible to construct RGB objects using any positive "
    "integral or decimal values in range [0, N] where N > 0 is "
    "the maximum allowed RGB value.")
{
    // Compiler also accepts negative values, but in this context
    // those are considered to be undefined behavior.
    GIVEN("some examples of valid RGB objects with integral values, "
        "where the last component (B) has the maximum allowed RGB "
        "value")
    {
        rgb<uint8_t>  rgb_8bit{ 100, 150, 255 };        // N = 2^8  - 1
        rgb<uint16_t> rgb_10bit{ 100, 150, 1023 };      // N = 2^10 - 1
        rgb<uint32_t> rgb_18bit{ 100, 150, 262143 };    // N = 2^18 - 1

        THEN("the developer can decide what accuracy is used to "
            "represent integer values as decimals ")
        {
            rgb<float>  rgb_8bit_{ rgb_8bit };
            rgb<double> rgb_10bit_{ rgb_10bit };
            rgb<float> rgb_18bit_{ rgb_18bit };

            AND_THEN(
                "the developer needs to apply the maximum allowed "
                "RGB value herself when normalizing RGB values to "
                "the range [0.0, 1.0]")
            {
                rgb_8bit_ /= (1 << 8) - 1;
                rgb_10bit_ /= (1 << 10) - 1;
                rgb_18bit_ /= (1 << 18) - 1;

                CHECK(rgb_8bit_.b == Approx(1.0));
                CHECK(rgb_10bit_.b == Approx(1.0));
                CHECK(rgb_18bit_.b == Approx(1.0));
            }
        }
    }
}

SCENARIO(
    "It is possible to convert between sRGB and XYZ colorspaces, "
    "where sRGB values either integers in range [0, 255] or decimals "
    "in range [0.0, 255.0], and normalized for XYZ conversion to [0, 1] "
    "and XYZ values are decimals in range [0.0 1.0].")
{
    GIVEN("a valid sRGB object with non-zero integer values")
    {
        rgb<uint8_t> original{ 100, 150, 90 };
        WHEN("the sRGB object is converted to decimal sRGB object")
        {
            rgb<float> original_decimal{ original };

            CHECK(original_decimal.r == Approx(100.0));
            CHECK(original_decimal.g == Approx(150.0));
            CHECK(original_decimal.b == Approx(90.0));

            AND_WHEN("XYZ object is constructed from the sRGB object "
                "that is normalized to [0.0, 1.0]")
            {
                xyz<float> xyz{ original_decimal / 255.0 };

                THEN("XYZ values should match precalculated values")
                {
                    CHECK(xyz.x == Approx(0.43577));
                    CHECK(xyz.y == Approx(0.52955));
                    CHECK(xyz.z == Approx(0.41310));
                }

                AND_WHEN("a new sRGB object is constructed from the XYZ object"
                    "and converted to [0, 255]")
                {
                    rgb<float> converted_decimal{ xyz };
                    rgb<uint8_t> converted{ converted_decimal * 255.0 };

                    THEN("original sRGB object should equal to the new sRGB object")
                    {
                        CHECK(original == converted);
                    }
                }
            }
        }
    }
}

SCENARIO(
    "It is possible to convert between sRGB and Lab colorspaces, "
    "where sRGB values either integers in range [0, 1023] or decimals "
    "in range [0.0, 1023.0], and normalized for XYZ conversion to [0, 1] "
    "and XYZ values are decimals in range [0.0 1.0].")
{
    GIVEN("a valid sRGB object with non-zero integer values")
    {
        rgb<uint16_t> original{ 400, 600, 360 };
        WHEN("the sRGB object is converted to decimal sRGB object")
        {
            rgb<float> original_decimal{ original };

            CHECK(original_decimal.r == Approx(400.0));
            CHECK(original_decimal.g == Approx(600.0));
            CHECK(original_decimal.b == Approx(360.0));

            AND_WHEN("Lab object is constructed from the sRGB object "
                "that is normalized to [0.0, 1.0]")
            {
                lab<float> lab{ original_decimal / 1023.0 };

                THEN("Lab values should match precalculated values")
                {
                    CHECK(lab.l == Approx(77.75688));
                    CHECK(lab.a == Approx(-18.95491));
                    CHECK(lab.b == Approx(17.00504));
                }

                AND_WHEN("a new sRGB object is constructed from the XYZ object"
                    "and converted to [0, 1023]")
                {
                    rgb<float> converted_decimal{ lab };
                    rgb<uint16_t> converted{ converted_decimal * 1023.0 };

                    THEN("original sRGB object should equal to the new sRGB object")
                    {
                        CHECK(original == converted);
                    }
                }
            }
        }
    }
}

SCENARIO(
    "It is possible to convert between sRGB and BT2020 colorspaces, "
    "where sRGB and BT2020 are defined under different illuminant (whitepoint). "
    "In addition both RGB color spaces are in range [0, 255]. "
    "When converting between different RGB colorspaces no normalization is needed.")
{
    GIVEN("a valid sRGB object with non-zero integer values")
    {
        rgb<int16_t> original{ 100, 150, 90 };
        WHEN("BT2100 object is constructed from the sRGB object")
        {
            rgb<int16_t, rgb_cs::BT2100, wp::D50> converted{ original };
            THEN("RGB values match to precalculated values")
            {
                rgb<int16_t, rgb_cs::BT2100, wp::D50> rgb_bt2100(128, 144, 72);
                CHECK(converted == rgb_bt2100);
            }
        }
    }
}

SCENARIO(
    "Use ColorCheckerClassic Lab D50 value conversion to sRGB as a "
    "regression test to validate internal calculations. "
    "Test converts Lab D50 to both linear sRGB and sRGB first to decimal "
    "values in range [0.0, 1.0] and then to integral value range [0, 255]. "
    "Finally the converted values are compared againts reference values.")
{
    std::vector<rgb<int16_t>> linear_srgb_ref = {
        {  44,  20,  13 },
        { 143,  71,  54 },
        {  27,  48,  84 },
        {  27,  38,  13 },
        {  58,  54, 109 },
        {  29, 129, 105 },
        { 190,  52,   8 },
        {  15,  26,  99 },
        { 143,  20,  29 },
        {  28,  11,  35 },
        {  85, 127,  11 },
        { 197,  91,   5 },
        {   5,  12,  72 },
        {  12,  74,  16 },
        { 114,   9,  10 },
        { 215, 146,   1 },
        { 133,  20,  74 },
        { -10,  60,  96 }, // cyan patch outside srgb gamut
        { 224, 226, 213 },
        { 149, 151, 149 },
        {  91,  94,  93 },
        {  48,  49,  48 },
        {  22,  23,  23 },
        {   8,   8,   8 },
    };

    std::vector<rgb<uint16_t>> srgb_ref = {
        { 116,  79,  65 },
        { 197, 144, 127 },
        {  91, 120, 155 },
        {  91, 108,  64 },
        { 131, 127, 175 },
        {  95, 189, 172 },
        { 224, 124,  48 },
        {  69,  90, 167 },
        { 197,  80,  95 },
        {  93,  58, 104 },
        { 156, 187,  58 },
        { 227, 161,  39 },
        {  40,  62, 145 },
        {  61, 147,  70 },
        { 178,  54,  57 },
        { 236, 199,  15 },
        { 191,  79, 146 },
        {   0, 133, 165 }, // cyan patch clipped to 0 in unsigned int
        { 241, 242, 235 },
        { 201, 202, 201 },
        { 161, 163, 163 },
        { 121, 121, 121 },
        {  83,  84,  85 },
        {  50,  50,  50 },
    };

    GIVEN("a valid Lab object with D50 white point")
    {
        std::vector<lab<float, wp::D50>> lab_ref = {
            { 37.54f,  14.37f,  14.92f },
            { 64.66f,  19.27f,  17.50f },
            { 49.32f,  -3.82f, -22.54f },
            { 43.46f, -12.74f,  22.72f },
            { 54.94f,   9.61f, -24.79f },
            { 70.48f, -32.26f,  -0.37f },
            { 62.73f,  35.83f,  56.50f },
            { 39.43f,  10.75f, -45.17f },
            { 50.57f,  48.64f,  16.67f },
            { 30.10f,  22.54f, -20.87f },
            { 71.77f, -24.13f,  58.19f },
            { 71.51f,  18.24f,  67.37f },
            { 28.37f,  15.42f, -49.80f },
            { 54.38f, -39.72f,  32.27f },
            { 42.43f,  51.05f,  28.62f },
            { 81.80f,   2.67f,  80.41f },
            { 50.63f,  51.28f, -14.12f },
            { 49.57f, -29.71f, -28.32f },
            { 95.19f,  -1.03f,  2.93f },
            { 81.29f,  -0.57f,  0.44f },
            { 66.89f,  -0.75f, -0.06f },
            { 50.76f,  -0.13f,  0.14f },
            { 35.63f,  -0.46f, -0.48f },
            { 20.64f,   0.07f, -0.46f },
        };

        WHEN("values are converted to linear sRGB in range [0,255]")
        {
            std::vector<rgb<int16_t>> linear_srgb;
            for (auto& item : lab_ref)
            {
                linear_srgb.push_back(rgb<float>(item) * 255.0);
            }
            THEN("linear sRGB should match precalculated values")
            {
                for (size_t i = 0; i < lab_ref.size(); ++i)
                {
                    CHECK(linear_srgb[i] == linear_srgb_ref[i]);
                }
            }
            AND_WHEN("lab values are converted to sRGB in range [0,255]")
            {
                std::vector<rgb<uint16_t>> srgb;
                for (auto& item : lab_ref)
                {
                    // Convert and apply gamma
                    srgb.push_back(gamma::set(rgb<float>(item)) * 255.0);
                }
                THEN("sRGB should match precalculated values")
                {
                    for (size_t i = 0; i < lab_ref.size(); ++i)
                    {
                        CHECK(srgb[i] == srgb_ref[i]);
                    }
                }
            }
        }
    }
}

SCENARIO(
    "It is possible to convert between linear RGB and nonlinear RGB values. "
    "All RGB values should be normalized to range [0, 1.0]. Gamma is not "
    "defined outside this range (however there will be output for these values). "
    "Gamma can only be calculated for floating point types, integral types are not "
    "supported."
)
{
    GIVEN("valid sRGB object in range [0.0, 255.0]")
    {
        rgb<double> rgb{ 100, 150, 120 };
        //rgb<uint16_t> b{ 10, 20, 30 };

        WHEN("linear sRGB object is normalized to range [0,1]")
        {
            // Normalize
            auto max_value = 255;
            rgb /= max_value;

            THEN("object is converted to non-linear")
            {
                gamma::set(rgb);
                //gamma::set(b); // Does not compile, integer type is not supported

                AND_THEN(
                    "object is converted back to linear and normalized to original range")
                {
                    gamma::unset(rgb);
                    //gamma::unset(b); // Does not compile, integer type is not supported

                    // Convert back to original range
                    rgb *= max_value;

                    AND_THEN("RGB values match to original values")
                    {
                        CHECK(rgb.r == Approx(100.0));
                        CHECK(rgb.g == Approx(150.0));
                        CHECK(rgb.b == Approx(120.0));
                    }
                }
            }
        }
    }
}

SCENARIO(
    "It is possible to convert between RGB and HCL values.")
{
    GIVEN("valid RGB object in range [0.0, 255.0]")
    {
        rgb<double, rgb_cs::Sensor, wp::None> rgb{ 0, 0.1, 0.2 };

        WHEN("rgb object is converted to hcl")
        {
            hcl<double> hcl(rgb);

            THEN("hcl values match to predefined ones")
            {
                CHECK(hcl.h == Approx(210.0));
                CHECK(hcl.c == Approx(0.1732));
                CHECK(hcl.l == Approx(0.1000));
            }
        }
    }

    GIVEN("valid HCL object")
    {
        hcl<double> hcl{ 100.0, 0.2, 0.3 };

        WHEN("hcl object is converted to rgb")
        {
            rgb<double, rgb_cs::Sensor, wp::None> rgb(hcl);

            THEN("rgb values match to predefined ones")
            {
                CHECK(rgb.r == Approx(0.27685));
                CHECK(rgb.g == Approx(0.42529));
                CHECK(rgb.b == Approx(0.19786));
            }
        }
    }
}

SCENARIO(
    "It is possible that the color of two objects look the same to one person but when we "
    "measure the color of the same objects, we can find slight differences between them. "
    "There is different equations to calculate color difference. Here we Test 3 different "
    "equations, which are deltaAB, delta E2000 and delta dc2000. For deltaAB and deltaE2000 "
    "we need difference in lightness, difference in red and green color spaces. "
    "and for calculating delta dc2000 we discard the lightness parameter ")
{
    GIVEN("some valid Lab objects with double values")
    {
        auto ref1_vector = std::vector<lab<double>>{
            {50, 2.6772, -79.7751},
            {50, 3.1571, -77.2803},
            {50, 2.8361, -74.0200},
            {50, -1.3802, -84.2814},
            {50, -1.1848, -84.8006},
            {50, -0.9009, -85.5211},
            {50, 0, 0},
            {50, -1, 2},
            {50, 2.4900, -0.0010},
            {50, 2.4900, -0.0010},
            {50, 2.4900, -0.0010},
            {50, 2.4900, -0.0010},
            {50, -0.0010, 2.4900},
            {50, -0.0010, 2.4900},
            {50, -0.0010, 2.4900},
            {50, 2.5000, 0},
            {50, 2.5000, 0},
            {50, 2.5000, 0},
            {50, 2.5000, 0},
            {50, 2.5000, 0},
            {50, 2.5000, 0},
            {50, 2.5000, 0},
            {50, 2.5000, 0},
            {50, 2.5000, 0},
            {60.2574, -34.0099, 36.2677},
            {63.0109, -31.0961, -5.8663},
            {61.2901, 3.7196, -5.3901},
            {35.0831, -44.1164, 3.7933},
            {22.7233, 20.0904, -46.6940},
            {36.4612, 47.8580, 18.3852},
            {90.8027, -2.0831, 1.4410},
            {90.9257, -0.5406, -0.9208},
            {6.7747, -0.2908, -2.4247},
            {2.0776, 0.0795, -1.1350}
        };

        auto ref2_vector = std::vector<lab<double>>{
            {50, 0, -82.7485},
            {50, 0, -82.7485},
            {50, 0, -82.7485},
            {50, 0, -82.7485},
            {50, 0, -82.7485},
            {50, 0, -82.7485},
            {50, -1, 2},
            {50, 0, 0},
            {50, -2.4900, 0.0009},
            {50, -2.4900, 0.0010},
            {50, -2.4900, 0.0011},
            {50, -2.4900, 0.0012},
            {50, 0.0009, -2.4900},
            {50, 0.0010, -2.4900},
            {50, 0.0011, -2.4900},
            {50, 0, -2.5},
            {73, 25, -18},
            {61, -5, 29},
            {56, -27, -3},
            {58, 24, 15},
            {50, 3.1736, 0.5854},
            {50, 3.2972, 0},
            {50, 1.8634, 0.5757},
            {50, 3.2592, 0.3350},
            {60.4626, -34.1751, 39.4387},
            {62.8187, -29.7946, -4.0864},
            {61.4294, 2.2480, -4.9620},
            {35.0232, -40.0716, 1.5901},
            {23.0331, 14.9730, -42.5619},
            {36.2715, 50.5065, 21.2231},
            {91.1528, -1.6435, 0.0447},
            {88.6381, -0.8985, -0.7239},
            {5.8714, -0.0985, -2.2286},
            {0.9033, -0.0636, -0.5514},
        };

        auto de2000_vector = std::vector<double>{
           2.04245,
           2.8615,
           3.4412,
           1.0000,
           1.0000,
           1.0000,
           2.3669,
           2.3669,
           7.1792,
           7.1792,
           7.2195,
           7.2195,
           4.8045,
           4.8045,
           4.7461,
           4.3065,
           27.1492,
           22.8977,
           31.9030,
           19.4535,
           1.0000,
           1.0000,
           1.0000,
           1.0000,
           1.2644,
           1.2630,
           1.8731,
           1.8645,
           2.0373,
           1.4146,
           1.4441,
           1.5381,
           0.6377,
           0.9082,

        };

        THEN("Color difference calculation by different equations."
        )
        {
            const static double e = 0.0001;
            std::vector<double> de2000_error;
            std::vector<double> dc2000_error;
            for (size_t i = 0; i < ref1_vector.size(); ++i)
            {
                de2000_error.push_back(color_diff<color_diff_type::DeltaE2000>::calculate(ref1_vector[i], ref2_vector[i]));
                dc2000_error.push_back(color_diff<color_diff_type::DeltaC2000>::calculate(ref1_vector[i], ref2_vector[i]));
            }

            AND_THEN("calculation for de2000")
            {
                for (size_t i = 0; i < de2000_error.size(); ++i)
                {
                    CHECK(de2000_error[i] == Approx(de2000_vector[i]).epsilon(e));
                }
            }
            AND_THEN("calculation for dc2000")
            {
                CHECK(dc2000_error[0] == Approx(2.0425).epsilon(e));
                CHECK(dc2000_error[1] == Approx(2.8615).epsilon(e));
                CHECK(dc2000_error[31] == Approx(0.5510).epsilon(e));
                CHECK(dc2000_error[32] == Approx(0.3280).epsilon(e));
                CHECK(dc2000_error[29] == Approx(1.4056).epsilon(e));
                CHECK(dc2000_error[17] == Approx(20.4310).epsilon(e));
            }
            AND_THEN("calculation for dab")
            {
                auto dAB_error = color_diff<color_diff_type::DeltaAB>::calculate(ref1_vector[0], ref2_vector[0]);
                CHECK(dAB_error == Approx(4.0011).epsilon(e));
            }
        }
    }
}