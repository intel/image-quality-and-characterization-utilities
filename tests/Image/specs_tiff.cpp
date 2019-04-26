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

#include "Teisko/Image/TIFF.hpp"
#include "catch.hpp"
#include <memory>
#include <random>
#include <algorithm>
#include <vector>
#include <ostream>
#include <fstream>
#include <map>

#include "_data/icc_profiles.hpp"
using namespace Teisko;

SCENARIO("HLG Icc profile can be embedded into 10-bit image")
{
    GIVEN("Opened stream and tiff class")
    {
        tiff_file tiff;
        std::ifstream file;
        file.open("C:\\temp\\tiff\\bar.tif", std::ios::binary);
        THEN("The file can be manipulated")
        {
            tiff.read(file);
            tiff.add_tag(tiff_tag<UNDEFINED>(tag_icc_profile, hlg_icc, hlg_icc_len));
            std::ofstream file2;
            file2.open("C:\\temp\\tiff\\test_with_hlg_profile.tif", std::ios::binary);
            tiff.write(file2);
        }
    }
}

SCENARIO("We can instantiate different tiff tags with different interfaces")
{
    // Single value
    CHECK_NOTHROW(tiff_tag<BYTE>(tag_new_subfile_type, 1));
    // Single string
    CHECK_NOTHROW(tiff_tag<ASCII>(tag_software, "IQStudio Unit Tests"));
    // Multiple values
    CHECK_NOTHROW(tiff_tag<SHORT>(tag_bits_per_sample, { 8, 8, 8 }));
    // Single 'Long'
    CHECK_NOTHROW(tiff_tag<LONG>(tag_new_subfile_type, 0));
    // Single Rational as an initializer list
    CHECK_NOTHROW(tiff_tag<RATIONAL>(tag_date_and_time, { 1, 2 }));
    // arbitrary data as an array and length
    CHECK_NOTHROW(tiff_tag<UNDEFINED>(tag_icc_profile, hlg_icc, hlg_icc_len));
}

SCENARIO("Data written as short can be retrieved as other types")
{
    tiff_tag_base foo = tiff_tag<SHORT>(tag_image_width, { 123, 321 });
    tiff_tag_base bar = tiff_tag<ASCII>(tag_image_width, "123 321");
    tiff_tag_base rational = tiff_tag<RATIONAL>(tag_image_width, { 1, 2, 3, 5 });
    auto array_as_string = foo.to_string();
    auto string_as_string = bar.to_string();
    CHECK(array_as_string == string_as_string);
    CHECK(rational.to_string() == "1/2 3/5");
    auto food = tiff_tag<SBYTE>(tag_image_height, { -'a', 'b', 'c' }).to_string();
}

SCENARIO("We can deserialize all types of data from a vector")
{
    std::vector<uint8_t> test1{ 0, 1, 3, 0, 2, 0, 0, 0, 9, 0, 8, 0 };
    std::vector<uint8_t> test2{ 1, 0, 0, 3, 0, 0, 0, 2, 0, 9, 0, 8 };

    auto little_endian_tag = tiff_tag_base(test1, 0, false);
    auto big_endian_tag = tiff_tag_base(test2, 0, true);
    CHECK_THROWS(tiff_tag_base(test2, 0, false));
    CHECK_THROWS(tiff_tag_base(test1, 0, true));
    CHECK(little_endian_tag.to_string() == big_endian_tag.to_string());
    auto vec = little_endian_tag.value();
    auto ref_vec = std::vector<int>{9, 8};
    CHECK(vec == ref_vec);
}

SCENARIO("Tiff class can read and write images")
{
    GIVEN("Opened stream and tiff class")
    {
        tiff_file tiff;
        std::ifstream file;
        file.open("C:\\temp\\tiff\\ColorChecker_AdobeRGB_from_Lab_D50_AfterNov2014.tif",
            std::ios::binary);
        tiff.read(file);
        THEN("We can add icc tag")
        {
            tiff.add_tag(
                tiff_tag<UNDEFINED>(tag_icc_profile, AdobeRGB1998_icc, AdobeRGB1998_icc_len));
            AND_THEN("The tiff can be written to another file")
            {
                std::ofstream out_file;
                out_file.
                    open("C:\\temp\\tiff\\ColorChecker_AdobeRGB_from_Lab_D50_with_icc_profile.tif",
                    std::ios::binary);
                tiff.write(out_file);
            }
        }
    }
}

SCENARIO("A Basic Single Image Tiff Writer Work Flow")
{
    GIVEN("Image dimensions for a tiff writer class")
    {
        tiff_file tiff;
        const int width = 5;
        const int height = 3;
        const int input_bits_per_sample = 8;
        WHEN("A single (gray scale) image is populated with small values")
        {
            std::vector<int> values{ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };
            image<int> gray_scale_image(height, width, values.data());
            THEN("The image can be written to string stream using the tiff class")
            {
                std::ostringstream stream;
                REQUIRE_NOTHROW(
                    tiff.write_image(stream, gray_scale_image, input_bits_per_sample,
                    tag_photometric_interpretation_bayer));
                auto position = (long)stream.tellp();
                CHECK(position == 174);
                AND_THEN("We can verify that the some tag was written correctly")
                {
                    CHECK(tiff.tag_as_string(tag_photometric_interpretation) == "1");
                }
            }
        }
    }
}

SCENARIO("Example 10-bit RGB Image Tiff Writing Work Flow")
{
    const std::string description = "Custom Description";
    GIVEN("Image dimensions for a tiff writer class")
    {
        tiff_file tiff;

        const int width = 256 + 1024 + 256;     // Square + two margins on left & right
        const int height = 1024;
        const int input_bits_per_sample = 10;   // any bit depth from 1 to 32 is accepted

        THEN("We can embed a custom ASCII string to the tag section")
        {
            CHECK_NOTHROW(tiff.add_tag(tiff_tag<ASCII>(tag_image_description, description)));
            CHECK(tiff.tag_as_string(tag_image_description) == description);
            AND_WHEN("The image is populated with values between 0 and 1023")
            {
                image<float> rgb_interleaved(height, width * 3); // Everything interleaved

                auto R = rgb_interleaved.subview(1, 3, 0, 0);
                auto G = rgb_interleaved.subview(1, 3, 0, 1);
                auto B = rgb_interleaved.subview(1, 3, 0, 2);
                for (int j = 0; j < height; j++)
                {
                    for (int i = 0; i < width; i++)
                    {
                        unsigned int x = (unsigned int)(i - 256);
                        if (x < 1024)       // Color gradient in the center square
                        {
                            R.at(j, i) = (float)j;           // Red has vertical gradient
                            G.at(j, i) = (float)x;           // Green has horizontal gradient
                            B.at(j, i) = (float)(2048 - j - x) / 2; // Blue has diagonal gradient
                        }
                        else                // Monochromatic gradient at left and right borders
                        {
                            R.at(j, i) = (float)j;
                            G.at(j, i) = (float)j;
                            B.at(j, i) = (float)j;
                        }
                    }
                }
                THEN("The image can be written to file as interleaved")
                {
                    std::ofstream file;
                    file.open("bar.tif", std::ios::binary);
                    REQUIRE_NOTHROW(
                        tiff.write_image(file, rgb_interleaved, input_bits_per_sample,
                        tag_photometric_interpretation_rgb));
                }
                AND_THEN("The image can be written as three separate planes with HLG profile")
                {
                    // First we have to apply HLG transfer function
                    short *g = hlg_gamma;
                    rgb_interleaved.for_each([g](float &inout) {
                        inout = (float)g[(short)(inout + 0.5f)];
                    });

                    tiff.add_tag(tiff_tag<UNDEFINED>(tag_icc_profile, hlg_icc, hlg_icc_len));
                    std::ofstream file;
                    file.open("bar2.tif", std::ios::binary);
                    image<float> planes[3] = { R, G, B };
                    REQUIRE_NOTHROW(
                        tiff.write_image(file, planes, input_bits_per_sample,
                        tag_photometric_interpretation_rgb));
                }
            }
        }
    }
}

SCENARIO("Sample flow for YUV tiff writing")
{
    GIVEN("Three color planes with half width chroma subsampling")
    {
        const int chroma_width = 100;
        const int chroma_height = 80;
        const int bits_per_pixel = 11;
        image<int> planes[3] = {
            image<int>(chroma_height * 2, chroma_width * 2),
            image<int>(chroma_height, chroma_width),
            image<int>(chroma_height, chroma_width)
        };
        for (auto &p : planes)
            p.fill(1230);
        THEN("The image can be written as YUV")
        {
            tiff_file tiff;
            std::ofstream file;
            file.open("bar3.tif", std::ios::binary);
            tiff.write_image(file, planes, bits_per_pixel, tag_photometric_interpretation_yuv);
        }
    }
}

SCENARIO("Tiff tag supports ASCII strings")
{
    GIVEN("A small string of length 3")
    {
        const tiff_tags tag_type = (tiff_tags)0x1231;
        auto small_string = "abc";
        tiff_tag<ASCII> small(tag_type, small_string);
        THEN("The tag type matches the given one")
        {
            CHECK(tag_type == small.tag());
        }
        AND_THEN("Stored string matches the given one")
        {
            CHECK(small_string == small.to_string());
        }
    }
}

SCENARIO("Tiff tag class exposes serialization methods")
{
    GIVEN("A 'small' string to be written to an output stream")
    {
        tiff_tag<ASCII> small_string(tag_image_width, "123");
        std::ostringstream ss;
        const bool is_bigendian = false;
        WHEN("The tag is serialized")
        {
            small_string.write_pod(ss, is_bigendian);
            THEN("The output stream contains the predefined 12 bytes")
            {
                std::vector<uint8_t> reference{
                    0, 1,   // 256 == 0x0100 as little endian
                    2, 0,   // ASCII == 2
                    4, 0, 0, 0,  // length = 4
                    '1', '2', '3', '\0'     // String 123 with ASCIIZ
                };
                std::string ref_string(reference.begin(), reference.end());

                CHECK(ss.str() == ref_string);
                AND_THEN("There is no additional data to be serialized")
                {
                    small_string.write_additional_data(ss, is_bigendian);
                    CHECK(ss.str() == ref_string);
                }
            }
        }
    }

    GIVEN("A 'long' string to be written to an output stream")
    {
        tiff_tag<ASCII> image_height(tag_image_height, "1234");
        std::ostringstream ss;
        const bool is_bigendian = false;
        WHEN("The tag with a larger string is serialized")
        {
            image_height.write_pod(ss, is_bigendian);
            image_height.write_additional_data(ss, is_bigendian);
            THEN("The output stream matches the reference stream")
            {
                std::vector<uint8_t> reference{
                    1, 1,                   // 257 == 0x0101 as little endian
                    2, 0,                   // ASCII == 2
                    5, 0, 0, 0,             // length == 5
                    0, 0, 0, 0,             // offset == 0
                    '1', '2', '3', '4',     // String
                    0,                      // ASCIIZ
                    0                       // Padding
                };
                std::string ref_string(reference.begin(), reference.end());
                CHECK(ss.str() == ref_string);
            }
        }
    }
}
