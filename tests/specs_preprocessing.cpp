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
#include "catch.hpp"
#include <memory>
#include <random>
#include <algorithm>
#include <vector>

using namespace Teisko;

namespace libprepro_algorithm_tests
{
    SCENARIO("Smoke Test of Initializing Black Level Model, Adding Vectors And Disposing")
    {
        CHECK_NOTHROW(black_level_model{}
            .set(1.0f, 0.03f, std::vector<float>{ 1.f, 2.f, 3.f })
            .set(1.0f, 0.01f, std::vector<float>{ 2.f, 3.f, 4.f }));
    }

    SCENARIO("How chromaticity is calculated")
    {
        GIVEN("One chromaticity calculator and four sensitivity values")
        {
            chromaticity_factory_f chroma;
            std::vector<double> sensitivities = { 1.0, 2.0, 3.0, 4.0 };

            WHEN("The factory is initialized with values corresponding to RGGB layout")
            {
                auto info = bayer_info_s(bayer_pattern_e::rggb);
                for (auto i : { 0, 1, 2, 3 })
                    chroma[info[i]] += sensitivities[i];

                chromaticity result = chroma;
                THEN("The chromaticity is R=1.0/2.5,B=4/2.5, I=0")
                {
                    CHECK(result._r_per_g == Approx(0.4).epsilon(1e-9));
                    CHECK(result._b_per_g == Approx(1.6).epsilon(1e-9));
                    CHECK(result._i_per_g == 0.0);
                }
            }

            AND_WHEN("The same values are interpreted as GBRI layout")
            {
                auto info = bayer_info_s(bayer_pattern_e::gbri);
                for (auto i : { 0, 1, 2, 3 })
                    chroma[info[i]] += sensitivities[i];

                chromaticity result = chroma;
                THEN("The chromaticity is B=2, R=3 and I=4")
                {
                    CHECK(result._b_per_g == 2.0);
                    CHECK(result._r_per_g == 3.0);
                    CHECK(result._i_per_g == 4.0);
                }
            }
        }
    }

    SCENARIO("The chromaticity of the image is calculated as the ratio of averaged R,B and I channels "
             "to the averaged G channels at the maximum sensitivity.")
    {
        GIVEN("A 2x2 raw image with four flat channels")
        {
            const int width = 70;
            const int height = 50;
            image<float> raw_2x2_image(height, width);
            auto sensitivities = std::vector<float>{ 2.0f, 3.0f, 4.0f, 6.0f };
            for (int ch = 0; ch < 4; ch++)
                raw_2x2_image.subview(2, 2, ch / 2, ch % 2).fill(sensitivities[ch]);

            WHEN("The chromaticity of the image is calculated as bayer order 'IBRG'")
            {
                bayer_image_s<float> img(raw_2x2_image, bayer_pattern_e::ibrg);
                auto chroma = get_flatfield_white_point(img);

                THEN("The chromaticity points for I,B,R should be 2/6, 3/6 and 4/6 respectively")
                {
                    CHECK(chroma._i_per_g == Approx(2.f / 6.f).epsilon(1e-7));
                    CHECK(chroma._b_per_g == Approx(3.f / 6.f).epsilon(1e-7));
                    CHECK(chroma._r_per_g == Approx(4.f / 6.f).epsilon(1e-7));
                }
            }

            AND_WHEN("The chromaticity of the image is calculated as bayer order 'GBRG'")
            {
                bayer_image_s<float> img(raw_2x2_image, bayer_pattern_e::gbrg);
                auto chroma = get_flatfield_white_point(img);
                THEN("The average of G is 4.0 and the chromaticity points for I,B and R should be 0, 3/4 and 4/4 respectively")
                {
                    CHECK(chroma._i_per_g == Approx(0.f / 4.f).epsilon(1e-7));
                    CHECK(chroma._b_per_g == Approx(3.f / 4.f).epsilon(1e-7));
                    CHECK(chroma._r_per_g == Approx(4.f / 4.f).epsilon(1e-7));
                }
            }
        }
    }

    SCENARIO("Libprepro can remove black level from images")
    {
        GIVEN("A 4x4 image and a preprocessor")
        {
            int img[4][4] =
            {
                { 2, 3, 4, 5 },
                { 5, 6, 7, 8 },
                { 0, 2, 4, 6 },
                { 7, 5, 3, 1 }
            };
            bayer_image_s<int> four_by_four(4, 4, &img[0][0], bayer_pattern_e::bggr);
            auto prepro = preprocessor_s{};

            WHEN("The Black level of 0,1,2,3 is removed the image")
            {
                std::vector<int> black_level{ 0, 1, 2, 3 };
                remove_black_level(four_by_four, black_level);
                THEN("The image should match the precalculated reference matrix with negative values clipped to zero")
                {
                    int A = black_level[0], B = black_level[1];
                    int C = black_level[2], D = black_level[3];
                    int img2[4][4] =
                    {
                        { 2 - A, 3 - B, 4 - A, 5 - B },
                        { 5 - C, 6 - D, 7 - C, 8 - D },
                        { 0 - A, 2 - B, 4 - A, 6 - B },
                        { 7 - C, 5 - D, 3 - C, 1 - D }
                    };
                    for (int i = 0; i < 4; i++)
                    {
                        for (int j = 0; j < 4; j++)
                        {
                            int ref = img2[i][j] < 0 ? 0 : img2[i][j];
                            CHECK(ref == img[i][j]);
                        }
                    }
                }
            }
        }
    }

    /// \brief make_simple_rgb_ir_struct    Generates rgb contamination removal grid with constant data
    /// The R,G,B are amounts of IR contamination in each channel -- each spatially flat
    /// The grid_indices will be filled by an indirection to the red, blue and green grids
    rgb_ir_struct make_simple_rgb_ir_struct(bayer_info_s &info, float r, float g, float b)
    {
        rgb_ir_struct result(8, 8);
        result._grid_data.resize(8 * 8 * 3);
        for (int i = 0; i < 8*8; i++)
        {
            result._grid_data[i] = r;
            result._grid_data[i + 64] = g;
            result._grid_data[i + 128] = b;
        }
        result._grid_indices = std::vector<char>(16, -1);

        /// loop over channels in the pattern while populating the 4x4 grid
        for (auto i : info)
        {
            auto color = info[i];
            if (color != color_info_e::ir)
            {
                auto w = info.get_width();
                REQUIRE(((w == 2) || (w == 4)));
                auto x = i & (w - 1);
                auto y = i >> (w / 2);
                result._grid_indices[x + y * 4] = (char)color;
            }
        }
        return result;
    }

    SCENARIO("Libprepro can remove RGB IR contamination from images")
    {
        GIVEN("A flat contamination model for R, G and B channels in IBGR layout")
        {
            auto i = 0.f;   // I channels is not contaminated
            auto b = 0.1f;  // how b channel is contaminated by I channel
            auto r = 0.15f; // how r channel is contaminated by I channel
            auto g = 0.2f;  // how g channel is contaminated by I channel

            auto info = bayer_info_s{ bayer_pattern_e::ibrg };
            auto model = rgb_ir_contamination{};
            auto grid = make_simple_rgb_ir_struct(info, r, g, b);
            model.add_data_grid(grid);

            WHEN("An a 4x4 image is contaminated by adding the I channel with a weight to a nearby color channel")
            {
                float img[4][4] =
                {
                    { 2, 3, 4, 5 },  // 2  4  == IR channel
                    { 5, 6, 7, 8 },
                    { 0, 2, 4, 6 },  // 0  4  == IR channel
                    { 7, 5, 3, 1 }
                };

                // TODO -- it would be nice to make this happen as
                // I =  getchannel(img, info, IR);  R = getchannel(img, info, RED); ...
                // R = R + I * i; G = G + I * g; B = B + I * b;

                float contaminated[4][4] =
                {
                    { 2 + i * 2, 3 + b * 2, 4 + i * 4, 5 + b * 4},
                    { 5 + r * 2, 6 + g * 2, 7 + r * 4, 8 + g * 4 },
                    { 0 + i * 0, 2 + b * 0, 4 + i * 4, 6 + b * 4 },
                    { 7 + r * 0, 5 + g * 0, 3 + r * 4, 1 + g * 4 }
                };

                auto teisko_img = bayer_image_s<float>(4, 4, &contaminated[0][0], bayer_pattern_e::ibrg);

                THEN("The removal succeeds to recover original image")
                {
                    CHECK_NOTHROW(model.remove_ir_contamination(teisko_img));
                    for (int y = 0; y < 4; y++)
                    {
                        for (int x = 0; x < 4; x++)
                        {
                            CHECK(contaminated[y][x] == Approx(img[y][x]).epsilon(1e-7));
                        }
                    }
                }
            };
        }
    }

    /// Fills each channel in 'img' by the contaminated values
    void make_contamination(bayer_image_s<int> &img, int r, int g, int b, int i, rgb_ir_struct &data)
    {
        auto grid_size = data._grid_width * data._grid_height;

        REQUIRE(data._grid_data.size() == 3 * grid_size);

        r += (int)(i * data._grid_data[0 * grid_size]);
        g += (int)(i * data._grid_data[1 * grid_size]);
        b += (int)(i * data._grid_data[2 * grid_size]);

        for (auto channel_idx : img)
        {
            int color;
            switch (img._layout[channel_idx])
            {
            case color_info_e::blue: color = b; break;
            case color_info_e::red: color = r; break;
            case color_info_e::green: color = g; break;
            default:
                color = i;
            }
            img[channel_idx].fill(color);
        }

    }

    SCENARIO("Libprepro can remove RGB IR contamination from all types of 4x4 RGB-IR images")
    {
        const int red = 12345;
        const int green = 65243;
        const int blue = 31311;
        const int ir = 10000;   // this is divisible evenly with 0.1,0.2 and 0.3

        GIVEN("A flat bayer image generated for all 4x4 rgb ir patterns")
        {
            const int height = 60;
            const int width = 80;
            REQUIRE(height % 4 == 0);
            REQUIRE(width % 4 == 0);
            std::vector<bayer_image_s<int>> images;
            std::vector<rgb_ir_struct> contaminations;
            for (auto pattern : bayer_info_s::get_ir_4x4_patterns())
            {
                images.emplace_back(height, width, pattern);
                auto &img = images.back();
                contaminations.push_back(make_simple_rgb_ir_struct(img._layout, 0.3f, 0.2f, 0.1f));
                make_contamination(img, red, green, blue, ir, contaminations.back());
                // Here we want to show that img has been changed and it's non-zero
                CHECK(img._img.at(0, 0) != 0);
                CHECK(img._img.at(0, 0) != red);
                CHECK(img._img.at(0, 0) != green);
                CHECK(img._img.at(0, 0) != blue);
            }

            WHEN("The images are decontaminated with the corresponding grids")
            {
                REQUIRE(images.size() == 8);
                REQUIRE(contaminations.size() == 8);
                for (int i = 0; i < 8; i++)
                {
                    rgb_ir_contamination model = {};
                    model.add_data_grid(contaminations[i]);
                    model.remove_ir_contamination(images[i]);
                }
                THEN("Each of the images contain only the original pixels")
                {
                    for (auto &image : images)
                    {
                        for (auto ch : image)
                        {
                            int color;
                            switch (image._layout[ch])
                            {
                            case color_info_e::blue: color = blue; break;
                            case color_info_e::red: color = red; break;
                            case color_info_e::green: color = green; break;
                            case color_info_e::ir: color = ir; break;
                            }
                            auto pixels = image[ch].to_vector();
                            auto ref_pixels = std::vector<int>(width * height / 16, color);
                            CHECK(pixels == ref_pixels);
                        }
                    }
                }
            }
        }
    }


    SCENARIO("Bayer image channels can be accessed by index")
    {
        GIVEN("A bayer image of 4x4 pixels as grbg order")
        {
            std::vector<int> img{
                2, 3, 4, 5,
                5, 6, 7, 8,
                0, 2, 4, 6,
                7, 5, 3, 1
            };
            auto bayer_image = bayer_image_s<int>(4, 4, img.data(), bayer_pattern_e::grbg);

            WHEN("The zeroeth channel is extracted")
            {
                auto pixel_vec = bayer_image[0].to_vector();

                THEN("The extracted vector matches pre-coded sub channel")
                {
                    auto ref_vector = std::vector<int>{ 2, 4, 0, 4};
                    CHECK(ref_vector == pixel_vec);
                }
            }

            WHEN("The [first] red channel is extracted")
            {
                auto pixel_vec = bayer_image[color_info_e::red].to_vector();

                THEN("The extracted vector matches pre-coded sub channel 1")
                {
                    auto ref_vector = std::vector<int>{ 3, 5, 2, 6};
                    CHECK(ref_vector == pixel_vec);
                }
            }

            THEN("Trying to extract non-existent channel throws")
            {
                CHECK_THROWS(bayer_image[color_info_e::ir]);
            }
        }
    }

    SCENARIO("Bayer_image<T> takes deep copies on assignments and copy constructors")
    {
        GIVEN("A bayer image of 4x4 pixels as grbg order")
        {
            std::vector<int> img{
                2, 3, 4, 5,
                5, 6, 7, 8,
                0, 2, 4, 6,
                7, 5, 3, 1
            };
            auto bayer_image = bayer_image_s<int>(4, 4, img.data(), bayer_pattern_e::grbg);

            THEN("The bayer_image points to the vector `img`")
            {
                CHECK(&bayer_image._img(0, 0) == img.data());
            }

            WHEN("New instances are created by constructors")
            {
                bayer_image_s<int> copy_1 = bayer_image;
                bayer_image_s<int> copy_2(bayer_image);     // Same thing

                THEN("The all three bayer images have unique data")
                {
                    auto ptr_orig = &bayer_image._img(0, 0);
                    auto ptr_copy_1 = &copy_1._img(0, 0);
                    auto ptr_copy_2 = &copy_2._img(0, 0);

                    CHECK(ptr_orig != ptr_copy_1);
                    CHECK(ptr_copy_1 != ptr_copy_2);
                    CHECK(ptr_orig != ptr_copy_2);
                }
            }

            WHEN("A bayer_image is assigned to an existing variable")
            {
                bayer_image = bayer_image_s<int>(4, 4, img.data(), bayer_pattern_e::grbg);

                THEN("The bayer_image has taken a copy and no longer points to the vector `img`")
                {
                    CHECK(&bayer_image._img(0, 0) != img.data());
                }
            }
        }
    }

    SCENARIO("Bayer_image<T> assignments and constructors allows conversion between pixel types")
    {
        GIVEN("Original 32 bit image as bayer_image")
        {
            // the data is selected so, that 65537 can be converted to a float exactly
            std::vector<int> data{ 65537, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16 };
            auto bayer_image = bayer_image_s<int>(4, 4, data.data(), bayer_pattern_e::grbg);

            THEN("The image can be converted to other types")
            {
                bayer_image_s<float> f_image(bayer_image);
                CHECK(f_image._img(0, 0) == 65537.0f);
            }

            THEN("The higher precision image can be converted to lower precision discarding the top bits")
            {
                bayer_image_s<char> eight_bit_image(bayer_image);
                auto first_pixel_value = eight_bit_image._img(0, 0);
                auto reference_value = (char)(65537 & 0xff);
                CHECK(reference_value == 1);
                CHECK(first_pixel_value == reference_value);
            }
        }
    }

    SCENARIO("Preprocessor reduces IR sensors to the nearest non-IR sensor")
    {
        /// checks the differences in the top 2x2 corners of
        /// two sensors (one of them being 2x2)
        auto calculate_differences = [](bayer_pattern_e ir, bayer_pattern_e regular)
        {
            int diffs = 0;
            bayer_info_s ir_info(ir);
            bayer_info_s rg_info(regular);
            REQUIRE(rg_info.is_2x2_sensor());
            int ir_width = ir_info.get_width();
            for (int y = 0; y < 2; y++)
            {
                for (int x = 0; x < 2; x++)
                {
                    diffs += (rg_info[y * 2 + x] != ir_info[y*ir_width + x]);
                }
            }
            return diffs;
        };

        GIVEN("A preprocessor")
        {
            preprocessor_s prepro;

            WHEN("RGIB sensor is converted to regular 2x2 sensor")
            {
                auto converted = prepro.reduce_ir_layout(bayer_pattern_e::rgib);
                THEN("the result RGGB differs by one color channel")
                {
                    CHECK(converted == bayer_pattern_e::rggb);

                    auto original_layout = bayer_info_s(bayer_pattern_e::rgib);
                    auto converted_layout = bayer_info_s(converted);
                    CHECK(original_layout[0] == converted_layout[0]);
                    CHECK(original_layout[1] == converted_layout[1]);
                    CHECK(original_layout[2] != converted_layout[2]);
                    CHECK(original_layout[3] == converted_layout[3]);
                }
            }

            WHEN("Each 2x2 and 4x4 ir sensors are reduced to the closest non-IR sensor")
            {
                std::vector<bayer_pattern_e> ir_sensors;
                std::vector<bayer_pattern_e> results;

                for (auto p : bayer_info_s::get_ir_2x2_patterns())
                {
                    ir_sensors.push_back(p);
                    results.push_back(prepro.reduce_ir_layout(p));
                }

                for (auto p : bayer_info_s::get_ir_4x4_patterns())
                {
                    ir_sensors.push_back(p);
                    results.push_back(prepro.reduce_ir_layout(p));
                }

                THEN("The closest sensor found differs by one color channel")
                {
                    size_t i = 0;
                    for (auto &sensor: ir_sensors)
                    {
                        CHECK(1 == calculate_differences(sensor, results[i++]));
                    }
                }
            }
        }
    }

    SCENARIO("Calling reconstruction for non IR sensor succeeds without changing a thing")
    {
        GIVEN("A preprocessor class and image data")
        {
            auto prepro = preprocessor_s{};
            std::vector<int> img
            {
                2, 3, 4, 5,
                5, 6, 7, 8,
                0, 2, 4, 6,
                7, 5, 3, 1
            };
            auto bayer_image = bayer_image_s<int>(4, 4, &img[0], bayer_pattern_e::rggb);
            WHEN("The 4x4 image is reconstructed as if it's of RGGB")
            {
                auto copy = bayer_image;
                CHECK_NOTHROW(prepro.reconstruct_bayer_cfa(bayer_image));
                THEN("The reconstructed image matches the original")
                {
                    CHECK(bayer_image._img.to_vector() == copy._img.to_vector());
                }
            }
            WHEN("The 4x4 image is reconstructed as if it's of RGB-IR type")
            {
                auto copy = bayer_image;
                bayer_image._layout = bayer_info_s(bayer_pattern_e::irbg);
                CHECK_NOTHROW(prepro.reconstruct_bayer_cfa(bayer_image));
                THEN("The bayer image has been changed")
                {
                    CHECK(bayer_image._img.to_vector() != copy._img.to_vector());
                }
            }
        }
    }

    std::vector<bayer_pattern_e> get_all_patterns()
    {
        auto append = [](std::vector<bayer_pattern_e> &a, std::vector<bayer_pattern_e> b)
        {
            a.insert(std::end(a), std::begin(b), std::end(b));
        };
        auto result = bayer_info_s::get_regular_2x2_patterns();
        append(result, bayer_info_s::get_ir_2x2_patterns());
        append(result, bayer_info_s::get_dp_4x2_patterns());
        append(result, bayer_info_s::get_ir_4x4_patterns());
        return result;
    }

    SCENARIO("Preprocess function fn0 normalizes all bayer patterns to 2x2 layout while doing a saturation")
    {
        using type = int;
        const type max_level = 6;
        std::vector<type> img_data
        {
            9, 3, 4, 5,
            5, 6, 7, 8,
            0, 2, 4, 6,
            7, 5, 3, 1
        };
        image<int> img(4, 4, img_data.data());
        GIVEN("A preprocessor initialized with the saturation level and a bayer_image_s of every kind of pattern")
        {
            auto prepro = preprocessor_s{};
            prepro.saturation = static_cast<float>(max_level);

            auto original_patterns = get_all_patterns();
            std::vector<bayer_image_s<type>> images;
            for (auto p : original_patterns)
                images.push_back(bayer_image_s<type>(img, p));

            WHEN("The images are reduced with the preprocess_fn0 function")
            {
                for (auto &i : images)
                {
                    CHECK_NOTHROW(prepro.preprocess_fn0(i));
                }

                THEN("The image pattern is one of the first 4 basic patterns")
                {
                    for (auto &i : images)
                    {
                        CHECK(true == i._layout.is_2x2_sensor());
                    }

                    AND_THEN("no pixel in the image is above the given threshold")
                    {
                        for (auto &i : images)
                        {
                            CHECK(i._img.foreach(maximum_f<type>{}) <= max_level);
                        }
                        // AND_THEN we check that DP sensors have been reduced in size
                        // - we use booleans here to ensure that both branches are visited
                        CHECK(original_patterns.size() == images.size());
                        bool there_are_dp_sensors = false;
                        bool there_are_non_dp_sensors = false;
                        size_t i = 0;
                        for (auto &p : original_patterns)
                        {
                            auto image_width = images[i++]._img._width;
                            if (bayer_info_s(p).is_dp_4x2_sensor())
                            {
                                there_are_dp_sensors = true;
                                CHECK(image_width == 2);
                            }
                            else
                            {
                                there_are_non_dp_sensors = true;
                                CHECK(image_width == 4);
                            }
                        }
                        CHECK(there_are_dp_sensors == true);
                        CHECK(there_are_non_dp_sensors == true);
                    }
                }
            }
        }
    }

    SCENARIO("Preprocess function fn1 removes black level while strecthing data to maximum value")
    {
        // Use constant values for each color channel
        auto ch1 = 14.0f;
        auto ch2 = 100.0f;
        auto ch3 = 1023.0f;
        auto ch4 = 500.0f;
        std::vector<float> img_data
        {
            ch1, ch2, ch1, ch2,
            ch3, ch4, ch3, ch4,
            ch1, ch2, ch1, ch2,
            ch3, ch4, ch3, ch4,
        };
        auto bayer_image = bayer_image_s<float>(4, 4, &img_data[0], bayer_pattern_e::rggb);

        GIVEN("An image, a preprocessor and black level")
        {
            std::vector<float> black_level{ 16.0f, 16.0f, 16.0f, 16.0f };
            preprocessor_s prepro;
            prepro.saturation = 1023.0f;
            prepro.bl_model.set(1.0, 300.0f, black_level);

            std::vector<float> m; // multiplier
            for (size_t i = 0; i < black_level.size(); ++i)
                m.emplace_back(prepro.saturation / (prepro.saturation - black_level[i]));

            WHEN("The image is preprocessed using preprocess_fn1")
            {
                prepro.preprocess_fn1(bayer_image);

                THEN("The image should match the precalculated reference matrix with negative values clipped to zero")
                {
                    auto A = black_level[0], B = black_level[1];
                    auto C = black_level[2], D = black_level[3];
                    float ref_img_data[4][4] =
                    {
                        { (ch1 - A) * m[0], (ch2 - B) * m[1], (ch1 - A) * m[0], (ch2 - B) * m[1] },
                        { (ch3 - C) * m[2], (ch4 - D) * m[3], (ch3 - C) * m[2], (ch4 - D) * m[3] },
                        { (ch1 - A) * m[0], (ch2 - B) * m[1], (ch1 - A) * m[0], (ch2 - B) * m[1] },
                        { (ch3 - C) * m[2], (ch4 - D) * m[3], (ch3 - C) * m[2], (ch4 - D) * m[3] }
                    };

                    for (int i = 0; i < 4; i++)
                    {
                        for (int j = 0; j < 4; j++)
                        {
                            auto ref = ref_img_data[i][j] < 0.0f ? 0.0f : ref_img_data[i][j];

                            auto is_equal = [](float a, float b)
                            {
                                return abs(a - b) < 0.0001f;
                            };

                            CHECK(is_equal(ref, bayer_image._img.at(i, j)));
                        }
                    }
                }
            }
        }

    }

    /// Makes a test image with constant R,G,B and I as defined by the vector of values
    template <typename T>
    bayer_image_s<T> make_test_image(int height, int width, bayer_pattern_e pat, std::vector<T> values)
    {
        auto img = bayer_image_s<T>(height, width, pat);
        for (auto channel_idx : img)
        {
            auto color_idx = static_cast<int>(img._layout[channel_idx]);
            if (color_idx >= 4)
                throw std::runtime_error("Index out of range");
            img[channel_idx].fill(values[color_idx]);
        }
        return img;
    }

    SCENARIO("Preprocessor can linearize/compress images with a piecewise linear curve")
    {
        GIVEN("A linearization_model with 3 knee-points and a 4x4 rggb image")
        {
            auto linearizer = linearization_s();
            linearizer.grid_indices = std::vector<int8_t>(16, 0);
            linearizer.knee_points.emplace_back();
            auto &kneepts = linearizer.knee_points.back();

            kneepts.emplace_back(1, 4);    // maps 0..1->4, 2->8, 3->12, 4->16
            kneepts.emplace_back(5, 20);   // maps 5->20, 6->21
            kneepts.emplace_back(7, 22);   // maps 7...maxint to 22

            auto raw_data = std::vector<uint16_t>(
            {
                0,0,1,2,
                3,4,5,6,
                7,8,9,10,
                11,12,130,65535
            });
            auto bayer = bayer_image_s<uint16_t>(4, 4, raw_data.data(), bayer_pattern_e::rggb);
            WHEN("The bayer image is linearized")
            {
                linearizer.linearize(bayer);
                THEN("The original image has been expanded by the transfer function")
                {
                    auto ref_data = std::vector<uint16_t>(
                    {
                        4, 4, 4, 8,
                        12, 16, 20, 21,
                        22, 22, 22, 22,
                        22, 22, 22, 22
                    });
                    CHECK(ref_data == raw_data);
                }
            }
        }
    }


    SCENARIO("Preprocessor can reconstruct all 2x2 and 4x4 IR sensors changing the pattern to the most matching non-ir 2x2 sensor, "
        "where there are no samples from the IR pixels and all channels contain only the corresponding R,G or B values")
    {
        preprocessor_s prepro{};
        std::vector<int> flat_image_sensitivities = { 131, 222, -123, 144 };

        std::vector<bayer_image_s<int>> images;

        GIVEN("All ir patterns converted to corresponding images")
        {
            for (auto p : bayer_info_s::get_ir_2x2_patterns())
                images.push_back(make_test_image(32, 48, p, flat_image_sensitivities));

            for (auto p : bayer_info_s::get_ir_4x4_patterns())
                images.push_back(make_test_image(32, 48, p, flat_image_sensitivities));

            WHEN("The images are reconstructed")
            {
                for (auto &image : images)
                {
                    CHECK_NOTHROW(prepro.reconstruct_bayer_cfa(image));
                }
                ///  if we would just re-interpret the 4x4 IR sensors as 2x2 sensors,
                ///  then sub-sampling the "R" or "B" colors would not be constant
                ///  but a 50:50 mixture of both
                ///  B g R g B g R g...   -- every 'i' and 'g' channel would be ok
                ///  g i g i g i g i...      but a 'B' or 'G' channel would not conceptually exist in 2x2 grid
                THEN("The resulting images do not contain IR pixels or R,G,B pixels on wrong positions")
                {
                    for (auto &image : images)
                    {
                        CHECK(image._layout.is_2x2_sensor());
                        for (auto i : image)
                        {
                            auto pixels = image[i].to_vector();
                            auto color = static_cast<int>(image._layout[i]);
                            auto ref_vec = std::vector<int>(pixels.size(), flat_image_sensitivities[color]);
                            CHECK(ref_vec == pixels);
                        }
                    }
                }
            }
        }
    }

    SCENARIO("Preprocessor returns medians of selected regions of images")
    {
        GIVEN("An image of 80x60 with 4x2 layout, each COLOR having unique constant value")
        {
            const int width = 80;
            const int height = 60;
            float red = 10.0f;
            float green = 20.0f;
            float blue = 30.0f;
            std::vector<float> channel_sensitivities = { red, green, blue, 0.0f };
            auto img = make_test_image(height, width, bayer_pattern_e::bggr_4x2, channel_sensitivities);
            WHEN("We can extract the median of each channel")
            {
                std::vector<roi_point> top_left;
                std::vector<roi_point> bot_right;
                top_left.emplace_back(-10, 40);     // x is outside of image, y is inside the image
                bot_right.emplace_back(10, 999);    // bottom right coordinate is strictly larger than the top left corner
                auto result = get_patch_white_point(img, top_left, bot_right);
                THEN("The result matches the color order `B, B, G, G, G, G, R, R`")
                {
                    auto ref_medians = std::vector<float> { blue, blue, green, green, green, green, red, red };
                    CHECK(result == ref_medians);
                }
            }
        }
    }

    /// Makes a bayer image with constant r,g,g,b, (i) values
    /// demosaics the image for R,G,B (and optionally I)
    /// checks that each R,G,B, I channel is flat and contains just the original pixel value
    template <typename T> void run_demosaic_test(bayer_pattern_e pattern)
    {
        int width = 12;
        int height = 8;

        std::vector<T> flat_image_sensitivities = { 15, 20, 30, 40 };   // matches r,g,b and I color channels
        auto img = make_test_image(height, width, pattern, flat_image_sensitivities);
        auto channels = std::vector<color_info_e>({ color_info_e::red, color_info_e::green, color_info_e::blue });
        if (img._layout.is_ir_sensor())
            channels.emplace_back(color_info_e::ir);

        // Must reduce DP before calling demosaic
        // Also must change the reference value to compensate the fact that
        // DP - reduction ADDs the left and right channel instead of averaging
        if (img._layout.is_dp_4x2_sensor())
        {
            preprocessor_s().preprocess_reduce_dp_sensor(img);
            for (auto &value : flat_image_sensitivities)
                value *= 2;
            width /= 2;
        }

        for (auto &color_code : channels)
        {
            int color_code_index = static_cast<int>(color_code);
            auto demosaiced_channel = demosaic_bilinear(img, color_code, REPLICATE_EVEN);
            auto reference = std::vector<T>(width * height, flat_image_sensitivities[color_code_index]);
            CHECK(demosaiced_channel.to_vector() == reference);
        }
    }

    SCENARIO("Smoke test: demosaic recognizes all bayer patterns. "
        "We generate a W*H image of all patterns with constant R,G,B and possibly I. "
        "The demosaiced image is of original size but contains just this single R,G,B or I value. ")
    {
        GIVEN("All known sensor patterns")
        {
            auto patterns = get_all_patterns();
            THEN("The demosaiced channels are flat and contain just the single channel values R,G,B or I")
            {
                for (auto &p : patterns)
                {
                    run_demosaic_test<float>(p);
                    run_demosaic_test<uint8_t>(p);
                    run_demosaic_test<uint16_t>(p);
                }
            }
        }
    }

    SCENARIO("Preprocessor can demosaic bilinearly all four basic bayer patterns")
    {
        auto middle = [](image<uint16_t> &x)
        {
            return x.region(4, 4, 1, 1).to_vector();
        };

        GIVEN("An 6x6 image")
        {
            uint16_t orig[6 * 6] = {
                86, 34, 18, 48, 77, 93,
                57, 50, 23, 48, 38, 94,
                54, 39, 41, 33, 23, 56,
                14, 7, 4, 89, 39, 5,
                84, 23, 89, 36, 9, 23,
                61, 12, 93, 11, 13, 34
            };
            WHEN("The image is demosaiced for color order 'rggb'")
            {
                bayer_image_s<uint16_t> img(6, 6, orig, bayer_pattern_e::rggb);

                auto r = demosaic_bilinear(img, color_info_e::red, ZERO_PAD);
                auto g = demosaic_bilinear(img, color_info_e::green, ZERO_PAD);
                auto b = demosaic_bilinear(img, color_info_e::blue, ZERO_PAD);
                THEN("The image matches the matlab generated reference values using ZERO padding")
                {
                    std::vector<uint16_t> reference_red_on_rggb = {
                        86, 52, 18, 48, 77, 39,
                        70, 50, 30, 40, 50, 25,
                        54, 48, 41, 32, 23, 12,
                        69, 67, 65, 41, 16, 8,
                        84, 87, 89, 49, 9, 5,
                        42, 43, 45, 25, 5, 2
                    };

                    std::vector<uint16_t>  reference_green_on_rggb = {
                        23, 34, 26, 48, 45, 93,
                        57, 38, 23, 36, 38, 47,
                        28, 39, 25, 33, 42, 56,
                        14, 20, 4, 28, 39, 30,
                        25, 23, 39, 36, 28, 23,
                        61, 44, 93, 36, 13, 9
                    };

                    std::vector<uint16_t>  reference_blue_on_rggb = {
                        13, 25, 25, 24, 36, 47,
                        25, 50, 49, 48, 71, 94,
                        14, 29, 49, 69, 59, 50,
                        4, 7, 48, 89, 47, 5,
                        5, 10, 30, 50, 35, 20,
                        6, 12, 12, 11, 23, 34
                    };
                    CHECK(reference_red_on_rggb == r.to_vector());
                    CHECK(reference_green_on_rggb == g.to_vector());
                    CHECK(reference_blue_on_rggb == b.to_vector());
                }
            }

            WHEN("The image is demosaiced for color order 'grbg'")
            {
                bayer_image_s<uint16_t> img(6, 6, orig, bayer_pattern_e::grbg);
                auto r2 = demosaic_bilinear(img, color_info_e::red);
                auto g2 = demosaic_bilinear(img, color_info_e::green);
                auto b2 = demosaic_bilinear(img, color_info_e::blue);
                THEN("The image (center) matches the matlab generated reference values")
                {
                    std::vector<uint16_t> reference_red_on_grbg = {
                        37, 39, 41, 58,
                        39, 36, 33, 45,
                        31, 33, 35, 37,
                        23, 30, 36, 30
                    };
                    std::vector<uint16_t> reference_green_on_grbg = {
                        50, 39, 48, 61,
                        38, 41, 50, 23,
                        7, 57, 89, 32,
                        48, 89, 50, 9
                    };
                    std::vector<uint16_t> reference_blue_on_grbg = {
                        40, 23, 31, 38,
                        25, 14, 26, 39,
                        9, 4, 22, 39,
                        43, 49, 37, 26
                    };

                    CHECK(reference_red_on_grbg == middle(r2));
                    CHECK(reference_green_on_grbg == middle(g2));
                    CHECK(reference_blue_on_grbg == middle(b2));
                }
            }

            WHEN("The image is demosaiced for color order 'gbrg'")
            {
                bayer_image_s<uint16_t> img(6, 6, orig, bayer_pattern_e::gbrg);

                auto r3 = demosaic_bilinear(img, color_info_e::red);
                auto g3 = demosaic_bilinear(img, color_info_e::green);
                auto b3 = demosaic_bilinear(img, color_info_e::blue);
                THEN("The image (center) matches the matlab generated reference values")
                {
                    std::vector<uint16_t> reference_red_on_gbrg = {
                        40, 23, 31, 38,
                        25, 14, 26, 39,
                        9, 4, 22, 39,
                        43, 49, 37, 26
                    };
                    std::vector<uint16_t> reference_green_on_gbrg = {
                        50, 39, 48, 61,
                        38, 41, 50, 23,
                        7, 57, 89, 32,
                        48, 89, 50, 9
                    };
                    std::vector<uint16_t> reference_blue_on_gbrg = {
                        37, 39, 41, 58,
                        39, 36, 33, 45,
                        31, 33, 35, 37,
                        23, 30, 36, 30
                    };
                    CHECK(reference_red_on_gbrg == middle(r3));
                    CHECK(reference_green_on_gbrg == middle(g3));
                    CHECK(reference_blue_on_gbrg == middle(b3));
                }
            }

            WHEN("The image is demosaiced for color order 'bggr'")
            {
                bayer_image_s<uint16_t> img(6, 6, orig, bayer_pattern_e::bggr);

                auto r4 = demosaic_bilinear(img, color_info_e::red);
                auto g4 = demosaic_bilinear(img, color_info_e::green);
                auto b4 = demosaic_bilinear(img, color_info_e::blue);
                THEN("The image (center) matches the matlab generated reference values")
                {
                    std::vector<uint16_t> reference_red_on_bggr = {
                        50, 49, 48, 71,
                        29, 49, 69, 59,
                        7, 48, 89, 47,
                        10, 30, 50, 35
                    };
                    std::vector<uint16_t> reference_green_on_bggr = {
                        38, 23, 36, 38,
                        39, 25, 33, 42,
                        20, 4, 28, 39,
                        23, 39, 36, 28
                    };
                    std::vector<uint16_t> reference_blue_on_bggr = {
                        50, 30, 40, 50,
                        48, 41, 32, 23,
                        67, 65, 41, 16,
                        87, 89, 49, 9
                    };
                    CHECK(reference_red_on_bggr == middle(r4));
                    CHECK(reference_green_on_bggr == middle(g4));
                    CHECK(reference_blue_on_bggr == middle(b4));
                }
            }
        }
    }
}

#if 0
SCENARIO("End to End gamut processing from ccc image")
{
    GIVEN("A prep file from cr folder")
    {
        for (auto s : { "C:\\work\\ccc.raw" })
        {
            auto img = image<uint16_t>(3120, 4208).read(s);
            auto bayer = bayer_image_s<uint16_t>(img, bayer_pattern_e::rggb);

            auto foo = calculate_image_gamut(bayer);
            CHECK(foo.size() > 0);
        }
    }
}
#endif


#if 0
/// TODO: Fix the loading routine to not enforce exact size
//        read image db environment variable
SCENARIO("End to End Chromaticity processing from Flat Field image")
{
    GIVEN("A flat field file and a preprocessor")
    {
        preprocessor_s prepro;
        //auto str = "C:\\temp\\name#UnitTestData~rev#0~type#ff~lsrc#1~ciex#456~ciey#414~cct#2780~ag#256~exp#35000~msid#001~iid#a_ff.raw";
        auto str = "C:\\temp\\koe.raw"; // the actual image has additional data
        auto width = 1920;
        auto height = 1080;
        bayer_image_s<uint16_t> img(height, width, bayer_pattern_e::grbg);
        img._img.read(str);

        prepro.bl_model.set(1.0, 300.0f, std::vector<float>{ 16, 16, 16, 16 });
        WHEN("The file is preprocessed")
        {
            prepro.preprocess_fn3a(img);        // use default ag, exp

            // img._img.write("C:\\work\\outq.raw");
            auto result = prepro.get_flat_field_white_point(img);
            THEN("The chromaticity matches the matlab model")
            {
                REQUIRE(result.roi.size() > 0);
                CHECK(result._chromaticity._b_per_g == Approx(0.411917953720527).epsilon(1e-3));
                CHECK(result._chromaticity._r_per_g == Approx(1.029704404115448).epsilon(1e-3));
            }

        }
    }
}
#endif