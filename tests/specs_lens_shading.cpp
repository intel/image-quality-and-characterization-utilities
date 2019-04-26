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


#include "Teisko/LensShading.hpp"
#include "catch.hpp"
#include <memory>
#include <random>
#include <algorithm>
#include <vector>

using namespace Teisko;

namespace teisko_liblensshading_tests
{
    template <typename T>
    struct test_image_maker
    {
        roi_point image_size;
        std::vector<T> image_data;      // We allocate the image data just once and recycle it
        test_image_maker(int width, int height)
            : image_size(width, height)
            , image_data(width * height) { }

        void resize(roi_point new_size)
        {
            image_size = new_size;
            image_data.resize(image_size._x * image_size._y);
        }

        // Creates a W*H image with channels filled with initial values (default = zero)
        bayer_image_s<T> make_flat_image(bayer_pattern_e pattern, std::vector<T> sensitivities = {})
        {
            auto result = bayer_image_s<T>(image_size._y, image_size._x, image_data.data(), pattern);
            for (auto i : result)
                result[i].fill(i < sensitivities.size() ? sensitivities[i] : 0);
            return result;
        }

        // Creates a W*H image with highest sensitivity in the middle
        // sensitivities -- vector of initial values
        // falloff -- amount of shading per channel at the corners (0.0 ... 1.0)
        bayer_image_s<T> make_curved_image(bayer_pattern_e pattern,
            std::vector<T> sensitivities = {}, std::vector<double> falloff = {})
        {
            auto result = bayer_image_s<T>(image_size._y, image_size._x, image_data.data(), pattern);
            for (auto i : result)
            {
                auto channel = result[i];
                channel.fill(i < sensitivities.size() ? sensitivities[i] : 0);
                if (i >= sensitivities.size() || i >= falloff.size() || falloff[i] < 0.0 || falloff[i] >= 1.0)
                    continue;

                // reduce luminance when the parameters make sense
                double w2 = channel._width * 0.5;
                double h2 = channel._height * 0.5;
                double level = (1.0 - falloff[i] * falloff[i]) / (w2*w2 + h2*h2);
                for (int y = 0; y < channel._height; y++)
                {
                    for (int x = 0; x < channel._width; x++)
                    {
                        double d2 = level * ((y - h2)*(y - h2) + (x - w2)*(x - w2));
                        double correction = sqrt(1.0 - d2);
                        channel(y, x) = reduce_to<T>(channel(y, x) * correction);
                    }
                }
            }
            return result;
        }
    };

    SCENARIO("Lens Shading Grids can be created by dimensions")
    {
        GIVEN("Width, Height and Channel count")
        {
            const int width = 2;
            const int height = 5;
            const int channels = 1;
            WHEN("The lens shading grid with the given dimensions is created")
            {
                auto lsg = lensshading_grid<double>(width, height, channels);
                THEN("The object contains a grid data of 10 ( = 2 * 5 * 1) zeros")
                {
                    auto ref_vec = std::vector<double>(10);
                    CHECK(ref_vec == lsg.grid);
                }
            }
        }
    }

    SCENARIO("Lens Shading grids can be converted between floating point and fixed point formats")
    {
        GIVEN("A 1x3 matrix of lens shading with data")
        {
            auto lsg = lensshading_grid<double>(3, 1, 1);
            lsg.grid = { 1.0, 1.5, 2.0 };
            WHEN("The grid is converted to fixed point with 11 bits of fractional values")
            {
                auto fixed = lensshading_grid<uint16_t>(lsg, 11);
                THEN("The fixed point vector contains values of 2048, 3072 and 4096")
                {
                    auto ref = std::vector<uint16_t>{ 2048, 3072, 4096 };
                    CHECK(ref == fixed.grid);
                }
                AND_WHEN("The fixed point vector is converted to a (single precision) floating point vector")
                {
                    auto float_vec = lensshading_grid<float>(fixed);
                    THEN("The floating point vector matches the original values")
                    {
                        auto ref3 = std::vector<float>{ 1.0f, 1.5f, 2.0f};
                        CHECK(ref3 == float_vec.grid);
                    }
                }
            }
            WHEN("The grid is converted to fixed point with 12 bits of fractional values")
            {
                auto fixed = lensshading_grid<uint16_t>(lsg, 12);
                THEN("The fixed point vector contains values of 4096, 6144 and 8192")
                {
                    auto ref = std::vector<uint16_t>{ 4096, 6144, 8192 };
                    CHECK(ref == fixed.grid);
                }
                AND_WHEN("The fixed point vector is converted to a (double precision) floating point vector")
                {
                    auto float_vec = lensshading_grid<double>(fixed);
                    THEN("The floating point vector matches the original values")
                    {
                        CHECK(float_vec.grid == lsg.grid);
                    }
                }
            }
        }
    }

    SCENARIO("Lens shading model can characterize an image")
    {
        auto imagemaker = test_image_maker<uint16_t>(40, 60);

        GIVEN("A lens shading model to produce grids of size 5x3")
        {
            auto const grid_width = 5;
            auto const grid_height = 3;
            auto model = lensshading_calculator{ grid_width, grid_height };
            WHEN("An (constant) image of non-zero dimensions is characterized")
            {
                auto image = imagemaker.make_flat_image(bayer_pattern_e::bggr_4x2);
                auto grid = model.calculate_grid(image);
                THEN("The resulting grid is of given size 5x3 and has 8 channels of all ones")
                {
                    auto ref = std::vector<double>(5 * 3 * 8, 1.0);
                    CHECK(grid.width == grid_width);
                    CHECK(grid.height == grid_height);
                    CHECK(grid.grid == ref);
                }
                AND_THEN("An image of other pattern fails to characterize")
                {
                    auto image_with_different_pattern = imagemaker.make_flat_image(bayer_pattern_e::rggb_4x2);
                    CHECK_THROWS(model.calculate_grid(image_with_different_pattern));
                }
                AND_THEN("An image with different size fails to characterize")
                {
                    imagemaker.resize({ 50, 30 });
                    auto image_with_different_size = imagemaker.make_flat_image(bayer_pattern_e::bggr_4x2);
                    CHECK_THROWS(model.calculate_grid(image_with_different_size));
                }
            }
        }
    }

    SCENARIO("Setting up Lens shading model with grid index vector")
    {
        auto imagemaker = test_image_maker<uint16_t>(40, 60);

        GIVEN("A lens shading model to produce grids of size 5x4 and an 16-entry grid index vector")
        {
            auto const grid_width = 5;
            auto const grid_height = 3;
            auto const grid_channels = 4;       // number of unique channels in the vec2x2
            auto vec_2x2 = std::vector<char>({
                +0, +1, -1, -1,
                +2, +3, -1, -1,
                -1, -1, -1, -1,
                -1, -1, -1, -1 });
            auto model = lensshading_calculator{ grid_width, grid_height, vec_2x2 };
            WHEN("An GBRG image with constant channel intensities of 1,2,4 and 3 is characterized")
            {
                const uint16_t gb = 1;
                const uint16_t b = 2;
                const uint16_t r = 4;
                const uint16_t gr = 3;

                auto image = imagemaker.make_flat_image(bayer_pattern_e::gbrg, { gb, b, r, gr });
                double avg_green = (gr + gb) * 0.5;

                auto grid = model.calculate_grid(image);
                THEN("The result grid has 3 channels worth of ones and chromaticity calculated from the channel averages")
                {
                    auto ref = std::vector<double>(grid_width * grid_height * grid_channels, 1.0);
                    auto ref_chroma = chromaticity(r / avg_green, b / avg_green, 0.0);
                    CHECK(grid.width == grid_width);
                    CHECK(grid.height == grid_height);
                    CHECK(grid.channels == grid_channels);
                    CHECK(ref == grid.grid);
                    CHECK(ref_chroma == grid.chroma);
                    AND_WHEN("The grid is quantized to uint16_t 5.11 format")
                    {
                        auto grid_as_fixed_point = lensshading_grid<uint16_t>(grid, 11);
                        THEN("The grid dimensions match and the grid data is all 2048")
                        {
                            auto ref_2 = std::vector<uint16_t>(grid_width * grid_height * grid_channels, 2048);
                            CHECK(grid_as_fixed_point.width == grid.width);
                            CHECK(grid_as_fixed_point.height == grid.height);
                            CHECK(grid_as_fixed_point.channels == grid.channels);
                            CHECK(grid_as_fixed_point.grid == ref_2);
                            CHECK(grid_as_fixed_point.chroma == grid.chroma);
                        }
                    }
                }
            }
        }
    }

    SCENARIO("A synthetic image characterization matches the IQS/matlab generated reference values")
    {
        GIVEN("A large synthetic non-flat image characterized in 63x47 grid")
        {
            const int grid_width = 63;
            const int grid_height = 47;

            const int image_width = 4208;
            const int image_height = 3120;

            auto imagemaker = test_image_maker<uint16_t>(image_width, image_height);
            auto sensitivities = std::vector<uint16_t>({ 200, 300, 440, 550 });
            auto shading_factors = std::vector<double>({ 0.2, 0.3, 0.4, 0.5 });
            auto image = imagemaker.make_curved_image(bayer_pattern_e::rggb, sensitivities, shading_factors);

            auto output = lensshading_calculator{ grid_width, grid_height }.calculate_grid(image);
            THEN("The center pixel values are all ones and the corner values match recorded gains within 1/1000")
            {
                const double e = 1e-3;
                double recorded_gains[4] = { 7580.0 / 2048.0, 5873.0 / 2048.0, 4707.0 / 2048.0, 3896.0 / 2048.0 };
                for (int i = 0; i < 4; i++)
                {
                    auto channel = output[i]; // I'th channel as image
                    CHECK(channel(0, 0) == Approx(recorded_gains[i]).epsilon(e));
                    CHECK(channel(roi_point(grid_width / 2, grid_height / 2)) == 1.0);
                }
                CHECK(output.chroma._r_per_g == Approx(0.5404965).epsilon(e));
                CHECK(output.chroma._b_per_g == Approx(1.486724).epsilon(e));
                CHECK(output.chroma._i_per_g == 0.0);
            }
        }
    }

    SCENARIO("An image multiplied by its characterization produces a flat image")
    {
        auto apply_gain_table = [](lensshading_grid<double> &grid, bayer_image_s<uint16_t> &img)
        {
            for (int i = 0; i < 4; i++)
            {
                auto gain = grid[i];
                img[i].foreach([](uint16_t &dst, double &gain)
                {
                    dst = static_cast<uint16_t>(dst * gain + 0.5);
                }, gain);
            }
        };

        GIVEN("A non flat image with four channels with known maximum sensitivities")
        {
            const int grid_width = 63;
            const int grid_height = 49;

            auto imagemaker = test_image_maker<uint16_t>(grid_width * 2, grid_height * 2);
            auto sensitivities = std::vector<uint16_t>({ 200, 300, 440, 550 });
            auto shading_factors = std::vector<double>({ 0.2, 0.3, 0.4, 0.5 });
            auto image = imagemaker.make_curved_image(bayer_pattern_e::rggb, sensitivities, shading_factors);

            WHEN("The image is characterized to the same dimensions as image channel")
            {
                auto model = lensshading_calculator{ grid_width, grid_height };
                auto output = model.calculate_grid(image);

                THEN("The maximum gain values of each characterization channels are approximate reciprocals of the shading factors")
                {
                    double e = 1e-2;
                    REQUIRE(output.channels == 4);
                    for (int i = 0; i < 4; i++)
                    {
                        auto max_pos = output[i].max_element();
                        auto max_gain = output[i](max_pos);
                        CHECK(max_gain == Approx(1.0 / shading_factors[i]).epsilon(e));
                    }

                    AND_WHEN("The image is multiplied by the characterization")
                    {
                        apply_gain_table(output, image);
                        THEN("The corrected image channels are close to the original maximum sensitivities values")
                        {
                            for (int i = 0; i < 4; i++)
                            {
                                uint16_t val = sensitivities[i];
                                uint16_t max_diff = 0;
                                image[i].foreach([val, &max_diff](uint16_t &src)
                                {
                                    uint16_t diff = src > val ? src - val : val - src;
                                    max_diff = std::max(max_diff, diff);
                                });
                                CHECK(max_diff < 5);
                            }
                        }
                    }
                }
            }
        }
    }

    SCENARIO("Lens shading model can integrate grids over repetitions and msids")
    {
        auto grid_with_gain = [](int gain)
        {
            auto grid = lensshading_grid<double>(1, 1, 1);
            grid.grid[0] = gain;
            return grid;
        };

        GIVEN("A lens shading model configured to produce 1x1 grids")
        {
            lsc_model model(1, 1, {});
            WHEN("The model is provided with grids of gains 2,3 and 7 each with different light source")
            {
                int light_source = 10;
                for (int gain : { 2, 3, 7 })
                {
                    auto grid = grid_with_gain(gain);
                    model.add_item(grid, 0, light_source++);
                }

                THEN("The model provides a map of three grids with light sources 10,11 and 12")
                {
                    auto output = model.get_static_tables();
                    CHECK(output.size() == 3);
                    CHECK(output.count(10) == 1);
                    CHECK(output.count(11) == 1);
                    CHECK(output.count(12) == 1);
                    // And of course the count of some other light source, such as 13 is 0
                    CHECK(output.count(13) == 0);
                }
            }
        }

        GIVEN("A lens shading model with two grids associated to msid 1 and one grid associated with another msid")
        {
            int msid_1 = 1;
            int msid_2 = 2;
            int the_only_lightsource = 77;

            lsc_model model(1, 1, {});

            model.add_item(grid_with_gain(2), msid_1, the_only_lightsource);
            model.add_item(grid_with_gain(4), msid_1, the_only_lightsource);
            model.add_item(grid_with_gain(5), msid_2, the_only_lightsource);

            WHEN("The Static Tables are retrieved")
            {
                auto output = model.get_static_tables();
                THEN("The model provides a single grid with the combined gain the average of 5 and 3 (being the average of 2 and 4)")
                {
                    REQUIRE(output.count(the_only_lightsource) == 1);
                    auto grid = output[the_only_lightsource];
                    uint16_t reference = 4 * 2048;
                    CHECK(grid.grid[0] == reference);
                }
            }
        }
    }

    SCENARIO("Lens shading model produces ratio tables when NVM grids are provided")
    {
        // Creates test NVM datas with N light sources
        // - each grid is of size W*H * 1 channel and contains the 'gain' replicated W*H times
        // e.g. (2,2, { 1.0, 2.0 }) produces
        //   grid for pair_index = 0, with data = { 1.0, 1.0, 1.0, 1.0 }
        //   grid for pair_index = 1, with data = { 2.0, 2.0, 2.0, 2.0 }
        // both grids wrapped as a vector
        auto make_nvm_test_data = [](int width, int height, std::vector<double> gains)
        {
            auto result = lsc_nvm(1);
            int pair = 0;
            for (auto g : gains)
            {
                // make a new grid
                auto lightsource = lensshading_grid<double>(width, height, 1);
                lightsource.grid = std::vector<double>(width * height, g);
                lightsource.pair_index = pair++;
                result[0].push_back(lightsource);
            }
            return result;
        };

        GIVEN("A lens shading model to produce 2x2 grids and NVM grids with three light sources")
        {
            int msid = 0;
            int light_source = 99;
            int grid_width = 2;
            int grid_height = 2;
            lsc_model model(grid_width, grid_height, {});
            model.nvm_grids = make_nvm_test_data(grid_width, grid_height, { 1.0, 2.0, 4.0 });
            WHEN("The lens shading model is added with a 2x2(x1) static grid with content { 1.0, 2.0, 3.0, 3.5 }")
            {
                auto lsc_grid_for_msid_0 = lensshading_grid<double>(grid_width, grid_height, 1);
                lsc_grid_for_msid_0.grid = { 1.0, 2.0, 3.0, 3.5 };

                model.add_item(lsc_grid_for_msid_0, msid, light_source);

                AND_WHEN("The model is queried for static and ratio tables")
                {
                    auto statics = model.get_static_tables();
                    auto ratios = model.get_ratio_tables();
                    THEN("Both tables contain one light source 99 and the model has produced the ratio table for gain 4.0")
                    {
                        CHECK(statics.size() == 1);
                        CHECK(ratios.size() == 1);
                        CHECK(statics.count(light_source) == 1);
                        CHECK(ratios.count(light_source) == 1);

                        CHECK(ratios[light_source].pair_index == 2);

                        auto ref = std::vector<uint16_t>(
                        {
                            static_cast<uint16_t>(2048 * 1.0 / 4.0 + 0.5),
                            static_cast<uint16_t>(2048 * 2.0 / 4.0 + 0.5),
                            static_cast<uint16_t>(2048 * 3.0 / 4.0 + 0.5),
                            static_cast<uint16_t>(2048 * 3.5 / 4.0 + 0.5)
                        });
                        CHECK(ref == ratios[light_source].grid);
                    }
                }
            }
        }
    }
}