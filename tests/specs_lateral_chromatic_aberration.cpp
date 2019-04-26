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


#include "Teisko/LateralChromaticAberration.hpp"
#include "Teisko/Preprocessing.hpp"

#include "catch.hpp"
#include <memory>
#include <random>
#include <algorithm>
#include <vector>
using namespace Teisko;

namespace teisko_liblateralchromaticaberration_tests
{
    // Auxiliary class to assist in creating flat images with spherical or square holes
    class lca_imagemaker
    {
        // The features are drawn to full resolution color planes
        image<uint16_t> red_plane;
        image<uint16_t> green_plane;
        image<uint16_t> blue_plane;
        image<uint16_t> ir_plane;
        std::vector<uint16_t> intensities;  // Used later for anti-aliasing a circle
    public:
        // Creates initial flat planar image with four color channels with constant intensities
        lca_imagemaker(int width, int height, std::vector<uint16_t> color_intensities)
            : red_plane(image<uint16_t>(height, width).fill(color_intensities.size() > 0 ? color_intensities[0] : 0))
            , green_plane(image<uint16_t>(height, width).fill(color_intensities.size() > 1 ? color_intensities[1] : 0))
            , blue_plane(image<uint16_t>(height, width).fill(color_intensities.size() > 2 ? color_intensities[2] : 0))
            , ir_plane(image<uint16_t>(height, width).fill(color_intensities.size() > 3? color_intensities[3] : 0))
            , intensities(color_intensities)
        {}

        enum feature_type { circle, square };

        lca_imagemaker& draw_feature(color_info_e plane, std::vector<point_xy> centers, double radius, feature_type type = circle)
        {
            double ramp = 2.0;  // size of the anti aliased edge of circle
            double r2 = radius * radius;                          // outer radius of circle
            double r1 = (radius - ramp) * (radius - ramp);        // inner radius (inside this the circle is all zero)
            auto img = get_plane(plane);
            for (auto &center : centers)
            {
                auto top_left = clip(roi_point(center.x - radius, center.y -radius));
                auto bot_right = clip(roi_point(center.x + radius, center.y + radius));
                for (int row = top_left._y; row < bot_right._y; row++)
                {
                    for (int col = top_left._x; col < bot_right._x; col++)
                    {
                        auto point = point_xy(col, row);
                        if (type == circle)
                        {
                            double d_squared = norm2(point - center);
                            if (d_squared > r2)
                                continue;
                            if (d_squared > r1)
                            {
                                // 0 at inner radius, 1 at outer radius
                                auto dist = (std::sqrt(d_squared) - (radius - ramp)) / ramp;
                                auto pix = static_cast<uint16_t>(std::lround(intensities[static_cast<int>(plane)] * dist));
                                img(row, col) = pix;
                                continue;
                            }
                        }
                       img(row, col) = 0;
                    }
                }
            }
            return *this;
        }

        lca_imagemaker& draw_feature(std::vector<point_xy> greens, std::vector<point_xy> reds, std::vector<point_xy> blues)
        {
            REQUIRE(greens.size() == reds.size());
            REQUIRE(greens.size() == blues.size());
            std::vector<point_xy> item(1);      // just one point at a time
            for (size_t i = 0; i < greens.size(); i++)
            {
                item.back() = greens[i];
                draw_feature(color_info_e::green, item, 20);
                item.back() = greens[i] + reds[i];
                draw_feature(color_info_e::red, item, 20);
                item.back() = greens[i] + blues[i];
                draw_feature(color_info_e::blue, item, 20);
            }
            return *this;
        }

        // Constructs a bayer image by sampling from the 3 or 4 planes according to given sensor type
        bayer_image_s<uint16_t> mosaic(bayer_pattern_e pattern)
        {
            auto is_dp_4x2 = bayer_info_s(pattern).is_dp_4x2_sensor();
            auto img = bayer_image_s<uint16_t>(
                red_plane._height,
                red_plane._width * (is_dp_4x2 ? 2 : 1), pattern);

            // size of underlying sensor: 2x2 or 4x4
            // For 4x2 sensor the data is in 2x2 format, but is just duplicated horizontally
            auto dims = is_dp_4x2 ? roi_point(2,2) : get_dim(img._layout);
            if (dims._x == 0)
                throw std::runtime_error("Illegal sensor width");

            for (auto &&index : img)
            {
                auto src_index = is_dp_4x2 ? index >> 1 : index;

                auto channel = img[index];      // Nth sub channel in bayer image
                channel.init_from(              // Initialize each sub channel
                    get_plane(img._layout[index])   // By R,G,B or Ir plane
                    .subview(dims, roi_point(src_index % dims._x, src_index / dims._x)));
            }
            return img;
        }

    private:
        image<uint16_t>& get_plane(color_info_e plane)
        {
            switch (plane)
            {
            case color_info_e::red: return red_plane;
            case color_info_e::green: return green_plane;
            case color_info_e::blue: return blue_plane;
            case color_info_e::ir: return ir_plane;
            default:
                throw std::runtime_error("Illegal plane index");
            }
        }

        // Ensure a (corner of bounding box) is fully within the image
        roi_point clip(roi_point p)
        {
            return{
                std::max(0, std::min(static_cast<int>(red_plane._width), p._x)),
                std::max(0, std::min(static_cast<int>(red_plane._height), p._y))
            };
        }
    };

    enum lca_grid_indices : size_t { blue_x = 0, blue_y, red_x, red_y };

    SCENARIO("Smoke test of lca characterization class")
    {
        GIVEN("A LCA characterization class")
        {
            lca_characterization lca;
            THEN("The class without given any images returns zero cell size, zero optical center and zero grid")
            {
                CHECK(lca.center() == roi_point(0, 0));
                CHECK(lca.cell_size() == roi_point(0, 0));

                auto grid_empty = lca.collect();
                CHECK(grid_empty.channels == 0);
                CHECK(grid_empty.width == 0);
                CHECK(grid_empty.height == 0);
            }
        }
    }

    SCENARIO("Basic flow of lca characterization class")
    {
        GIVEN("A LCA characterization class and value for shift to be injected")
        {
            lca_characterization lca;
            double r = 20.0;        // radius of feature
            point_xy red_shift{ -1.1, 1.33 };
            point_xy blue_shift{ 4.2, -1.86 };

            AND_WHEN("We add an empty Bayer image of size 1000x600 to the model with one feature in the middle with given shift")
            {
                auto green = point_xy(300, 200);    // location of green feature
                auto img = lca_imagemaker(1000, 600, { 100, 200, 300, 0 })
                    .draw_feature(color_info_e::red, { green + red_shift }, r, lca_imagemaker::circle)
                    .draw_feature(color_info_e::green, { green }, r, lca_imagemaker::circle)
                    .draw_feature(color_info_e::blue, { green + blue_shift }, r, lca_imagemaker::circle)
                    .mosaic(bayer_pattern_e::bggr);

                auto features = lca += img;
                THEN("The returned vector contains one feature")
                {
                    CHECK(features.size() == 1);

                    AND_THEN("The model returns the optical center, the cell size and the grid no larger than 64x64")
                    {
                        auto center = lca.center();
                        auto cell_size = lca.cell_size();
                        auto grid = lca.collect();
                        CHECK(center == roi_point(1000 / 2, 600 / 2));
                        CHECK(cell_size._x > 0);
                        CHECK(cell_size._y > 0);
                        CHECK(grid.channels == 4);
                        CHECK(grid.width > 0);
                        CHECK(grid.height > 0);
                        AND_THEN("The returned grid is close to the given distortion")
                        {
                            double max_diff = 0.1;      // one tenth of a pixel
                            auto blue_diff_x = grid[blue_x].to_vector();
                            auto blue_diff_y = grid[blue_y].to_vector();
                            auto red_diff_x = grid[red_x].to_vector();
                            auto red_diff_y = grid[red_y].to_vector();
                            for (size_t i = 0; i < 1; i++)// blue_diff_x.size(); i++)
                            {
                                CHECK(std::fabs(blue_diff_x[i] - blue_shift.x) <= max_diff);
                                CHECK(std::fabs(blue_diff_y[i] - blue_shift.y) <= max_diff);
                                CHECK(std::fabs(red_diff_x[i] - red_shift.x) <= max_diff);
                                CHECK(std::fabs(red_diff_y[i] - red_shift.y) <= max_diff);
                            }
                        }
                    }
                }
            }
        }
    }


    SCENARIO("Basic flow of lca characterization class with four features")
    {
        GIVEN("A LCA characterization class and values for shift to be injected")
        {
            lca_characterization lca;
            std::vector<point_xy> greens({ { 200, 100 }, { 800, 100 }, { 500, 300 }, { 200, 500 }, { 800, 500 } });
            std::vector<point_xy> reds({ { -1.1, 0.9 }, { 2.0, 0.4 }, { 0, 0 }, { -2.5, -4.1 }, { 0.5, 3.0 } });
            std::vector<point_xy> blues({ { -2.1, -0.5 }, { 1.1, 1.5 }, { 0, 0 }, { 1.4, 0.5 }, { 1.1, 1.5 } });
            std::vector<uint16_t> max_pixel_intensities({ 100, 200, 300, 0 });
            const int width = 1000;
            const int height = 600;

            AND_WHEN("We add an empty Bayer image of size 1000x600 to the model with five features")
            {
                auto img = lca_imagemaker(width, height, max_pixel_intensities)
                    .draw_feature(greens, reds, blues)
                    .mosaic(bayer_pattern_e::rgib);

                lca += img;

                THEN("The model returns the optical center, the cell size and the grid no larger than 64x64")
                {
                    auto center = lca.center();
                    auto cell_size = lca.cell_size();
                    auto grid = lca.collect();
                    CHECK(center == roi_point(1000 / 2, 600 / 2));
                    CHECK(cell_size._x > 0);
                    CHECK(cell_size._y > 0);
                    CHECK(grid.channels == 4);
                    CHECK(grid.width > 0);
                    CHECK(grid.height > 0);
                    AND_THEN("The few values in the grid close to injected features match approximately the injected chromatic aberration")
                    {
                        for (size_t n = 0; n < greens.size(); n++)
                        {
                            auto idx = roi_point(greens[n].x * (1.0 / cell_size._x), greens[n].y * (1.0 / cell_size._y));
                            REQUIRE(idx._x > 1);
                            REQUIRE(idx._x + 2 < grid.width);
                            REQUIRE(idx._y > 1);
                            REQUIRE(idx._y + 2 < grid.height);
                            auto const max_difference = 0.2;
                            CHECK(std::fabs(grid[0](idx) - blues[n].x) < max_difference);
                            CHECK(std::fabs(grid[1](idx) - blues[n].y) < max_difference);
                            CHECK(std::fabs(grid[2](idx) - reds[n].x) < max_difference);
                            CHECK(std::fabs(grid[3](idx) - reds[n].y) < max_difference);
                        }
                    }
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

    SCENARIO("LCA characterization class gives (almost) identical results with all bayer patterns")
    {
        GIVEN("An image with bunch of features and a reference characterization in RGGB pattern")
        {
            const int width = 1000;
            const int height = 600;
            std::vector<point_xy> greens({ { 200, 100 }, { 800, 100 }, { 500, 300 }, { 200, 500 }, { 800, 500 } });
            std::vector<point_xy> reds({ { -1.1, 0.9 }, { 2.0, 0.4 }, { 0, 0 }, { -2.5, -4.1 }, { 0.5, 3.0 } });
            std::vector<point_xy> blues({ { -2.1, -0.5 }, { 1.1, 1.5 }, { 0, 0 }, { 1.4, 0.5 }, { 1.1, 1.5 } });
            std::vector<uint16_t> max_pixel_intensities({ 100, 200, 300, 0 });

            auto maker = lca_imagemaker(width, height, max_pixel_intensities)
                .draw_feature(greens, reds, blues);

            auto bayer_image = maker.mosaic(bayer_pattern_e::rggb);
            lca_characterization lca;
            lca += bayer_image;

            auto grid = lca.collect();
            CHECK(grid.channels == 4);
            CHECK(grid.width > 0);
            CHECK(grid.height > 0);

            WHEN("The same features are mosaiced in different pattern, the characterizations are almost equal")
            {
                for (auto &p : get_all_patterns())
                {
                    auto bayer_other = maker.mosaic(p);
                    lca_characterization lca_other;
                    lca_other += bayer_other;
                    auto grid_other = lca_other.collect();
                    const auto max_rms = 0.1;
                    REQUIRE(grid_other.channels == 4);
                    REQUIRE(grid_other.width == grid.width);
                    REQUIRE(grid_other.height == grid.height);
                    double rms = 0;
                    for (size_t i = 0; i < grid_other.grid.size(); i++)
                    {
                        double diff = grid_other.grid[i] - grid.grid[i];
                        rms += diff * diff;
                    }

                    //grid_other[0].write("C:\\temp\\bx.raw");
                    //grid_other[1].write("C:\\temp\\by.raw");
                    //grid_other[2].write("C:\\temp\\rx.raw");
                    //grid_other[3].write("C:\\temp\\ry.raw");
                    CHECK((rms / grid_other.grid.size()) < max_rms);
                }
            }
        }
    }

// #define DEBUG_FILE_OPERATIONS_SVEMD
#ifdef DEBUG_FILE_OPERATIONS_SVEMD
    SCENARIO("LCA Can characterize sve md images")
    {
        auto dims = roi_point(8064, 3024);
        lca_characterization lca;

        GIVEN("A file")
        {
            auto img = image<uint16_t>(dims).read("C:\\work\\lca_svemd\\name#IMX362~type#dot~lux#700~cct#6500~ag#256~exp#60000~cmnt#50cm~rep#1.raw");
            auto bayer_img = bayer_image_s<uint16_t>(img, bayer_pattern_e::rggb_4x2);

            auto prepro = preprocessor_s{};
            prepro.sve_matrix = 11915;
            prepro.saturation = 65535.0f;
            prepro.sve_md_demosaic_bilinear(bayer_img);

            lca += bayer_img;
            auto result = lca.collect();
            // lca.dump("C:\\temp\\nonwfov.raw");
            result[0].write("C:\\temp\\nonwfovbx.raw");
            result[1].write("C:\\temp\\nonwfovby.raw");
            result[2].write("C:\\temp\\nonwfovrx.raw");
            result[3].write("C:\\temp\\nonwfovry.raw");
        }
    }

#endif

//#define DEBUG_FILE_OPERATIONS_SINGLE
#ifdef DEBUG_FILE_OPERATIONS_SINGLE
    SCENARIO("LCA Can locate Regular features")
    {
        auto dims = roi_point(4208, 3120);
        lca_characterization lca;

        GIVEN("A file")
        {
            auto img = image<uint16_t>(dims).read("C:\\work\\dot1.raw");
            auto bayer_img = bayer_image_s<uint16_t>(img, bayer_pattern_e::rggb);
            lca += bayer_img;
            auto result = lca.collect();
            // lca.dump("C:\\temp\\nonwfov.raw");
            result[0].write("C:\\temp\\nonwfovbx.raw");
            result[1].write("C:\\temp\\nonwfovby.raw");
            result[2].write("C:\\temp\\nonwfovrx.raw");
            result[3].write("C:\\temp\\nonwfovry.raw");
        }
    }

#endif

//#define DEBUG_FILE_OPERATIONS_WFOV
#ifdef DEBUG_FILE_OPERATIONS_WFOV
    SCENARIO("LCA Can locate WFoV features")
    {
        auto dims = roi_point(4608, 3456);
        lca_characterization lca;

        GIVEN("A file")
        {
            auto img = image<uint16_t>(dims).read("C:\\work\\WFOV\\DB_WFOV\\lens_0\\name#OV16860~type#studio~ag#256~exp#80000~lsrc#19~cmnt#dot~rep#01.raw");
            auto bayer_img = bayer_image_s<uint16_t>(img, bayer_pattern_e::bggr);
            lca += bayer_img;
            img.read("C:\\work\\WFOV\\DB_WFOV\\lens_0\\name#OV16860~type#studio~ag#256~exp#80000~lsrc#19~cmnt#dot~rep#02.raw"); lca += bayer_img;
            img.read("C:\\work\\WFOV\\DB_WFOV\\lens_0\\name#OV16860~type#studio~ag#256~exp#80000~lsrc#19~cmnt#dot~rep#03.raw"); lca += bayer_img;
            img.read("C:\\work\\WFOV\\DB_WFOV\\lens_0\\name#OV16860~type#studio~ag#256~exp#80000~lsrc#19~cmnt#dot~rep#04.raw"); lca += bayer_img;
            img.read("C:\\work\\WFOV\\DB_WFOV\\lens_0\\name#OV16860~type#studio~ag#256~exp#80000~lsrc#19~cmnt#dot~rep#05.raw"); lca += bayer_img;
            img.read("C:\\work\\WFOV\\DB_WFOV\\lens_0\\name#OV16860~type#studio~ag#256~exp#80000~lsrc#19~cmnt#dot~rep#06.raw"); lca += bayer_img;
            img.read("C:\\work\\WFOV\\DB_WFOV\\lens_0\\name#OV16860~type#studio~ag#256~exp#80000~lsrc#19~cmnt#dot~rep#07.raw"); lca += bayer_img;
            img.read("C:\\work\\WFOV\\DB_WFOV\\lens_0\\name#OV16860~type#studio~ag#256~exp#80000~lsrc#19~cmnt#dot~rep#08.raw"); lca += bayer_img;
            img.read("C:\\work\\WFOV\\DB_WFOV\\lens_0\\name#OV16860~type#studio~ag#256~exp#80000~lsrc#19~cmnt#dot~rep#09.raw"); lca += bayer_img;
            img.read("C:\\work\\WFOV\\DB_WFOV\\lens_0\\name#OV16860~type#studio~ag#256~exp#80000~lsrc#19~cmnt#dot~rep#10.raw"); lca += bayer_img;
            img.read("C:\\work\\WFOV\\DB_WFOV\\lens_0\\name#OV16860~type#studio~ag#256~exp#80000~lsrc#19~cmnt#dot~rep#11.raw"); lca += bayer_img;
            img.read("C:\\work\\WFOV\\DB_WFOV\\lens_0\\name#OV16860~type#studio~ag#256~exp#80000~lsrc#19~cmnt#dot~rep#12.raw"); lca += bayer_img;
            img.read("C:\\work\\WFOV\\DB_WFOV\\lens_0\\name#OV16860~type#studio~ag#256~exp#80000~lsrc#19~cmnt#dot~rep#13.raw"); lca += bayer_img;
            img.read("C:\\work\\WFOV\\DB_WFOV\\lens_0\\name#OV16860~type#studio~ag#256~exp#80000~lsrc#19~cmnt#dot~rep#14.raw"); lca += bayer_img;
            img.read("C:\\work\\WFOV\\DB_WFOV\\lens_0\\name#OV16860~type#studio~ag#256~exp#80000~lsrc#19~cmnt#dot~rep#15.raw"); lca += bayer_img;
            img.read("C:\\work\\WFOV\\DB_WFOV\\lens_0\\name#OV16860~type#studio~ag#256~exp#80000~lsrc#19~cmnt#dot~rep#16.raw"); lca += bayer_img;
            auto result = lca.collect();
            result[0].write("C:\\temp\\result_bx.raw");
            result[1].write("C:\\temp\\result_by.raw");
            result[2].write("C:\\temp\\result_rx.raw");
            result[3].write("C:\\temp\\result_ry.raw");
        }
    }
#endif
}