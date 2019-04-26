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


#include "Teisko/MacbethDetector.hpp"
#include "Teisko/Image/RGB.hpp"
#include "catch.hpp"
#include <memory>
#include <random>
#include <algorithm>
#include <vector>
#include <string>
#include <fstream>

//#define INTERACTIVE_VISUALIZATION
//#define USE_VISUALIZATION_AND_JPEG
#ifdef USE_VISUALIZATION_AND_JPEG
#include "c:\work\tinyjpeg\CImg.h"              // Visualization of images + loading of .jpg
#include "c:\work\tinyjpeg\tinyjpeg.hpp"        // ImageDB contains jpeg files not readable by CImg.h
#endif

using namespace Teisko;

struct macbeth_generator
{
    rgb_image_s<uint16_t> canvas;                           // collection of three rgb color planes
    std::vector<uint16_t> get_patch_values(int idx = -1)    // gray scale, vs red, green, blue
    {
        const uint8_t srgb[24][3] = {
            /*  sRGB                 R    G    B */
            /*  dark skin     */ { 116,  79,  65 },
            /*  light skin    */ { 197, 144, 127 },
            /*  blue sky      */ {  91, 120, 155 },
            /*  foliage       */ {  91, 108,  64 },
            /*  blue flower   */ { 131, 127, 175 },
            /*  bluish green  */ {  95, 189, 172 },
            /*  orange        */ { 224, 124,  48 },
            /*  purplish blue */ {  69,  90, 167 },
            /*  moderate red  */ { 197,  80,  95 },
            /*  purple        */ {  93,  58, 104 },
            /*  yellow green  */ { 156, 187,  58 },
            /*  orange yellow */ { 227, 161,  39 },
            /*  blue          */ {  40,  62, 145 },
            /*  green         */ {  61, 147,  70 },
            /*  red           */ { 178,  54,  57 },
            /*  yellow        */ { 236, 199,  15 },
            /*  magenta       */ { 191,  79, 146 },
            /*  cyan          */ {   0, 133, 165 },
            /*  white         */ { 241, 242, 235 },
            /*  neutral 8     */ { 201, 202, 201 },
            /*  neutral 6.5   */ { 161, 163, 163 },
            /*  neutral 5     */ { 121, 121, 121 },
            /*  neutral 3.5   */ {  83,  84,  85 },
            /*  black         */ {  50,  50,  50 },
        };
        auto patch_values = std::vector<uint16_t>(24);
        for (int i = 0; i < 24; i++)
        {
            if (idx < 0 || idx > 2)
                patch_values[i] = srgb[i][0] + srgb[i][1] + srgb[i][2];
            else
                patch_values[i] = srgb[i][idx];
        }
        return patch_values;
    }

    // Generates srgb macbeth chart on RGB image of given height
    // - the macbeth patch size is currently constant -- 50 pixels
    macbeth_generator(int canvas_width, int canvas_height) : canvas(roi_point(canvas_width, canvas_height))
    {
        // Draw chart
        auto patch_size = roi_point(50, 50);
        auto patch_distance = roi_point(58, 58);
        auto patch_offset = roi_point(300, 200);

        for (int ch : {0, 1, 2})
        {
            auto patch_values = get_patch_values(ch);
            for (int row = 0; row < 4; row++)
            {
                for (int col = 0; col < 6; col++)
                {
                    auto p = roi_point(col, row);
                    canvas[ch].region(patch_size, patch_offset + patch_distance * p).fill(patch_values[row * 6 + col]);
                }
            }
        }
    }

    // Constructs a bayer image by sampling from the 3 or 4 planes according to given sensor type
    bayer_image_s<uint16_t> mosaic(bayer_pattern_e pattern)
    {
        auto is_dp_4x2 = bayer_info_s(pattern).is_dp_4x2_sensor();
        auto size = canvas[0].size();
        auto img = bayer_image_s<uint16_t>(size._y, size._x * (is_dp_4x2 ? 2 : 1), pattern);

        // size of underlying sensor: 2x2 or 4x4
        // For 4x2 sensor the data is in 2x2 format, but is just duplicated horizontally
        auto dims = is_dp_4x2 ? roi_point(2, 2) : get_dim(img._layout);
        if (dims._x == 0)
            throw std::runtime_error("Illegal sensor width");

        for (auto &&index : img)
        {
            auto src_index = is_dp_4x2 ? index >> 1 : index;
            auto color_code = static_cast<int>(img._layout[index]); // color as 0,1 or 2
            auto channel = img[index];      // Nth sub channel in bayer image
            if (color_code > 2)
                channel.fill(0);            // no IR
            else
                channel.init_from(canvas[color_code].subview(dims, roi_point(src_index % dims._x, src_index / dims._x)));
        }
        return img;
    }

    // Constructs a bayer image by sampling from the 3 or 4 planes according to given sensor type
    image<uint16_t> to_gray_scale()
    {
        auto size = canvas[0].size();
        auto img = image<uint16_t>(size._y, size._x).init_from(canvas[0]);
        for (auto i : { 1, 2 })
            img.foreach([](uint16_t &sum, uint16_t &src) {sum += src; }, canvas[i]);

        return img;
    }
};


namespace teisko_libmacbeth_tests
{
    SCENARIO("Macbeth detector founds a chart from artificial RGB image")
    {
        auto model = macbeth_generator(640, 480);
        GIVEN("An artificial RGB macbeth chart of size 640x480")
        {
            rgb_image_s<uint16_t> rgb_image = model.canvas;     // type explicitly visible

            WHEN("The chart is detected using the macbeth chart finder")
            {
                auto detection = macbeth_chart().find(rgb_image);
                THEN("The detection is valid and it returns 3*24 patch values matching exactly to the model")
                {
                    REQUIRE(detection.is_valid() == true);

                    auto patches_red = get_patch_trimmed_mean<uint16_t>(rgb_image[rgb_color_e::red], detection.polygons);
                    auto patches_green = get_patch_trimmed_mean<uint16_t>(rgb_image[rgb_color_e::green], detection.polygons);
                    auto patches_blue = get_patch_trimmed_mean<uint16_t>(rgb_image[rgb_color_e::blue], detection.polygons);

                    CHECK(patches_red == model.get_patch_values(0));
                    CHECK(patches_green == model.get_patch_values(1));
                    CHECK(patches_blue == model.get_patch_values(2));
                }
            }
        }
    }

    SCENARIO("Macbeth detector founds a chart from artificial bayer image")
    {
        auto model = macbeth_generator(640, 480);
        GIVEN("An artificial macbeth chart of size 640x480 seen through bayer RGGB filter")
        {
            bayer_image_s<uint16_t> bayer_rggb = model.mosaic(bayer_pattern_e::rggb);

            WHEN("The chart is detected using the macbeth chart finder")
            {
                auto detection = macbeth_chart().find(bayer_rggb);
                THEN("The detection is valid and it returns 24 gray scale patch values matching exactly to the model")
                {
                    auto scale = point_xy(0.5, 0.5);        // each color channel is half the size of bayer image
                    REQUIRE(detection.is_valid() == true);
                    auto patches_red = get_patch_trimmed_mean<uint16_t>(bayer_rggb[0], detection.polygons, scale);
                    auto patches_g_red = get_patch_trimmed_mean<uint16_t>(bayer_rggb[1], detection.polygons, scale);
                    auto patches_g_blue = get_patch_trimmed_mean<uint16_t>(bayer_rggb[2], detection.polygons, scale);
                    auto patches_blue = get_patch_trimmed_mean<uint16_t>(bayer_rggb[3], detection.polygons, scale);

                    CHECK(patches_red == model.get_patch_values(0));
                    CHECK(patches_g_red == model.get_patch_values(1));
                    CHECK(patches_g_blue == model.get_patch_values(1));
                    CHECK(patches_blue == model.get_patch_values(2));
                }
            }
        }
    }

    SCENARIO("Macbeth detector founds a chart from artificial gray scale image. "
        "In this scenario we also show that the chart locator fails to find the chart (by design) "
        "from those otherwise perfectly discovered charts of 6x4 patches, where there is no "
        "descending progression of intensities in the (logical) bottom row.")
    {
        auto model = macbeth_generator(640, 480);
        GIVEN("An artificial grayscale macbeth chart of size 640x480, which we rotate any arbitrary angle")
        {
            image<uint16_t> gray_scale_image = model.to_gray_scale();
            gray_scale_image = rotate_image(gray_scale_image, 212.0);

            WHEN("The chart is detected using the macbeth chart finder")
            {
                auto detection = macbeth_chart().find(gray_scale_image);
                THEN("The detection is valid and it returns 24 gray scale patch values matching exactly to the model")
                {
                    REQUIRE(detection.is_valid() == true);

                    auto patches_gray = get_patch_trimmed_mean<uint16_t>(std::move(gray_scale_image), detection.polygons);

                    CHECK(patches_gray == model.get_patch_values(-1));

                    AND_WHEN("The image is corrupted so that the bottom 6 patches are not in descending order")
                    {
                        uint16_t look_for = patches_gray[23];
                        uint16_t change_to = patches_gray[20];
                        gray_scale_image.foreach([look_for, change_to](uint16_t &a) { if (a == look_for) a = change_to; });
                        THEN("The image is no longer detected with macbeth chart")
                        {
                            auto detect2 = macbeth_chart().find(gray_scale_image);
                            CHECK(detect2.is_valid() == false);
                        }
                    }
                    AND_WHEN("We mirror [or flip] the image")
                    {
                        auto mirrored = gray_scale_image.mirror();
                        THEN("The image is no longer detected with macbeth chart")
                        {
                            auto detect3 = macbeth_chart().find(mirrored);
                            CHECK(detect3.is_valid() == false);
                        }
                    }
                }
            }
        }
    }
}

#ifdef USE_VISUALIZATION_AND_JPEG
namespace teisko_libmacbeth_tests_with_files
{
    image<uint8_t> combine_rgb(image<uint8_t> &r, image<uint8_t> &g, image<uint8_t> &b)
    {
        r.foreach([](uint8_t &red, uint8_t &green, uint8_t &blue)
        {
            red = std::max({ red, green, blue });
        }, g, b);
        return r.convert_to();
    }

    std::vector<uint8_t> read_file(std::string filename)
    {
        std::ifstream input;
        input.open(filename, std::ios::binary);
        input.seekg(0, std::ios::end);
        std::streampos file_len = input.tellg();
        input.seekg(0, std::ios::beg);
        std::vector<uint8_t> data(static_cast<size_t>(file_len));
        input.read(reinterpret_cast<char*>(data.data()), file_len);
        return data;
    }

#ifdef USE_CIMG_FOR_JPG
    image<uint16_t> load_jpeg(std::string fname)
    {
        cimg_library::CImg<unsigned char> jpeg(fname.c_str());
        auto *ptr = jpeg.data();

        auto width = jpeg.width();
        auto height = jpeg.height();
        auto ptr = jpeg.data();

        auto red = image<uint8_t>(height, width, reinterpret_cast<uint8_t*>(ptr));
        auto green = image<uint8_t>(height, width, reinterpret_cast<uint8_t*>(ptr) + width * height);
        auto blue = image<uint8_t>(height, width, reinterpret_cast<uint8_t*>(ptr)+ 2 * width * height);
        return combine_rgb(red, green, blue);
    }
#else
    image<uint8_t> load_jpeg(std::string fname)
    {
        auto vec = read_file(fname);
        Jpeg::Decoder decoder(vec.data(), vec.size());
        if (decoder.GetResult() != Jpeg::Decoder::OK)
            throw std::runtime_error("Can't load");

        auto width = decoder.GetWidth();
        auto height = decoder.GetHeight();
        auto ptr = decoder.GetImage();
        auto img = image<uint8_t>(height, width * 3, ptr);
        auto red = img.subview(1, 3, 0, 0);
        auto green = img.subview(1, 3, 0, 1);
        auto blue = img.subview(1, 3, 0, 2);
        auto result = combine_rgb(red, green, blue);
        return result;
    }
#endif
    SCENARIO("Using macbeth detector to detect images from raw files")
    {
        GIVEN("A file")
        {
            struct rawfile {
                roi_point dims;
                std::string name;
                rawfile(int w, int h, std::string fname) : dims(w, h), name(fname) {  }
            };
            std::vector<rawfile> files;
            files.emplace_back(5634, 3752, "C:\\git_repo\\koe1.raw");
            std::ifstream file("H:\\raws.txt");
            if (file.is_open()) {
                std::string line;
                while (getline(file, line)) {
                    std::string prefix = "H:\\ccc_raw_";
                    if (line.compare(0, prefix.length(), prefix) == 0)
                    {
                        auto rest = line.substr(prefix.length());
                        auto w = std::stoi(rest);
                        auto xpos = rest.find("x");
                        auto h = std::stoi(rest.substr(xpos + 1));
                        files.emplace_back(w, h, line);
                    }
                }
                file.close();
            }

            WHEN("The file is decoded")
            {
                //for (size_t i : {2937, 2952, 931, 1143, 1228, 1424, 1509, 1730, 1774, 1866, 1935, 2145, 2631 })
                for (size_t i = 0; i < 1; i++)// files.size(); i++)
                {
                    auto &file = files[i];
                    auto raw16 = image<uint16_t>(file.dims).read(file.name);
                    auto maxval = raw16(raw16.max_element());
                    int shift = 0;
                    while (maxval > 255)
                    {
                        shift++;
                        maxval /= 2;
                    }
                    auto r_planar = raw16.transform<uint8_t>([shift](uint16_t &src){ return static_cast<uint8_t>(src >> shift); });

                    auto detection = macbeth_chart().find(raw16);
                    CHECK(detection.is_valid() == true);
                    if (!detection.is_valid())
                    {
                        std::string msg = "Chart not found from " + std::to_string(i) + " :" + file.name;
                        WARN(msg);
                    }

#ifdef INTERACTIVE_VISUALIZATION
                    // Configure this to show all or just the failed ones in bulk test
                    if (!detection.is_valid())
                    {
                        // Visual feedback and browse through images by
                        // toggling arrow keys (forward, backward, fast forward, fast backward)
                        cimg_library::CImg<uint8_t> loaded_image(r_planar._begin, r_planar._width, r_planar._height);

                        int poly_idx = 0;
                        for (auto &p : detection.polygons)
                        {
                            uint8_t color[] = { poly_idx++ == 20 ? 255 : 155 };
                            cimg_library::CImg<int> points(4, 2);
                            for (int i = 0; i < 4; i++)
                            {
                                points(i, 0) = static_cast<int>(p[i].x);
                                points(i, 1) = static_cast<int>(p[i].y);
                            }
                            loaded_image.draw_polygon(points, color);
                        }

                        auto title = std::to_string(i) + "/" + std::to_string(files.size()) + " " + file.name;

                        cimg_library::CImgDisplay main_disp(loaded_image, title.c_str());
                        while (!main_disp.is_closed())
                        {
                            if (i > 0 && main_disp.is_key(cimg_library::cimg::keyARROWLEFT))
                            {
                                i -= 2;
                                break;
                            }
                            if (main_disp.is_key(cimg_library::cimg::keyARROWRIGHT))
                            {
                                break;
                            }
                            if (main_disp.is_key(cimg_library::cimg::keyARROWDOWN))
                            {
                                i += 20;
                                break;
                            }
                            if (main_disp.is_key(cimg_library::cimg::keyARROWUP))
                            {
                                i -= 20;
                                if (i > files.size())
                                {
                                    i = 0; --i;
                                }
                                break;
                            }
                        }
#endif
                    }
                }
            }
        }
    }

    std::vector<point_xy> calculate_centroids(std::vector<std::vector<point_xy>> &polygons)
    {
        auto output = std::vector<point_xy>();
        for (size_t i = 0; i < polygons.size(); i++)
        {
            auto p = polygons[i];
            if (p.size() != 4)
                return{};
            double area = 0.0;
            point_xy centroid(0, 0);
            for (int j = 0; j < 4; j++)
            {
                int jplus = (j + 1) & 3;  // next polygon vertex, wrapping from 3 to 0
                double cross = p[j].x * p[jplus].y - p[jplus].x*p[j].y;
                area += cross;
                centroid = centroid + (p[j] + p[jplus]) * cross;
            }
            area = area == 0.0 ? 0.0 : 1.0 / (6.0 * 0.5 * area);
            output.push_back(centroid * area);
        }
        return output;
    }

    SCENARIO("Detecting WFOV test image produces almost identical results to opencv version")
    {
        GIVEN("Path to image of given size")
        {
            std::string name = "C:\\git_repo\\koe1.raw";
            roi_point dims(5634, 3752);
            auto raw_image = image<uint16_t>(dims).read(name);
            auto bayer_image = bayer_image_s<uint16_t>(raw_image, bayer_pattern_e::bggr);
            WHEN("A detection is created from the image")
            {
                auto detection = macbeth_chart().find(bayer_image);
                THEN("The detection coordinates matches the precalculated set")
                {
                    double max_diff = dims._y * 0.005;      // relatively small, absolutely quite a lot -- 18 pixels
                    auto reference = std::vector<point_xy>({
                        { 2167, 1470 }, { 2406, 1459 }, { 2681, 1458 }, { 2971, 1470 }, { 3243, 1494 }, { 3482, 1525 },
                        { 2151, 1731 }, { 2390, 1734 }, { 2667, 1741 }, { 2960, 1753 }, { 3237, 1768 }, { 3479, 1785 },
                        { 2151, 1993 }, { 2386, 2011 }, { 2657, 2027 }, { 2945, 2039 }, { 3219, 2045 }, { 3460, 2047 },
                        { 2166, 2234 }, { 2393, 2264 }, { 2654, 2288 }, { 2928, 2301 }, { 3191, 2303 }, { 3426, 2294 }
                    });

                    for (size_t i = 0; i < reference.size(); i++)
                    {
                        auto points = calculate_centroids(detection.polygons);
                        CHECK(std::fabs(points[i].x - reference[i].x) < max_diff);
                        CHECK(std::fabs(points[i].y - reference[i].y) < max_diff);
                    }

                    AND_THEN("The data retrieved matches the collected rgb values")
                    {
                        // These values are checked by matlab:
                        // fh = fopen('c:\git_repo\koe1.raw', 'r'); bw = reshape(fread(fh,'ushort'),[5634 3752])'; fclose(fh);
                        // c0=bw(1:2:end,1:2:end);
                        // med=[];for p=1:24;k=c0(tl(p,2):br(p,2),tl(p,1):br(p,1)); med=[med trimmean(k(:),50)]; end; m=[m med'];
                        // .. repeat for the other three channels c0=bw(1:2:end, 2:2:end);
                        //  - the c++ implementation uses polygon tests
                        auto reference_medians = std::vector<std::vector<double>>({
                            { 19.3603, 4.4294, 23.4496, 19.3985 },
                            { 72.5982, 20.2721, 83.1954, 72.5365 },
                            { 45.6382, 21.2892, 27.5053, 45.7661 },
                            { 34.0573, 7.0154, 23.4367, 34.0479 },
                            { 52.1888, 25.2718, 40.6174, 52.3506 },
                            { 97.2171, 32.5541, 45.8439, 97.3343 },
                            { 52.2490, 6.2349, 78.5581, 52.0471 },
                            { 27.0655, 19.9731, 15.5734, 27.2082 },
                            { 27.7309, 7.4742, 60.9050, 27.5516 },
                            { 12.8238, 6.8958, 14.8261, 12.8188 },
                            { 87.0847, 11.0000, 57.1042, 87.0435 },
                            { 75.3358, 7.8987, 89.6190, 75.1516 },
                            { 10.7821, 10.6009, 5.1197, 10.8681 },
                            { 43.9204, 8.8258, 20.4208, 43.9648 },
                            { 15.8563, 3.5062, 42.9638, 15.7443 },
                            { 97.5846, 10.2226, 95.0672, 97.3734 },
                            { 27.5504, 12.7779, 51.3337, 27.4753 },
                            { 35.0751, 19.5430, 10.1648, 35.1507 },
                            { 139.450, 44.3984, 106.8198, 139.2354 },
                            { 93.1904, 30.2685, 70.9819, 93.1731 },
                            { 58.4140, 19.0949, 44.4164, 58.4959 },
                            { 28.6083, 9.3325, 21.7963, 28.6627 },
                            { 13.5084, 4.4198, 10.1604, 13.5745 },
                            { 4.0915, 1.3914, 3.1781, 4.2099 }
                        });

                        for (int ch = 0; ch < 4; ch++)
                        {
                            double max_diff = 0.50;
                            auto scale = point_xy(0.5, 0.5);        // due to getting subchannels
                            auto vec = get_patch_trimmed_mean(bayer_image[ch], detection.polygons, scale);
                            REQUIRE(vec.size() == 24);
                            for (size_t j = 0; j < 24; j++)
                            {
                                CHECK(std::fabs(reference_medians[j][ch] - vec[j]) < max_diff);
                            }
                        }
                     }
                }
            }
        }
    }
}
#endif
