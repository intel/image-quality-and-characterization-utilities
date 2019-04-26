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


#pragma once
#include "Teisko/Image/API.hpp"
#include "Teisko/Image/RGB.hpp"
#include "Teisko/Image/Algorithms.hpp"
#include "Teisko/Preprocessing.hpp"
#include "Teisko/Algorithm/VectorMedian.hpp"
#include "Teisko/Algorithm/TrimmedMean.hpp"
#include "Teisko/Algorithm/NelderMead.hpp"

#include <cstdint>
#include <vector>
#include <map>
#include <array>
#include <list>
#include <string>
#include <algorithm>
#include <numeric>

namespace Teisko
{
    // Utility class to accumulate variance
    //  - this is stable enough when mean is close to zero
    //  - improvement is to accumulate (item - K), where K is close to mean
    //    or within the range of the samples (i.e. 0)
    struct variance_s
    {
        double sum = 0;
        double sum2 = 0;
        size_t n = 0;
        variance_s() : sum(0), sum2(0), n(0) { }
        variance_s& operator+=(const double sample)
        {
            sum += sample;
            sum2 += sample * sample;
            n++;
            return *this;
        }
        // returns true std (or unscaled variance if normalization or sqrt is not needed)
        double operator() (bool scaled = true)
        {
            if (scaled)
                return n <= 1 ? 0 : std::sqrt((sum2 - (sum*sum) / (double)n) / (double)(n - 1));
            return sum2*n - sum*sum;
        }
    };

    /// Fills the inside of triangle given in points a,b and c with function `void func(T &pixel);`
    /// To the rectangular canvas located between `offset` ... `offset + canvas.size()`
    template <typename T, typename fill_function>
    inline void triangle_fill(image<T> &canvas, roi_point offset, point_xy a, point_xy b, point_xy c, fill_function func)
    {
        // the canvas is logically contained between `offset` to `offset + canvas.size()`
        for (decltype(canvas._height) y = 0; y < canvas._height; y++)
        {
            auto p = point_xy(offset._x, y + offset._y);
            // delta between each x-increment
            auto d_a = b.y - c.y;
            auto d_b = c.y - a.y;
            auto d_c = a.y - b.y;
            // barycentric weight at each line
            auto w_a = d_a * (p.x - b.x) + (c.x - b.x) * (p.y - b.y);
            auto w_b = d_b * (p.x - c.x) + (a.x - c.x) * (p.y - c.y);
            auto w_c = d_c * (p.x - a.x) + (b.x - a.x) * (p.y - a.y);
            for (decltype(canvas._width) x = 0; x < canvas._width; x++, w_a += d_a, w_b += d_b, w_c += d_c)
            {
                if (w_a >= 0 && w_b >= 0 && w_c >= 0)
                    func(canvas(y, x));
            }
        }
    }

    // Special case of polygon fill -- handles quads (by splitting them to two triangles)
    // - calculates the signed distance of all points in ROI to all polygon edges
    // - for convex quads/polygons, a point is inside if the signed distance >= 0 for all edges
    // - for non-convex quad, the point is inside, if the point is inside either of the triangles
    //   - the triangle must be split along the vertex forming non-right turning angle
    template <typename T>
    inline void quad_fill(image<T> &canvas, std::vector<point_xy> &quad, roi_point offset = { 0, 0 })
    {
        // test for an angle between poins a-b-p being a monotonically right turning
        // - a quad (or any other polygon) can be fan-triangulated from the unique concave vertex
        auto turns_right = [](point_xy a, point_xy b, point_xy p)
        {
            auto d = b - a;
            p = p - a;
            return d.x * p.y > d.y * p.x;
        };

        auto sz = quad.size();
        if (sz != 4)
            throw std::runtime_error("This function only handles quads");

        // Two possible triangulations -- split along the concave (or 180 degree) angle
        // if that doesn't exists, either triangulation is OK
        // The other triangulation happens between points 0,1,2/0,2,3
        // The other with points 0,1,3/1,2,3
        auto i = (!turns_right(quad[3], quad[0], quad[1]) || !turns_right(quad[1], quad[2], quad[3]))
            ? 0 : 1;
        triangle_fill(canvas, offset, quad[0], quad[1], quad[2 + i], [](T &x){ x = 1; });
        triangle_fill(canvas, offset, quad[0 + i], quad[2], quad[3], [](T &x){ x = 1; });
    };

    // Quantizes the point x,y to an integral multiple of multiplier
    inline roi_point lower_multiple(roi_point x, roi_point multiplier)
    {
        return (multiplier._x < 2 || multiplier._y < 2)
            ? x : (x / multiplier) * multiplier;
    }

    // Quantizes the point x,y to next integral multiple of multiplier
    inline roi_point upper_multiple(roi_point x, roi_point multiplier)
    {
        return (multiplier._x < 2 || multiplier._y < 2)
            ? x : lower_multiple(x + multiplier - 1, multiplier);
    }

    // Clamps a point p between [0,0], [size._x, size._y] inclusive
    inline roi_point clamp_to(roi_point p, roi_point size)
    {
        return roi_point(
            std::max(std::min(p._x, size._x), 0),
            std::max(std::min(p._y, size._y), 0)
            );
    }

    // Returns a bw image of W*H being a multiple of sensor_dimensions (if given)
    // - where each pixel contained within the quad is 1, otherwise 0
    template <typename T = uint8_t>
    image<T> get_patch_stencil(
        roi_point size,
        std::vector<point_xy> &quad,
        roi_point &offset,
        roi_point sensor_dimensions = roi_point(0,0))
    {
        if (quad.size() != 4)
            throw std::runtime_error("Expecting a quad as polygon");

        auto min_x = std::floor(std::min({ quad[0].x, quad[1].x, quad[2].x, quad[3].x }));
        auto max_x = std::ceil(std::max({ quad[0].x, quad[1].x, quad[2].x, quad[3].x }));
        auto min_y = std::floor(std::min({ quad[0].y, quad[1].y, quad[2].y, quad[3].y }));
        auto max_y = std::ceil(std::max({ quad[0].y, quad[1].y, quad[2].y, quad[3].y }));

        auto top_left = clamp_to(lower_multiple(roi_point(min_x, min_y), sensor_dimensions), size);
        auto bot_right = clamp_to(upper_multiple(roi_point(max_x, max_y), sensor_dimensions), size);

        auto output = image<T>(bot_right - top_left).fill(0);
        quad_fill(output, quad, top_left);

        offset = top_left;
        return output;
    }

    /// Get patch pixel data under the stencil
    template <typename T>
    void get_pixel_data(image<T> &pixels, image<uint8_t> &stencil, std::vector<T> &result)
    {
        auto max_dims = stencil.size();
        result.clear();
        result.reserve(max_dims._x * max_dims._y);
        stencil.foreach([&result](uint8_t &mask, T &src)
        {
            if (mask) result.emplace_back(src);
        }, pixels);
    }

    /// Get patch data from an image based on a stencil (mask)
    /// The image data top left corner must be pre-aligned to stencil top left corner
    ///  - the image data must be at least as large as the stencil
    ///  - returns the standard deviation of the pixels included by the stencil
    template <typename T>
    double get_patch_deviation(image<T> &&channel, image<uint8_t> &stencil)
    {
        variance_s my_functor;
        stencil.foreach([&my_functor](uint8_t &mask, T &src)
        {
            if (mask) my_functor += static_cast<double>(src);
        }, channel);
        return my_functor();
    }

    /// Get patch data from an image based on a stencil (mask)
    /// The image data top left corner must be pre-aligned to stencil top left corner
    ///  - the image data must be at least as large as the stencil
    ///  - returns the trimmed mean of those pixels with the corresponding non-zero pixel in stencil
    template <typename T>
    double get_patch_trimmed_mean(image<T> &&channel, image<uint8_t> &stencil, double save_proportion = 0.5)
    {
        std::vector<T> data;
        get_pixel_data(channel, stencil, data);
        return data.size() == 0 ? 0.0 : trimmed_mean(data, save_proportion);
    }

    /// Get data from a channel based on a vector of polygons
    ///  - the polygons can have an optional scaling factor
    template <typename U = double, typename image_type>
    std::vector<U> get_patch_trimmed_mean(
        image_type &&channel,
        std::vector<std::vector<point_xy>> &polygons,
        point_xy scale = point_xy(1.0, 1.0))
    {
        auto result = std::vector<U>();
        std::transform(polygons.begin(), polygons.end(), std::back_inserter(result),
            [&scale, &channel](std::vector<point_xy> quad)
        {
            for (auto &point : quad)
                point = point * scale;
            auto offset = roi_point(0, 0);
            auto stencil = get_patch_stencil(channel.size(), quad, offset);
            return reduce_to<U>(get_patch_trimmed_mean(channel.region(stencil.size(), offset), stencil));
        });
        return result;
    }


    /// macbeth_chart is a class holding and finding polygonal patches
    /// from images and extracting patch values from images based on the polygons
    class macbeth_chart
    {
    public:
        std::vector<std::vector<point_xy>> polygons;    /// Result of macbeth detection

        /// Default ctor -- set patch fill ratio controlling the size of generated polygons
        ///  - 0.0 = point, 1.0 = patches share corners, 0.65 = default
        macbeth_chart(double patch_fill_ratio = 0.65) : fill_ratio(patch_fill_ratio)
        {
            if (patch_fill_ratio < 0.0 || patch_fill_ratio > 1.0)
                throw std::runtime_error("Macbeth chart fill ratio should be between 0 and 1");
       }

        /// Returns true, if the chart contains 24 polygons
        bool is_valid() {
            const int chart_size = 4 * 6;
            return polygons.size() == chart_size;
        }

        /// Locates a Macbeth chart from a bayer_image
        macbeth_chart& find(bayer_image_s<uint16_t> &img)
        {
            auto green = demosaic_bilinear(img, color_info_e::green, MIRROR_EVEN);
            auto red = demosaic_bilinear(img, color_info_e::red, MIRROR_EVEN);
            auto blue = demosaic_bilinear(img, color_info_e::blue, MIRROR_EVEN);
            return find(green, red, blue);
        }

        /// Locates a Macbeth chart from an RGB image
        template <typename T>
        macbeth_chart& find(rgb_image_s<T> &img)
        {
            auto green = img[rgb_color_e::green].template convert_to<uint16_t>();
            auto red = img[rgb_color_e::red].template convert_to<uint16_t>();
            auto blue = img[rgb_color_e::blue].template convert_to<uint16_t>();
            return find(green, red, blue);
        }

        /// Locates a Macbeth chart from a single gray scale (or R,G,B) channel
        template <typename T>
        macbeth_chart& find(image<T> &img)
        {
            detect_macbeth_chart(img.template convert_to<uint16_t>());
            return *this;
        }

    private:
        double fill_ratio;
        /// Sums up R,G,B channels to locate the macbeth chart
        /// With this prototype being private we save the burden of checking that channel dimensions match
        macbeth_chart& find(image<uint16_t> &green, image<uint16_t> &red, image<uint16_t> &blue)
        {
            green.foreach([](uint16_t &g, uint16_t &r, uint16_t &b)
            {
                g = static_cast<uint16_t>((g + r + b) / 3);
            }, red, blue);
            return find(green);
        }

        // Properties of located features
        // - detecting squares from other candidates is done primarily by the feature
        //   properties domain
        struct feature_s
        {
            // Primary statistics
            std::vector<point_xy> path;     // List of pixels at the outer boundary of object
            std::vector<point_xy> corners;  // Four outermost pixels from the centroid
            point_xy centroid;              // Centroid of a feature

            // Additional statistics mainly used for screening
            double area;                    // Area of feature in pixels
            double corrected_path_length;   // Length of path after barrel correction
            double minimum_distance;        // Distance to closest feature
            double fft_rms;                 // Similarity of path to square
            int group;                      // Assigns a feature to 'group' by proximity

            feature_s(region_props &props, std::vector<roi_point> &perimeter)
                : centroid(props.centroid_x, props.centroid_y)
                , area(props.area)
                , corrected_path_length(0.0)
                , minimum_distance(0.0)
                , fft_rms(0.0)
                , group(0)
            {
                path.reserve(perimeter.size());
                for (auto &p : perimeter)
                    path.emplace_back(p._x, p._y);
                extract_corners();
            }

            void extract_corners()
            {
                corners.clear();
                corners.reserve(4);
                // locates the furthest point of the feature from the centroid
                auto furthest_1 = furthest_point(path, centroid);
                // and then the opposite point of the previous
                auto furthest_2 = furthest_point(path, furthest_1);
                // and locate two extreme points to the line between the two calculated extrema
                auto parallels = max_distance_from_line(furthest_1, furthest_2, path);

                corners.emplace_back(furthest_1);
                corners.push_back(parallels[0]);
                corners.emplace_back(furthest_2);
                corners.push_back(parallels[1]);
            }

            double calculate_squareness()
            {
                // Matches the centroid distance curve of the feature perimeter to precalculated
                // model of a square in FFT domain, returning the root mean square error
                fft_rms = -1.0;
                auto sz = path.size();
                if (sz < 32)
                    return fft_rms;

                // make a lookup table for the sines and cosines
                double angle = 6.283185307179586 / sz;
                std::vector<point_xy> trigs(sz);
                for (size_t i = 0; i < sz; i++)
                    trigs[i] = { std::cos(angle * i), std::sin(angle * i) };

                const int fft_points = 16;
                const double sqFFT[fft_points - 1] = {  // first term is 1.0 after normalization
                    0.0010871038, 0.0012745257, 0.0019206153, 0.0760958988, 0.0007661686,
                    0.0000199350, 0.0005197369, 0.0158663285, 0.0003992480, 0.0000085283,
                    0.0003177932, 0.0070672313, 0.0002918863, 0.0000040173, 0.0002226015
                };

                std::vector<double> fft_magnitude;  // evaluate FFT at 1 + 15 points
                for (size_t i = 0; i < fft_points; i++)
                {
                    size_t idx = 0;
                    point_xy fft_item;
                    for (auto &p : path)
                    {
                        auto term = norm(p - centroid) * trigs[idx];
                        fft_item += term;
                        idx += i;
                        if (idx >= sz)
                            idx -= sz;
                    }
                    fft_magnitude.emplace_back(norm(fft_item));
                }

                double rms = 0.0;
                double scale = 1.0 / fft_magnitude[0];
                for (int i = 0; i < fft_points - 1; i++)
                {
                    auto diff = sqFFT[i] - fft_magnitude[i + 1] * scale;
                    rms += diff * diff;
                }
                fft_rms = std::sqrt(rms * (2.0 / (2 * fft_points - 1)));  // 31 items contribute to the "mean"

                return fft_rms;
            }

            // Matches a feature to all other features in a container
            // (dismissing those that have identical coordinates)
            double find_minimum_distance(std::list<feature_s> &contours)
            {
                minimum_distance = -1.0;
                for (auto &feature : contours)
                {
                    auto dist2 = norm2(centroid - feature.centroid);
                    if (dist2 == 0.0)
                        continue;       // comparing to self results in zero distance
                    if (minimum_distance < 0 || dist2 < minimum_distance)
                        minimum_distance = dist2;
                }
                if (minimum_distance >= 0)
                    minimum_distance = std::sqrt(minimum_distance);

                return minimum_distance;
            }

            // Return true if two features are close enough and have about same size (+50%)
            bool is_connected(const feature_s &other) const
            {
                double min_area = area;
                double max_area = other.area;
                if (min_area > max_area)
                   std::swap(min_area, max_area);
                if (min_area * 1.5 < max_area)
                    return false;
                auto max_dist = std::min(norm2(corners[0] - corners[2]), norm2(corners[1] - corners[3]));
                if (norm2(centroid - other.centroid) > 3 * max_dist)
                    return false;
                for (auto c0: corners)
                    for (auto c1 : other.corners)
                    {
                        if (norm2(c0 - c1) * 2 < max_dist)
                            return true;
                    }

                return false;

                // return min_area * 1.5 > max_area && norm2(centroid - other.centroid) < min_area * 3.0;
            }

            // Finds the two most distant points from a line defined by two end points
            //  - one has positive signed distance and the other has negative distance
            std::array<point_xy, 2> max_distance_from_line(point_xy &left, point_xy &right,
                std::vector<point_xy> &vec)
            {
                std::array<point_xy, 2> result = { left, right };
                point_xy &min = result[0];
                point_xy &max = result[1];
                point_xy delta = left - right;
                double dxy = left.y * delta.x - left.x * delta.y;
                double max_dist = 0;
                double min_dist = 0;
                for (auto &p : vec)
                {
                    double signed_distance = dxy + p.x * delta.y - p.y * delta.x;
                    if (signed_distance > max_dist)
                    {
                        max_dist = signed_distance;
                        max = p;
                    }
                    else if (signed_distance < min_dist)
                    {
                        min_dist = signed_distance;
                        min = p;
                    }
                }
                return result;
            }

            // Given a set of points, returns the furthest point (or self)
            point_xy furthest_point(std::vector<point_xy> &vec, point_xy &self)
            {
                point_xy best_so_far = self;
                double best_distance = 0.0;
                for (auto &p : vec)
                {
                    double d = norm2(p - self);
                    if (d > best_distance)
                    {
                        best_distance = d;
                        best_so_far = p;
                    }
                }
                return best_so_far;
            }
        };

        // Quantifies image by threshold after edge detection
        // - horizontal + vertical filtering with 9 tap kernel
        image<uint16_t> gradient_threshold(image<uint16_t> &img, double threshold = 0.3)
        {
            // from matlab where W = normalized magnitude between 0..1 from hypot(img conv h, img conv h')
            // W = W.^(1./rolloffFactor);   %% rolloffFactor = 3
            // W = (1 - W). / (1 + W);
            // W(W < weightCutoff) = floorOfW;
            // result = W < threshold
            // we don't calculate sqrt of horizontal and vertical gradients, so we have to scale by
            // exponent of six...
            double y = std::pow(std::abs(threshold - 1) / (threshold + 1), 6);

            auto temp_dst = image<int>(img.size());
            auto temp_src = img.make_borders(4, 4, 4, 4, REPLICATE);
            auto center = temp_src.region(temp_dst.size(), roi_point(4));
            int stride = temp_src._skip_y;
            int max_level = 0;
            temp_dst.foreach([stride, &max_level](int &dst, uint16_t &src)
            {
                uint16_t *s = &src;
                int v =
                    (s[stride * 4] - s[-stride * 4]) * 1228 +
                    (s[stride * 3] - s[-stride * 3]) * 2210 +
                    (s[stride * 2] - s[-stride * 2]) * 2752 +
                    (s[stride * 1] - s[-stride * 1]) * 2002;        // sum = 8192
                int h =
                    (s[4] - s[-4]) * 1228 +
                    (s[3] - s[-3]) * 2210 +
                    (s[2] - s[-2]) * 2752 +
                    (s[1] - s[-1]) * 2002;
                v /= 16384;   // max sum of 'v' is +-65535 * 8192
                h /= 16384;   // -- we have to divide by 2 in order of not to overflow 'dst'
                dst = v*v + h*h;
                if (dst > max_level)
                    max_level = dst;
            }, center);
            max_level = std::lround(y * max_level);
            return temp_dst.transform<uint16_t>([max_level](int &src) -> uint16_t { return src < max_level; });
        }

        image<uint16_t> image_open(image<uint16_t> &img)
        {
            return img.filter<3, 3>([](support<3, 3, uint16_t> &data) { return data.maximum(); });
        }

        image<uint16_t> image_close(image<uint16_t> &img)
        {
            return img.filter<3, 3>([](support<3, 3, uint16_t> &data) { return data.minimum(); });
        }

        bool is_full_row_missing(std::vector<int> &indices, size_t start_index)
        {
            for (size_t i = 0; i < 6; i++)
                if (start_index + i >= indices.size() || indices[start_index + i] >= 0)
                    return false;
            return true;
        }

        void shift_indices(std::vector<int> &indices, size_t start_index)
        {
            std::vector<int> tmp = indices;
            auto tmp_size = tmp.size();
            for (size_t i = 0; i < tmp_size; i++)
                indices[i] = tmp[(i + start_index) % tmp_size];
        }

        bool is_increasing_sequence(std::vector<double> &data)
        {
            if (data.size() != 24)
                return false;
            for (int i = 0; i < 5; i++)
                if (data[22-i] <= data[23-i])
                    return false;
            return true;
        }

        // Checks that the detection in "polygon" / "point"
        // has been able to locate an increasing sequence of intensities
        // from patches 23,22,21,20,19 and 18
        bool validate_grid(image<uint16_t> &item, double scale)
        {
            auto result = get_patch_trimmed_mean(item, polygons, point_xy(scale, scale));
            if (is_increasing_sequence(result))
                return true;
            // try the other way round
            std::reverse(polygons.begin(), polygons.end());
            std::reverse(result.begin(), result.end());
            return is_increasing_sequence(result);
        }

        // Blur image by block filter
        template <typename T>
        void image_blur(image<T> &img, int kernel_size)
        {
            if (kernel_size < 2)
                return;
            if (kernel_size == 2)
            {
                img = img.template filter<2, 2>([](support<2, 2, T> &support)
                {
                    return round_shift<T, 2>(support[0][0] + support[0][1] + support[1][0] + support[1][1]);
                }, REPLICATE);
            }
            else
            {
                filter_separable(img, std::vector<double>(kernel_size, 1.0 / kernel_size), REPLICATE);
            }
        }

        /// Median filter an image
        template <typename T>
        void median_filter_3x3(image<T> &img)
        {
            img = img.template filter<3, 3>([](support<3, 3, T> &support) { return support.median(); }, REPLICATE);
        }

        // actual implementation --
        // receives one channel (or a gray scale channel)
        // - prescaling_factor = 2 for bayer_image_s and 1 for RGB images
        void detect_macbeth_chart(image<uint16_t> gray_scale)
        {
            // The image is first rescaled to ~1000 x 750 range
            auto scale = resize_image_to_1MP_range(gray_scale);
            std::list<std::list<feature_s>> contours_candidates;  // make a few attempts at locating contours

            auto bw_img_without_filtering = adaptive_threshold<uint16_t>(gray_scale);
            contours_candidates.emplace_back(find_square_candidates(bw_img_without_filtering));

            // Image is filtered depending on the size with median and/or box filter (up to 2x)?
            // - the filtering should be more adaptive
            // - images with sharp edges / no gap between patches suffer from filtering
            // - images with high noise require more filtering

            image_blur(gray_scale, (gray_scale._height + 256) / 512);

            if (gray_scale._height > 300)
                median_filter_3x3(gray_scale);

            // This seems to work quite fine for dark images
            auto bw_img = adaptive_threshold<uint16_t>(gray_scale);
            contours_candidates.emplace_back(find_square_candidates(bw_img));

            // Occasionally the image should be morphologically opened/closed
            // Closing by 5x5 support, then opening by 5x5 support is found empirically
            // to work in few cases in the offline image database of 4000+ images
            auto bw_copy = image_close(bw_img);
            bw_copy = image_close(bw_copy);
            bw_copy = image_open(bw_copy);
            bw_copy = image_open(bw_copy);
            contours_candidates.emplace_back(find_square_candidates(bw_copy));

            // and this seems to work quite fine for noisy images
            // We call image_open, because gradient threshold will typically produce quite wide borders
            // - we widen the squares inside these borders
            bw_img = gradient_threshold(gray_scale);
            bw_img = image_open(bw_img);
            contours_candidates.emplace_back(find_square_candidates(bw_img));

            sort_list_of_list_by_size(contours_candidates);
            for (auto &contours : contours_candidates)
            {
                point_xy center(bw_img._width * 0.5, bw_img._height * 0.5);

                auto contours_all = contours;

                // Iterate over all groups of adjacent features (starting from the largest group)
                auto groups = find_adjacent_features(contours);

                // Estimate global barrel correction from one connected group (at a time)
                // When the group is of proper size, we can try to use that as a global solution
                for (auto &g : groups)
                {
                    if (g.size() <= 2 || g.size() > 24)
                        continue;

                    double alpha = g.size() < 5 ? 0.0 : estimate_barrel_distortion(g, center);
                    auto alpha_per_r2 = alpha / norm2(center);

                    // take a copy of all contours, or the current connected group if it is large enough
                    contours = g.size() >= 8 && contours_all.size() > 24 ? g : contours_all;

                    // Apply the correction...
                    for (auto &feature : contours)
                    {
                        barrel_correction(feature.centroid, center, alpha_per_r2);
                        barrel_correction(feature.path, center, alpha_per_r2);
                    }

                    if (!remove_non_squares(contours))
                        continue;

                    if (!remove_by_mismatched_pathlengths(contours))
                        continue;

                    auto median_centroid_distance = remove_by_centroid_distance(contours);

                    if (median_centroid_distance < 0)
                        continue;

                    auto angle = find_representative_angle(contours);

                    std::vector<point_xy> centroids;
                    for (auto &feature : contours)
                        centroids.emplace_back(feature.centroid * scale);

                    auto scaled_center = center * scale;
                    auto indices = fit_centroids_to_grid(centroids, angle, median_centroid_distance * scale, scaled_center);
                    if (indices.size() == 0)
                        continue;

                    // Try at most two places...
                    for (int i = 0; i < 2; i++)
                    {
                        auto points = fill_missing_items(centroids, indices);
                        if (points.size() == 24)
                        {
                            calculate_polygons(points, scaled_center, alpha / norm2(scaled_center));
                            if (validate_grid(gray_scale, 1.0 / scale))
                                return;
                        }
                        if (is_full_row_missing(indices, 0))
                            shift_indices(indices, 6);
                        else if (is_full_row_missing(indices, 18))
                            shift_indices(indices, 18);
                        else break;
                    }
                    polygons.clear();
                }
            }
            // return without finding a chart
        }

        // downscales image data to about 1000 x 750 size -- returns the integer scaling factor
        // target size is predetermined by author of the algorithm
        int resize_image_to_1MP_range(image<uint16_t> &img)
        {
            int scale = 1;
            auto size = img.size();
            if ((size._y > 750) || (size._x > 1000))
            {
                scale = std::lround((size._y / 750.0 + size._x / 1000.0) * 0.5);
                if (scale <= 1)
                    return 1;
                img = resize_image(img, point_xy(1.0 / scale, 1.0 / scale));
                return scale;
            }
            return scale;
        }

        // modifies image data from grayscale to bw with local thresholding:
        // paints flat areas with white, where flat is described as
        // maximum of 3x3 neighborhood is not larger than X * minimum of 3x3 neighborhood
        template <typename U, typename T>
        image<U> adaptive_threshold(image<T> &img)
        {
            // erode
            uint64_t sum_min = 0;
            auto min3x3 = img.template filter<3, 3>([&sum_min](support<3, 3, T> &pixels)
            {
                auto mn = pixels.minimum();
                sum_min += static_cast<uint64_t>(mn);
                return mn;
            });

            // dilate
            uint64_t sum_max = 0;
            auto max3x3 = img.template filter<3, 3>([&sum_max](support<3, 3, T> &pixels)
            {
                auto mx = pixels.maximum();
                sum_max += static_cast<uint64_t>(mx);
                return mx;
            });

            // we could also use ratio = 1.333;
            // initial guess of threshold is the ratio of averaged mx and mn (e.g. 1.02)
            auto ratio = static_cast<float>((sum_max * 1.0) / (sum_min * 1.0));

            // then we refine the ratio by averaging only those mn/mx that are classified as edges
            sum_min = 0;
            sum_max = 0;
            min3x3.foreach([&sum_min, &sum_max, &ratio](T &mn, T &mx)
            {
                if (mn * ratio < mx)
                {
                    sum_min += mn;
                    sum_max += mx;
                }
            }, max3x3);

            // and then we threshold with the refined ratio
            ratio = static_cast<float>((sum_max * 1.0) / (sum_min * 1.0));

            auto output = image<U>(img.size());
            output.foreach([ratio](U &dst, T &mn, T &mx)
            {
                dst = mn * ratio > mx;
            }, min3x3, max3x3);

            return output;
        }

        // Given a bw image, returns descriptor for all connected (white) areas
        // of correct size and reasonable squareness (area to perimeter ratio)
        std::list<feature_s> find_square_candidates(image<uint16_t> &bw_img)
        {
            // Prune out items based on area or squareness area/perimeter
            // Actual formula:  0.65 <= 4piA / p^2 <= 1.05
            //  - then we invert this to get limits for the min/max perimeter length
            double max_area = (bw_img._width * bw_img._height) / 24.0;
            double min_area = std::min(100.0, max_area / 100.0);
            double max_perimeter = std::sqrt(19.3329 * max_area) * 1.5;
            double min_perimeter = std::sqrt(11.9680 * min_area) * (1.0 / 1.5);

            bwlabel(bw_img);
            std::vector<roi_point> path;
            path.reserve(10000);
            std::list<feature_s> passed;
            uint16_t next_label = 1;
            for (auto j = 0; j < bw_img._height; j++)
            {
                for (auto i = 0; i < bw_img._width; i++)
                {
                    roi_point p(i, j);
                    if (bw_img(p) == next_label)
                    {
                        ++next_label;
                        path.clear();
                        auto props = get_regionprops<uint16_t, true>(bw_img, p, &path);
                        if (props.area >= min_area && props.area <= max_area &&
                            props.perimeter >= min_perimeter && props.perimeter <= max_perimeter)
                        {
                            const double magic = 4.0 * 3.141592653589793;
                            double squareness = magic * props.area / (props.perimeter * props.perimeter);
                            if (0.65 <= squareness && squareness <= 1.05)
                                passed.emplace_back(props, path);
                        }
                    }
                }
            }
            return passed;
        }

        // Given multiple lists, sort them the longest list first
        void sort_list_of_list_by_size(std::list<std::list<feature_s>> &list_of_list)
        {
            // Then sort by number of features in each group
            list_of_list.sort([](const std::list<feature_s> &a, const std::list<feature_s> &b)
            {
                return a.size() > b.size();
            });
        }

        // Moves the content in feature list 'feats' to lists of adjacent items
        // Complexity: O(N^2)
        std::list<std::list<feature_s>> find_adjacent_features(std::list<feature_s> &feats)
        {
            std::list<std::list<feature_s>> results;
            while (feats.size())
            {
                // move one (first) item from feature list to a new group
                results.push_back(std::list<feature_s>());
                auto &group = results.back();
                group.splice(group.begin(), feats, feats.begin());

                // Move new members from features to the group
                for (auto it = group.begin(); it != group.end(); ++it)
                {
                    for (auto it2 = feats.begin(); it2 != feats.end();)
                    {
                        auto next = std::next(it2);
                        if (it->is_connected(*it2))
                            group.splice(group.end(), feats, it2);
                        it2 = next;
                    }
                }
            }
            sort_list_of_list_by_size(results);
            return results;
        }

        // barrel_correct a point from screen space to cartesian
        static void barrel_correction(point_xy &point, point_xy center, double alpha_per_r2)
        {
            point -= center;
            point = center + point / (1 + alpha_per_r2 * norm2(point));
        }

        // barrel correct a vector of points from "raw" to Cartesian
        static void barrel_correction(std::vector<point_xy> &a, point_xy center, double alpha_per_r2)
        {
            if (alpha_per_r2 != 0.0)
            {
                for (auto &p : a)
                    barrel_correction(p, center, alpha_per_r2);
            }
        }

        // distorts a point (x,y) from undistorted domain back to input screen space
        void inverse_barrel(point_xy &p, point_xy center, double alpha_per_r2)
        {
            auto original_point = p - center;
            auto converging_point = original_point;
            const double accepted_difference = 0.05;        // pixels
            for (int iterations = 0; iterations < 15; iterations++)
            {
                auto old_point = converging_point;
                double inverse_distortion = 1 + alpha_per_r2 * norm2(old_point);
                converging_point = original_point * inverse_distortion;

                if (norm2(converging_point - old_point) < (accepted_difference * accepted_difference))
                    break;
            }
            p = converging_point + center;
        }

        // distort the vector to input screen space
        void inverse_barrel(std::vector<point_xy> &vec, point_xy center, double alpha_per_r2)
        {
            if (alpha_per_r2 != 0)
            {
                for (auto &p : vec)
                    inverse_barrel(p, center, alpha_per_r2);
            }
        }

        // calculates a squareness estimate of group of 4-sided polygons
        // the two diagonals in a perfect square have equal lengths
        // - minimizing variance of the diagonal difference provides a good estimate of overall squareness
        double squareness_estimate(std::vector<point_xy> &corners)
        {
            variance_s squareness;
            auto sz = corners.size();
            if (sz == 0 || sz % 4)
                return 0.0;

            for (decltype(sz) i = 0; i < sz; i += 4)
                squareness += norm(corners[i] - corners[i + 2]) - norm(corners[i + 1] - corners[i + 3]);

            return squareness();
        }


        std::vector<point_xy> extract_corners(std::list<feature_s> &contours)
        {
            std::vector<point_xy> output;
            output.reserve(contours.size() * 4);
            for (auto &c : contours)
            {
                for (auto &corner : c.corners)
                    output.emplace_back(corner);
            }
            return output;
        }

        // returns estimated alpha for R = r * ( 1+ alpha*r^2)
        // requires image center (given as image size)
        // The contours are scanned to detect four extreme points
        //  - the best alpha minimizes the STD (or variance) of the angles between the corners
        //  - ideally all corners are 90 degrees - thus the variance and mean are both zero
        double estimate_barrel_distortion(std::list<feature_s> &contours, point_xy center)
        {
            // Checks if all elements in a 2x2 neighborhood are missing
            auto is_empty_quad = [](std::vector<int> &indices, int offset)
            {
                return
                    indices[offset + 0] < 0 && indices[offset + 1] < 0 &&
                    indices[offset + 6] < 0 && indices[offset + 7] < 0;
            };

            auto angle = find_representative_angle(contours);

            std::vector<point_xy> centroids;
            for (auto &feature : contours)
                centroids.emplace_back(feature.centroid);

            auto median_distance = get_centroid_median_distance(contours);

            auto indices = fit_centroids_to_grid(centroids, angle, median_distance, center);
            if (indices.size() != 24)
                return 0.0;

            // check that there's representation in all 4 corners
            if (is_empty_quad(indices, 0) || is_empty_quad(indices, 4) ||
                is_empty_quad(indices, 12) || is_empty_quad(indices, 16))
                return 0.0;

            auto corners = extract_corners(contours);
            if (corners.empty())
                return 0.0;

            std::vector<double> measure;
            std::vector<double> alphas;

            // Loop for a range of doubles in [-3.2 : 0.1 : 3.2] inclusive
            // this seems to be a practical range of allowed distortion
            for (int k = -32; k <= 32; k++)
            {
                auto alpha = k * 0.1;
                auto distorted = corners;
                barrel_correction(distorted, center, alpha / norm2(center));
                measure.push_back(squareness_estimate(distorted));
                alphas.push_back(alpha);
            }
            auto it = std::min_element(measure.begin(), measure.end());
            auto idx = std::distance(measure.begin(), it);
            return alphas[idx];
        }

        // Using FFT method, removes all non squares from the set of features
        // Returns false for empty set
        bool remove_non_squares(std::list<feature_s> &contours)
        {
            contours.remove_if([](feature_s &f) {
                const double threshold = 0.02;      // empirical value
                auto rms = f.calculate_squareness();
                return rms < 0 || rms > threshold;
            });
            return !contours.empty();
        }

        // Calculates the path length of barrel corrected perimeter
        double calculate_undistorted_pathlength(feature_s &feature)
        {
            auto &path = feature.path;
            auto count = path.size();
            double sum = 0;

            for (decltype(count) i = 1; i < count; i++)
                sum += norm(path[i] - path[i - 1]);

            return sum;
        }

        // Finds the median length of the feature
        // - when there are more than 24 items, the item lenghts are binned to three bin histogram
        //   then we locate those features, that count less than 24
        // - this is mostly needed for images with both a 6x4 macbeth and some 10x10 chart
        double locate_typical_feature_length(std::vector<double> &sorted_lengths)
        {
            auto nonzero_paths = sorted_lengths.size();
            double typical_length = 0;

            if (nonzero_paths <= 24)
            {
                auto clip = (nonzero_paths + 2) >> 2;
                auto first = sorted_lengths.begin() + clip - 1;
                auto last = sorted_lengths.end() - clip;
                auto items = (double)std::distance(first, last);
                typical_length = std::accumulate(first, last, 0.0) / items;
            }
            else
            {
                // make a histogram of three bins, locating the maximum count that is less than 25
                // and the corresponding bin value
                // the original version in matlab locates the maximum size -- which is probably untested
                double diff = sorted_lengths.back() - sorted_lengths.front();
                double edge1 = sorted_lengths.front() + (1.0 / 3.0) * diff;
                double edge2 = sorted_lengths.front() + (2.0 / 3.0) * diff;
                int bins[3] = { 0, 0, 0 };
                double avg[3] = { 0, 0, 0 };
                for (auto &x : sorted_lengths)
                {
                    auto idx = x < edge1 ? 0 : ((x > edge2) ? 2 : 1);
                    bins[idx]++;
                    avg[idx] += x;
                }

                int max_count = 0;
                for (int i = 0; i < 3; i++)
                {
                    if (bins[i] <= 24 && bins[i] > max_count)
                    {
                        max_count = bins[i];
                        typical_length = avg[i] / bins[i];
                    }
                }
            }
            return typical_length;
        }

        // Locates the median of (minimum) distance across a group of contours
        double get_centroid_median_distance(std::list<feature_s> &contours)
        {
            std::vector<double> minimum_centroid_distance;

            for (auto &feature : contours)
                minimum_centroid_distance.push_back(feature.find_minimum_distance(contours));

            return vector_median(minimum_centroid_distance);
        }

        // removes features too far from other features
        //  - this will remove isolated squares, however that should
        //  not be a problem with the local thresholding algorithm
        double remove_by_centroid_distance(std::list<feature_s> &contours)
        {
            // multiplier = 1.1 * sqrt(2) == 1.55536
            auto median_centroid_distance = get_centroid_median_distance(contours);
            if (median_centroid_distance < 0)
                return median_centroid_distance;

            auto max_centroid_distance = 1.5556 * median_centroid_distance;
            contours.remove_if([=](const feature_s &feature)
            {
                return feature.minimum_distance > max_centroid_distance;
            });

            if (contours.empty())
                return -1.0;

            return median_centroid_distance;
        }

        // removes features with path length outside the expected mean
        bool remove_by_mismatched_pathlengths(std::list<feature_s> &contours)
        {
            std::vector<double> path_lengths;
            // Store the path length both in the feature and in a vector to be sorted
            for (auto &feature : contours)
            {
                double length = calculate_undistorted_pathlength(feature);
                feature.corrected_path_length = length;
                path_lengths.push_back(length);
            }

            auto count = path_lengths.size();
            if (count <= 1)
                return false;

            std::sort(path_lengths.begin(), path_lengths.end());

            auto typical_length = locate_typical_feature_length(path_lengths);
            auto max_difference = 0.25 * typical_length;

            contours.remove_if([=](const feature_s &feature)
            {
                return std::abs(feature.corrected_path_length - typical_length) > max_difference;
            });
            return !contours.empty();
        }

        // Finds rotation angle of the feature set by binning all edges of the feature set
        // to a histogram. Gives angles between +-90 degrees occasionally rotating almost perfectly
        // axis aligned chart by 90 degrees -- this is then compensated in grid fitting
        double find_representative_angle(std::list<feature_s> &contours)
        {
            // We expect the angles to deviate by multiples of 90 degrees.
            // By binning the angles by 60 degrees we ensure that there will be enough samples
            // in some of the slots
            std::vector<std::vector<double>> angle_histogram(3);        // 3 bins
            const auto bin1_angle = 1.047197551196598;  // pi/3
            const auto bin2_angle = 2.094395102393195;  // 2*pi/3
            auto corners = extract_corners(contours);
            auto sz = corners.size();
            std::vector<double> angles_all;
            for (size_t i = 0; i < sz; i++)
            {
                // Treat the 'corner' data as a 2d-array of corners[contours.size()][4];
                // Then calculate direction of all edges in all polygons
                auto base_of_i = i & ~3;          // index of corner 0 or any polygon
                auto next_idx = (i + 1) & 3;      // next index modulo 4
                auto diff = corners[i] - corners[base_of_i + next_idx];
                auto angle = std::atan2(diff.y, diff.x);  // 0 < angle <= pi
                angles_all.push_back(angle);
                if (angle < 0)
                    angle += (bin1_angle + bin2_angle);
                int slot = angle < bin1_angle ? 0 : (angle > bin2_angle ? 2 : 1);
                angle_histogram[slot].push_back(angle);
            }
            // Select the largest bin as representative
            std::vector<double> best_hist = angle_histogram[0];
            for (auto &h : angle_histogram)
                if (h.size() > best_hist.size())
                    best_hist = h;

            return trimmed_mean(best_hist, 0.50);
        }

        // Bins vectors by x/y coordinate to integral bins by locating a discontinuity
        // of at least median_dist/2 in the sorted array of x or y coordinates
        std::vector<int> histogram_dynamic(std::vector<point_xy> &vec, bool is_x, double median_dist)
        {
            if (vec.size() < 2)
                return{ 0 };

            std::vector<std::pair<double, size_t>> items;
            for (size_t i = 0; i < vec.size(); i++)
                items.emplace_back(is_x ? vec[i].x : vec[i].y, i);

            std::sort(items.begin(), items.end(),
                [](const std::pair<double, size_t> &a, const std::pair<double, size_t> &b) { return a.first < b.first; });

            // Place all successive elements to a group, if they fall within median_dist/2
            // from the start of the group
            //  - assign each output vector with an initial group number from 0..G
            //  - this group number will be adjusted if we later find out, that the difference between
            //    some groups are much larger (about 2x)
            auto result = std::vector<int>(vec.size());
            int group_number = 0;
            double max_dist = median_dist * 0.5;
            std::vector<double> middle_of_groups;
            std::vector<double> group(1, items[0].first);
            result[items[0].second] = 0;        // redundant ...
            for (size_t i = 1; i < items.size(); i++)
            {
                if (items[i].first - group.front() > max_dist)
                {
                    middle_of_groups.emplace_back(vector_median(group));
                    group_number++;
                    group.clear();
                }
                group.emplace_back(items[i].first);
                result[items[i].second] = group_number;
            }
            middle_of_groups.emplace_back(vector_median(group));
            if (middle_of_groups.size() <= 1 || middle_of_groups.size() > 6)
                return{};

            // Calculate the typical true distances within the groups
            auto group_diffs = std::vector<double>();
            for (size_t i = 1; i < middle_of_groups.size(); i++)
                group_diffs.emplace_back(middle_of_groups[i] - middle_of_groups[i - 1]);
            auto median_of_medians = vector_median(group_diffs);

            // Then find out if one or more distances are much larger than the others
            // - currently we only try to check if delta[i] ~ 2x median_delta
            // - theoretically we could be prepared for 3x difference or 3/2 difference
            auto remap = std::vector<int>(1);
            for (size_t i = 1; i < middle_of_groups.size(); i++)
            {
                if (middle_of_groups[i] - middle_of_groups[i - 1] > 1.5 * median_of_medians)
                    remap.emplace_back(remap.back() + 2);
                else
                    remap.emplace_back(remap.back() + 1);
            }

            for (auto &index : result)
                index = remap[index];

            return result;
        }

        // places N centroids in a grid of 6x4 points
        //  - the unoccupied slots have index -1
        //  - returns empty vector on failure
        std::vector<int> fit_centroids_to_grid(std::vector<point_xy> &centroids, double angle, double centroid_distance, point_xy center)
        {
            auto corrected = rotate_vector(centroids, -angle, center);
            auto hist_x = histogram_dynamic(corrected, true, centroid_distance);
            auto hist_y = histogram_dynamic(corrected, false, centroid_distance);

            if (hist_x.size() == 0 || hist_y.size() != hist_x.size())
                return{};

            auto result = std::vector<int>(24, -1);
            // In case we have toggled the histograms, we also need to invert one of the axis
            // (toggling means transposing, toggle + invert means rotation)
            auto must_swap = *std::max_element(hist_x.begin(), hist_x.end()) < *std::max_element(hist_y.begin(), hist_y.end());
            if (must_swap)
                std::swap(hist_x, hist_y);

            for (size_t i = 0; i < hist_x.size(); i++)
            {
                int y = must_swap ? 3 - hist_y[i] : hist_y[i];
                int x = hist_x[i];
                if (y < 0 || y > 3 || x < 0 || x > 5)
                    return{};       // illegal indexing
                int idx = x + y * 6;
                if (result[idx] >= 0)
                    return{};       // pigeon slot full
                result[idx] = static_cast<int>(i);
            }
            return result;
        }

        // locates missing items from an axis-aligned grid
        std::vector<point_xy> fill_missing_items(std::vector<point_xy> &centroids, std::vector<int> &indices)
        {
            std::vector<point_xy> grid(24);

            // calculate some statistics
            int hist_x[6] = { 0 };
            int hist_y[4] = { 0 };
            std::vector<roi_point> must_reconstruct;
            for (int j = 0; j < 4; j++)
            {
                for (int i = 0; i < 6; i++)
                {
                    int idx = indices[i + j * 6];
                    if (idx < 0)
                        must_reconstruct.emplace_back(i, j);
                    else
                    {
                        hist_x[i]++;
                        hist_y[j]++;
                        grid[i + j * 6] = centroids[idx];
                    }
                }
            }

            // Incrementally locate the item that can be interpolated from
            // largest amount of parametric equations
            while (must_reconstruct.size() != 0)
            {
                // build a vector of tuples of { quality, index }
                auto best_quality = 0;
                auto best_index = -1;

                for (size_t i = 0; i < must_reconstruct.size(); i++)
                {
                    auto &a = must_reconstruct[i];
                    int xa = hist_x[a._x];
                    int ya = hist_y[a._y];
                    int quality = (xa < 2 ? 0 : xa) + (ya < 2 ? 0 : ya);
                    if (quality > best_quality)
                    {
                        best_quality = quality;
                        best_index = static_cast<int>(i);
                    }
                }
                if (best_index < 0)
                    return{};

                // swap with last element and shrink vector by one
                auto best_item = must_reconstruct[best_index];
                std::swap(must_reconstruct[best_index], must_reconstruct.back());
                must_reconstruct.pop_back();

                point_xy interpolated(0, 0);
                int count = 0;
                if (hist_x[best_item._x] >= 2)
                {
                    parametric_model vertical;
                    for (int y = 0; y < 4; y++)
                        vertical.add_pt(grid[best_item._x + y * 6], y);

                    ++count;
                    interpolated += vertical(best_item._y);
                }

                if (hist_y[best_item._y] >= 2)
                {
                    parametric_model horizontal;
                    for (int x = 0; x < 6; x++)
                        horizontal.add_pt(grid[best_item._y * 6 + x], x);
                    ++count;
                    interpolated += horizontal(best_item._x);
                }

                interpolated *= (count == 2 ? 0.5 : 1.0);

                grid[best_item._x + best_item._y * 6] = interpolated;
                hist_x[best_item._x]++;
                hist_y[best_item._y]++;
            }
            return grid;
        }

        // Reference:
        // mathportal.org/calculators/statistic-calculator/correlation-and-regression-calculator.php
        static point_xy linear_regression(std::vector<point_xy> set)
        {
            double x = 0, y = 0, xy = 0, xx = 0;
            for (auto &p : set)
            {
                x += p.x;
                y += p.y;
                xy += p.x * p.y;
                xx += p.x * p.x;
            }
            double n = (double)set.size();
            double d = n * xx - x*x;
            if (n > 1 && d != 0)
                return point_xy((n*xy - x*y) / d, (y*xx - x*xy) / d);
            return point_xy(std::numeric_limits<double>::infinity(), x / n);
        }

        // Constructs a parametric model (x,y) = t * (dx,dy) + (x0,y0) from binned centroids
        // - allows estimation of missing centroids and generation of polygon corner points
        struct parametric_model
        {
            std::vector<point_xy> set_x;        // samples for x(t==0), x(t==1), ... x(t==5)
            std::vector<point_xy> set_y;        // samples for y(t==0), y(t==1), ... y(t==5)
            std::vector<point_xy> equations;    // line equations for {x,y} = t * delta + origin
            void add_pt(point_xy data, int n)
            {
                if (!(data == point_xy(0, 0)))  // do not add points to be interpolated
                {
                    set_x.emplace_back(n, data.x);
                    set_y.emplace_back(n, data.y);
                }
            }
            // evaluate the model at point n
            point_xy operator()(double n)
            {
                // cache the equations on first usage -- assumes that no further points are added later
                if (equations.size() != 2)
                {
                    equations.push_back(linear_regression(set_x));
                    equations.push_back(linear_regression(set_y));
                }
                return point_xy(equations[0].x * n + equations[0].y, equations[1].x * n + equations[1].y);
            }
        };

        // Extrapolate / interpolate a polygon around each centroid
        // make a parametric model for each horizontal line going approximately through the centroids
        // --*---*---*---*---*---*--  line 0
        // --*---*---*---*---*---*--  line 1
        // --*---*---*---*---*---*--  line 2
        // --*---*---*---*---*---*--  line 3
        // Sample those at int n= 0..5 +- 0.3, forming 12 vertical parametric models
        //   - sample those 12 parametric models at y={0,1,2,3} +- s
        //     to acquire 6x4 polygons
        //   - inverse barrel correct the polygons to remap back to (distorted) screen coordinates
        void calculate_polygons(std::vector<point_xy> &centers, point_xy scaled_center, double alpha_per_r2)
        {
            const double s = fill_ratio * 0.5;
            const int grid_height = 4;
            const int grid_width = 6;
            parametric_model horizontal_lines[grid_height];
            for (int y = 0; y < grid_height; y++)
                for (int x = 0; x < grid_width; x++)
                    horizontal_lines[y].add_pt(centers[x + y * grid_width], x);

            parametric_model vertical_lines[grid_width * 2];
            for (int y = 0; y < grid_height; y++)
            {
                for (int x = 0; x < grid_width; x++)
                {
                    vertical_lines[x * 2 + 0].add_pt(horizontal_lines[y](x - s), y);
                    vertical_lines[x * 2 + 1].add_pt(horizontal_lines[y](x + s), y);
                }
            }

            polygons.clear();
            polygons.reserve(grid_width * grid_height);

            for (int i = 0; i < grid_height; i++)
            {
                for (int j = 0; j < grid_width; j++)
                {
                    auto temp_poly = std::vector<point_xy>();

                    temp_poly.emplace_back(vertical_lines[j * 2 + 0](i - s));
                    temp_poly.emplace_back(vertical_lines[j * 2 + 1](i - s));
                    temp_poly.emplace_back(vertical_lines[j * 2 + 1](i + s));
                    temp_poly.emplace_back(vertical_lines[j * 2 + 0](i + s));
                    inverse_barrel(temp_poly, scaled_center, alpha_per_r2);

                    polygons.emplace_back();
                    auto &polygon = polygons.back();
                    for (auto &p : temp_poly)
                        polygon.emplace_back(p.x, p.y);
                }
            }
        }

        // generates the (undistorted) coordinates for x = {0..5}, y= {0..3} from homography matrix
        static point_xy point_reconstruction_function(point_xy p, double *a)
        {
            p = p * 2 + point_xy(-5, -3);
            auto z = a[6] * p.x + a[7] * p.y + a[8];
            auto x = a[0] * p.x + a[1] * p.y + a[2];
            auto y = a[3] * p.x + a[4] * p.y + a[5];
            return point_xy(x / z, y / z);
        }

        // Models the rotation/perspective transform of the macbeth chart as a function of x={0-5}, y={0-3} to screen
        // - suffers from remaining distortion after the barrel correction
        // - attempts to make a combined corrector with alpha/optical center estimator tried but failed
        std::vector<double> calculate_homography(std::vector<point_xy> &grid, point_xy mins, double median_distance)
        {
            struct homography_solver
            {
                std::vector<point_xy> data;

                homography_solver(std::vector<point_xy> &grid) : data(grid) { }
                double operator() (double *a)
                {
                    double sum = 0.0;
                    int idx = 0;
                    auto const zero = point_xy(0, 0);
                    for (int i = 0; i < 4; i++)
                    {
                        for (int j = 0; j < 6; j++)
                        {
                            if (data[idx] == zero)
                                continue;
                            auto projected = point_reconstruction_function(point_xy(j,i), a);
                            sum += norm2(data[idx] - projected);
                            idx++;
                        }
                    }
                    return sum;
                }
            } my_solver(grid);
            auto initial_values = std::vector<double>({
                median_distance / 2, 0, mins.x + median_distance * 2.5,
                0.0, median_distance / 2, mins.y + median_distance * 1.5,
                0, 0, 1.0 });
            auto result = nelder_mead_simplex(initial_values, my_solver, { 1500 });
            return result.second;
        }

        // Direct generation of polygons from homography -- not used currently
        void calculate_polygons_from_homography(point_xy scaled_center, double alpha_per_r2, std::vector<double> homography)
        {
            const double x = 0.5 * 0.5;
            auto patch = std::vector<point_xy>({
                { -x, -x }, { 0, -x }, { x, -x },  // top row
                { x, 0 }, { x, x },               // right column
                { 0, x }, { -x, x }, { -x, 0 }    // bottom row and left
            });

            std::vector<point_xy> polygon;
            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 6; j++)
                {
                    polygon = patch;
                    for (auto &p : polygon)
                        p = point_reconstruction_function(point_xy(j, i) + p, homography.data());

                    inverse_barrel(polygon, scaled_center, alpha_per_r2);

                    polygons.emplace_back();
                    auto &last = polygons.back();
                    // convert the polygon to integers
                    for (auto &p : polygon)
                        last.emplace_back(p.x, p.y);
                }
            }
        }
    };
}