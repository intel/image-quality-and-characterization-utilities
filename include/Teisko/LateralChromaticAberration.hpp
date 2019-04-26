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
#define NOMINMAX
#include "Teisko/Algorithm/DelaunayTriangulation.hpp"
#include "Teisko/Algorithm/LinearSpace.hpp"
#include "Teisko/Image/API.hpp"       // bayer_image_s
#include "Teisko/LensShading.hpp"     // use lensshading_grid as container
#include <cstdint>
#include <vector>

namespace Teisko
{
    // Domain specific algorithms for Lateral Chromatic Aberration Characterization

    using lca_grid = lensshading_grid<double>;

    /// Structure to hold + calculate center of green feature and
    /// relative displacements of red/blue feature centroids
    struct lca_feature_s
    {
        point_xy green;        // centroid of green feature
        point_xy diff_blue;    // calculated blue difference at point
        point_xy diff_red;     // calculated red difference at point
        double area = 0;       // area of feature

        lca_feature_s() = default;
        lca_feature_s(double a, point_xy r, point_xy g, point_xy b)
            : green(g), diff_blue(b - g), diff_red(r - g), area(a) { }

        // returns reference to the payload by index
        // - helps current implementation of filtering outliers
        // - int x should be in range 0..3
        double& operator[](int x)
        {
            double *ptrs[4] = { &diff_blue.x, &diff_blue.y, &diff_red.x, &diff_red.y };
            return *(ptrs[x & 3]);
        }

        // Pairwise scaling of feature
        // Used in barycentric interpolation: a * feature_a + b * feature_b + c * feature_c
        lca_feature_s operator *(const double x) {
            auto copy = *this;
            copy.area *= x;
            copy.green *= x;
            copy.diff_blue *= x;
            copy.diff_red *= x;
            return copy;
        }

        // Pairwise sum of features
        lca_feature_s operator +(const lca_feature_s &other) {
            auto copy = *this;
            copy.area += other.area;
            copy.green += other.green;
            copy.diff_blue += other.diff_blue;
            copy.diff_red += other.diff_red;
            return copy;
        }
    };

    /// How to use:
    /// Step 1)     lca_characterization implementation;
    /// Step 2)     for (auto &i : my_images)  implementation += i;
    ///                 - throws if image sizes mismatch
    /// Step 3)     auto output = implementation.collect();

    class lca_characterization
    {

    public:
        roi_point center() { return _center; }
        roi_point cell_size() { return _cell_size; }
        lca_characterization() : _init_done(false) {}

        /// Add image to model -- returns the newly added features
        std::vector<lca_feature_s> operator += (bayer_image_s<uint16_t> &image)
        {
            if (image._layout.is_dp_4x2_sensor())
            {
                preprocessor_s preprocessor;
                preprocessor.saturation = 65535;
                preprocessor.preprocess_reduce_dp_sensor(image);
            }

            validate_image_size(image._img.size());
            return add_image(image);
        }

        /// Collect results
        ///  - should we try to remove outliers by proximity before averaging?
        ///    `remove_outliers_by_median(9)`
        lca_grid collect()
        {
            // remove_outliers_by_median(9);
            auto output = average_results();
            remove_outliers_by_median_on_regular_grid(output);
            return output;
        }

    private:
        bool _init_done = false;        // Guard for allowing image size initialize only once
        roi_point _center;              // Optical center -- half of image dimensions
        roi_point _cell_size;           // Power of Two cell size producing max 64x64 grid
        roi_point _dims;                // Dimensions of the images (after 4x2->2x2 collapsing)

        // Called to assign image size related parameters on first run
        // and then to check that succeeding images are of the same size
        void validate_image_size(roi_point size)
        {
            if (_init_done)
            {
                if (_dims == size)
                    return;

                throw std::runtime_error("Image size mismatch");
            }
            const int max_grid_dims = 64;
            auto max_items = (size - 1) / (max_grid_dims - 1);
            _init_done = true;
            _dims = size;
            _center = _dims / 2;
            _cell_size = {
                static_cast<int>(next_power_of_two(max_items._x)),
                static_cast<int>(next_power_of_two(max_items._y))
            };
        }

        std::vector<std::vector<lca_feature_s>> _features;

        // Check if a feature is valid in LCA context
        //  - by area
        //  - by perimeter to area ratio (i.e. being simple enough)
        //  - the maximum area/perimeter^2 ratio is for circle -- so we don't want to limit those
        //    in the other end there are very thin ellipses or line segments or very complex shapes
        bool is_valid(region_props &feature, double max_area)
        {
            const double minimum_area = max_area * (0.5e-4);
            const double maximum_area = max_area * (8e-3);
            const double magic = 4.0 * 3.141592653589793;  // 4*pi -- that's area/perimeter^2 limit for circle
            if (feature.area < minimum_area || feature.area > maximum_area)
                return false;
            auto ratio = magic * feature.area / (feature.perimeter * feature.perimeter);
            return ratio > 0.5;             // empirical constant 0.5 for regular enough shape
        }

        // Quantifies image by threshold after edge detection
        // - horizontal + vertical filtering with 9 tap kernel
        void gradient_threshold(image<uint16_t> &img, double threshold = 0.4)
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
            img.foreach([max_level](uint16_t &dst, int &src) { dst = src >= max_level; }, temp_dst);
        }

        std::vector<region_props> get_valid_features(image<uint16_t> &channel)
        {
            std::vector<double> gaussian_7x7 = { 0.0702, 0.1311, 0.1907, 0.2161, 0.1907, 0.1311, 0.0702 };

            auto resampled = channel.convert_to();      // take a copy
            filter_separable(resampled, gaussian_7x7, REPLICATE);

            gradient_threshold(resampled, 0.4);

            auto count = bwlabel(resampled);
            int labels = 0;

            std::vector<region_props> features;
            features.reserve(count);

            auto dimensions = resampled.size();
            double max_area = dimensions._x * dimensions._y;
            for (int y = 0; y < dimensions._y; y++)
            {
                for (int x = 0; x < dimensions._x; x++)
                {
                    if (resampled(y, x) == labels + 1)
                    {
                        ++labels;
                        auto feature = get_regionprops(resampled, roi_point(x, y));
                        if (is_valid(feature, max_area))
                            features.push_back(feature);
                    }
                }
            }
            std::sort(features.begin(), features.end(), [](const region_props &a, const region_props &b)
            {
                return a.centroid_x < b.centroid_x;
            });
            return features;
        }

        // Given coarsely found set of green features locate the corresponding features
        // from each color plane using Otsu thresholding
        //  - discard all newly found features touching the bounding box
        //    or with blue to green / red to green magnitude larger than 7 pixels
        std::vector<lca_feature_s> find_matching_features(
            bayer_image_s<uint16_t> &bayer,
            image<uint16_t> &green,
            std::vector<region_props> &green_features)
        {
            auto clamp_point = [](roi_point p, roi_point size)
            {
                return roi_point(std::max(std::min(size._x - 1, p._x), 0),
                    std::max(std::min(size._y - 1, p._y), 0));
            };

            auto max_abs_distance = [](point_xy a, point_xy b)
            {
                auto diff = a - b;
                return std::max(std::abs(diff.x), std::abs(diff.y));
            };

            auto red = demosaic_bilinear(bayer, color_info_e::red);
            auto blue = demosaic_bilinear(bayer, color_info_e::blue);

            auto size = green.size();
            auto margin = roi_point(15);

            quantizer q_red, q_blue, q_green;

            double max_area = size._x * size._y;
            auto current_set_of_features = std::vector<lca_feature_s>();

            for (auto &g : green_features)
            {
                auto top_left = clamp_point(g.top_left - margin, size);
                auto roi_size = clamp_point(g.bot_right + margin, size) - top_left + 1;
                auto centroid_coarse = roi_point(g.centroid_x, g.centroid_y) - top_left;

                // Quantize each matching region (with reusing the histogram / context)
                auto bitonal_g = q_green(green.region(roi_size, top_left));
                auto bitonal_b = q_blue(blue.region(roi_size, top_left));
                auto bitonal_r = q_red(red.region(roi_size, top_left));

                // Require that the center of the original feature has been quantized to zero (a hole)
                if (bitonal_g(centroid_coarse) != 0 ||
                    bitonal_b(centroid_coarse) != 0 ||
                    bitonal_r(centroid_coarse) != 0)
                    continue;

                auto props_green = get_regionprops(bitonal_g, centroid_coarse);
                auto props_blue = get_regionprops(bitonal_b, centroid_coarse);
                auto props_red = get_regionprops(bitonal_r, centroid_coarse);

                // re-check that the region properties are within expected in all color channels
                if (!is_valid(props_green, max_area) ||
                    !is_valid(props_blue, max_area) ||
                    !is_valid(props_red, max_area))
                    continue;

                // Check that the areas are within +-25% of the median area
                // and that the centroids are located withing the maximum distance
                double area[3] = { props_green.area, props_red.area, props_blue.area };
                std::sort(area, area + 3);


                auto centroid_g = point_xy(props_green.centroid_x + top_left._x, props_green.centroid_y + top_left._y);
                auto centroid_r = point_xy(props_red.centroid_x + top_left._x, props_red.centroid_y + top_left._y);
                auto centroid_b = point_xy(props_blue.centroid_x + top_left._x, props_blue.centroid_y + top_left._y);

                const double max_distance = 7.0;
                const double max_area_ratio = 1.25;     // allows one area to be 25% larger than the other

                if (area[0] * max_area_ratio < area[1] || area[1] * max_area_ratio < area[2] ||
                    max_abs_distance(centroid_g, centroid_r) > max_distance ||
                    max_abs_distance(centroid_g, centroid_b) > max_distance)
                    continue;

                // else add the item to found items
                current_set_of_features.emplace_back(area[1], centroid_r, centroid_g, centroid_b);
            }
            return current_set_of_features;
        }

        // Locates the green features and their difference to red and blue features
        std::vector<lca_feature_s> add_image(bayer_image_s<uint16_t> &img)
        {
            auto resampled_green = demosaic_bilinear(img, color_info_e::green);
            auto green_features = get_valid_features(resampled_green);

            auto new_features = find_matching_features(img, resampled_green, green_features);

            _features.push_back(new_features);

            return new_features;
        }

        uint32_t next_power_of_two(uint32_t a)
        {
            if (a == 0)
                return 1;
            // Two examples:       a1 = 10000000         a2 = 100000011  (a1 is already power of two, a2 isn't)
            --a;                //      01111111              100000010   after subtraction
            a |= a >> 1;        //      01111111              110000011   after a |= a >> 1
            a |= a >> 2;        //      01111111              111100011   after next stage
            a |= a >> 4;        //      01111111              111111111   after oring in groups of four
            a |= a >> 8;        //  each iteration closes 1,2,4,8,..16 intermediate zeros/holes in the variable
            a |= a >> 16;       //  at this point all zeros have vanished, a is of form 0000...111111
            return ++a;         //      10000000             1000000000   Only one bit set!
        }

        struct heap_element_s
        {
            double distance;
            lca_feature_s *p;

            bool operator <(heap_element_s &other) { return distance < other.distance; }
            heap_element_s() : distance(0.0), p(nullptr) { }
            heap_element_s(lca_feature_s *ptr, point_xy point)
                : distance(norm2(ptr->green - point))
                , p(ptr) { }
        };

        // With starting point of left_idx == right_idx (within vec)
        // adds the N closest neighbors of the middle elements to a heap
        // - adjusts the "left_idx" and "right_idx"
        void init_heap(std::vector<heap_element_s> &heap, std::vector<lca_feature_s*> &vec,
            size_t &left_idx, size_t &right_idx, point_xy point)
        {
            auto n = heap.size();
            auto len = vec.size();
            for (size_t k = 0; k < n;)
            {
                if (left_idx - 1 < len)
                {
                    heap[k++] = heap_element_s(vec[--left_idx], point);
                    if (k == n)
                        break;
                }
                if (right_idx + 1 < len)
                    heap[k++] = heap_element_s(vec[++right_idx], point);
            }
            std::make_heap(heap.begin(), heap.end());
        }

        // `idx` is the last item added to heap (or points outside the vector)
        void trial_add(std::vector<heap_element_s> &heap, std::vector<lca_feature_s*> &vec, size_t &idx, point_xy p)
        {
            auto len = vec.size();
            // we scan from "middle" index to left or right
            // until we face an element whose X distance is greater than Kth smallest distance
            //  - after this point there are no more better candidates on that direction

            if (idx >= len || (vec[idx]->green.x - p.x) * (vec[idx]->green.x - p.x) > heap.front().distance)
            {
                idx = len;
                return;
            }

            heap_element_s trial(vec[idx], p);
            if (trial.distance < heap.front().distance)
            {
                std::pop_heap(heap.begin(), heap.end());   // moves the largest to the end
                heap.back() = trial;                       // replaces the largest element
                std::push_heap(heap.begin(), heap.end());  // rearranges the heap
            }
        }

        // Removes outliers from the full set of data
        // Steps 1) sort all items by x-coordinate
        //       2) for each value, locate the N nearest neighbours
        //           - because the set is sorted, we can stop moving left/right
        //             after there can not come any improvement
        //          3) calculate the median of absolute difference to median of set
        //             - reject the sample, if its more than 3 sigma away from the median
        //       4) After all samples are processed, fill the data
        // Also considered: put all samples in 2D grid (multimap, hash table)
        //  - then scan in a spiral path until no advances can be made
        void remove_outliers_by_median(size_t n)
        {
            // Sort all features in all images by x-coordinate
            std::vector<lca_feature_s*> vec;
            for (auto &image : _features)
                for (auto &f : image)
                    vec.emplace_back(&f);

            std::sort(vec.begin(), vec.end(), [](lca_feature_s *a, lca_feature_s *b)
            {
                return a->green.x < b->green.x;
            });

            auto len = vec.size();
            if (len < n)
                return;

            struct replacements
            {
                lca_feature_s* what;
                int parameter;
                double value;
                replacements(lca_feature_s* p, int channel, double val) : what(p), parameter(channel), value(val) { }
            };
            std::vector<replacements> outliers;

            // Create a heap to store pointers and distances of N closest items to each point in the set
            std::vector<heap_element_s> heap(n - 1);
            std::vector<double> parameter(n);
            std::vector<double> distance(n);

            for (size_t i = 0; i < len; i++)
            {
                auto left_idx = i;
                auto right_idx = i;

                init_heap(heap, vec, left_idx, right_idx, vec[i]->green);

                // Add items from left and right
                while (left_idx < len || right_idx < len)
                {
                    trial_add(heap, vec, --left_idx, vec[i]->green);
                    trial_add(heap, vec, ++right_idx, vec[i]->green);
                }

                double three_sigma = 4.447806655516805;

                for (int p = 0; p < 4; p++)
                {
                    parameter[0] = (*vec[i])[p];
                    for (size_t k = 1; k < n; k++)
                        parameter[k] = (*heap[k-1].p)[p];

                    std::sort(parameter.begin(), parameter.end());

                    for (size_t k = 0; k < n; k++)
                        distance[k] = std::abs(parameter[k] - parameter[3]);

                    std::sort(distance.begin(), distance.end());
                    if (std::abs((*vec[i])[p] - parameter[n / 2]) > three_sigma * distance[n / 2])
                        outliers.emplace_back(vec[i], p, parameter[n/2]);
                }
            }

            // and then we replace
            for (auto &rep : outliers)
                (*rep.what)[rep.parameter] = rep.value;
        }

        // Takes 3x3 neighbors replacing each bx, by, rx, ry
        // by the local median if the point is judged an outlier
        //   - when 'x' needs to be replaced, we replace also the corresponding 'y'
        //     and vice versa
        void remove_outliers_by_median_on_regular_grid(lca_grid &grid)
        {
            for (int ch = 0; ch < grid.channels; ch += 2)
            {
                auto copy_x = grid[ch].make_borders(1, 1, 1, 1, REPLICATE);
                auto copy_y = grid[ch+1].make_borders(1, 1, 1, 1, REPLICATE);
                int stride = copy_x._skip_y;
                int offset = grid.width * grid.height;
                grid[ch].foreach([offset, stride](double &dst, double &src_x, double &src_y)
                {
                    double three_sigma = 4.447806655516805;
                    double *dx = &dst;
                    double *dy = &dst + offset;
                    double *sx = &src_x;
                    double *sy = &src_y;
                    point_xy values[9];
                    for (int i = 0; i < 3; i++)
                    {
                        for (int j = 0; j < 3; j++)
                        {
                            values[i * 3 + j].x = sx[i * stride + j];
                            values[i * 3 + j].y = sy[i * stride + j];
                        }
                    }
                    // Sort values by x
                    std::sort(values, values + 9, [](point_xy a, point_xy b) { return a.x < b.x; });
                    double median_x = values[4].x;
                    double diffs[9];
                    for (int i = 0; i < 9; i++)
                        diffs[i] = std::abs(values[i].x - median_x);
                    std::sort(diffs, diffs + 9);
                    double median_abs_diff_x = diffs[4];
                    if (std::abs(*dx - median_x) > three_sigma * median_abs_diff_x)
                    {
                        for (int i = 0; i < 9; i++)
                            if (i != 4 && values[i].x == *dx)
                            {
                                values[i].x = median_x;
                                break;
                            }
                        *dx = median_x;
                        *dy = values[4].y;
                    }

                    // Sort values by y
                    std::sort(values, values + 9, [](point_xy a, point_xy b) { return a.y < b.y; });
                    double median_y = values[4].y;
                    for (int i = 0; i < 9; i++)
                        diffs[i] = std::abs(values[i].y - median_y);
                    std::sort(diffs, diffs + 9);
                    double median_abs_diff_y = diffs[4];
                    if (std::abs(*dy - median_y) > three_sigma * median_abs_diff_y)
                    {
                        for (int i = 0; i < 9; i++)
                            if (i != 4 && values[i].y == *dy)
                            {
                                values[i].y = median_y;
                                break;
                            }
                        *dy = median_y;
                        *dx = values[4].x;
                    }
                }, copy_x, copy_y);
            }
        }

        //  Interpolates each separate characterization (one per image) on a regular grid
        //  Combines the interpolations for final grid
        lca_grid average_results()
        {
            if (_dims._x == 0 || _dims._y == 0 || _features.size() == 0 || _cell_size._x == 0 || _cell_size._y == 0)
                return{};

            // Calculates the grid size so that it has a pow of two spacing and max 64 items per dimension
            auto grid_dims = 1 + (_dims + _cell_size - 2) / _cell_size;
            auto elements = grid_dims._x * grid_dims._y;
            auto grid = std::vector<lca_feature_s>(elements);
            auto best_distance = std::vector<double>(elements, std::numeric_limits<double>::infinity());
            auto best_area = std::vector<double>(elements, 0.0);

            // resample separate files/images with scattered_interpolator
            // - accumulate the interpolated sample points to grid
            for (auto &i : _features)
            {
                if (i.size() == 0)
                    continue;

                std::vector<point_xy> points;
                points.reserve(i.size());
                for (auto &p : i)
                    points.emplace_back(p.green);

                auto x_vec = linear_space(0.0, (double)_cell_size._x * (grid_dims._x - 1), grid_dims._x);
                auto y_vec = linear_space(0.0, (double)_cell_size._y * (grid_dims._y - 1), grid_dims._y);

                // Interpolate the data
                auto weights = scattered_interpolation(points, x_vec, y_vec);
                auto interpolated_grid = barycentric_interpolation(i, weights);

                for (int j = 0; j < elements; j++)
                {
                    auto &w = weights[j];
                    auto &lca_data = interpolated_grid[j];
                    // Select the largest feature by area to represent the inside convex hull sample
                    // - alternative: accumulate total area, then normalize
                    if (w.v.a != w.v.b)
                    {
                        if (lca_data.area > best_area[j])
                        {
                            best_area[j] = lca_data.area;
                            grid[j] = lca_data;
                        }
                    }
                    // else select the closest point to represent outside the convex hull sample
                    else if (best_area[j] == 0)
                    {
                        int y = j / grid_dims._x;
                        int x = j % grid_dims._x;
                        auto distance = norm2(points[w.v.a] - point_xy(x_vec[x], y_vec[y]));
                        if (distance < best_distance[j])
                        {
                            best_distance[j] = distance;
                            grid[j] = lca_data;
                        }
                    }
                }
            }
            // convert the array of struct to struct of arrays
            auto result = lca_grid(grid_dims._x, grid_dims._y, 4);
            for (int ch = 0; ch < 4; ch++)
            {
                auto *p = &result.grid[ch * elements];
                for (int j = 0; j < elements; j++)
                    p[j] = grid[j][ch];
            }
            return result;
        }
    };
}