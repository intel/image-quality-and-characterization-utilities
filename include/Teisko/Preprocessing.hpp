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
#include "Teisko/Image/Polyscale.hpp"
#include "Teisko/Image/Algorithms.hpp"
#include "Teisko/BayerImage.hpp"
#include "Teisko/BayerInfo.hpp"
#include "Teisko/Chromaticity.hpp"
#include "Teisko/Algorithm/Interpolate.hpp"     // needed by black_level_model
#include "Teisko/Algorithm/ConvexHull.hpp"
#include "Teisko/Algorithm/VectorMedian.hpp"

#include <cstdint>
#include <memory>
#include <vector>
#include <limits>               // numeric_limits<T>::min() used in int->uint conversion
#include <algorithm>            // std::sort, std::find
#include <map>

namespace Teisko
{
    // Shift an integral value by N bits while rounding
    // - used inside filter kernels to allow same handling for integer and fp types
    template <typename T, int shift, typename U>
    typename std::enable_if<std::is_integral<T>::value && std::is_integral<U>::value, T>::type
        round_shift(U value) {
        return static_cast<T>((value + (1 << (shift - 1))) >> shift);
    }

    // "Shift right" a floating point value -- without explicit rounding
    // - used inside filter kernels
    template <typename T, int shift>
    typename std::enable_if<std::is_floating_point<T>::value, T>::type
        round_shift(T value) {
        return value * static_cast<T>((T)1 / (1 << shift));
    }

    // Functions to perform averaging in a neighborhood
    //  - "I" = identity -- do nothing
    //  - "H" = horizontal average
    //  - "V" = vertical average
    //  - "X" = diagonal average
    //  - "O" = orthogonal average
    //  - "O2" = orthogonal average 2 pixels away
    enum class kernel_functions_e { I, H, V, X, O, O2 };
    template<typename T, kernel_functions_e type>
    struct kernel_function {
        static void func(T *, int, int skipx = 1);
    };

    // Partial specialization of "I" = identity -- do nothing
    template<typename T>
    struct kernel_function<T, kernel_functions_e::I> {
        static void func(T *, int, int) {}
    };

    // Partial specialization of "H" = horizontal average
    template<typename T>
    struct kernel_function<T, kernel_functions_e::H> {
        static void func(T *m, int, int skipx = 1) {
            m[0] = round_shift<T, 1>(m[-skipx] + m[skipx]);
        }
    };

    // Partial specialization of "V" = vertical average
    template<typename T>
    struct kernel_function<T, kernel_functions_e::V> {
        static void func(T *m, int skipy, int)
        {
            m[0] = round_shift<T, 1>(m[-skipy] + m[skipy]);
        }
    };

    // Partial specialization of "X" = diagonal average
    template<typename T>
    struct kernel_function<T, kernel_functions_e::X> {
        static void func(T *m, int skipy, int skipx = 1)
        {
            m[0] = round_shift<T, 2>(m[-skipy - skipx] + m[-skipy + skipx] + m[skipy - skipx] + m[skipy + skipx]);
        }
    };

    // Partial specialization of "O" = Orthogonal average of H and V neighbors
    template<typename T>
    struct kernel_function<T, kernel_functions_e::O> {
        static void func(T *m, int skipy, int skipx = 1)
        {
            m[0] = round_shift<T, 2>(m[-skipy] + m[skipy] + m[-skipx] + m[skipx]);
        }
    };

    // Partial specialization of "O2" = Orthogonal average of H and V neighbors 2 pixels away
    template<typename T>
    struct kernel_function<T, kernel_functions_e::O2> {
        static void func(T *m, int skipy, int skipx = 1)
        {
            m[0] = round_shift<T, 2>(m[-2*skipy] + m[2*skipy] + m[-2*skipx] + m[2*skipx]);
        }
    };

    /// Reconstructs Red/Blue components (if needed) in 4x4 IR image so that
    /// there are two G components and one requested color component in each quad
    /// - returns the color order in the all the quads after reconstruction
    template <typename T>
    std::vector<color_info_e> reconstruct_cfa_bilinear(image<T> &inout, bayer_info_s &info, color_info_e channel, uint32_t rep_scheme = ZERO_PAD)
    {
        auto pattern = info.get_color_pattern();
        if (info.is_4x4_ir_sensor())
        {
            if (channel == color_info_e::red || channel == color_info_e::blue)
            {
                // reconstruct red on blue and blue on red
                auto other_channel = channel != color_info_e::red ? color_info_e::red : color_info_e::blue;
                auto center = inout.region(inout.size() - 8, roi_point(4));
                int skipx = inout._skip_x * 2;
                int skipy = inout._skip_y * 2;
                for (int ch = 0; ch < 16; ch++)
                {
                    if (pattern[ch] == other_channel)
                    {
                        pattern[ch] = channel;      // fix the channel
                        center.subview(4, 4, ch / 4, ch % 4).foreach([skipx, skipy](T &middle)
                        {
                            auto *inout = &middle;
                            middle = round_shift<T, 2>(inout[-skipx] + inout[skipx] + inout[-skipy] + inout[skipy]);
                        });
                    }
                }
                // After finishing the middle, copy the border again for pass #2
                inout.copy_borders(4, 4, 4, 4, rep_scheme);
            }
            // get the top 2x2 portion of the color codes containing one 'channel'
            pattern = std::vector<color_info_e>({ pattern[0], pattern[1], pattern[4], pattern[5] });
        }
        return pattern;
    }

    // Top level function to perform I,H,V,X or O function to all four subpixels in a quad
    template <typename T, kernel_functions_e tl, kernel_functions_e tr, kernel_functions_e bl, kernel_functions_e br>
    struct quad_filter
    {
        int skipy;
        int skipx;
        quad_filter(int stride_y, int stride_x = 1) : skipy(stride_y), skipx(stride_x) { }

        void operator()(T &center)
        {
            kernel_function<T, tl>::func(&center, skipy, skipx);              // Apply the 'top_left' function to top left position
            kernel_function<T, tr>::func(&center + skipx, skipy, skipx);          // Apply the 'top_right' function to top right position
            kernel_function<T, bl>::func(&center + skipy, skipy, skipx);       // bottom left
            kernel_function<T, br>::func(&center + skipy + skipx, skipy, skipx);   // bottom right
        }
    };

    /// \brief  Demosaics R,G,B from bayer image using simple bilinear interpolation
    /// Processes one quad at a time with following interpolation patterns
    /// [1 +] for Green        [1 H] for Green on RGBI
    /// [+ 1] on rggb/rggi     [V X]     or NonGreen on any 2x2 pattern
    /// 4x4 RGBI patterns need to first recover the missing R/B values on every 2x2 quad
    /// rep_scheme would be typically zero pad or replicate_even (which copies full 4x4 patterns)
    template <typename T>
    image<T> demosaic_bilinear(bayer_image_s<T> &input, color_info_e channel, uint32_t rep_scheme = ZERO_PAD)
    {
        if (input._layout.is_dp_4x2_sensor())
            throw std::runtime_error("4x2 bayer images must be reduced to 2x2 format first");

        auto dims = get_dim(input._layout);
        auto copy = input._img.make_borders(dims, rep_scheme);
        auto pattern = reconstruct_cfa_bilinear(copy, input._layout, channel, rep_scheme);
        auto size = input._img.size();

        auto quad_view = copy.region(size, dims).subview(2, 2);
        int stride = copy._skip_y;

        auto const I = kernel_functions_e::I;
        auto const H = kernel_functions_e::H;
        auto const V = kernel_functions_e::V;
        // There is X defined in global namespace in Color.hpp at line 1418
        auto const X1 = kernel_functions_e::X;
        auto const O = kernel_functions_e::O;

        if (channel == pattern[0] && channel == pattern[3])
        {
            quad_view.foreach(quad_filter<T, I, O, O, I>(stride));
        }
        else if (channel == pattern[1] && channel == pattern[2])
        {
            quad_view.foreach(quad_filter<T, O, I, I, O>(stride));
        }
        else if (channel == pattern[0])
        {
            quad_view.foreach(quad_filter<T, I, H, V, X1>(stride));
        }
        else if (channel == pattern[1])
        {
            quad_view.foreach(quad_filter<T, H, I, X1, V>(stride));
        }
        else if (channel == pattern[2])
        {
            quad_view.foreach(quad_filter<T, V, X1, I, H>(stride));
        }
        else if (channel == pattern[3])
        {
            quad_view.foreach(quad_filter<T, X1, V, H, I>(stride));
        }
        else
            throw std::runtime_error("Unrecognized bayer pattern");

        return copy.region(size, dims);
    }

    // Top level function to perform I,X or O(2) function to all 4x4 subpixels in SVE pattern
    // There are 16 functions, but only 8 of those are independent
    template <typename T,
        kernel_functions_e k0, kernel_functions_e k1, kernel_functions_e k2, kernel_functions_e k3,
        kernel_functions_e k4, kernel_functions_e k5, kernel_functions_e k6, kernel_functions_e k7>
    struct sve_filter
    {
        int skipy;
        int skipx;
        sve_filter(int stride_y, int stride_x = 1) : skipy(stride_y), skipx(stride_x) { }
        sve_filter(roi_point strides) : skipy(strides._y), skipx(strides._x) { }

        void operator()(T &center)
        {
            kernel_function<T, k0>::func(&center + 0 * skipx + 0 * skipy, skipy, skipx);
            kernel_function<T, k1>::func(&center + 1 * skipx + 0 * skipy, skipy, skipx);
            kernel_function<T, k2>::func(&center + 2 * skipx + 0 * skipy, skipy, skipx);
            kernel_function<T, k3>::func(&center + 3 * skipx + 0 * skipy, skipy, skipx);
            kernel_function<T, k4>::func(&center + 0 * skipx + 1 * skipy, skipy, skipx);
            kernel_function<T, k5>::func(&center + 1 * skipx + 1 * skipy, skipy, skipx);
            kernel_function<T, k6>::func(&center + 2 * skipx + 1 * skipy, skipy, skipx);
            kernel_function<T, k7>::func(&center + 3 * skipx + 1 * skipy, skipy, skipx);
            kernel_function<T, k2>::func(&center + 0 * skipx + 2 * skipy, skipy, skipx);
            kernel_function<T, k3>::func(&center + 1 * skipx + 2 * skipy, skipy, skipx);
            kernel_function<T, k0>::func(&center + 2 * skipx + 2 * skipy, skipy, skipx);
            kernel_function<T, k1>::func(&center + 3 * skipx + 2 * skipy, skipy, skipx);
            kernel_function<T, k6>::func(&center + 0 * skipx + 3 * skipy, skipy, skipx);
            kernel_function<T, k7>::func(&center + 1 * skipx + 3 * skipy, skipy, skipx);
            kernel_function<T, k4>::func(&center + 2 * skipx + 3 * skipy, skipy, skipx);
            kernel_function<T, k5>::func(&center + 3 * skipx + 3 * skipy, skipy, skipx);
        }
    };

    /// \brief  Reconstructs short pixels for long (or other way round) for 2x2 RGGB, BGGR, GRBG and GBRG sensors
    /// Processes 4x4 patterns at a time applying either unity, orthogonal filter with skip factor of 2
    /// or diagonal (X) filter for green component
    ///  - There are 4 fundamental SVE patterns supported per quad -- top, down, left or right are reconstructed
    ///  - [* *] [- *]   - When the top 4x2 pixels use e.g. TOP + RIGHT, the bottom pixels use RIGHT + TOP order
    ///  - [- -] [- *]   - one can select between TR, RT, BR, RB, TL, LT, BL, LB
    ///  The Green Pixel must be located at the intersection of these T/B and L/R combinations limiting the total
    ///  amount of combination to those eight
    /// \param inout    bayer image to be reconstructed
    /// \param pattern  bit mask for channels to be reconstructed -- one = reconstruct, zero = don't reconstruct
    /// \returns true   if the sensor order was successfully reconstructed using the given pattern
    template <typename T>
    bool sve_demosaic_bilinear(bayer_image_s<T> &inout, uint16_t pattern)
    {
        auto is_sve_pattern = [](uint16_t pattern, int quad1_mask, int quad2_mask)
        {
            return pattern == (quad1_mask | (quad2_mask << 2) | (quad2_mask << 8) | (quad1_mask << 10));
        };

        if (pattern == 0 || pattern == 0xffff)
            return true;

        if (!inout._layout.is_2x2_sensor())
            return false;

        auto dims = roi_point(4, 4);
        auto copy = inout._img.make_borders(dims, REPLICATE_EVEN);
        auto size = inout._img.size();

        auto quad_view = copy.region(size, dims).subview(dims);
        auto skipxy = roi_point(copy._skip_x, copy._skip_y);

        auto green_first = inout._layout[0] == color_info_e::green;
        auto color_first = !green_first;
        enum { b = 0x30, l = 0x11, r = 0x22, t = 0x03 };  // bits set under 4x2 pattern

        auto const I = kernel_functions_e::I;
        auto const G = kernel_functions_e::X;       // type A of pixels = Green
        auto const C = kernel_functions_e::O2;      // type B of pixels = Color

        if (color_first && is_sve_pattern(pattern, t, r))                       // CG .G
            quad_view.foreach(sve_filter<T, C, G, I, G, I, I, I, C>(skipxy));  // .. .C

        else if (color_first && is_sve_pattern(pattern, r, t))                  // .G CG
            quad_view.foreach(sve_filter<T, I, G, C, G, I, C, I, I>(skipxy));  // .C ..

        else if (green_first && is_sve_pattern(pattern, b, r))                  // .. .C
            quad_view.foreach(sve_filter<T, I, I, I, C, C, G, I, G>(skipxy));  // CG .G

        else if (green_first && is_sve_pattern(pattern, r, b))                  // .C ..
            quad_view.foreach(sve_filter<T, I, C, I, I, I, G, C, G>(skipxy));  // .G CG

        else if (green_first && is_sve_pattern(pattern, t, l))                  // GC G.
            quad_view.foreach(sve_filter<T, G, C, G, I, I, I, C, I>(skipxy));  // .. C.

        else if (green_first && is_sve_pattern(pattern, l, t))                  // G. GC
            quad_view.foreach(sve_filter<T, G, I, G, C, C, I, I, I>(skipxy));  // C. ..

        else if (color_first && is_sve_pattern(pattern, b, l))                  // .. C.
            quad_view.foreach(sve_filter<T, I, I, C, I, G, C, G, I>(skipxy));  // GC G.

        else if (color_first && is_sve_pattern(pattern, l, b))                  // C. ..
            quad_view.foreach(sve_filter<T, C, I, I, I, G, I, G, C>(skipxy));  // G. GC
        else
            return false;

        inout._img.init_from(copy);
        return true;
    }

    /// Converts rgb to grayscale (in all channels)
    template <typename T>
    void convert_to_grayscale(rgb_image_s<T> &img)
    {
        auto r = img[rgb_color_e::red];
        auto g = img[rgb_color_e::green];
        auto b = img[rgb_color_e::blue];
        r.foreach([](T &r, T &g, T &b)
        {
            auto value = reduce_to<T>(0.299 * r + 0.587 * g + 0.114 * b);
            r = value;
            g = value;
            b = value;
        }, g, b);
    }

    /// \brief  Calculates  chromaticity from flat field bayer image
    /// \param  img         image to calculate
    /// \param  pattern     String encoding of Bayer sensor order e.g. "RGGB"
    /// \returns    Chromaticity R/G B/G and I/G
    template <typename T>
    chromaticity get_flatfield_white_point(bayer_image_s<T> &b_img)
    {
        using pixel_type = T;
        chromaticity_factory_f statistics{};
        const int kernel_size = 7;

        for (auto &i : b_img)
        {
            // Gets one color channel
            auto subchannel = b_img[i]; // get_sub_channel(i, img, pattern);
            // Break that channel to chunks of 7x7
            //  - restricting the size of image to full multiples of 7x7
            auto subsubchan = subchannel.subview(kernel_size, kernel_size)
                .region(subchannel._height / kernel_size, subchannel._width / kernel_size);

            /// functor to calculate a single median in 7x7 grid -- and take the max
            struct
            {
                pixel_type max;
                int stride_x;
                int stride_y;
                pixel_type arr[49];
                // 49 pixels in, 1 pixel out
                void operator() (pixel_type &top_left_corner)
                {
                    pixel_type *src = &top_left_corner;
                    pixel_type *dst = arr;
                    for (int row = 0; row < 7; row++)
                    {
                        for (int j = 0; j < 7; j++)
                            *dst++ = src[j * stride_x];
                        src += stride_y;
                    }
                    std::sort(arr, arr + 49);
                    max = std::max(arr[24], max);
                }
            } magic_max_median_7x7 = { 0, subchannel._skip_x, subchannel._skip_y, { 0 } };

            // Then we evaluate the max
            subsubchan.foreach(magic_max_median_7x7);
            // Update the chromaticity averaging functor
            statistics[b_img._layout[i]] += (double)magic_max_median_7x7.max;
        }
        // And finally get the stats
        return statistics;
    };

    /// \brief  Domain specific Black Level Interpolator with two dimensions
    struct black_level_model : private interpolate_x_dim<float>
    {
        using type = float;
        /// Add data at a grid point
        black_level_model& set(type ag, type exp, std::vector<type> data)
        {
            interpolate_x_dim<type>::set({ ag, exp }, data);
            return *this;
        }

        /// Interpolates the model data by at point
        std::vector<type> get(type ag, type exp)
        {
            return interpolate_x_dim<type>::get({ ag, exp });
        }
    };

    /// Removes black level from teisko_image<T> with associated bayer_info
    template <typename T, typename U>
    static void remove_black_level(bayer_image_s<T> &img, std::vector<U> &black_level)
    {
        if (black_level.size() == 0)
            return;
        if (img._layout.get_channels() != black_level.size())
        {
            throw std::runtime_error("Black Level vector size doesn't match image type");
        }

        for (auto i : img)
        {
            T item = (T)black_level[i];
            img[i].foreach([item](T &pixel) {
                pixel = (pixel > item) ? pixel - item : 0;
            });
        }
    }

    /// Single entry of RGB IR contamination correction grid
    /// - also a single entry for LSC correction grid
    struct rgb_ir_struct
    {
        chromaticity _chromaticity;         ///< r_per_g, b_per_g, i_per_g of the contamination
        uint32_t _grid_width;               ///< Width of correction grid (per channel)
        uint32_t _grid_height;              ///< Height of correction grid (per channel)
        std::vector<char> _grid_indices;    ///< Index of Correction Channel -- 4x4 layout
        std::vector<float> _grid_data;      ///< W*H *N_Channels worth of correction grid

        /// \brief rgb_ir_contamination Creates a new RGB IR contamination data set
        rgb_ir_struct(uint32_t width = 0, uint32_t height = 0)
            : _grid_width(width)
            , _grid_height(height)
        {}

        rgb_ir_struct(double r, double b, double i, uint32_t w, uint32_t h,
            std::vector<char> indices, std::vector<float> &data)
            : _chromaticity(r, b, i), _grid_width(w), _grid_height(h)
            , _grid_indices(indices), _grid_data(data)
        {}
    };

    /// Container and algorithm for RGB IR Contamination Removal
    ///  - dependencies to:  image<T>
    ///                      polynomial (up)scaling / modeling
    ///                      rgb_ir_struct, chromaticity
    ///                      get_sensor_pattern ... and other sensor info
    ///                      average / chromaticity_factory_f
    class rgb_ir_contamination
    {
    public:
        /// \brief Adds one rgb IR contamination data point to model
        /// if it's unique enough and consistent with existing data points
        /// \param  data_point  Correction grid associated with chromaticity point
        /// \returns  true if the data point passes consistency tests
        bool add_data_grid(rgb_ir_struct &data_point)
        {
            // Do not allow zero grid size
            auto vec_size = data_point._grid_data.size();
            auto grid_size = data_point._grid_width * data_point._grid_height;
            if (vec_size == 0 || grid_size == 0 || vec_size % grid_size != 0)
                return false;

            _channel_count = (uint32_t)(vec_size / grid_size);

            // Check that the added data point dimensions match those already added
            if (_rgb_ir.size() != 0)
            {
                auto first_grid = _rgb_ir.front();
                if (data_point._grid_width != first_grid._grid_width ||
                    data_point._grid_height != first_grid._grid_height ||
                    vec_size != first_grid._grid_data.size() ||
                    data_point._grid_indices != first_grid._grid_indices)
                    return false;
            }

            // Ignore duplicate entries by chromaticity distance
            for (auto &rgb : _rgb_ir)
            {
                const double epsilon = 1e-10;
                if (rgb._chromaticity.dist2(data_point._chromaticity) < epsilon)
                    return true;
            }
            _rgb_ir.push_back(data_point);

            return true;
        }

        /// \brief Removes the IR contamination from image
        /// \param img      Image to modify with black level removed
        template <typename T>
        void remove_ir_contamination(bayer_image_s<T> &img)
        {
            //  Step 0:  Determine the IR channels for the layout
            //  Step 1:  Calculate image chromaticities
            //  Step 2:  Interpolate between different RGB-IR contamination grids based on image chromaticity
            //  Step 3:  Upscale IR contamination grid to image size using 4th order polynomial fit
            //  Step 4:  Remove the effect of IR channel from each non-IR channel

            auto &layout = img._layout;
            if (!layout.is_ir_sensor())
                return;

            std::vector<int> ir_ch_indices;

            if (layout.is_4x4_ir_sensor())
            {
                // The first IR channel will be either at 0,1, 5 or 6
                // This is then copied over the other 3 2x2 sub-blocks
                auto ir_pos = (char)layout.locate_color(bayer_info_s::color_info_e::ir);
                ir_ch_indices = {
                    ir_pos, ir_pos, ir_pos + 2, ir_pos + 2,
                    ir_pos, ir_pos, ir_pos + 2, ir_pos + 2,
                    ir_pos + 8, ir_pos + 8, ir_pos + 10, ir_pos + 10,
                    ir_pos + 8, ir_pos + 8, ir_pos + 10, ir_pos + 10
                };
            }
            else if (layout.is_2x2_ir_sensor())
            {
                // locate the position of the I channel repeating it to all four locations
                auto ir_pos = (char)layout.locate_color(bayer_info_s::color_info_e::ir);
                ir_ch_indices = { ir_pos, ir_pos, ir_pos, ir_pos };
            }
            else
            {
                throw std::runtime_error("Unexpected IR sensor configuration");
            }

            // Calculate Image Chromaticity
            chromaticity white_point = get_flatfield_white_point(img);

            // Interpolate the grids in Chromaticity space
            auto interpolated_grid = interpolate(white_point);

            if (interpolated_grid._grid_data.size() == 0)
                throw std::runtime_error("Interpolation by white point failed");

            // Upscale IR to bayer image dimensions
            auto upscaled_grid = upscale(interpolated_grid, get_channel_dim(img));

            // get the top 2x2, 4x2 or 4x4 portion of the grid indices
            auto correction_grid_remapped = remap_4x4_vector(upscaled_grid._grid_indices, img._layout);

            //  Step 4:  Remove IR contamination from the image data (J = I - wX),
            //                                                  where J = restored pixel value
            //                                                        I = original pixel value
            //                                                        X = IR pixel value
            //                                                        w = associated IR weight from grid
            for (auto &i : img)
            {
                auto correction_grid_ind = correction_grid_remapped[i];
                if (correction_grid_ind < 0 || correction_grid_ind >(int)_channel_count)
                    continue; // No IR contamination is removed from this channel

                auto color_chan = img[i];   // channel needing contamination removal
                auto ir_channel = img[ir_ch_indices[i]];    // IR channel
                auto ir_contamination = get_channel(upscaled_grid, correction_grid_ind);

                color_chan.foreach([](T &dst, T &ir, float &weight)
                {
                    auto ir_reduction = static_cast<float>(ir) * weight;
                    auto pixel = static_cast<float>(dst);
                    dst = (pixel < ir_reduction) ? (T)0 : reduce_to<T>(pixel - ir_reduction);
                }, ir_channel, ir_contamination);
            }
        }

        rgb_ir_contamination() : _channel_count(0) { }

        bool is_initialized() {
            return _channel_count > 0 && _rgb_ir.size() > 0;
        }

        std::vector<rgb_ir_struct> _rgb_ir;      ///< Stored data

    private:
        uint32_t _channel_count;                 ///< For internal use -- updated from first ADD

        // Interprets contamination removal grid as Teisko image
        image<float> get_channel(rgb_ir_struct &orig, uint32_t channel)
        {
            uint32_t channel_size = orig._grid_height * orig._grid_width;
            return image<float>(orig._grid_height, orig._grid_width,
                orig._grid_data.data() + channel * channel_size);
        }

        // Resamples given contamination removal grid to full image size
        rgb_ir_struct upscale(rgb_ir_struct &orig, roi_point size)
        {
            rgb_ir_struct result(size._x, size._y);
            result._chromaticity = orig._chromaticity;
            result._grid_indices = orig._grid_indices;
            result._grid_data.resize(size._x * size._y * _channel_count);
            for (uint32_t i = 0; i < _channel_count; i++)
            {
                // Create scaler and input/output parameters
                poly_4_scaler_2d<float> polynomial;
                auto in_grid = get_channel(orig, i);
                auto out_grid = get_channel(result, i);
                polynomial.scale(in_grid, out_grid);
            }
            return result;
        }

        // Returns RGB grid interpolated from all input grids by distance to given white point
        rgb_ir_struct interpolate(chromaticity white_point)
        {
            size_t sz = _rgb_ir.size();
            if (sz == 0)
                return rgb_ir_struct{};

            // Use the first item in the vector to copy dimensions
            auto grid = _rgb_ir.front();
            size_t grid_size = grid._grid_data.size();

            rgb_ir_struct result(grid._grid_width, grid._grid_height);
            result._grid_indices = grid._grid_indices;
            result._chromaticity = white_point;
            result._grid_data.resize(grid_size);

            // Calculate the squared distance from each item to white_point
            std::vector<double> distance_vector;
            for (auto &dp : _rgb_ir)
            {
                distance_vector.emplace_back(dp._chromaticity.dist2(white_point));
            }

            for (size_t i = 0; i < sz; ++i)
            {
                double distance_ratio = 0;
                const double distance_threshold = 0.000001;
                for (size_t j = 0; j < sz; ++j)
                {
                    if (i != j)
                    {
                        // If one of the points has very short distance, weight calculations are not needed.
                        if (distance_vector[j] >= distance_threshold)
                            distance_ratio = distance_ratio + distance_vector[i] / distance_vector[j];
                        else
                            distance_ratio = distance_ratio + distance_vector[i] / distance_threshold;
                    }
                }
                double weight = 1.0 / (1.0 + distance_ratio);

                // Accumulate the input data with a weight to final grid
                auto output = result._grid_data.data();
                auto input = _rgb_ir[i]._grid_data.data();
                for (size_t j = 0; j < grid_size; ++j)
                {
                    output[j] += (float)(input[j] * weight);
                }
            }
            return result;
        }
    };

    // Clamps ROI between {0,0} and {width-1, height-1}
    inline void fix_roi_to_image_boundaries(roi_point &roi, roi_point &max_dims)
    {
        if (max_dims._x == 0 || max_dims._y == 0)
        {
            roi = { 0, 0 };
        }
        else
        {
            roi._x = std::min(std::max(roi._x, 0), max_dims._x - 1);
            roi._y = std::min(std::max(roi._y, 0), max_dims._y - 1);
        }
    }

    /// Returns median of each channel for each ROI
    /// Roi top_left/bot_right are scaled to full image dimensions (e.g. 0..4207, 0..3120)
    /// Returns 4 values for 2x2 sensors, 8 values for 4x2 sensor etc. for each ROI rectangle
    template <typename T>
    std::vector<T> get_patch_white_point(bayer_image_s<T> &img,
        std::vector<roi_point> top_left,
        std::vector<roi_point> bot_right)
    {
        if (top_left.size() != bot_right.size())
            throw std::runtime_error("input vectors don't match in size");
        std::vector<T> result;
        result.reserve(img._layout.get_channels() * top_left.size());

        auto ch_dims = get_dim(img._layout);
        auto half_dims = ch_dims / 2;
        auto channel_max_dim = get_dim(img) / ch_dims;

        for (size_t i = 0; i < top_left.size(); i++)
        {
            // scale the ROI to channel resolution
            auto tl = (top_left[i] + half_dims) / ch_dims;
            auto br = (bot_right[i] + half_dims) / ch_dims;

            fix_roi_to_image_boundaries(tl, channel_max_dim);
            fix_roi_to_image_boundaries(br, channel_max_dim);

            auto roi_size = (br - tl) + 1;
            if (roi_size._x <= 0 || roi_size._y <= 0)
                throw std::runtime_error("ROI is zero or negative");

            for (auto &j : img)
            {
                auto chan = img[j].region(roi_size._y, roi_size._x, tl._y, tl._x)
                    .to_vector();

                result.emplace_back(vector_median(chan));
            }
        }
        return result;
    }

    struct linearization_s
    {
        std::vector<int8_t> grid_indices;   // 4x4 grid of channels
        std::vector<std::vector<roi_point>> knee_points;

        linearization_s() = default;

        void linearize(bayer_image_s<uint16_t> &img, int max_val = 65535)
        {
            if (grid_indices.size() != 16 || knee_points.size() == 0)
                return;

            for (auto ch : img) // loop over channel indices
            {
                auto lut = make_lut(img, ch, max_val);
                auto size = lut.size();
                if (size == 0)
                    continue;
                auto chan = img[ch];

                uint16_t *lut_ptr = &lut[0];
                uint16_t max_x = static_cast<uint16_t>(size - 1);
                chan.foreach([lut_ptr, max_x](uint16_t &d) { d = lut_ptr[d > max_x ? max_x : d]; });
            }
        }

        // Convert the knee points to LUT
        std::vector<uint16_t> make_lut(bayer_image_s<uint16_t> &img, int channel, int max_output)
        {
            auto clip = [](int val, int max_val)
            {
                return val <= 0 ? 0 : (val >= max_val ? max_val : val);
            };
            if (grid_indices.size() != 16)
                return{};

            auto grid = remap_4x4_vector(grid_indices, img._layout);
            if (channel < 0 || channel >= static_cast<int>(grid.size()) ||
                grid[channel] < 0 || static_cast<size_t>(grid[channel]) > knee_points.size())
                return{};

            auto &kneepoints = knee_points[grid[channel]];
            auto ks = kneepoints.size();
            if (ks < 2)
                return{};

            auto max_input = 65535;
            auto last_x = clip(kneepoints.back()._x, max_input);
            auto first_y = clip(kneepoints.front()._y, max_output);

            auto lut = std::vector<uint16_t>(last_x + 1, static_cast<uint16_t>(first_y));

            // Interpolate between the rest of knee points
            for (size_t j = 1; j < ks; j++)
            {
                int x0 = kneepoints[j - 1]._x;
                int span = kneepoints[j]._x - x0;
                int y0 = kneepoints[j - 1]._y;
                int slope = kneepoints[j]._y - y0;
                if ((x0 < 0) || ((x0 + span) > (last_x + 1)))
                    throw std::runtime_error("Invalid kneepoints");

                for (int i = 0; i < span; i++)
                {
                    lut[x0 + i] = static_cast<uint16_t>(clip(y0 + (i * slope + span / 2) / span, max_output));
                }
            }
            lut.back() = static_cast<uint16_t>(clip(kneepoints.back()._y, max_output));
            return lut;
        }
    };

    // Commonly used preprocessing operations for
    //  - black level removal
    //  - rgb ir contamination removal
    //  - clamping
    //  - bayer pattern restoration (by nearest neighbor interpolation)
    //  - coming up: linearization
    //  - coming up: simple demosaicing
    struct preprocessor_s
    {
        rgb_ir_contamination rgb_ir_model{};
        black_level_model bl_model{};
        float saturation{ 1023.0f };    // for 10-bit sensor
        uint16_t sve_matrix{ 0 };       // longest exposure

        /// Reconstructs MD, SVE-MD and SVE images
        bool sve_md_demosaic_bilinear(bayer_image_s<uint16_t> &img)
        {
            preprocess_reduce_dp_sensor(img);
            return sve_demosaic_bilinear(img, sve_matrix ^ 0xffff);        // reconstruct the _short_ pixels
        }

        /// Calculates the Nth percentile of green channel from uint16_t typed image.
        /// Generalization to float type is a bit harder and can wait
        uint16_t get_green_percentile(bayer_image_s<uint16_t> &img, double percentile)
        {
            sve_demosaic_nearest_neighbor(img);
            // make enough storage for long and short histograms
            const int quantization_levels = 65536;
            std::vector<uint32_t> histogram(quantization_levels);
            uint32_t *hist = histogram.data();
            uint32_t max_count = 0;
            for (auto i : img)
            {
                if (img._layout[i] == color_info_e::green)
                {
                    auto channel = img[i];
                    max_count += channel._width * channel._height;
                    channel.foreach([hist](uint16_t &src) {
                        hist[src]++;
                    });
                }
            }

            // Find the Nth item based on percentile -- clamp between 1..max_count
            // to guarantee that we find minimum and maximum
            auto find_count = std::min(max_count, std::max(1u, (uint32_t)(max_count * (percentile / 100.0) + 0.5)));
            uint32_t sum = 0;
            for (size_t i = 0; i < quantization_levels; i++)
            {
                sum += hist[i];
                if (sum >= find_count)
                    return static_cast<uint16_t>(i);
            }
            return std::numeric_limits<uint16_t>::max();
        }

        // Fixes the black level data to proper format (0->N, 4->8 or 4->16)
        // or throws, if no conversion is available
        std::vector<float> get_black_level(bayer_info_s &layout, float ag, float exp)
        {
            auto bl_vec = bl_model.get(ag, exp);

            if (bl_vec.size() == layout.get_channels())
                return bl_vec;

            if (bl_vec.size() == 0)
                return std::vector<float>(layout.get_channels());

            if (bl_vec.size() == 4 && layout.is_dp_4x2_sensor())
                return{
                bl_vec[0], bl_vec[0], bl_vec[1], bl_vec[1],
                bl_vec[2], bl_vec[2], bl_vec[3], bl_vec[3]
            };

            if (bl_vec.size() == 4 && layout.is_4x4_ir_sensor())
                return{
                bl_vec[0], bl_vec[1], bl_vec[0], bl_vec[1],
                bl_vec[2], bl_vec[3], bl_vec[2], bl_vec[3],
                bl_vec[0], bl_vec[1], bl_vec[0], bl_vec[1],
                bl_vec[2], bl_vec[3], bl_vec[2], bl_vec[3]
            };
            throw std::runtime_error("Can't convert black level model to given channel size");
        }

        /// Removes black level and stretches the value between 0 and saturation
        /// Multiplier of 2 is used for DP 4x2 sensors where the stretching is done
        /// before summing up left and right channels
        template<typename T>
        void stretch_to_max_level(bayer_image_s<T> &img, std::vector<float> &bl, float sat, int multiplier)
        {
            auto bl_size = img._layout.get_channels();
            if (bl.size() != bl_size)
                throw std::runtime_error("Black Level vector size doesn't match image type");

            for (auto i : img)
            {
                float offset = bl[i];
                float stretch = sat / (sat - multiplier * offset);
                img[i].foreach([offset, stretch, sat](T &inout)
                {
                    auto result = ((float)inout - offset) * stretch;
                    if (result < 0)
                        result = 0;
                    if (result > sat)
                        result = sat;
                    inout = (T)result;
                });
            }
        }

        /// Calculates the chromaticity of preprocessed image in a ROI
        /// first each component in 2x2, 4x2 or 4x4 grid is filtered though a median
        template <typename T>
        chromaticity get_patch_white_point(bayer_image_s<T> &img,
            roi_point top_left,
            roi_point bot_right)
        {
            chromaticity_factory_f statistics{};
            auto tl_vec = std::vector<roi_point>(1, top_left);
            auto br_vec = std::vector<roi_point>(1, bot_right);
            auto channels = Teisko::get_patch_white_point(img, tl_vec, br_vec);

            if (channels.size() != (size_t)img._layout.get_channels())
                throw std::runtime_error("Conflict in channel sizes");  // should not happen

            for (auto i : img._layout)
                statistics[img._layout[i]] += channels[i];

            return statistics;
        }

        template <typename T>
        void preprocess_reduce_dp_sensor(bayer_image_s<T> &dp_image)
        {
            T sat = (T)saturation;
            if (!dp_image._layout.is_dp_4x2_sensor())
                return;
            auto left = dp_image._img.subview(1, 2, 0, 0);
            auto right = dp_image._img.subview(1, 2, 0, 1);
            left.foreach([sat](T &left, T &right)
            {
                left += right;
                if (left > sat)
                    left = sat;
            }, right);
            dp_image._img = left;   // replace the image by its subview
            dp_image._layout = bayer_info_s(reduce_dp_layout(dp_image._layout));
        }

        /// Normalizes bayer pattern to regular rggb, bggr, grbg abd gbrg formats
        /// by resampling I samples and summing up Left/Right from 4x2 DP sensors
        template <typename T>
        void preprocess_fn0(bayer_image_s<T> &image)
        {
            T sat = (T)saturation;
            auto &img = image._img;
            auto &layout = image._layout;

            if (layout.is_dp_4x2_sensor())
            {
                preprocess_reduce_dp_sensor(image);
                return;
            }

            // reconstruct IR (nothing done for non ir -sensors)
            reconstruct_bayer_cfa(image);

            // saturate to max level
            img.foreach([sat](T &pixel)
            {
                if (pixel > sat)
                    pixel = sat;
            });
        }

        /// intra channel interpolating based on the SVE matrix
        /// first locate the horizontal, then vertical and finally diagonal neighbor
        /// sve_pattern must be pre-shifted to contain the data at 0,2,8th and 10th bit
        template <typename T>
        void fix_color_channel(image<T> &&channel, uint16_t sve)
        {
            int masks[4] = {
                (sve & (1 << 0)), (sve & (1 << 2)),
                (sve & (1 << 8)), (sve & (1 << 10))
            };

            for (auto i : { 0u, 1u, 2u, 3u })
            {
                if (masks[i] != 0) continue;    // do nothing for long channel
                for (auto j : { 1u, 2u, 3u })
                {
                    // for channel i = 0 --> the search channels are 0 ^ {1,2,3} = 1,2,3
                    //             i = 1 --> 1 ^ {1,2,3}  = 0,3,2
                    //             i = 2 --> 2 ^ {1,2,3}  = 3,0,1
                    //             i = 3 --> 3 ^ {1,2,3}  = 2,1,0
                    auto search_channel = i ^ j;
                    if (masks[search_channel] != 0)
                    {
                        //auto long_chan = ;
                        channel.subview(2, 2, i / 2, i % 2)
                            .init_from(channel.subview(2, 2, search_channel / 2, search_channel % 2));
                        break;
                    }
                }
            }
        }

        // Returns left or right side of 4x2 DP image
        template <typename T>
        bayer_image_s<T> get_dp_image_view(bayer_image_s<T> &image, int side)
        {
            auto layout = reduce_dp_layout(image._layout);
            auto left_or_right = image._img.subview(1, 2, 0, side);
            return bayer_image_s<T>(left_or_right, layout);
        }

        /// reconstruct short SVE pixels from neighbouring pixels of same kind
        /// - can handle 2x2 regular layouts (and 4x2 layouts by recursion)
        /// - we currently support just a few models for recovery:
        ///   - gb is recovered from gr or vice versa
        ///     red can have 0..2 short subchannels
        ///     blue can have 0..2 short subchannels
        ///     - the short channel is recovered from horizontal, vertical or diagonal neighbor
        template <typename T>
        bool sve_demosaic_nearest_neighbor(bayer_image_s<T> &image)
        {
            if (sve_matrix == 0 || sve_matrix == 0xffff)
                return true;
            if (image._layout.is_dp_4x2_sensor())
            {
                for (auto side : { 0, 1 })
                {
                    auto view = get_dp_image_view(image, side);
                    if (sve_demosaic_nearest_neighbor(view) == false)
                        return false;
                }
                return true;
            }
            if (!image._layout.is_2x2_sensor())
                return false;
            /// we use deep bit magic to find out if all SVE bits are zero or one for a given channel
            /// when the channel indices are 0 1  and the bits are located at 0 1
            ///                              2 3                              4 5
            /// then the horizontal neighbor in all positions is i ^ 1
            /// then the vertical neighbor is i ^ 2 or i ^ 4 at the right table
            /// then the diagonal neighbor is i ^ 3 or i ^ 5 at the right table

            auto green_index = image._layout[0] == color_info_e::green ? 0 : 1;
            uint16_t pixel_mask = (1 << 0) | (1 << 2) | (1 << 8) | (1 << 10);
            uint16_t green_mask0 = pixel_mask << green_index;           // shifted by 0 or 1
            uint16_t green_mask1 = pixel_mask << (green_index ^ 5);     // shifted by 5 or 4

            // Check if we can recover red & blue channels
            uint16_t color_mask_top = pixel_mask << (green_index ^ 1);     // top color exposure times
            uint16_t color_mask_bot = pixel_mask << (green_index ^ 4);     // bottom color exposure times
            if ((sve_matrix & color_mask_top) == 0 || (sve_matrix & color_mask_bot) == 0)
                return false;

            // First we recover the green channel from the other green channel
            if (((sve_matrix & green_mask0) == 0) && ((sve_matrix & green_mask1) == green_mask1))
            {
                // top row is short, bottom row is long -- copy bottom row to top row
                image[green_index].init_from(image[green_index ^ 3]);
            }
            else if (((sve_matrix & green_mask1) == 0) && ((sve_matrix & green_mask0) == green_mask0))
            {
                // bottom row is short, top row is long -- copy top to bottom
                image[green_index ^ 3].init_from(image[green_index]);
            }
            else
            {
                // we can try to recover the green channels independently just like the non-green channels
                if (((sve_matrix & green_mask0) == 0) || ((sve_matrix & green_mask1) == 0))
                    return false;
                fix_color_channel(image[green_index], sve_matrix >> green_index);
                fix_color_channel(image[green_index ^ 3], sve_matrix >> (green_index ^ 5));
            }

            fix_color_channel(image[green_index ^ 1], sve_matrix >> (green_index ^ 1));
            fix_color_channel(image[green_index ^ 2], sve_matrix >> (green_index ^ 4));
            return true;
        }

        // Removes the black level global from the data and stretches the data between 0 ... maximum
        template <typename T>
        void preprocess_fn1(bayer_image_s<T> &image, float ag = 1.0f, float exp = 300.0f)
        {
            int multiplier = image._layout.is_dp_4x2_sensor() ? 2 : 1;
            auto bl_vec = get_black_level(image._layout, ag, exp);
            stretch_to_max_level(image, bl_vec, saturation, multiplier);
        }

        //     fn3a(image, meta,
        //          colorOrder, blackLevel, saturationLevel, irCalibration, bayerBlackLevelMultiplier)
        //
        //  - the black levels are adjusted and the intensities stretched between 0 ... saturation_level
        //  - in case of RGB IR sensor, the IR contamination is removes, and the I channel
        //    is resampled to nearest Green channel, changing the bayer layout to regular 2x2 patterns
        //  - The reference image is converted to 1/x image working as a gain
        //  - The original image is multiplied with the inverted image
        //
        //  And finally
        //  - DP 4x2 sensor layout images are collapsed to 2x2 layout
        //  - it's checked that no pixel value exceed the saturation level
        template <typename T>
        void preprocess_fn3a(bayer_image_s<T> &ff_img, float ag = 1.0f, float exp = 300.0f)
        {
            // We need to calculate the gain tables etc at least as floats -- take a copy
            bayer_image_s<float> img = ff_img;

            int multiplier = ff_img._layout.is_dp_4x2_sensor() ? 2 : 1;
            auto bl_vec = get_black_level(ff_img._layout, ag, exp);
            stretch_to_max_level(img, bl_vec, saturation, multiplier);

            // IR contamination removal needed after stretching -- will alter the colorOrder!
            rgb_ir_model.remove_ir_contamination(img);
            reconstruct_bayer_cfa(img);
            sve_demosaic_nearest_neighbor(img);

            // Multiply by self
            calculate_and_apply_lsc_gain_table(img, img);

            // Clip to range / reduce 4x2 DP sensor while clipping to range
            preprocess_fn0(img);

            ff_img = img;   // reduce float -> int
        }

        //     fn3b(image, meta, ffImage, ffMeta, ...
        //          colorOrder, blackLevel, saturationLevel, irCalibration, bayerBlackLevelMultiplier)
        //     stretches the dynamic range of both images (after black level reduceion)
        //     multiplies "image" by the inverse of ff_Image
        //     returns regular 2x2 sensor pattern by removing IR contamination and sums up DP pixels
        template <typename T>
        void preprocess_fn3b(
            bayer_image_s<T> &ccc_img,
            bayer_image_s<T> &ff_image,
            float ag = 1.0f, float exp = 300.0f)
        {
            // We need to calculate the gain tables etc at least as floats -- take a copy
            bayer_image_s<float> img = ccc_img;
            bayer_image_s<float> img_pair = ff_image;
            auto &layout = ccc_img._layout;

            int multiplier = layout.is_dp_4x2_sensor() ? 2 : 1;
            auto bl = get_black_level(layout, ag, exp);
            stretch_to_max_level(img, bl, saturation, multiplier);
            stretch_to_max_level(img_pair, bl, saturation, multiplier);

            // IR contamination removal needed after stretching -- will alter the colorOrder!
            rgb_ir_model.remove_ir_contamination(img);
            reconstruct_bayer_cfa(img);
            sve_demosaic_nearest_neighbor(img);

            rgb_ir_model.remove_ir_contamination(img_pair);
            reconstruct_bayer_cfa(img_pair);
            sve_demosaic_nearest_neighbor(img_pair);

            calculate_and_apply_lsc_gain_table(img, img_pair);

            // Sums Left and Right subchannels of DP 4x2 sensor
            preprocess_fn0(img);

            ccc_img = img;
        }
        bayer_pattern_e reduce_dp_layout(bayer_pattern_e pattern)
        {
            switch (pattern)
            {
            case bayer_pattern_e::bggr_4x2: return bayer_pattern_e::bggr;
            case bayer_pattern_e::rggb_4x2: return bayer_pattern_e::rggb;
            case bayer_pattern_e::gbrg_4x2: return bayer_pattern_e::gbrg;
            case bayer_pattern_e::grbg_4x2: return bayer_pattern_e::grbg;
            default:
                throw std::runtime_error("Can't reduce non dp image");
            }
        }

        // IR sensor reduceion to 2x2 symbolic interpretation:
        // all 2x2 sensors are of form
        //   RG  -> RG   i.e. I is changed to G
        //   IB     GB
        // all 4x4 sensors (or their top-left 2x2 subsets) are of form
        //   BG     BG   i.e. I is changed to either B, or R
        //   GI ->  GR        which ever has only zero occurrences
        static bayer_pattern_e reduce_ir_layout(bayer_pattern_e pattern)
        {
            // Replace the "I" by "G" in the enumeration
            switch (pattern)
            {
            case bayer_pattern_e::bgir:
            case bayer_pattern_e::bigr:
                return bayer_pattern_e::bggr;

            case bayer_pattern_e::rgib:
            case bayer_pattern_e::rigb:
                return bayer_pattern_e::rggb;

            case bayer_pattern_e::gbri:
            case bayer_pattern_e::ibrg:
                return bayer_pattern_e::gbrg;

            case bayer_pattern_e::grbi:
            case bayer_pattern_e::irbg:
                return bayer_pattern_e::grbg;

            case bayer_pattern_e::bgrg_gigi_rgbg_gigi: return bayer_pattern_e::bggr;
            case bayer_pattern_e::grgb_igig_gbgr_igig: return bayer_pattern_e::grbg;
            case bayer_pattern_e::rgbg_gigi_bgrg_gigi: return bayer_pattern_e::rggb;
            case bayer_pattern_e::gbgr_igig_grgb_igig: return bayer_pattern_e::gbrg;
            case bayer_pattern_e::gigi_rgbg_gigi_bgrg: return bayer_pattern_e::gbrg;
            case bayer_pattern_e::igig_gbgr_igig_grgb: return bayer_pattern_e::rggb;
            case bayer_pattern_e::gigi_bgrg_gigi_rgbg: return bayer_pattern_e::grbg;
            case bayer_pattern_e::igig_grgb_igig_gbgr: return bayer_pattern_e::bggr;

            default:
                throw std::runtime_error("Unknown IR sensor");
            }
        }

        /// inplace copy nearest G channel over I channel
        /// reduces the bayer sensor layout to one of four basic layouts
        template <typename T>
        void reconstruct_bayer_cfa(bayer_image_s<T> &bayer_image)
        {
            auto &layout = bayer_image._layout;
            if (!layout.is_ir_sensor())
                return;

            if (layout.is_2x2_ir_sensor())
            {
                // Reconstruct the I channel from G channel
                bayer_image[bayer_info_s::color_info_e::ir].init_from(bayer_image[bayer_info_s::color_info_e::green]);
            }
            else if (layout.is_4x4_ir_sensor())
            {
                //   The first quad determines the destination order of all 4 sub-quads:
                //   R G -> R G | b G   'b' and 'i' are in wrong order --> 'i'
                //   G I    G B | G i   is copied from the non-green color in the first quad
                //          b G   r G    in the third and fourth quad the 'b' is copied to 'i'
                //          G i   G i    and then the 'r' is replaced over the 'b'
                //

                // Quadrants 1 and 2 -- locates the 'I' from quadrant 1
                auto ir_index = layout.locate_color(bayer_info_s::color_info_e::ir);  // 0,1, 4 or 5
                auto color_idx = ir_index ^ 5;  // locates the color opposite of 'I' by toggling 0<->5 and 1<->4

                // copy the color from quadrant 2 over 'I' (in both quadrants)
                bayer_image[ir_index].init_from(bayer_image[color_idx + 2]);
                bayer_image[ir_index + 2].init_from(bayer_image[color_idx + 2]);
                // -- then copy the color form quadrant 1 over color in quad 2
                bayer_image[color_idx + 2].init_from(bayer_image[color_idx]);

                // Copy the color from quadrant 3 over 'I' in quadrants 3 and 4
                bayer_image[ir_index + 8].init_from(bayer_image[color_idx + 8]);
                bayer_image[ir_index + 10].init_from(bayer_image[color_idx + 8]);
                // -- and replace the quadrant 3 color by quadrant 4 color
                bayer_image[color_idx + 8].init_from(bayer_image[color_idx + 10]);
            }
            layout = bayer_info_s(reduce_ir_layout(layout));
        }

        /// Given two images, the first is multiplied by the reciprocal of the second (channel-wise)
        /// Note that both images can point to the same object
        /// - A copy is made in the inner loop where "dst_chan" and "chan" point to different objects
        /// \TODO - if the reference image happens to be different size, it's scaled
        template<typename T>
        void calculate_and_apply_lsc_gain_table(bayer_image_s<T> &img, bayer_image_s<T> &img_ff)
        {
            if (!img.matches_by_size_and_type(img_ff))
                throw std::runtime_error("TODO: LSC Gain must have images of equal dimensions");

            /// idiomatic way to loop over all channels in image -- i = 0..N-1
            for (auto &i : img)
            {
                auto chan = img_ff[i];
                auto dst_chan = img[i];

                auto medi = median7x7(chan, SYMMETRIC);
                gaussian35x35(medi, REPLICATE);

                auto maxitem = (float)(T)medi.foreach(maximum_f<T>{});

                dst_chan.foreach([maxitem](T &dst, T&src)
                {
                    float gain = maxitem / (float)src;
                    if (gain > 10.0f)
                        gain = 16.0f;
                    dst = reduce_to<T>(gain * dst);
                }, medi);
            }
        }
    };


    inline image<uint16_t> median_subsample_7x7(image<uint16_t> &big)
    {
        auto block_size = roi_point(7,7);
        auto output = image<uint16_t>(big.size() / block_size);
        int skipx = big._skip_x;
        int skipy = big._skip_y;
        output.foreach([skipx, skipy](uint16_t &dst, uint16_t &src)
        {
            uint16_t arr[49];
            uint16_t *s = &src;
            uint16_t *a = &arr[0];
            for (int y = 0; y < 7; y++, s += skipy)
            {
                for (int x = 0; x < 7; x++)
                    *a++ = s[x * skipx];
            }
            std::sort(arr, arr + 49);
            dst = arr[24];
        }, big.subview(block_size));
        return output;
    }

    /// Calculates the convex hull of R/G, B/G chromaticity coordinates
    /// from a (preprocessed) image
    /// The gamut points are quantized to 8q8 format -- and as such they
    /// are domain specific to Chromaticity Response calculation
    /// A better API would externalize the r,g,b -> r/g b/g conversion
    inline std::vector<point_xy> calculate_image_gamut(bayer_image_s<uint16_t> &img)
    {
        auto red_idx = img._layout.locate_color(color_info_e::red);
        auto blue_idx = img._layout.locate_color(color_info_e::blue);
        auto green_idx = img._layout.locate_color(color_info_e::green);

        // Step 1 - filter / subsample image in 7x7 blocks
        // Get the first red, green and blue channel, then average the 2 green channels
        auto red = img[red_idx];
        auto blue = img[blue_idx];
        auto green = img[green_idx];
        // when on 2x2 sensors, the opposite index == 3 - index (or index ^ 3)
        // when the opposite index of first green != green
        // this function can calculate the gamut anyway (e.g. for RGBIR 2x2, or rgbir 4x4)
        // using just 3 channels out of 16.
        if (img._layout[3 - green_idx] == color_info_e::green)
        {
            green.foreach([](uint16_t &green1, uint16_t &green2)
            {
                green1 = static_cast<uint16_t>((green1 + green2 + 1) >> 1);
            }, img[3 - green_idx]);
        }

        // ignore elements smaller than this intensity
        // previously this was fixed at 51.2f (with 10 bit sensors)
        auto min_intensity = std::max(51, green(green.max_element()) / 20);

        // reduce convex hull complexity by sub sampling by 7x7 (taking median of each support)
        auto red_small = median_subsample_7x7(red);
        auto green_small = median_subsample_7x7(green);
        auto blue_small = median_subsample_7x7(blue);

        // Step 2 - convert reduced sized planes to r/g, b/g
        struct fxy
        {
            float x;
            float y;
            bool operator ==(const fxy &other) { return x == other.x && y == other.y; }
            fxy() : x(0.0f), y(0.0f) { }
        };
        auto size_small = red_small.size();

        std::vector<fxy> xy_points(size_small._x * size_small._y);  // allocate maximum amount of points
        auto *ptr = &xy_points[0];                                  // propagate raw pointer to the lambda

        red_small.foreach([&ptr, min_intensity](uint16_t &red, uint16_t &green, uint16_t &blue)
        {
            if (green > min_intensity && (red > min_intensity || blue > min_intensity))
            {
                auto g_round = green / 2;
                // calculate r/g, b/g in steps of 1/256
                ptr->x = (1.0f / 256.0f) * ((red * 256 + g_round) / green);
                ptr->y = (1.0f / 256.0f) * ((blue * 256 + g_round) / green);
                ptr++;
            }
        }, green_small, blue_small);

        xy_points.resize(ptr - xy_points.data());       // remove excess allocation

        // Step 3 - calculate convex hull from the points
        convex_hull::quick_hull(xy_points);

        // Step 4 - prune points that don't contribute much
        //        - keep minimum of 8 points, max 16 points
        //        - and features of size 1e-4 (to original polygon)
        convex_hull::prune_hull(xy_points, 8, 16, 1e-4);

        auto result = std::vector<point_xy>(xy_points.size());
        for (size_t i = 0; i < xy_points.size(); i++)
        {
            result[i] = point_xy(xy_points[i].x, xy_points[i].y);
        }

        return result;
    }
};
