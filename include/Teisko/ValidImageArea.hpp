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
#include "Teisko/BayerImage.hpp"
#include "Teisko/Preprocessing.hpp"
#include <cstdint>
#include <vector>
#include <map>
#include <array>
#include <string>
#include <algorithm>  // std::find

namespace Teisko
{
    /// ValidImageArea removes pixels from an initially "white" image of size W*H
    /// according to new images added to the model
    struct ValidImageArea
    {
    public:
        roi_point grid_size;
        std::vector<uint16_t> data;

        Teisko::image<uint16_t> get_grid() { return image<uint16_t>(grid_size, data.data()); }

        ValidImageArea(int width = 60, int height = 60)
            : grid_size{ width, height }, data(width * height, 1) { }

        void Add(Teisko::bayer_image_s<uint16_t> &image) {

            auto green_channel = demosaic_bilinear(image, color_info_e::green, REPLICATE_EVEN).convert_to<double>();

            auto resized = green_channel.convert_to();                  // make a copy
            // decimate by 2 until the image is small enough
            while (resized._width > 200 || resized._height > 200)
                resized = resize_image(resized, point_xy(0.5, 0.5), REPLICATE);

            // Fit the resized model to a polynomial, and evaluate the quality of it
            auto scaler = poly_4_scaler_2d<double>{};
            auto poly_fit = Teisko::image<double>(resized.size());
            scaler.scale(resized, poly_fit);

            auto fit_negative_values = exclude_negative(resized, poly_fit);

            if (fit_negative_values)
                scaler.scale(resized, poly_fit);

            if (is_bad_global_model(resized, poly_fit))
                return;

            // apply the model to original image sized image, then evaluate the goodness block per block
            auto orig_sized_model = Teisko::image<double>(green_channel.size());
            scaler.interp(orig_sized_model);

            auto grid = get_grid();
            auto kernel_size = green_channel.size() / grid_size;
            auto last_coordinates = green_channel.size() - kernel_size;
            auto grid_pos_y = calculate_grid_positions(grid._height, 0, last_coordinates._y);
            auto grid_pos_x = calculate_grid_positions(grid._width, 0, last_coordinates._x);

            for (auto &y : grid._height)
            {
                for (auto &x : grid._width)
                {
                    auto offset = roi_point(grid_pos_x[x], grid_pos_y[y]);
                    auto block_input = green_channel.region(kernel_size, offset);
                    auto block_output = orig_sized_model.region(kernel_size, offset);
                    if (is_bad_block(block_input, block_output))
                        grid(y, x) = 0;
                }
            }
        }

    private:
        // Given input image and same sized output image from polyfitting
        // mark those input pixels as Nan, which have a corresponding negative output in the model
        bool exclude_negative(Teisko::image<double> &input, Teisko::image<double> &model)
        {
            auto has_corrected_values = false;
            input.foreach([&has_corrected_values](double &i, double &o)
            {
                if (o < 0.0)
                {
                    has_corrected_values = true;
                    i = std::numeric_limits<double>::quiet_NaN();
                }
            }, model);
            return has_corrected_values;
        }

        // evaluates how the polynomial model fits the global data
        bool is_bad_global_model(Teisko::image<double> &input, Teisko::image<double> &model)
        {
            return is_bad_model(input, model, true);
        }

        // evaluates how the polynomial model fits a local block
        bool is_bad_block(Teisko::image<double> &input, Teisko::image<double> &model)
        {
            return is_bad_model(model, input, false);
        }

        // calculates a correspondence metric of input matching a model
        bool is_bad_model(Teisko::image<double> &input, Teisko::image<double> &model, bool calculate_mean)
        {
            double sum = 0.0;
            if (calculate_mean)
                model.foreach([&sum](double &x) { sum += x; });

            struct {
                double mean;
                double ss_res; // Sum of squares of residuals
                double ss_tot; // Total sum of squares
                double pow2(double a) { return a*a; }
                double operator()() { return ss_tot == 0.0 ? 0.0 : 1 - ss_res / ss_tot; }
                void operator()(double &x, double &y)
                {
                    auto d0 = pow2(x - y);
                    auto d1 = pow2(x - mean);
                    if (!std::isnan(d0))
                    {
                        ss_res += d0;
                        ss_tot += d1;
                    }
                }
            } functor_r_square{ sum / (model._width * model._height), 0, 0 };

            return model.foreach(functor_r_square, input)() < 0.95;
        }

        // Divides evenly items from 'first' to 'last' (both inclusive) into a grid
        std::vector<int> calculate_grid_positions(int elements, int first, int last)
        {
            std::vector<int> grid;
            grid.reserve(elements);
            double temp = (elements <= 1) ? 0.0 : (double)(last - first) / (double)(elements - 1);

            for (int i = 0; i < elements; i++)
                grid.emplace_back((int)(round(temp * i) + first));
            return grid;
        }
    };
}