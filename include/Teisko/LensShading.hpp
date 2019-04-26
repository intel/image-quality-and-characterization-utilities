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
#include "Teisko/Image/Algorithms.hpp"
#include "Teisko/Preprocessing.hpp"
#include <cstdint>
#include <vector>

namespace Teisko
{
    template <typename T>
    struct lensshading_grid
    {
        int width;                      // output grid width
        int height;                     // output grid height
        int channels;                   // number of unique channels
        int scale;                      // scaling factor as 2^(scale)
        int pair_index;                 // Best matching pair for ratio tables, 0 otherwise
        std::vector<T> grid;            // Gain table -- Channels * width * height elements
        std::vector<char> ch_indices;   // 16-entry vector for channel indices
        chromaticity chroma;            // calculated chromaticity

        // Returns image<T>:: access class for Nth lens shading grid
        image<T> operator[] (size_t channel)
        {
            return image<T>(roi_point{ width, height }, grid.data() + channel * (height * width));
        }

        // returns true if two grids match in dimensions
        // - this is prerequisite e.g. for combining grids
        bool have_same_dimensions(const lensshading_grid &other)
        {
            return
                width == other.width && height == other.height &&
                channels == other.channels && scale == other.scale &&
                grid.size() == other.grid.size();
        }

        // Accumulates two grids for averaging
        lensshading_grid& operator += (const lensshading_grid &other)
        {
            if (have_same_dimensions(other))
            {
                chroma += other.chroma;
                for (size_t i = 0; i < grid.size(); i++)
                    grid[i] += other.grid[i];
            }
            return *this;
        }

        // scales grid and chromaticity after averaging
        lensshading_grid& operator /= (int count)
        {
            if (count > 1)
            {
                chroma /= count;
                for (auto &x : grid)
                    x /= count;
            }
            return *this;
        }

        // Divides (static) grid by (NVM) grid of same size
        // - sets the pair index according to the divider
        lensshading_grid& operator /= (const lensshading_grid &other)
        {
            // Validates that the dimensions of this and the other match
            if (have_same_dimensions(other))
            {
                for (size_t i = 0; i < grid.size(); i++)
                    grid[i] /= other.grid[i];
                pair_index = other.pair_index;
            }
            return *this;
        }

        // Creates a lens shading grid from given dimensions
        lensshading_grid(int w = 0, int h = 0, int chans = 0, int fractional_bits = 0)
            : width(w), height(h), channels(chans), scale(fractional_bits), pair_index(0)
            , grid(w * h * chans), ch_indices(0)
        { }

        // Up/Down -scales the grid between fixed-point and floating point representations
        template <typename U>
        lensshading_grid(const lensshading_grid<U> &other, int fractional_bits = 0)
            : width(other.width), height(other.height), channels(other.channels), scale(fractional_bits)
            , pair_index(other.pair_index), grid(other.grid.size()), ch_indices(other.ch_indices)
            , chroma(other.chroma)
        {
            // factor = std::pow(2.0, scale - other.scale);
            auto factor = scale >= other.scale
                ? 1.0 * (1 << (scale - other.scale))
                : 1.0 / (1 << (other.scale - scale));

            for (size_t i = 0; i < grid.size(); i++)
                grid[i] = reduce_to<T>((double)other.grid[i] * factor);
        }

        // Spatial up/down scaling the grid using polynomial fitting
        // Used e.g. in adding LSC NVM data to model
        template <typename U>
        lensshading_grid(const lensshading_grid<U> &input, int dst_width, int dst_height, int fractional_bits = 0)
            : width(dst_width), height(dst_height), channels(input.channels), scale(fractional_bits)
            , pair_index(input.pair_index), grid(dst_width * dst_height * channels)
            , ch_indices(input.ch_indices), chroma(input.chroma)
        {
            poly_4_scaler_2d<T> polynomial_scaler;
            auto tmp_copy = lensshading_grid<T>(input, fractional_bits);
            for (int i = 0; i < channels; i++)
            {
                auto input_channel = tmp_copy[i];
                auto output_channel = (*this)[i];
                polynomial_scaler.scale(input_channel, output_channel);
            }
        }
    };

    using lsc_nvm = std::vector<std::vector<lensshading_grid<double>>>;

    /// Lens Shading Calculator class produces e.g. 63x47 sized gain table from a flat field image.
    /// Multiplying the original image by the stretched gain table should produce an approximately
    /// flat image -- i.e. with constant intensity in all channels.
    class lensshading_calculator
    {
    public:
        // Initialize a model to be used for all succeeding image characterizations
        //  - grid_indices can be empty, in which case all channels in the image are characterized separately
        //  - otherwise the grid index vector can contain numbers from 0..N without gaps
        //    where all subchannels containing same number are averaged
        //       e.g. [0 0 1 1; 2 2 3 3; -1, -1, ... ] can be used to average L+R pixels on DP 4x2 sensor
        //  - don't care indices are marked by value -1
        //     - when the image sensor type is WxH, only the top left WxH portion of the indices is considered
        lensshading_calculator(int grid_width = 63, int grid_height = 47, std::vector<char> grid_indices = {})
            : grid_dims(grid_width, grid_height), indices(grid_indices)
        {
            if (grid_indices.size() != 0 && grid_indices.size() != 16)
                throw std::runtime_error("Expecting zero length or 16 items in grid indices");
        }

        float get_aspect_ratio()
        {
            return channel_dims._y == 0 ? 0.0f : (float)channel_dims._x / (float)channel_dims._y;
        }

        // Calculates the shading grid and chromaticity for a preprocessed image
        // - image should have black level and rgb-ir contamination removed
        template <typename T>
        lensshading_grid<double> calculate_grid(bayer_image_s<T> &image)
        {
            const double grid_mean_threshold = 0.5;    // Rejects 25% + 25% of outliers
            model_setup(get_channel_dim(image), image._layout);

            auto result = lensshading_grid<double>(grid_dims._x, grid_dims._y, (int)weights.size());
            auto temp_grid = Teisko::image<double>(grid_dims);

            chromaticity_factory_f chroma_stats;
            for (auto i : image)
            {
                auto tgt_index = grid_indices_reshaped[i];
                if (tgt_index < 0 || static_cast<size_t>(tgt_index) > weights.size())
                    continue;
                auto weight = weights[tgt_index];

                // First average the channel with a box filter of K*K
                // intermediate image pixel type = uint32_t
                auto averaged = image[i].template convert_to<uint32_t>();
                box_filter(averaged, kernel_size, REPLICATE);

                // Calculate the trimmed mean around each grid point
                for (int y = 0; y < grid_dims._y; y++)
                {
                    roi_point offset(0, y_positions[y]);
                    for (int x = 0; x < grid_dims._x; x++)
                    {
                        offset._x = x_positions[x];
                        temp_grid(y, x) = mean_calculator(
                            averaged.region(roi_point(kernel_size), offset), grid_mean_threshold);
                    }
                }

                // Locates the centroid of all items in grid confirming to the highest value
                // Invert/normalize the grid to the maximum value
                // then accumulate it to result -- for regular sensors weights are all ones
                // - for 4x4 RGB-IR sensors or some other cases we can average e.g. I channels directly

                auto max_position = temp_grid.max_element();
                auto max_value = temp_grid(max_position);
                result[tgt_index].foreach([max_value, weight](double &dst, double &val)
                {
                    // Following border value (16) was selected experimentally
                    // to make sure output of WFOV image is decent
                    const double max_gain = 16.0;
                    const double unity = 1.0;
                    double gain;
                    if (val <= 0.0 || max_value <= 0.0)
                    {
                        // 0/0 will produce unity gain, non-zero / 0 produces "inf"
                        gain = max_value <= 0.0 ? unity : max_gain;
                    }
                    else
                    {
                        gain = max_value >= max_gain * val ? max_gain : max_value / val;
                    }
                    dst += gain * weight;
                }, temp_grid);

                // Locates the trimmed mean around the fine tuned maximum sensitivity point
                max_position._x = x_positions[max_position._x];
                max_position._y = y_positions[max_position._y];
                max_value = find_max_sensitivity(max_position + kernel_size / 2, averaged);
                chroma_stats[image._layout[i]] += max_value;
            }

            result.ch_indices = indices;
            result.chroma = chroma_stats;

            return result;
        }

    private:
        /// The operation model for this class is
        /// 1) create the model
        ///   -- initializing grid_dimensions and the grid index vector
        roi_point grid_dims;            // width and height of lens shading grid
        std::vector<char> indices;      // indices to be stored to each lsc grid

        /// 2) Process images one by one
        ///   -- the first image will initialize the context
        ///      setting the parameters derived from the image dimensions

        static const int minimum_awb_kernel = 10;
        std::vector<char> grid_indices_reshaped;    // grid indices linearized to 2x2, 4x2 or 4x4 format
        std::vector<double> weights;                // Weight for each unique channel
        int kernel_size = 0;
        bayer_pattern_e pattern = bayer_pattern_e::bggr;
        roi_point channel_dims;                     // width and height of a single channel in image
        roi_point max_search_area;                  // width and height of awb sensitivity search window
        std::vector<int> x_positions;               // centers of kernel positions
        std::vector<int> y_positions;               // same here -- derived from first input image

        trimmed_mean_f<uint32_t> mean_calculator;     // can calculate trimmed mean of vectors of uint32_t

        // Uniformly samples the vector (first:last) rounding to nearest integer
        // constraints: all numbers to be positive or zero
        // may overflow -- this should never be an issue for images of reasonable size
        template <typename T>
        std::vector<T> lin_space(int elements, T first, T last)
        {
            std::vector<T> grid(elements--, first);
            for (int i = 1; i <= elements; i++)
                grid[i] = static_cast<T>(first + ((last - first) * i + elements / 2) / elements);
            return grid;
        }

        // returns element wise smallest coordinate
        roi_point min(roi_point a, roi_point b)
        {
            return roi_point(std::min(a._x, b._x), std::min(a._y, b._y));
        }

        // returns element wise largest coordinate
        roi_point max(roi_point a, roi_point b)
        {
            return roi_point(std::max(a._x, b._x), std::max(a._y, b._y));
        }

        // Record/cache the first valid image dimension and layout
        // and setup the internal values with the dimension
        //  - throw errors if there is a mismatch or illegal argument
        // Called from calculate_grid()
        void model_setup(roi_point dims, bayer_info_s &layout)
        {
            if (dims._x <= 0 || dims._y <= 0)
                throw std::runtime_error("Image dimensions must be positive");

            if (channel_dims == dims && layout == pattern)
                return;

            if (!(channel_dims == roi_point()))
                throw std::runtime_error("Image dimensions and type must match between iterations");

            channel_dims = dims;
            pattern = layout;

            // median filter AND box filter size -- kernel size must be ODD for symmetry
            // The channel width is divided by experimental value of 64, which is close
            // to the actual fixed grid width 63. This value is probably good for channel widths
            // close to 4k, meaning that near by final grid values are not completely independent
            const int experimental_divisor = 64;
            kernel_size = (dims._x / experimental_divisor) | 1;

            x_positions = lin_space<int>(grid_dims._x, 0, dims._x - kernel_size);
            y_positions = lin_space<int>(grid_dims._y, 0, dims._y - kernel_size);

            // calculate search area that spans at least to next grid center (but still at least 10 pixels)
            max_search_area = max(roi_point(minimum_awb_kernel, minimum_awb_kernel), (dims + (grid_dims / 2)) / grid_dims);

            // Make a grid index map corresponding to number of color channels in bayer pattern
            if (indices.size() != 16)
            {
                for (uint32_t i = 0; i < layout.get_channels(); i++)
                {
                    grid_indices_reshaped.push_back((char)i);
                    weights.push_back(1.0);
                }
            }
            else
            {
                grid_indices_reshaped = remap_4x4_vector(indices, layout);
                weights = get_weights_from_indices(indices);
            }
        }

        // Returns the maximum sensitivity from the `channel` given a coarse estimation of search center
        template <typename T>
        double find_max_sensitivity(roi_point max_pos, image<T> &channel)
        {
            const double sensitivity_threshold = 0.4;  // Rejects 20% + 20% of outliers
            // search window is limited by max area and the distance of the position to image boundary
            auto search_window = min(max_search_area, min(max_pos, channel_dims - 1 - max_pos));
            auto avg_kernel_size = std::min({ minimum_awb_kernel, search_window._x, search_window._y });

            // box filter around the search area (producing a fresh copy)
            // locate the maximum and remap that once again to `channel`
            max_pos = max_pos - search_window;
            auto averaged = channel.region(search_window * 2 + 1, max_pos).template convert_to<double>();
            box_filter(averaged, avg_kernel_size * 2 + 1, REPLICATE);

            max_pos = max_pos + averaged.max_element(); // position according to filtering

            // Find out if we can fit 41x41 trimmed mean averaging kernel around this new position
            auto distance_to_border = min(max_pos, channel_dims - 1 - max_pos);
            auto half_kernel_size = std::min({ 20, distance_to_border._x, distance_to_border._y });

            auto mean = mean_calculator(channel
                .region(roi_point(half_kernel_size * 2 + 1), max_pos - half_kernel_size),
                sensitivity_threshold);
            return mean;
        }

        // Returns the vector of weights for channels 0..N, N < 16
        //  - typical case for 2x2 vector is to have a weight of 1.0 for all channels
        //  - in 4x4 RGB-IR sensors some channels are averaged in order to produce
        //    less than 16 channels for characterization
        // Returns empty vector in case of malformed input
        inline std::vector<double> get_weights_from_indices(std::vector<char> &vec)
        {
            if (vec.size() != 16)
                return{};

            // first make a histogram of entries
            auto result = std::vector<double>(16);
            for (auto &i : vec)
            {
                if (i > 15)
                    return{};
                if (i >= 0)
                    result[i]++;
            }
            // remove trailing zeros from the histogram
            while (!result.empty())
            {
                if (result.back() == 0.0)
                    result.pop_back();
                else
                    break;
            }
            // invert the histogram (checking that all elements are non-zero)
            for (auto &x : result)
            {
                if (x == 0.0)
                    return{};
                x = 1.0 / x;
            }
            return result;
        }
    };

    /// Domain model for LSC calculation
    ///   - Incremental calculation of images and storage of grids
    ///     preprocesses images for
    ///       - black level
    ///       - todo: linearization
    ///       - SVE pattern removal
    ///       - rgb ir (de)contamination
    ///   - Get Static Tables integrates the current state so far
    ///   - Get Ratio Tables integrates the current state so far
    class lsc_model
    {
    public:
        lsc_model(int width, int height, std::vector<char> indices)
            : calculator(width, height, indices) {}

        float get_aspect_ratio() {
            return calculator.get_aspect_ratio();
        }

        std::vector<std::vector<lensshading_grid<double>>> nvm_grids;

        // Adds lsc grid to internal state -- primarily for testing purposes
        void add_item(lensshading_grid<double> grid, int msid, int light_source)
        {
            state[light_source][msid].push_back(grid);
        }

        void add_item(bayer_image_s<uint16_t> image, int msid, int light_source)
        {
            auto grid = calculator.calculate_grid(image);
            add_item(grid, light_source, msid);
        }

        std::map<int, lensshading_grid<uint16_t>> get_ratio_tables(int exponent = 5)
        {
            return get_tables(exponent, true);
        }

        std::map<int, lensshading_grid<uint16_t>> get_static_tables(int exponent = 5)
        {
            return get_tables(exponent, false);
        }

    private:
        lensshading_calculator calculator;

        // state hierarchy
        //  state[1]["abcd"] =  { grid0, grid1, grid2 };
        //  state[1]["cdef"] =  { grid3, grid4 };
        //  state[6]["abcd"] =  { grid5 ];
        //  with two light sources {1, 6} and two msids { "abcd" and "cdef" }
        //  - when creating the final static grid for light source `1`
        //    we first average horizontally for grids012 and grids34, then vertically
        //    these two grids with equal weights
        std::map<int, std::map<int, std::vector<lensshading_grid<double>>>> state;

        // Given a vector of grids, calculate the average
        // - applies to chromaticity and gains
        lensshading_grid<double> average_of_grids(std::vector<lensshading_grid<double>> &vec)
        {
            if (vec.size() == 0)
                return{};

            auto result = vec[0];
            for (size_t i = 1; i < vec.size(); i++)
                result += vec[i];

            if (vec.size() > 1)
                result /= static_cast<int>(vec.size());

            return result;
        }

        // Divides static table by ratio table determined by the msid and light_source
        // returns false in case of error, in which case the `tables` is left in
        // unspecified state and must be rejected by caller
        bool divide_by_light_source(
            std::vector<lensshading_grid<double>> &tables,
            std::vector<int> &msids,
            size_t light_source)
        {
            size_t n = tables.size();
            if (n != msids.size() || nvm_grids.size() == 0)
                return false;

            for (size_t i = 0; i < n; i++)
            {
                int msid = msids[i];
                if (msid < 0 || static_cast<size_t>(msid) >= nvm_grids.size())
                    return false;
                tables[i] /= nvm_grids[msid][light_source];
            }
            return true;
        }

        //  divides `tables` by (best) nvm light source
        //  returns false in case of error leaving `tables` to unspecified state
        //  - implements a voting scheme for smallest total difference to unity
        //    in case that nvm data contains grids for multiple light sources
        bool convert_to_ratio_tables(
            std::vector<lensshading_grid<double>> &tables,
            std::vector<int> &msids)
        {
            if (nvm_grids.size() == 0)
                return false;

            size_t n_light_sources = nvm_grids[0].size();
            if (n_light_sources == 0)
                return false;

            // ensure that all msidX.nvms have same amount of light sources
            for (auto &nvm_grid : nvm_grids)
            {
                if (nvm_grid.size() != n_light_sources)
                    return false;
            }

            if (n_light_sources == 1)
                return divide_by_light_source(tables, msids, 0);

            // make N copies  -- divide each copy by the Nth light source
            std::vector<std::vector<lensshading_grid<double>>> copies(n_light_sources, tables);
            for (size_t i = 0; i < n_light_sources; i++)
            {
                if (divide_by_light_source(copies[i], msids, i) == false)
                    return false;
            }

            // Have each MSID vote for best light source
            auto votes = std::vector<int>(n_light_sources);
            for (size_t j = 0; j < tables.size(); j++)
            {
                std::vector<double> sad(n_light_sources);
                for (size_t i = 0; i < n_light_sources; i++)
                    sad[i] = calculate_sad_to_unity(copies[i][j]);
                auto it = std::min_element(std::begin(sad), std::end(sad));
                votes[std::distance(std::begin(sad), it)]++;
            }

            // find maximum vote count -- and [move] assign the best light source as result
            auto it = std::max_element(std::begin(votes), std::end(votes));
            auto winner = std::distance(std::begin(votes), it);
            tables = std::move(copies[winner]);     // actual ownership transfer
            return true;
        }

        /// Returns sum of absolute difference of gains to 1.0 for
        /// estimation of best match to NVM light source
        template <typename T>
        T calculate_sad_to_unity(lensshading_grid<T> &grid)
        {
            T sum = 0.0;
            for (auto &item : grid.grid)
                sum += std::abs(item - (T)1.0);
            return sum;
        }

        /// Produces final averaged static and ratio tables
        /// Both cases need to average the repetitions (within each lightsource and msid)
        /// and then by giving equal weight to each msid (within each lightsource)
        /// Produces ratio tables (or empty table) when is_ratio is true and nvm contains proper amount of nvm lsc grids
        /// Produces static tables when `is_ratio` is false
        std::map<int, lensshading_grid<uint16_t>> get_tables(int exponent, bool is_ratio)
        {
            std::map<int, lensshading_grid<uint16_t>> tables;

            for (auto light_source_map : state)
            {
                // Average all repetitions for each msid
                std::vector<lensshading_grid<double>> grids;
                std::vector<int> ids;
                for (auto msid : light_source_map.second)
                {
                    ids.push_back(msid.first);
                    grids.push_back(average_of_grids(msid.second));
                }

                // Divide the static tables by the nvm model -- don't produce
                // a ratio table if it fails -- for static tables don't try to convert
                if (is_ratio && convert_to_ratio_tables(grids, ids) == false)
                    continue;

                // Both ratio and static tables are averaged giving equal weight to each contributing MSID
                auto avg_final_grid = average_of_grids(grids);
                tables[light_source_map.first] =
                    lensshading_grid<uint16_t>(avg_final_grid, 16 - exponent);
            }
            return tables;
        }
    };
}