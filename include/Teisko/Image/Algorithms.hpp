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
#include "Teisko/Algorithm/PointXY.hpp"
#include "Teisko/Algorithm/Bit.hpp"
#include "Teisko/Algorithm/ReduceTo.hpp"
#include "Teisko/Image/API.hpp"

#include <cstdint>
#include <vector>
#include <limits>               // numeric_limits<T>::min() used in int->uint conversion

namespace Teisko
{
    // Image Processing Algorithms
    // Median filter an image of type T with rectangular support 7x7
    // border processing schemes support zeroing, replication and mirroring

    template <typename T> inline
    image<T> median7x7(image<T> &src, int mirror_scheme = ZERO_PAD);

    template <typename T> inline
    void gaussian35x35(image<T> &src, int rep_scheme = ZERO_PAD);

    // Copies a block of 8 rows, 2 columns from `origin`
    // deinterleaving them to bitplanes
    template <typename T> inline
    void update_masks(uint64_t(&bits)[sizeof(T) * 8], T *origin, int stride);

    // Runs N iterations for pixel data as bit planes in `bits`
    // to recover the K(th) largest item (as set by initial threshold)
    // The parameter `mask` must be initialized to contain a set bit for all those bits
    // of interest in the corresponding pixel data bits[0..N-1]
    template <int N> inline
    uint64_t median8x8_iteration(uint64_t(&bits)[N], uint64_t mask, uint64_t threshold)
    {
        uint64_t result = 0;
        int i = 0;
        do
        {
            uint64_t ones = mask & bits[i];
            uint64_t ones_size = popcount(ones);
            uint64_t mask_size = popcount(mask);
            auto zero_size = mask_size - ones_size;
            int new_bit = 0;
            if (zero_size < threshold)
            {
                new_bit = 1;
                threshold -= zero_size;
                mask = 0;
            }
            result = result * 2 + new_bit;
            mask ^= ones;
        } while (++i < N);
        return result;
    }

    /// Floats are of format [Sgn][exp][mantissa]
    /// The most negative value (i.e. smallest value) is 0xffffffff; (when NaN's are omitted)
    /// the most positive value is 0x7fffffff; again not counting NaNs
    /// We'll convert these to 32-bit integers (preserving the order) by adding top bit to positive values
    /// and ones complementing all negative values -- and reverse the process afterwards
    template <> inline
    image<float> median7x7(image<float> &src, int mirror_scheme)
    {
        auto copy = src.convert_to();
        uint32_t *begin = reinterpret_cast<uint32_t*>(copy._begin);
        auto size = src._width * src._height;
        for (decltype(size) i = 0; i < size; ++i)
        {
            uint32_t toggle_mask = begin[i] & 0x80000000 ? 0xFFFFFFFF : 0x80000000;
            begin[i] ^= toggle_mask;
        }
        auto input_as_uint32 = image<uint32_t>(copy._height, copy._width, begin);

        auto result_as_uint32 = median7x7(input_as_uint32, mirror_scheme);

        auto result_as_float = image<float>(copy._height, copy._width);
        for (decltype(copy._height) i = 0; i < result_as_uint32._height; i++)
        {
            uint32_t *dst_row = (uint32_t*)&result_as_float.at(i, 0);
            uint32_t *src_row = (uint32_t*)&result_as_uint32.at(i, 0);
            for (decltype(copy._width) j = 0; j < result_as_uint32._width; j++)
            {
                uint32_t word = *src_row++;
                uint32_t toggle_mask = word & 0x80000000 ? 0x80000000 : 0xFFFFFFFF;
                *dst_row++ = word ^ toggle_mask;
            }
        }
        return result_as_float;
    }

    template <typename T>
    image<T> median7x7(image<T> &src, int mirror_scheme)
    {
        static_assert(std::is_integral<T>::value, "Floating point types not supported");
        static_assert(sizeof(T) <= 4, "Operation not supported for types larger than 4 bytes per element");
        using type = T;
        // add 3 pixel replicated border everywhere and pad to multiple of two)
        // if image dimension is odd, add one more pixel to allow kernel
        // to process everything in chunks of 2x2 blocks

        auto img = src.make_borders(3, 3, 3 + (src._height & 1), 3 + (src._width & 1), mirror_scheme);

        // we must iterate a support of 8x2 at all possible locations
        auto height = img._height - 6;
        auto width = img._width;
        uint64_t bits[sizeof(type) * 8] = { 0 };
        T correction = std::numeric_limits<T>::min();       // correction for signed int

        for (decltype(height) row = 0; row + 1 < height; row += 2)
        {
            for (decltype(width) col = 0; col + 1 < width; col += 2)
            {
                update_masks(bits, &img.at(row, col), width);
                if (col >= 6)
                {
                    const uint64_t mask = 0xfefefefefefefe00ULL;
                    const uint64_t threshold = 25;   // middle 25th item out of 49
                    img.at(row, col - 6) = (T)median8x8_iteration(bits, mask >> 9, threshold) ^ correction;
                    img.at(row, col - 5) = (T)median8x8_iteration(bits, mask >> 1, threshold) ^ correction;
                    img.at(row + 1, col - 6) = (T)median8x8_iteration(bits, mask >> 8, threshold) ^ correction;
                    img.at(row + 1, col - 5) = (T)median8x8_iteration(bits, mask >> 0, threshold) ^ correction;
                }
            }
        }
        return img.region(src._height, src._width);
    }

    template <int N> inline __m128i extract_top_byte(__m128i (&data)[N]);

    template <> inline __m128i extract_top_byte<1>(__m128i (&data)[1]) { return data[0]; }
    template <> inline __m128i extract_top_byte<2>(__m128i (&data)[2])
    {
        auto topA = _mm_srli_epi16(data[0], 8);
        data[0] = _mm_slli_epi16(data[0], 8);
        auto topB = _mm_srli_epi16(data[1], 8);
        data[1] = _mm_slli_epi16(data[1], 8);
        return _mm_packus_epi16(topA, topB);
    }

    inline __m128i top_from_i32(__m128i &item)
    {
        auto res = _mm_srli_epi32(item, 24);
        item = _mm_slli_epi32(item, 8);
        return res;
    }

    template <> inline __m128i extract_top_byte<4>(__m128i (&data)[4])
    {
        auto d0 = top_from_i32(data[0]);
        auto d1 = top_from_i32(data[1]);
        auto d2 = top_from_i32(data[2]);
        auto d3 = top_from_i32(data[3]);
        return _mm_packus_epi16(_mm_packus_epi32(d0, d1), _mm_packus_epi32(d2, d3));
    }

    template <typename T> inline
    void update_masks(uint64_t (&bits)[sizeof(T) * 8], T *origin, int stride)
    {
        T correction_factor = std::numeric_limits<T>::min();
        // Load 8x2 items to temporary memory (interleaving as 2x8 units)
        T tmp_data[2][8];
        for (int r = 0; r < 8; r++)
        {
            tmp_data[0][r] = origin[0] ^ correction_factor;
            tmp_data[1][r] = origin[1] ^ correction_factor;
            origin += stride;
        }

        // reorder data (interpreted as uint8_t data[16][number_of_bytes_per_element]
        __m128i data[sizeof(T)];
        for (int i = 0; i < (int)sizeof(T); i++)
            data[i] = _mm_loadu_si128(((__m128i *)tmp_data) + i);

        for (int i = 0; i < (int)sizeof(T); i++)
        {
            __m128i hibyte = extract_top_byte(data);
            uint64_t *b = bits + i * 8;
            for (int j = 0; j < 8; j++)
            {
                uint64_t topbits = (uint16_t)_mm_movemask_epi8(hibyte);
                hibyte = _mm_add_epi8(hibyte, hibyte);
                b[j] = (b[j] >> 16) | (topbits << 48);
            }
        }
    }

    /// given an array of kernel coefficients of length k_len
    /// produces a matrix of W*W coefficients, with W = min(k_len, image_width)
    /// In case that the image/kernel has to be mirrored and the the kernel_width is still larger
    /// than the extended width, the last image item is repeated
    /// Image:                             A B C
    /// Image Mirrored         =     C B A|A B C|C B A
    /// Mirrored and extended  = C C C B A|A B C|C B A A A A A ...
    /// While there is no practical reason to have so large filter, this is implemented
    /// as a free prevention of undefined behavior
    template <typename T>
    std::vector<T> replicate_filter_kernel(const T *kernel, int k_len, int image_width, int rep_scheme)
    {
        // For odd size filter |shift_left| == shift_right
        // for even size filter |shift_left| == shift_right + 1
        const int w = k_len < image_width ? k_len : image_width;
        const int kernel_mid_point = (k_len - 1) >> 1;

        std::vector<T> result(w * w, 0);
        std::vector<int> idx_map(k_len + w);  // maps index 0...k_len - mid_point to indices within 0..w

        for (int i = 0; i < (k_len + w); i++)
        {
            int idx = i - kernel_mid_point;
            if (idx < 0)
            {
                if ((rep_scheme & LEFT) == (ZERO_PAD & LEFT))
                    idx = -1;
                else if ((rep_scheme & LEFT) == (REPLICATE & LEFT))
                    idx = 0;
                else if ((rep_scheme & LEFT) == (MIRROR_EVEN & LEFT))
                    idx = std::min(-1 - idx, w - 1);
                else
                    idx = std::min(-idx, w - 1);
            }
            else if (idx >= w)
            {
                if ((rep_scheme & RIGHT) == (ZERO_PAD & RIGHT))
                    idx = -1;
                else if ((rep_scheme & RIGHT) == (REPLICATE & RIGHT))
                    idx = w - 1;
                else if ((rep_scheme & RIGHT) == (MIRROR_EVEN & RIGHT))
                    idx = std::max(0, w - (idx - w));
                else
                    idx = std::max(0, w - (idx - w + 1));
            }
            idx_map[i] = idx;
        }

        // loop over the rows of the result --
        // here the kernel middle point is thought to align with 'i'
        for (int i = 0; i < w; i++)
        {
            T *row = &result[i * w];
            // Sum up the kernel coefficient k_len with the proper index
            for (int j = 0; j < k_len; j++)
            {
                int idx = idx_map[i + j];
                if (idx >= 0)
                    row[idx] += kernel[j];
            }
        }
        return result;
    }

    template <typename T, typename U>
    U rounded_dot_product(const double *coeffs, const T *src, const int skip_x, const int N)
    {
        double sum = 0.0;
        for (int i = 0; i < N; i++)
            sum += coeffs[i] * src[i * skip_x];
        return reduce_to<U>(sum);
    }

    /// Filter a horizontal line using pre-computed kernels of double[w][k]
    /// where 'k' is the (maximum) filter kernel length and 'w' is the width of the output image
    /// the index vector is of length (w) and contains the first column in source image
    /// to be multiply accumulated with the corresponding kernel
    ///  -- 0<= x <w,  output[x] = sum i=0..k-1 input(offset[x] + i) * kernel[x][i];
    template <typename T, typename U>
    void filter1d_horizontal(std::vector<double> &coeffs, std::vector<int> &offsets, image<T> &src, image<U> &&dst)
    {
        const int h = static_cast<int>(src._height);
        const int w = static_cast<int>(src._width);
        const int h_out = static_cast<int>(dst._height);
        const int w_out = static_cast<int>(dst._width);

        if (coeffs.size() == 0 || offsets.size() == 0 || coeffs.size() % offsets.size() != 0)
            throw std::runtime_error("Coeffs size must be k * output_width");

        auto k = static_cast<int>(coeffs.size() / offsets.size());
        if (k > w)
            throw std::runtime_error("Kernel width can't be larger than input size");

        if (offsets.size() != static_cast<size_t>(w_out))
            throw std::runtime_error("Offsets size must match destination width");

        if (h_out > h)
            throw std::runtime_error("Destination image is not large enough");

        for (int y = 0; y < h_out; y++)
        {
            T *src_row = &src.at(y, 0);
            U *dst_row = &dst.at(y, 0);
            double *coeff_row = coeffs.data();
            int *offset = offsets.data();
            int dst_skipx = dst._skip_x;
            int src_skipx = src._skip_x;
            for (int x = 0; x < w_out; x++, dst_row += dst_skipx, coeff_row += k)
                *dst_row = rounded_dot_product<T, U>(coeff_row, src_row + src_skipx * offset[x], src_skipx, k);
        }
    }

    /// Filter a horizontal line
    /// The replication scheme is built in the array of coefficients
    /// coeffs == square matrix of at most size k*k or width*width
    /// ceil(k/2) first rows of the matrix contain coefficients for first k/2 columns in the src image
    /// the middle row is used for all middle columns in src
    /// and the last fix(k/2) rows are used for the last k/2 columns in the src image
    template <typename T, typename U>
    void filter1d_horizontal(std::vector<double> &coeffs, int k, image<T> &src, image<U> &dst)
    {
        const int h = (int)src._height;
        const int w = (int)src._width;

        if (dst._height < src._height || dst._width < src._width)
            throw std::runtime_error("Destination image is not large enough");  // Should we resize it?

        if ((int)coeffs.size() != k*k || k > w)
            throw std::runtime_error("Expecting square matrix not wider than image width");

        int bot_elems = k >> 1;
        int top_elems = (k - 1) >> 1;
        int middle_elems = w - top_elems - bot_elems;

        for (int j = 0; j < h; j++)
        {
            // iterate halfway the matrix using the same k items from the source row
            // and always a new row from the matrix
            T *src_row = &src.at(j, 0);
            U *dst_row = &dst.at(j, 0);
            int skip_x = dst._skip_x;
            int sx = src._skip_x;
            double *coeff_row = coeffs.data();
            // Iterate first half of the coeffs using the same source, but different coefficients
            for (int i = 0; i < top_elems; i++, coeff_row += k, dst_row += skip_x)
                *dst_row = rounded_dot_product<T,U>(coeff_row, src_row, sx, k);

            // Iterate the main (middle) area of the source image using the original kernel
            for (int i = 0; i < middle_elems; i++, src_row += sx, dst_row += skip_x)
                *dst_row = rounded_dot_product<T, U>(coeff_row, src_row, sx, k);

            if (middle_elems)
            {
                coeff_row += k;
                src_row -= sx;
            }

            // Iterate the last half of the coeffs using the same source (and different coefficients)
            for (int i = 0; i < bot_elems; i++, coeff_row += k, dst_row += skip_x)
                *dst_row = rounded_dot_product<T,U>(coeff_row, src_row, sx, k);
        }
    }

    /// Gaussian filter for LSC preprocessing
    template <typename T>
    void filter_separable(image<T> &src, std::vector<double> kernels, int rep_scheme)
    {
        // make a new image of H*W, then transpose the iterators
        // to make it W*H
        auto size = static_cast<int>(kernels.size());

        if (kernels.size() == 0)
            return;

        if (kernels.size() == 1 && kernels[0] != 1.0)
        {
            double scale = kernels[0];
            src.foreach([scale](T &dst)
            {
                dst = reduce_to<T>(scale * dst);
            });
            return;
        }

        auto tmp_image = image<double>(src._width, src._height); // transposed dimensions
        auto transpose = tmp_image.transpose();

        auto coeffs = replicate_filter_kernel(kernels.data(), size, src._width, rep_scheme);

        filter1d_horizontal(coeffs, size, src, transpose);

        coeffs = replicate_filter_kernel(kernels.data(), size, tmp_image._width, rep_scheme >> 8);
        auto src_tp = src.transpose();

        filter1d_horizontal(coeffs, size, tmp_image, src_tp);
    }

    /// Gaussian filter for LSC preprocessing
    template <typename T>
    void gaussian35x35(image<T> &src, int rep_scheme)
    {
        const double gaussian_coeffs35_sigma12[35] = {
            0.014248645487202, 0.015978517790508, 0.017794404801498, 0.019679520280453, 0.021613724692408,
            0.023573756271676, 0.025533598726188, 0.027464983023134, 0.029338013695678, 0.031121903026171,
            0.032785789743931, 0.034299612985591, 0.035635007639430, 0.036766184208162, 0.037670755260384,
            0.038330471561864, 0.038731834108016, 0.038866553395413, 0.038731834108016, 0.038330471561864,
            0.037670755260384, 0.036766184208162, 0.035635007639430, 0.034299612985591, 0.032785789743931,
            0.031121903026171, 0.029338013695678, 0.027464983023134, 0.025533598726188, 0.023573756271676,
            0.021613724692408, 0.019679520280453, 0.017794404801498, 0.015978517790508, 0.014248645487202
        };

        auto kernels = std::vector<double>(gaussian_coeffs35_sigma12, gaussian_coeffs35_sigma12 + 35);
        filter_separable(src, kernels, rep_scheme);
    }

    // maps the negative indexes `i` to positive range between [0...max_i[
    // - returns negative for zero pad case
    // - throws for error
    inline int mirror(int i, uint32_t rep_scheme, int max_i)
    {
        if (i < 0)
        {
            switch (rep_scheme & LEFT)
            {
            case (ZERO_PAD & LEFT) :
                return -1;
            case (REPLICATE & LEFT) :
                return 0;
            case (MIRROR_ODD & LEFT) :
                i = -i;
                break;
            case (MIRROR_EVEN & LEFT) :
                i = -i - 1;
                break;
            default:
                i = max_i;      // throws e.g. for REPLICATE_EVEN
                break;
            }
        }
        else if (i >= max_i)
        {
            switch (rep_scheme & RIGHT)
            {
            case (ZERO_PAD & RIGHT) :
                return -1;
            case (REPLICATE & RIGHT) :
                return max_i - 1;
            case (MIRROR_ODD & RIGHT) :
                i = max_i - (i - max_i + 2);
                break;
            case (MIRROR_EVEN & RIGHT) :
                i = max_i - (i - max_i + 1);
                break;
            default:
                i = max_i;      // throws e.g. for REPLICATE_EVEN
                break;
            }
        }
        if (i >= max_i)
            throw std::runtime_error("Failed to map index to range");
        return i;
    }

    // Given weights for indices between [left_index ... left_index + weights.size() [
    //   - fold the negative indices / indices larger than input width (w_in) back to [0 ... w_in[
    //   - remove leading zeros from weights (if possible)
    //   - return the index of last non-zero in the vector (plus one)
    // throws for invalid rep_scheme
    //   - replicate_even must know the block size to be replicated -- this mode is not applicable anyway
    //   - mirroring schemes may run out of data -- it's hard to understand why these would be needed...
    //   - zero pad and replicate (last column) should always succeed
    // Returns a remapped vector of weights -- adjusts the starting indexes if needed
    inline std::vector<double> remap_to_range(
        std::vector<double> &weights,       // [w_out][k'] sized vector
        std::vector<int> &indices,          // [w_out] indices to input image
        int w_in,                           // width out input image -- for detection of out of image weights
        int min_x,                          // first index in each row to process
        int max_x,                          // last index in each row to process
        uint32_t rep_scheme)                // how to map out-of-image pixels
    {
        if (weights.size() == 0 || indices.size() == 0 || weights.size() % indices.size() != 0)
            throw std::runtime_error("weights size must be an integer multiple of index vector size");

        auto k = static_cast<int>(weights.size() / indices.size());
        if (min_x < 0 || min_x >= k || max_x < min_x || max_x >= k)
            throw std::runtime_error("Min_x, max_x must be contained within the `k` in weights.size() == indices.size() * k");

        auto k_prime = max_x - min_x + 1;            // Size of the non-zero portion of the kernels (weights)
        auto k_final = std::min(k_prime, w_in);      // There can't be more columns than the size of the input

        auto output = std::vector<double>(k_final * indices.size());    // double output[w_out][k_final];

        for (size_t x = 0; x < indices.size(); x++)
        {
            enum out_of_range_pixels_e { none = 0, left = 1, right = 2, both = left + right };
            auto &left_index = indices[x];
            left_index += min_x;
            // remap the indexes from "left_index + [min_x ... max_x]"  to [0 .. k_final[
            if (left_index < 0 || left_index + k_prime > w_in)
            {
                auto offset = left_index < 0 ? left_index : left_index + k_final - w_in;
                left_index -= offset;

                for (int i = 0; i < k_prime; i++)
                {
                    int i_prime = mirror(i + offset, rep_scheme, k_final);
                    if (i_prime >= 0)
                        output[x * k_final + i_prime] += weights[x * k + i + min_x];
                }
            }
            else
            {
                for (int i = 0; i < k_prime; i++)
                    output[x * k_final + i] = weights[x * k + i + min_x];
            }
        }

        return output;
    }

    /// Given input image width w_in, output image width w_out and scaling factor
    /// Generate a matrix of double[w_out][k'] coefficients and an associated
    /// vector of indexes[w_out], where index[i] means the offset in input image
    /// to be multiplied by the corresponding kernel coefficient coeff[i][0]
    /// Each row of coefficients are first normalized to unity gain, and then
    /// remapped to input image range [0..w_in[ using the given repetition scheme
    /// - supports:   mirror even / odd, zero pad and replicate (no replicate even)
    inline std::vector<double> make_bicubic_filter_kernel(
        int w_in, int w_out, double scale, std::vector<int> &indices, uint32_t rep_scheme)
    {
        auto bicubic = [](double x)
        {
            x = std::abs(x);
            auto x2 = x * x;
            auto x3 = x * x2;
            if (x <= 1.0)
                return 1.5 * x3 - 2.5 * x2 + 1;
            if (x <= 2.0)
                return -0.5 * x3 + 2.5*x2 - 4.0 * x + 2;
            return 0.0;     // out of range -2...2
        };

        if (w_out <= 0 || w_in <= 0)
            return{};

        if (scale == 0)
            scale = (double)w_out / (double)w_in;

        // when up-scaling, we use kernel width = 4 as in bicubic resampling
        // - when down-scaling we combine the resampling filter with low pass filtering
        auto kernel_scale = std::min(1.0, scale);
        auto kernel_width = scale >= 1.0 ? 4.0 : 4.0 / scale;

        auto k = static_cast<int>(std::ceil(kernel_width) + 2);

        std::vector<double> output(k * w_out);
        indices.resize(w_out);
        auto min_non_zero_index = k;      // keep track of full columns to prune out
        auto max_non_zero_index = 0;
        for (int x = 0; x < w_out; x++)
        {
            // center of kernel aligned to input space
            auto u = ((x + 1) / scale + 0.5*(1 - 1.0 / scale));

            // left most pixel in input image possibly involved in calculation (+1)
            auto left_idx = static_cast<int>(std::floor(u - kernel_width * 0.5));

            // Step 1)
            // Calculate the interpolation weights for [left_ind ... + k[
            auto sum = 0.0;
            for (int i = 0; i < k; i++)
            {
                auto x_fractional = u - (left_idx + i);
                auto weight = kernel_scale * bicubic(x_fractional * kernel_scale);
                if (weight != 0)
                {
                    min_non_zero_index = std::min(min_non_zero_index, i);
                    max_non_zero_index = std::max(max_non_zero_index, i);
                }
                output[x * k + i] = weight;
                sum += weight;
            }
            // Step 2) Normalize the current row of weights to sum == 1.0 (+- epsilon)
            auto const epsilon = 1e-7;
            if (sum != 0.0 && std::fabs(sum - 1.0) > epsilon)
            {
                sum = 1.0 / sum;
                for (int i = 0; i < k; i++)
                    output[x * k + i] *= sum;
            }

            indices[x] = left_idx - 1;          // convert to 0-based index
        }

        return remap_to_range(output, indices, w_in, min_non_zero_index, max_non_zero_index, rep_scheme);
    }

    /// Resize image with given ratio - final image size = ceil(input_size * scale)
    template <typename T>
    image<T> resize_image(image<T> &input, point_xy scale, uint32_t rep_scheme = REPLICATE)
    {
        auto input_size = input.size();
        auto dst_size = roi_point(std::ceil(input_size._x * scale.x), std::ceil(input_size._y * scale.y));
        if (dst_size._x == 0 || dst_size._y == 0)
            throw std::runtime_error("Zero output image size");

        auto dst_image = image<T>(dst_size);
        std::vector<int> indices;

        // there's an intermediate image of type double
        auto tmp_image = image<double>(dst_image._width, input._height); // transposed dimensions
        auto coeffs = make_bicubic_filter_kernel(input._width, dst_image._width, scale.x, indices, rep_scheme);
        filter1d_horizontal(coeffs, indices, input, tmp_image.transpose());

        coeffs = make_bicubic_filter_kernel(input._height, dst_image._height, scale.y, indices, rep_scheme >> 8);
        filter1d_horizontal(coeffs, indices, tmp_image, dst_image.transpose());
        return dst_image;
    }

    /// Resize image to given pixel size (and type)
    template <typename T, typename U>
    void resize_image(image<T> &input, image<U> &dst, uint32_t rep_scheme = REPLICATE)
    {
        auto input_size = input.size();
        auto dst_size = dst.size();
        if (dst_size._x == 0 || dst_size._y == 0 || input_size._x == 0 || input_size._y == 0)
            throw std::runtime_error("Zero image size");

        std::vector<int> indices;
        point_xy scale((double)dst_size._x / input_size._x, (double)dst_size._y / input_size._y);

        // there's an intermediate image of type double
        auto tmp_image = image<double>(dst._width, input._height);
        auto coeffs = make_bicubic_filter_kernel(input._width, dst._width, scale.x, indices, rep_scheme);
        filter1d_horizontal(coeffs, indices, input, tmp_image.transpose());

        coeffs = make_bicubic_filter_kernel(input._height, dst._height, scale.y, indices, rep_scheme >> 8);
        filter1d_horizontal(coeffs, indices, tmp_image, dst.transpose());
    }

    // Rotate image to output of given size
    template <typename T, typename U>
    void rotate_image(image<T> &input, image<U> &output, double degrees, U outside_value = 0)
    {
        using point_xy = point_xy;
        auto rotate_by_pt = [](point_xy p, point_xy rot)
        {
            return point_xy(p.x * rot.x - p.y * rot.y, p.y * rot.x + p.x * rot.y);
        };
        auto in_size = input.size();
        auto out_size = output.size();
        auto radians = degrees * (3.141592653589793 / 180.0);
        auto offset = point_xy(in_size._x - 1.0, in_size._y - 1.0) * 0.5;
        auto offset_out = point_xy(out_size._x - 1.0, out_size._y - 1.0) * 0.5;
        auto delta = point_xy(std::cos(radians), std::sin(radians));

        for (int y = 0; y < out_size._y; y++)
        {
            // Calculate sampling point in input space for the first pixel of each scanline
            auto start = offset + rotate_by_pt(point_xy(-offset_out.x, y - offset_out.y), delta);
            for (int x = 0; x < out_size._x; x++)
            {
                auto f = start + delta * x;      // floating point sampling point of the rest of the points in scanline
                // These four lines could be part of libimage api
                // - output = input(point_xy, interpolation_e method = interpolation_e::nearest, 0.0)
                auto p = roi_point(f.x, f.y);    // point rounded to nearest integral point
                output(y, x) = (p._x < 0 || p._y < 0 || p._x >= in_size._x || p._y >= in_size._y)
                    ? outside_value
                    : static_cast<U>(input(p));
            }
        }
    }

    template <typename T>
    image<T> rotate_image(image<T> &input, double degrees, T outside_value = 0)
    {
        using point_xy = point_xy;
        auto size = input.size();
        auto half_size = point_xy((size._x - 1) * 0.5, (size._y - 1) * 0.5);
        auto radians = degrees * (3.141592653589793 / 180.0);
        // Create a vector containing two extreme points to find out extremes of outside image
        auto vec = std::vector<point_xy>({ half_size, half_size * point_xy(1.0, -1.0)});
        auto rotated = rotate_vector(vec, radians);
        auto max_x = 1 + 2 * static_cast<int>(std::ceil(std::max(std::abs(rotated[0].x), std::abs(rotated[1].x))));
        auto max_y = 1 + 2 * static_cast<int>(std::ceil(std::max(std::abs(rotated[0].y), std::abs(rotated[1].y))));
        auto output = image<T>(max_y, max_x);

        rotate_image(input, output, degrees, outside_value);
        return output;
    }


    /// Trimmed mean finder -- reference implementation using std::sort just to make it work
    ///  - the vector storing the blocks will be reused for performance
    template <typename T>
    struct trimmed_mean_f
    {
        trimmed_mean_f() = default;
        std::vector<T> data;            // Re-use the data for successive calls

        /// \brief              Calculate mean of image region without outliers
        /// \param  block       Region to evaluate
        /// \param  trim        Total percentage of samples to exclude (between 0 and 1.0)
        /// \return             Average of sample with outliers removed
        double operator()(image<T> &&block, double trim)
        {
            size_t size = static_cast<size_t>(block._width * block._height);
            size_t first_idx = std::min(size, static_cast<size_t>(std::max(0.0, 0.5 + 0.5 * trim * size)));
            size_t last_idx = size - first_idx;
            if (first_idx >= last_idx)
                return 0.0;

            if (size > data.size())
                data = std::vector<T>(size);

            block.copy_to(data.data());     // serialize the data for sorting

            // No need to sort if we want mean of everything
            if (first_idx != 0)
                std::sort(data.begin(), data.begin() + size);

            double sum = 0.0;
            for (size_t i = first_idx; i < last_idx; i++)
                sum += data[i];
            return sum / static_cast<double>(last_idx - first_idx);
        }
    };

    // Calculates Integral of image in-place
    template <typename T>
    void summed_area_table(image<T> &input)
    {
        // accumulate first row
        for (decltype(input._width) i = 1; i < input._width; i++)
            input(0, i) += input(0, i - 1);

        // accumulate first column
        for (decltype(input._height) i = 1; i < input._height; i++)
            input(i, 0) += input(i - 1, 0);

        int skipx = -input._skip_x;
        int skipy = -input._skip_y;
        input.region(input.size() - 1, roi_point(1)).
            foreach([skipx, skipy](T &dst)
        {
            T *src = &dst;
            dst += src[skipx] + src[skipy] - src[skipx + skipy];
        });
    }

    /// Box Filter -- filters input image with average 1/(kernel_size * kernel_size)
    template <typename T>
    void box_filter(image<T> &input, int kernel_size, int repetition_scheme)
    {
        if (kernel_size <= 1)
            return;

        // copy the input image to a temporary buffer for SAT generation
        auto excess = std::max(kernel_size, 1) - 1;
        auto left = excess / 2;
        auto right = excess - left;
        auto temp_image = input.make_borders(left + 1, left + 1, right, right, repetition_scheme);

        temp_image.rows(1, 0).fill(0);
        temp_image.columns(1, 0).fill(0);

        summed_area_table(temp_image);
        auto temp_image_aligned = temp_image.region(input.size(), roi_point(kernel_size));

        auto skip_x = -kernel_size;
        auto skip_y = -kernel_size * temp_image_aligned._skip_y;
        double factor = 1.0 / (kernel_size * kernel_size);

        input.foreach([skip_x, skip_y, factor](T &dst, T &src)
        {
            T *s = &src;
            double tmp = factor * (s[0] - s[skip_x] - s[skip_y] + s[skip_x + skip_y]);
            dst = reduce_to<T>(tmp);
        }, temp_image_aligned);
    }

    // Basic properties of connected components
    struct region_props
    {
        double area = 0.0;
        double perimeter = 0.0;
        double centroid_x = 0.0;
        double centroid_y = 0.0;
        roi_point top_left;     // bounding box
        roi_point bot_right;    // bounding box
        region_props() = default;
    };

    // Short cut to updating minimum roi point
    inline roi_point& operator <<= (roi_point &a, const roi_point &b)
    {
        a._x = std::min(a._x, b._x);
        a._y = std::min(a._y, b._y);
        return a;
    }

    // Short cut to updating maximum roi point
    inline roi_point& operator >>= (roi_point &a, const roi_point &b)
    {
        a._x = std::max(a._x, b._x);
        a._y = std::max(a._y, b._y);
        return a;
    }

    template <bool add> void add_to_list(std::vector<roi_point> *visited_points, roi_point p);
    template <> inline void add_to_list<false>(std::vector<roi_point> *, roi_point) { };
    template <> inline void add_to_list<true>(std::vector<roi_point> *visited_points, roi_point p)
    {
        if (visited_points != nullptr)
            visited_points->push_back(p);
    };


    // Locates centroids, area, perimeter and squareness by tracing the perimeter of a feature
    // - treats the feature 8-connected
    // - automatically removes holes
    // todo: static_assert(connectivity == 4 || connectivity == 8, "Region props supports only 4 and 8-connected neighborhoods");
    template <typename T, bool record_visited_points = false>
    region_props get_regionprops(image<T> &im, roi_point point, std::vector<roi_point> *list_of_points = nullptr)
    {
        auto is_inside = [](roi_point &p, roi_point &dims)
        {
            return p._x >= 0 && p._x < dims._x && p._y >= 0 && p._y < dims._y;
        };

        int32_t segments_straight = 0;
        int32_t segments_diagonal = 0;
        int64_t sum_x = 0;
        int64_t sum_y = 0;
        int32_t area = 0;

        auto size = im.size();
        if (!is_inside(point, size))
            return{};

        T start = im(point);
        while (point._x > 0 && im(point._y, point._x - 1) == start)
            point._x--;

        auto origin = point;
        roi_point top_left = origin;
        roi_point bot_right = origin;

        // we have traversed to left side of the feature - start tracing the boundary
        // initial direction is "up"
        // we can iterate at least "distance_to_border" items without checking for boundary conditions
        // int distance_to_border = std::min({ point._x, point._y, size._x - point._x, size._y - point._y });
        int dir = 0;
        const roi_point north[4] = { { 0, -1 }, { 1, 0 }, { 0, 1 }, { -1, 0 } };
        const roi_point nw[4] = { { -1, -1 }, { 1, -1 }, { 1, 1 }, { -1, 1 } };

        add_to_list<record_visited_points>(list_of_points, point);

        do
        {
            top_left <<= point;     // update bounding box statistics: minimum of two points
            bot_right >>= point;    // maximum of two points

            // update statistics if we go up / down
            if (dir == 0)
            {
                area -= point._x;
                sum_y -= point._x * point._y;
                sum_x -= point._x * (point._x - 1); // actually / 2, but we delay it
            }
            else if (dir == 2)
            {
                area += (point._x + 1);
                sum_y += (point._x + 1) * point._y;
                sum_x += (point._x + 1) * point._x; //  actually / 2, but we delay it
            }

            auto p = point + nw[dir];
            // if we can go north west, we go and turn left
            if (is_inside(p, size) && im(p) == start)
            {
                point = p;
                dir = (dir - 1) & 3;
                ++segments_diagonal;
                add_to_list<record_visited_points>(list_of_points, point);
            }
            else
            {
                // if we can go north, we do
                auto n = point + north[dir];
                if (is_inside(n, size) && im(n) == start)
                {
                    point = n;
                    ++segments_straight;
                    add_to_list<record_visited_points>(list_of_points, point);
                }
                else
                {
                    // turn right on the spot
                    dir = (dir + 1) & 3;
                }
            }
        } while (dir || !(point == origin));

        const double sqrt_2 = 1.414213562373095;
        region_props result;
        if (area <= 0)
            return result;

        result.area = area;
        result.perimeter = sqrt_2 * segments_diagonal + segments_straight;
        result.centroid_x = (sum_x / 2) / result.area;
        result.centroid_y = sum_y / result.area;
        result.bot_right = bot_right;
        result.top_left = top_left;
        return result;
    }


    inline int calculate_otsu_threshold(std::vector<int> &histogram, int min_index, int max_index, size_t total_pixels)
    {
        if (histogram.size() <= static_cast<size_t>(max_index))
            throw std::runtime_error("Histogram vector too short");

        auto weight_A = 0.0;
        auto weight_B = static_cast<double>(total_pixels);
        auto sum_A = 0.0;
        auto sum_B = 0.0;
        auto max = 0.0;
        int idx = min_index;

        for (int i = min_index; i <= max_index; i++)
            sum_B += (double)i * histogram[i];

        for (int i = min_index; i <= max_index; i++)
        {
            if (histogram[i] == 0)
                continue;

            double data = histogram[i];
            weight_A += data;
            weight_B -= data;
            sum_A += i * data;
            sum_B -= i * data;

            // theoretically these weights should not be zero -- since they are cumulative histograms
            if (weight_A == 0 || weight_B == 0)
                continue;

            double diff = sum_A / weight_A - sum_B / weight_B;
            double variance = weight_A * weight_B * diff * diff;
            if (variance >= max)
            {
                max = variance;
                idx = i;
            }
            else
                break;      // One local maximum
        }
        return idx;
    }

    inline int calculate_otsu_threshold(std::vector<int> &histogram)
    {
        if (histogram.size() == 0)
            return 0;

        size_t sum = 0;
        for (auto &x : histogram)
            sum += x;

        return calculate_otsu_threshold(histogram, 0, static_cast<int>(histogram.size() - 1), sum);
    }

    // The quantizer manages the memory for histogram and quantized image data
    // reusing the image data between each call to the quantizer
    // usage:    `quantizer red; for (auto &x: views) { tmp = red(x); do_something_with(tmp); }`
    class quantizer
    {
    public:
        quantizer() : histogram(65536) { }
        image<uint16_t> operator() (image<uint16_t> view)
        {
            auto size = view.size();
            auto count = static_cast<size_t>(size._x * size._y);
            if (image_data.size() < count)
                image_data.resize(count);

            auto result = image<uint16_t>(size, image_data.data());
            if (count == 0)
                return result;

            // serialize the data to container while taking some statistics
            //  - histogram
            //  - minimum and maximum
            int *hist = histogram.data();

            uint16_t value_min = view(0, 0);
            uint16_t value_max = value_min;

            result.foreach([hist, &value_min, &value_max](uint16_t &dst, uint16_t &src)
            {
                dst = src;
                if (dst < value_min)
                    value_min = dst;
                if (dst > value_max)
                    value_max = dst;
                hist[src]++;
            }, view);

            auto threshold = static_cast<uint16_t>(calculate_otsu_threshold(histogram, value_min, value_max, count));
            for (unsigned int i = value_min; i <= value_max; i++)
                hist[i] = 0;

            result.foreach([threshold](uint16_t &inout) { inout = inout > threshold ? 1 : 0; });
            return result;
        }
    private:
        std::vector<int> histogram;
        std::vector<uint16_t> image_data;
    };

    /// Section Computational Geometry

    // Find (index to) Root of label with path compression:
    //  - all encountered indices will be replaced by direct links to root
    //   - reduces worst case complexity of ANY scheme to log M (M == number of disjoint sets)
    // Root is marked by a special tag
    //  - TOP 4 bits encode Rank
    //  - Bottom 28 bits of ROOT are freely reusable
    // No range check done
    //  - circular link will result in stack overflow
    //  - when join by rank is utilized, the maximum length of the chain is 4 (or 5)
    //    and (tail) recursion can be replaced by unrolling the loop
    static inline uint32_t parent(uint32_t label, uint32_t *all_labels)
    {
        uint32_t tmp = all_labels[label];
        if (tmp & 0xf0000000)
            return label;
        tmp = parent(tmp, all_labels);
        all_labels[label] = tmp;
        return tmp;
    }

    /// \brief      Label each connected region with unique label
    /// \param      image - image to convert
    /// \param      connectivity - 4 or 8
    /// \returns    Number of connected areas, negative for errors
    template <typename T> int32_t bwlabel(image<T> &img, int connectivity = 8)
    {
        const uint32_t rankmask = 0xF0000000;
        switch (connectivity)
        {
        case 4:
            connectivity = 0;
            break;
        case 8:
            connectivity = 1;
            break;
        default:
            // Illegal parameter
            return -1;
        }
        // Make a buffer for labels from original image, surrounded by 1 pixel wide zero border
        //  - no need to check availability of NW/N/NE neighbors
        auto bitonal = image<uint32_t>(img.size() + 2);
        bitonal.copy_borders(1, 1, 1, 1, ZERO_PAD);
        bitonal.region(img.size(), roi_point(1)).init_from(img);

        int32_t stride_left = bitonal._skip_y + connectivity;
        int32_t stride_right = bitonal._skip_y - connectivity;

        uint32_t *all_labels = &bitonal(0, 0);          // Store the labels to self
        uint32_t *ptr = &bitonal(1, 1);              // Start of image
        uint32_t *end = &bitonal(bitonal._height - 1, 0);  // End of image

        // Pass 1 -- Separate the whole image into spans of non-zeros
        // Unify current span and all spans 'A', 'B',  ... 'C' directly above the current span
        //       [8][A]   [B ]  [C    ][8]   -- when 8-connectivity is used, the search range
        //          [Current Span     ]      -- is extended by one pixel in each direction
        while (ptr != end)
        {
            // Skip zeros
            while (ptr != end && ptr[0] == 0)
                ptr++;
            auto *begin = ptr;
            while (ptr != end && ptr[0] != 0)
                ptr++;
            // Found a span of non-zeroes (scan the previous line +- 1 pixel for 8-connected areas)
            if (begin == ptr)
                continue;
            uint32_t *prev = begin - stride_left;
            uint32_t *pend = ptr - stride_right;
            uint32_t top_root = 0;     // memorizes the label with highest rank in previous line

            // Split previous line (and it's extensions) to spans with identical labels
            while (prev != pend)
            {
                uint32_t last_label = 0;
                while (prev != pend && (last_label = *prev) == 0)
                    prev++;

                if (last_label == 0)
                    break;

                last_label = parent(last_label & ~rankmask, all_labels);

                while (++prev != pend && *prev != 0);

                // If current span is not associated with a label, take the first from previous line
                if (top_root == 0)
                {
                    top_root = last_label;
                }
                else if (top_root != last_label)
                {
                    // Current span has different parent than in the previous line
                    // - reroot the smaller branch to larger branch
                    //   - if branches are of equal height (rank), select one arbitrarily and increase its rank
                    uint32_t top_rank = all_labels[top_root] & rankmask;
                    uint32_t tmp_rank = all_labels[last_label] & rankmask;

                    if (top_rank > tmp_rank)
                    {
                        all_labels[last_label] = top_root;
                    }
                    else if (tmp_rank > top_rank)
                    {
                        all_labels[top_root] = last_label;
                        top_root = last_label;
                    }
                    else
                    {
                        all_labels[last_label] = top_root;
                        all_labels[top_root] -= rankmask;  // Increase by 1 == decrease by -1
                    }
                }
            }
            // Must generate new root -- the first item in the new span is root with rank = 1
            // The 28 low bits in the root are redundant -- these could be used to store a single statistic (e.g. area)
            // or to point to a struct containing the statistics (including centroid, representative pixel, bounding box)
            if (top_root == 0)
            {
                top_root = (uint32_t)(begin - all_labels);      // Index to top_root
                *begin++ = top_root - rankmask;                // Set rank to 1 or -(-1)
            }
            // Fill the current span with index to parent
            while (begin != ptr)
                *begin++ = top_root;
        }

        // Pass 2 -- Rename temporary labels to final labels
        //  - if the parent would store the pixel count of the label, we could immediately filter results by size
        uint32_t last_label = 0;
        uint32_t last_parent = 0;
        struct {
            uint32_t *labels;
            uint32_t label_count;
            uint32_t rename(uint32_t label)
            {
                const uint32_t rank_mask = 0xF0000000;
                if ((label & rank_mask) == 0)
                {
                    // Replace indirections to root with root
                    uint32_t tmp = labels[label];
                    tmp = rename(tmp);
                    labels[label] = tmp;
                    label = tmp;
                }
                else
                {
                    if (label < rank_mask)
                    {
                        uint32_t tmp = ++label_count;
                        labels[label & ~rank_mask] = tmp;
                        label = tmp;
                    }
                }
                return label;
            }
        } tmp = { all_labels, rankmask };      // First unique label must have rankmask bits set

        for (int j = 0; j < (int)img._height; j++)
        {
            const int skipx = img._skip_x;
            T *dst = &img(j, 0);
            uint32_t *src = &bitonal(j + 1, 1);
            uint32_t *end_of_line = src + img._width;
            while (src != end_of_line)
            {
                // Rename all non-zero labels to root labels
                // - memorize the root of last non-zero label
                uint32_t label = *src++;
                if (label > 0)
                {
                    if (label != last_label)
                    {
                        last_label = label;
                        last_parent = tmp.rename(last_label) & ~rankmask;
                    }
                    label = last_parent;
                }
                *dst = (T)label;
                dst += skipx;
            }
        }

        return tmp.label_count & ~rankmask;
    }
}


