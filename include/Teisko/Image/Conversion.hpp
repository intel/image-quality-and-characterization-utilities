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
#include <cstdint>
#include <memory>
#include <vector>
#include <algorithm>
#include <iterator>

#if (defined(WIN32) || defined(_WIN32))
#include "intrin.h"
#else
// Linux
#include "x86intrin.h"
#endif

namespace Teisko
{
    template <typename T, T shift>
    typename std::enable_if<std::is_integral<T>::value, T>::type
        round_shift(int value) {
        return (T)((value + (1 << (shift - 1))) >> shift);
    }

    template <typename T, T shift>
    typename std::enable_if<std::is_floating_point<T>::value, T>::type
        round_shift(T value) {
        return value * ((T)1 / (1 << shift));
    }

    template <typename T>
    image<T> chroma_upscale(image<T> &input)
    {
        // Assumes Chroma is subsampled symmetrically in center of the quad
        //  - based on the fact that the reconstruction kernel is symmetric
        //  -- K =  [ 1 3 ]   [3 1 ]
        //          [ 3 9 ]   [9 3 ]
        //          ----------------
        //          [ 3 9 ]   [9 3]
        //          [ 1 3 ]   [3 1]

        // For Each U or V sample, we generate four outputs using the kernels above
        // The Anchor point is always the item with weight '9'
        //   - out of image pixels (in U/V) follow replication scheme
        auto temp = input.make_borders(1, 1, 1, 1, REPLICATE);
        auto output = image<T>(input.size() * 2);
        int skipy_src = temp._skip_y;
        int skipy_dst = output._skip_y;
        // We don't align the source buffer to the middle in `temp` -- but to the top left corner
        //  - this allows us to use positive offsets -- 0,1,2 horizontally
        //    and 0,skipy,skipy*2 vertically
        auto temp_top_left = temp.region(input.size());
        output.subview(2, 2).foreach([skipy_dst, skipy_src](T &dst, T &src)
        {
            T *d = &dst;
            T *s = &src;
            d[0] = round_shift<T, 4>(s[0] + s[1] * 3 + s[skipy_src] * 3 + s[skipy_src + 1] * 9);
            d[1] = round_shift<T, 4>(3 * s[1] + s[2] + s[skipy_src + 1] * 9 + s[skipy_src + 2] * 3);
            d[skipy_dst] = round_shift<T, 4>(s[skipy_src] * 3 + s[skipy_src + 1] * 9 + s[2 * skipy_src] + s[2 * skipy_src + 1] * 3);
            d[skipy_dst + 1] = round_shift<T, 4>(s[skipy_src + 1] * 9 + s[skipy_src + 2] * 3 + s[2 * skipy_src + 1] * 3 + s[2 * skipy_src + 2]);
        }, temp_top_left);
        return output;
    }

    ///   b) convert from uint16_t to int32_t
    ///      |a b|
    ///      |c d| --> |a+3c, e+3c, b+3d, f+3d|
    ///      |e f|
    ///   - strategy:   read 'a b', 'e f' in other register,  'c d' in other
    inline void sum_up_two_rows(__m128i &result, char *ptr, int stride)
    {
        __m128i mid_row = _mm_cvtsi32_si128(*reinterpret_cast<int*>(ptr + 2 * stride));
        __m128i top_row = _mm_cvtsi32_si128(*reinterpret_cast<int*>(ptr + 0 * stride));
        __m128i bot_row = _mm_cvtsi32_si128(*reinterpret_cast<int*>(ptr + 4 * stride));
        mid_row = _mm_unpacklo_epi16(mid_row, mid_row);     // c c d d 0 0 0 0
        top_row = _mm_unpacklo_epi16(top_row, bot_row);     // a e b f 0 0 0 0
        mid_row = _mm_cvtepu16_epi32(mid_row);              // c c d d as int32
        top_row = _mm_cvtepu16_epi32(top_row);              // a e b f as int32
        top_row = _mm_add_epi32(top_row, mid_row);
        mid_row = _mm_add_epi32(mid_row, mid_row);
        result = _mm_add_epi32(top_row, mid_row);
    }

    /// Generic upscaling using kernels of form [a + 3b + 3c + 9d + 8]>>16
    /// uses 32-bit arithmetic -- is suitable for bit depths 13..16
    /// Performance is about 2x as fast as the non-SIMD, but ~33% slower than the fully 16-bit version
    /// Uses one register as delay line containing top+3*mid pixels and bot+3*mid pixels as 32-bit integers
    inline image<uint16_t> chroma_upscale_asm_32(image<uint16_t> &input)
    {
        auto size = input.size();
        auto right_pad = size._x & 1;   // we consume 2 items at a time -- pad odd widths to even
        auto temp = input.make_borders(1, 1, 1, 1 + right_pad, REPLICATE);
        auto output = image<uint16_t>(input.size() * 2);
        auto stride_y = temp._skip_y;
        const __m128i eight = _mm_set1_epi32(8);        // for rounding
        const __m128i shuffle = _mm_set_epi8(15, 14, 7, 6, 11, 10, 3, 2, 13, 12, 5, 4, 9, 8, 1, 0);

        for (int j = 0; j < size._y; j++)
        {
            auto *ptr = reinterpret_cast<char*>(&temp(j, 0));               // beginning of the row
            auto *dst1 = &output(j * 2, 0);
            auto *dst2 = &output(j * 2 + 1, 0);
            __m128i prev = _mm_setzero_si128();
            __m128i next = _mm_setzero_si128();
            sum_up_two_rows(next, ptr, stride_y);
            auto width = size._x;
            while (width > 0)
            {
                width -= 2;
                ptr += 4;
                prev = next;
                sum_up_two_rows(next, ptr, stride_y);
                // prev | next
                // a c | e g  --> calculate [a c] + [3c 3e],  [3c 3e] + [e g]
                // b d | f h                [b d] + [3d 3f],  [3d 3f] + [f h]
                auto odd = _mm_sub_epi32(next, prev);
                auto mid = _mm_alignr_epi8(next, prev, 8);  // b d e g
                prev = _mm_add_epi32(prev, mid);    // left + mid
                mid = _mm_add_epi32(mid, mid);      // mid * 2
                prev = _mm_add_epi32(prev, mid);    // a + c + 2c   = output_even(without shift + round)
                odd = _mm_add_epi32(odd, prev);     // a + c + 2c + e - a = output_odd(without shift / round)
                prev = _mm_add_epi32(prev, eight);
                odd = _mm_add_epi32(odd, eight);
                prev = _mm_srli_epi32(prev, 4);
                odd = _mm_srli_epi32(odd, 4);
                prev = _mm_packus_epi32(prev, odd); // result: top0 bot0 top2 bot2 top1 bot1 top3 bot3 as uint16_t
                // must be converted to order of top0 top1 top2 top3 ... bot0 bot1 bot2 bot3
                prev = _mm_shuffle_epi8(prev, shuffle);
                if (width >= 0)
                {
                    _mm_storel_epi64(reinterpret_cast<__m128i*>(dst1), prev);
                    _mm_storeh_pi(reinterpret_cast<__m64*>(dst2), _mm_castsi128_ps(prev));
                    dst1 += 4;
                    dst2 += 4;
                }
            }
            if (width < 0)
            {
                const __m128i mask_one = _mm_set_epi32(0, 0, 0, -1);
                auto top_to_bot = _mm_shuffle_epi32(prev, 0xee);
                _mm_maskmoveu_si128(prev, mask_one, reinterpret_cast<char*>(dst1));
                _mm_maskmoveu_si128(top_to_bot, mask_one, reinterpret_cast<char*>(dst2));
            }
        }
        return output;
    }

    // Sum up incoming rows:    weighted sum of three rows at `ptr[stride * {0,1,2}]` into two outputs
    //    top_next = top_row + 3 * mid_row;
    //    bot_next = bot_row + 3 * mid_row;
    static inline void sum_up_two_rows(__m128i &top, __m128i &bot, uint16_t *ptr, int stride, __m128i multiplier)
    {
        auto mid = _mm_mullo_epi16(_mm_loadu_si128(reinterpret_cast<__m128i*>(ptr + stride)), multiplier);
        top = _mm_loadu_si128(reinterpret_cast<__m128i*>(ptr));
        bot = _mm_loadu_si128(reinterpret_cast<__m128i*>(ptr + 2 * stride));
        top = _mm_add_epi16(top, mid);
        bot = _mm_add_epi16(bot, mid);
    }

    ///   Processes 8 pixels at a time using SIMD
    ///   In regular domain we would use `[a b;g h]` to derive the top left output
    ///   - SIMD implementation uses vectors containing `[abcd... bcde...;ghij... hijk...]`
    ///   [a][b c d e][f]    Each row will be stored as ......ab | cdef
    ///   [g][h i j k][l]
    ///   [m][n o p q][r]
    ///   - other mechanisms considered
    ///   a) read [a b c d ...] in one registers, [b c d e ... ] in second, [c d e f ...] in third
    ///   b) pad the incoming items in column major (after `sum up two rows`)
    ///      0 2 4 6|8 a c e
    ///      1 3 5 7|9 b d f
    ///      - this has the possibility of processing just one item at a time
    inline image<uint16_t> chroma_upscale_asm(image<uint16_t> &input, int bits)
    {
        auto size = input.size();
        if (size == roi_point(0))
            return input;

        // The algorithm can only tolerate bit depths of 12 or less
        // since the intermediate value (1a + 3b + 3c + 9d + 8) must fit into uint16_t
        if (bits > 12)
        {
            // Instead revert to the more generic 32-bit assembler version
            return chroma_upscale_asm_32(input);
        }

        const __m128i three = _mm_set1_epi16(3);
        const __m128i zero = _mm_setzero_si128();
        const __m128i masks[5] = {
            _mm_set_epi32(0, 0, 0, 0),   // write nothing
            _mm_set_epi32(0, 0, 0, -1),   // writes 2 x uint16_t
            _mm_set_epi32(0, 0, -1, -1),   // writes 4 x uint16_t
            _mm_set_epi32(0, -1, -1, -1),   // writes 6 x uint16_t
            _mm_set_epi32(-1, -1, -1, -1)    // writes 8 x uint16_t
        };

        // make_borders should generate a contiguous buffer (at least skip_x = 1)
        // We add padding to the right border so that the inner kernel can always consume
        // integral number of items (with the last vector containing 1-7 extra elements)
        auto output = image<uint16_t>(size * 2);
        auto right_pad = (-(size._x + 2)) & 7;      //  6 5 4 3 2 1 0 7  when size % 8 == 0...7
        auto temp = input.make_borders(1, 1, 1, 1 + right_pad, REPLICATE);

        // We minimize the operations (additions + shifting)
        // by summing up the incoming 3 rows to two rows -->
        // top_row = highest_input_row + 3 * middle_input_row
        // bot_row = lowest_input_row + 3 * middle_input_row
        auto stride_y = temp._skip_y;
        for (int j = 0; j < size._y; j++)
        {
            auto *ptr = &temp(j, 0);     // beginning of the row
            auto *dst1 = reinterpret_cast<__m128i*>(&output(j * 2, 0));
            auto *dst2 = reinterpret_cast<__m128i*>(&output(j * 2 + 1, 0));
            __m128i top_prev = zero;     // result for top row 0..7
            __m128i bot_prev = zero;     // result for bot row 0..7
            __m128i top_prev2 = zero;    // result for top row 8..15
            __m128i bot_prev2 = zero;    // result for bot row 8..15
            __m128i top_next = zero;
            __m128i bot_next = zero;
            sum_up_two_rows(top_next, bot_next, ptr, stride_y, three);
            // In the (unlikely) case that input width <= 6, we don't have the next block
            auto width = size._x;
            if (width > 6)
                ptr += 8;

            while (width > 0)
            {
                top_prev = top_next;
                bot_prev = bot_next;
                sum_up_two_rows(top_next, bot_next, ptr, stride_y, three);
                ptr += 8;
                width -= 8;
                // top_prev|top_next  -- we need only the first 10 items of these 2 vectors
                // bot_prev|bot_next
                // results are:  top_left = top_prev + 3 * items1..8(top row)
                //               bot_left = bot_prev + 3 * items1..8(bottom row)
                // top_right = 3 * items1..8(top row) + items2..9(top row)
                // bot_right = 3 * items1..8(bottom row) + items2..9(bottom row)
                // ( and all of those have to be shifted and rounded
                //   - strategies are:  shift by 3, then avg_epu16(x, zero_register)
                //                      add 8, then shift by 4
                //                      single instruction mulhrs by 1/16 -- requires bitdepth <= 11
                //                      mulhi_epu16 by 0x4000
                auto top_mid = _mm_mullo_epi16(_mm_alignr_epi8(top_next, top_prev, 2), three);
                auto bot_mid = _mm_mullo_epi16(_mm_alignr_epi8(bot_next, bot_prev, 2), three);
                auto top_right = _mm_alignr_epi8(top_next, top_prev, 4);
                auto bot_right = _mm_alignr_epi8(bot_next, bot_prev, 4);
                top_prev = _mm_add_epi16(top_prev, top_mid);
                bot_prev = _mm_add_epi16(bot_prev, bot_mid);
                top_right = _mm_add_epi16(top_right, top_mid);
                bot_right = _mm_add_epi16(bot_right, bot_mid);
                // Shift and round -- choosing the shift + avg strategy
                top_prev = _mm_srli_epi16(top_prev, 3);
                bot_prev = _mm_srli_epi16(bot_prev, 3);
                top_right = _mm_srli_epi16(top_right, 3);
                bot_right = _mm_srli_epi16(bot_right, 3);
                top_prev = _mm_avg_epu16(top_prev, zero);
                bot_prev = _mm_avg_epu16(bot_prev, zero);
                top_right = _mm_avg_epu16(top_right, zero);
                bot_right = _mm_avg_epu16(bot_right, zero);
                // Finally interleave the even/odd to get items 0..7 and 8..15 on both rows
                top_prev2 = _mm_unpackhi_epi16(top_prev, top_right);
                top_prev = _mm_unpacklo_epi16(top_prev, top_right);
                bot_prev2 = _mm_unpackhi_epi16(bot_prev, bot_right);
                bot_prev = _mm_unpacklo_epi16(bot_prev, bot_right);

                if (width >= 0)
                {
                    _mm_storeu_si128(dst1, top_prev);
                    _mm_storeu_si128(dst1 + 1, top_prev2);
                    _mm_storeu_si128(dst2, bot_prev);
                    _mm_storeu_si128(dst2 + 1, bot_prev2);
                    dst1 += 2;
                    dst2 += 2;
                }
            }
            if (width < 0)
            {
                auto mask = masks[width & 3];
                if ((width + 8) >= 4)
                {
                    // Write the first 4 items without mask
                    _mm_storeu_si128(dst1, top_prev);
                    _mm_storeu_si128(dst2, bot_prev);
                    _mm_maskmoveu_si128(top_prev2, mask, reinterpret_cast<char*>(dst1 + 1));
                    _mm_maskmoveu_si128(bot_prev2, mask, reinterpret_cast<char*>(dst2 + 1));
                }
                else
                {
                    _mm_maskmoveu_si128(top_prev, mask, reinterpret_cast<char*>(dst1));
                    _mm_maskmoveu_si128(bot_prev, mask, reinterpret_cast<char*>(dst2));
                }
            }
        }
        return output;
    }

    // Calculates the matrix multiplication and clipping for yuv->rgb conversion
    // for `count * 8` number of successive items of 'uint16_t's
    //   r = Clamp(0, 1<<bits - 1, (y * m00 + (u-bias) * m01 + (v-bias) * m02 + 8192) >> 14)
    //   g = Clamp(0, 1<<bits - 1, (y * m10 + (u-bias) * m11 + (v-bias) * m12 + 8192) >> 14)
    //   b = Clamp(0, 1<<bits - 1, (y * m20 + (u-bias) * m21 + (v-bias) * m22 + 8192) >> 14)
    // coeff[0] = m00 1         m00..m22 have 'shift' number of fractional bits
    // coeff[1] = m01 m02
    // coeff[2] = m10 1
    // coeff[3] = m11 m12
    // coeff[4] = m20 1
    // coeff[5] = m21 m22
    // coeff[6] = 1023 x 8  == max_value
    // coeff[7] = 1 << 13   == 0.5  in s1q14 fixed point format
    // coeff[8] = 512 x 8   ==  bias
    template <int shift>
    inline void yuv_to_rgb(__m128i *yr, __m128i *ug, __m128i *vb, __m128i *coeff, unsigned int count)
    {
        // Interleave u,v -- this allows later uv uv uv uv tuples to be multiply-accumulated without precision loss
        // the "madd" aka "multiply accumulate" works as u*a + v*b, where u,a,v,b are all signed int16
        // The multiplying of y*c is coupled with multiplying two constants, 1 and 'half'
        //  - the equation (u*a + v*b + y*c + 1 * half) is carried with two instruction

        while (count-- > 0)
        {
            // load and re-arrange data for matrix multiplication
            auto y_r = _mm_loadu_si128(yr);
            auto u_g = _mm_sub_epi16(_mm_loadu_si128(ug), coeff[8]);    // subtract bias
            auto v_b = _mm_sub_epi16(_mm_loadu_si128(vb), coeff[8]);    // subtract [same] bias
            auto y_lo = _mm_unpacklo_epi16(y_r, coeff[7]);      // Y0 0.5 Y1 0.5 Y2 0.5 Y3 0.5
            auto y_hi = _mm_unpackhi_epi16(y_r, coeff[7]);      // Y4 0.5 Y5 0.5 Y6 0.5 Y7 0.5
            auto uv_lo = _mm_unpacklo_epi16(u_g, v_b);          // U0 V0 U1 V1 U2 V2 U3 V3
            auto uv_hi = _mm_unpackhi_epi16(u_g, v_b);          // U4 V4 U5 V5 U6 V6 U7 V7
            // multiply
            auto r0 = _mm_madd_epi16(y_lo, coeff[0]);           // Y0 * m00 + 0.5 ... Y3 * m00 + 0.5
            auto r1 = _mm_madd_epi16(y_hi, coeff[0]);           // Y4 * m00 + 0.5 ... Y7 * m00 + 0.5
            auto r01 = _mm_madd_epi16(uv_lo, coeff[1]);         // U0 * m01 + V0 * m02 ... U3 * m01 + V3 * m02
            auto r11 = _mm_madd_epi16(uv_hi, coeff[1]);         // U4 * m01 + V4 * m02 ... U7 * m01 + V7 * m02
            auto g0 = _mm_madd_epi16(y_lo, coeff[2]);           // Y0..Y3 * m10 + half
            auto g1 = _mm_madd_epi16(y_hi, coeff[2]);           // Y4..Y7 * m10 + half
            auto g01 = _mm_madd_epi16(uv_lo, coeff[3]);         // U0..3 * m11 + V0..3 * m12
            auto g11 = _mm_madd_epi16(uv_hi, coeff[3]);         // U4..7 * m11 + V4..7 * m12
            auto b0 = _mm_madd_epi16(y_lo, coeff[4]);           // Y0..Y3 * m20 + 0.5
            auto b1 = _mm_madd_epi16(y_hi, coeff[4]);           // Y4..Y7 * m20 + 0.5
            auto b01 = _mm_madd_epi16(uv_lo, coeff[5]);         // U0..3 * m21 + V0..3 * m22
            auto b11 = _mm_madd_epi16(uv_hi, coeff[5]);         // U4..7 * m21 + V4..7 * m22
            // sum up the products and the rounding term, then shift
            r0 = _mm_add_epi32(r0, r01);                        // sum of red components R0..3
            r1 = _mm_add_epi32(r1, r11);                        // sum of red components R4..7
            r0 = _mm_srai_epi32(r0, shift);                     // Red0..3 shifted
            r1 = _mm_srai_epi32(r1, shift);                     // Red4..7 shifted
            g0 = _mm_add_epi32(g0, g01);                        // Greens 0..7 summed
            g1 = _mm_add_epi32(g1, g11);                        //
            g0 = _mm_srai_epi32(g0, shift);                     // Greens 0..7 shifted
            g1 = _mm_srai_epi32(g1, shift);                     //
            b0 = _mm_add_epi32(b0, b01);                        // Blues 0..7 summed
            b1 = _mm_add_epi32(b1, b11);                        //
            b0 = _mm_srai_epi32(b0, shift);                     // Blues 0..7 shifted
            b1 = _mm_srai_epi32(b1, shift);
            // Convert intermediate 32-bit integers to final representation while clamping to 0..B
            // Saturating conversion from int -> ushort will clip negative values as a side effect
            r0 = _mm_packus_epi32(r0, r1);                      // Combines 8 x int32 -> 8 x uint16_t
            g0 = _mm_packus_epi32(g0, g1);                      // same here for green
            b0 = _mm_packus_epi32(b0, b1);                      // and red
            // Then we clip against the max value and store
            _mm_storeu_si128(yr++, _mm_min_epu16(r0, coeff[6]));
            _mm_storeu_si128(ug++, _mm_min_epu16(g0, coeff[6]));
            _mm_storeu_si128(vb++, _mm_min_epu16(b0, coeff[6]));
        }
    }

    /// Given coefficients with `shift` amount of fractional bits
    /// convert y,u,v inplace to r,g,b  clipped to `bits` number of bits
    template <int shift>
    inline void yuv_to_rgb(uint16_t &y_r, uint16_t &u_g, uint16_t &v_b, const int *coeffs, int bits)
    {
        const int half = 1 << (shift - 1);
        const int bias = 1 << (bits - 1);

        auto y = (int)y_r;
        auto u = (int)u_g - bias;
        auto v = (int)v_b - bias;

        auto r = (y * coeffs[0] + u * coeffs[1] + v * coeffs[2] + half) >> shift;
        auto g = (y * coeffs[3] + u * coeffs[4] + v * coeffs[5] + half) >> shift;
        auto b = (y * coeffs[6] + u * coeffs[7] + v * coeffs[8] + half) >> shift;

        y_r = static_cast<uint16_t>(std::max(0, std::min(r, bias * 2 - 1)));
        u_g = static_cast<uint16_t>(std::max(0, std::min(g, bias * 2 - 1)));
        v_b = static_cast<uint16_t>(std::max(0, std::min(b, bias * 2 - 1)));
    };

    /// Quantizes double precision matrix to given integral format
    template <typename T, int bits>
    typename std::enable_if<std::is_integral<T>::value, std::vector<T>>::type
    inline to_fixed_point(const double(&matrix)[9])
    {
        std::vector<T> coeffs;
        std::transform(matrix, matrix + 9, std::back_inserter(coeffs),
            [](double a) -> T { return static_cast<T>(round(a * (1 << bits))); });
        return coeffs;
    }

    /// Inplace color conversion -- YUV to RGB
    /// The output data is clipped to range 0..2^bits - 1
    /// This should work well for 8,10,12 bits of input
    inline void yuv_to_rgb_interleaved(image<uint16_t> &input, int bits, const double(&matrix)[9])
    {
        if (input._width % 3 != 0)
            throw std::runtime_error("Image is not in interleaved RGB or YUV format");

        auto coeffs = to_fixed_point<int, 14>(matrix);
        auto *data = coeffs.data();

        int skipx = input._skip_x;
        input.subview(1, 3).foreach([skipx, bits, data](uint16_t &dst)
        {
            uint16_t *ptr = &dst;
            yuv_to_rgb<14>(ptr[0], ptr[skipx], ptr[skipx * 2], data, bits);
        });
    }

    /// Inplace color conversion -- YUV to RGB in three planes
    /// The output data is clipped to range 0..2^bits - 1
    /// This should work well for 8,10,12 bits of input
    inline void yuv_to_rgb_planar(
        image<uint16_t> &y_r,
        image<uint16_t> &u_g,
        image<uint16_t> &v_b,
        int bits, const double(&matrix)[9])
    {
        if (!(y_r.size() == u_g.size() && y_r.size() == v_b.size()))
            throw std::runtime_error("Y, U and V planes mismatch in dimensions");

        auto coeffs = to_fixed_point<int, 14>(matrix);
        auto *data = coeffs.data();

        y_r.foreach([data, bits](uint16_t &y_r, uint16_t &u_g, uint16_t &v_b)
        {
            yuv_to_rgb<14>(y_r, u_g, v_b, data, bits);
        }, u_g, v_b);
    }

    /// Inplace color conversion -- YUV to RGB in three planes
    /// The output data is clipped to range 0..2^bits - 1
    /// This should work well for 8,10,12 bits of input
    inline void yuv_to_rgb_planar_asm(
        image<uint16_t> &y_r,
        image<uint16_t> &u_g,
        image<uint16_t> &v_b,
        int bits, const double(&matrix)[9])
    {
        if (!(y_r.size() == u_g.size() && y_r.size() == v_b.size()))
            throw std::runtime_error("Y, U and V planes mismatch in dimensions");

        auto coeffs = to_fixed_point<int, 14>(matrix);
        auto *data = coeffs.data();

        if (y_r.is_contiguous() && u_g.is_contiguous() && v_b.is_contiguous())
        {
            __m128i coeff[9] =
            {
                _mm_set1_epi32((coeffs[0] & 0xffff) | (1 << 16)),           // coeff[0] = m00 1
                _mm_set1_epi32((coeffs[1] & 0xffff) | (coeffs[2] << 16)),   // coeff[1] = m01 m02
                _mm_set1_epi32((coeffs[3] & 0xffff) | (1 << 16)),           // coeff[0] = m10 1
                _mm_set1_epi32((coeffs[4] & 0xffff) | (coeffs[5] << 16)),   // coeff[1] = m11 m12
                _mm_set1_epi32((coeffs[6] & 0xffff) | (1 << 16)),           // coeff[0] = m20 1
                _mm_set1_epi32((coeffs[7] & 0xffff) | (coeffs[8] << 16)),   // coeff[1] = m21 m22
                _mm_set1_epi16((1 << bits) - 1),                            // coeff[6] = 1023 x 8  == max_value
                _mm_set1_epi16(1 << (14 - 1)),                              // coeff[7] = 1 << 13   == 0.5  in s1q14 fixed point format
                _mm_set1_epi16(1 << (bits - 1))                             // coeff[8] = 512 x 8   ==  bias
            };

            // calculate as many elements as possible using assembler
            auto size = u_g._width * u_g._height;
            yuv_to_rgb<14>(
                reinterpret_cast<__m128i*>(y_r._begin),
                reinterpret_cast<__m128i*>(u_g._begin),
                reinterpret_cast<__m128i*>(v_b._begin),
                coeff, size / 8);

            // and then iterate the rest of the items using the reference method
            for (auto x = size & ~7; x < size; x++)
                yuv_to_rgb<14>(y_r._begin[x], u_g._begin[x], v_b._begin[x], data, bits);
        }
        else
        {
            // Iterate using the regular method
            y_r.foreach([data, bits](uint16_t &y_r, uint16_t &u_g, uint16_t &v_b)
            {
                yuv_to_rgb<14>(y_r, u_g, v_b, data, bits);
            }, u_g, v_b);
        }
    }
};
