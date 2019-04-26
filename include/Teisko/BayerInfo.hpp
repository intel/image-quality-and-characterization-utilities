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
#include "Teisko/Algorithm/Iterators.hpp"
#include <cstdint>
#include <vector>
#include <map>
#include <array>
#include <string>
#include <algorithm>  // std::find

namespace Teisko
{
    struct bayer_info_s
    {
        // Copy of liblab_sensor_layout_e
        enum class bayer_pattern_e
        {
            rggb = 0,
            grbg,
            gbrg,
            bggr,
            grbi,
            irbg,
            rgib,
            rigb,
            bgir,
            bigr,
            gbri,
            ibrg,
            grbg_4x2 = 128,             // 1st row: ggrr, 2nd row: bbgg (4x2)
            rggb_4x2,                   // 1st row: rrgg, 2nd row: ggbb (4x2)
            bggr_4x2,                   // 1st row: bbgg, 2nd row: ggrr (4x2)
            gbrg_4x2,                   // 1st row: ggbb, 2nd row: rrgg (4x2)
            bgrg_gigi_rgbg_gigi = 256,  // 1st row: bgrg, 2nd row: gigi, 3rd row: rgbg, 4th row: gigi (4x4)
            grgb_igig_gbgr_igig,        // 256 shifted horizontally by one column
            rgbg_gigi_bgrg_gigi,
            gbgr_igig_grgb_igig,
            gigi_rgbg_gigi_bgrg,
            igig_gbgr_igig_grgb,
            gigi_bgrg_gigi_rgbg,        // imx488
            igig_grgb_igig_gbgr
        };

        static inline std::vector<bayer_pattern_e> get_regular_2x2_patterns()
        {
            return{
                bayer_pattern_e::rggb, bayer_pattern_e::bggr,
                bayer_pattern_e::gbrg, bayer_pattern_e::grbg
            };
        }

        static inline std::vector<bayer_pattern_e> get_ir_2x2_patterns()
        {
            return{ bayer_pattern_e::grbi, bayer_pattern_e::irbg, bayer_pattern_e::rgib,
                bayer_pattern_e::rigb, bayer_pattern_e::bgir, bayer_pattern_e::bigr,
                bayer_pattern_e::gbri, bayer_pattern_e::ibrg };
        }

        static inline std::vector<bayer_pattern_e> get_dp_4x2_patterns()
        {
            return{
                bayer_pattern_e::grbg_4x2, bayer_pattern_e::rggb_4x2,
                bayer_pattern_e::bggr_4x2, bayer_pattern_e::gbrg_4x2
            };
        }

        static inline std::vector<bayer_pattern_e> get_ir_4x4_patterns()
        {
            return{ bayer_pattern_e::bgrg_gigi_rgbg_gigi, bayer_pattern_e::grgb_igig_gbgr_igig,
                bayer_pattern_e::rgbg_gigi_bgrg_gigi, bayer_pattern_e::gbgr_igig_grgb_igig,
                bayer_pattern_e::gigi_rgbg_gigi_bgrg, bayer_pattern_e::igig_gbgr_igig_grgb,
                bayer_pattern_e::gigi_bgrg_gigi_rgbg, bayer_pattern_e::igig_grgb_igig_gbgr
            };
        }

        enum class color_info_e { red, green, blue, ir };

        inline operator int() { return (int)_layout; }
        inline operator bayer_pattern_e() { return _layout; }

        /// returns the color type of 'i'th channel
        color_info_e operator[](uint32_t index)
        {
            if (index >= _channels)
                throw std::runtime_error("Index out of range int bayer_info_s::operator()(uint32_t)");
            return _pattern[index];
        }

        /// Ranged for necessity -- start of iterator
        channel_iterator begin() { return channel_iterator{ 0 }; }
        /// Ranged for necessity -- end of iterator
        channel_iterator end() { return channel_iterator{ _channels }; }

        /// Returns copy of color patterns
        std::vector<color_info_e> get_color_pattern() { return _pattern; }
        /// Return number of color components in bayer pattern
        uint32_t get_channels() { return _channels; }
        /// Return number of horizontal components in bayer pattern
        uint32_t get_width() { return _width; }
        /// Return number of vertical components in bayer pattern
        uint32_t get_height() { return _height; }

        /// Returns true if the sensor layout is one of rggb, bggr, grbg, gbrg
        bool is_2x2_sensor() { return (int)_layout <= (int)bayer_pattern_e::bggr; }
        /// Returns true if the 2x2 sensor has I component
        bool is_2x2_ir_sensor() { return _is_ir && _channels == 4; }
        /// Returns true if the sensor if one of four dual phase 4x2 sensors
        /// with color components split to left and right half pixels
        bool is_dp_4x2_sensor() { return _channels == 8; }
        /// Return true if the sensor is 4x4 and has I components
        bool is_4x4_ir_sensor() { return _is_ir && _channels == 16; }
        /// Return true for 2x2 and 4x4 sensors with Infra Red pixels
        bool is_ir_sensor() { return _is_ir; }

        /// Locate the first channel by color enumeration
        /// or -1, if not found
        int locate_color(color_info_e color, int offset = 0)
        {
            auto it = std::find(_pattern.begin() + offset, _pattern.end(), color);
            if (it == _pattern.end())
                return -1;
            return (int)std::distance(_pattern.begin(), it);
        }

        bayer_info_s(bayer_pattern_e info = bayer_pattern_e::rggb)
            : _layout(info)
            , _pattern(0)
            , _width(0)
            , _height(0)
            , _channels(0)
            , _is_ir(false)
        {
            const auto r = color_info_e::red;
            const auto g = color_info_e::green;
            const auto b = color_info_e::blue;
            const auto i = color_info_e::ir;
            // First make the pattern
            switch (_layout)
            {
            case bayer_pattern_e::rggb: _pattern = { r, g, g, b }; break;
            case bayer_pattern_e::grbg: _pattern = { g, r, b, g }; break;
            case bayer_pattern_e::gbrg: _pattern = { g, b, r, g }; break;
            case bayer_pattern_e::bggr: _pattern = { b, g, g, r }; break;
            case bayer_pattern_e::grbi: _pattern = { g, r, b, i }; break;
            case bayer_pattern_e::irbg: _pattern = { i, r, b, g }; break;
            case bayer_pattern_e::rgib: _pattern = { r, g, i, b }; break;
            case bayer_pattern_e::rigb: _pattern = { r, i, g, b }; break;
            case bayer_pattern_e::bgir: _pattern = { b, g, i, r }; break;
            case bayer_pattern_e::bigr: _pattern = { b, i, g, r }; break;
            case bayer_pattern_e::gbri: _pattern = { g, b, r, i }; break;
            case bayer_pattern_e::ibrg: _pattern = { i, b, r, g }; break;
            case bayer_pattern_e::grbg_4x2: _pattern = { g, g, r, r, b, b, g, g }; break;
            case bayer_pattern_e::rggb_4x2: _pattern = { r, r, g, g, g, g, b, b }; break;
            case bayer_pattern_e::bggr_4x2: _pattern = { b, b, g, g, g, g, r, r }; break;
            case bayer_pattern_e::gbrg_4x2: _pattern = { g, g, b, b, r, r, g, g }; break;
            case bayer_pattern_e::bgrg_gigi_rgbg_gigi: _pattern = { b, g, r, g, g, i, g, i, r, g, b, g, g, i, g, i }; break;
            case bayer_pattern_e::grgb_igig_gbgr_igig: _pattern = { g, r, g, b, i, g, i, g, g, b, g, r, i, g, i, g }; break;
            case bayer_pattern_e::rgbg_gigi_bgrg_gigi: _pattern = { r, g, b, g, g, i, g, i, b, g, r, g, g, i, g, i }; break;
            case bayer_pattern_e::gbgr_igig_grgb_igig: _pattern = { g, b, g, r, i, g, i, g, g, r, g, b, i, g, i, g }; break;
            case bayer_pattern_e::gigi_rgbg_gigi_bgrg: _pattern = { g, i, g, i, r, g, b, g, g, i, g, i, b, g, r, g }; break;
            case bayer_pattern_e::igig_gbgr_igig_grgb: _pattern = { i, g, i, g, g, b, g, r, i, g, i, g, g, r, g, b }; break;
            case bayer_pattern_e::gigi_bgrg_gigi_rgbg: _pattern = { g, i, g, i, b, g, r, g, g, i, g, i, r, g, b, g }; break;
            case bayer_pattern_e::igig_grgb_igig_gbgr: _pattern = { i, g, i, g, g, r, g, b, i, g, i, g, g, b, g, r }; break;
            default:
                break;
            }

            // Derive width, height and channel count (from pattern)

            switch (_pattern.size())
            {
            case 4:
                _width = 2;
                _height = 2;
                break;
            case 8:
                _width = 4;
                _height = 2;
                break;
            case 16:
                _width = 4;
                _height = 4;
                break;
            }

            _channels = _width * _height;

            // derive IR status from pattern
            for (auto &x : _pattern)
                if (x == i)
                    _is_ir = true;
        }

        private:

        bayer_pattern_e _layout;
        std::vector<color_info_e> _pattern;
        uint32_t _width;
        uint32_t _height;
        uint32_t _channels;
        bool _is_ir;
    };

    using bayer_pattern_e = bayer_info_s::bayer_pattern_e;
    using color_info_e = bayer_info_s::color_info_e;

    inline bool operator ==(bayer_info_s &a, bayer_info_s &b)
    {
        return (bayer_pattern_e)a == (bayer_pattern_e)b;
    }

    template <typename T>
    inline std::vector<T> remap_4x4_vector(std::vector<T> four_by_four_data, bayer_info_s &info)
    {
        if (four_by_four_data.size() != 16)
            throw std::runtime_error("The 4x4 data must have 16 items!");

        switch (info.get_channels())
        {
        case 4:
            return{ four_by_four_data[0], four_by_four_data[1], four_by_four_data[4], four_by_four_data[5] };
        case 8:
            four_by_four_data.resize(8);
            return four_by_four_data;
        case 16:
            return four_by_four_data;
        default:
            throw std::runtime_error("Unexpected sensor layout");
        }
    }
}