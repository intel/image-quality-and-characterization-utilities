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
#include "Teisko/Image/API.hpp"
#include "Teisko/BayerInfo.hpp"

namespace Teisko
{
    // These forward declarations are needed on GCC but not on MSVC++.
    // The operator [](uint32_t i) in the bayer_image_s class depends
    // and GCC is more stricter about this than MSVC++.
    template <typename T>
    image<T> get_sub_channel(uint32_t ch, image<T> &img, bayer_info_s &info);

    template <typename T>
    image<T> get_sub_channel(bayer_info_s::color_info_e color, image<T> &img,
        bayer_info_s &info);

    /// Container for image data and it's properties as used by preprocessing
    /// Q: class or struct? Should we allow bayer_image_s to be constructed
    ///    without valid bayer pattern? -- if not, there's an invariant
    ///    and this should be a class. See c++ CoreGuideLines C.2
    template <typename T>
    class bayer_image_s
    {
    public:
        // functions allowing bayer_image_s<T> to be range for looped
        // by channels:  -- use as `for (auto i: image)`
        channel_iterator begin() {
            return channel_iterator{ 0 };
        }
        channel_iterator end() {
            return channel_iterator{ _layout.get_channels() };
        }

        // Returns Ith sub-channels of image as view
        image<T> operator [](uint32_t i)
        {
            return get_sub_channel(i, _img, _layout);
        }

        /// Accessing by color is mostly usable for rgbi sensors
        /// due to API limitation of advancing the index to acquire next channel of kind
        image<T> operator [](color_info_e color)
        {
            return get_sub_channel(color, _img, _layout);
        }

        bool matches_by_size_and_type(bayer_image_s<T> &other)
        {
            return (_img._height == other._img._height) &
                (_img._width == other._img._width) &
                (_layout == other._layout);
        }

        // Creates a new self contained image
        bayer_image_s(int height, int width, bayer_pattern_e color_info)
            : _img(height, width), _layout(color_info)
        { }

        // Creates a new view to image
        bayer_image_s(int height, int width, T *ptr, bayer_pattern_e color_info)
            : _img(height, width, ptr), _layout(color_info)
        { }

        // Creates a new view from other view
        bayer_image_s(image<T> &view, bayer_pattern_e color_info)
            : _img(view), _layout(color_info)
        { }

        // Creates a new view from other view
        template <typename image_type>
        bayer_image_s(image_type &&view, bayer_pattern_e color_info)
            : _img(std::forward(view)), _layout(color_info)
        { }

        // Copy constructor takes deep copy of same type
        // we have a templated version, but that is not enough
        // this non-templated version is a must
        bayer_image_s(const bayer_image_s &other)
            : _img(other._img.convert_to()), _layout(other._layout)
        {
        }

        // Copy constructor takes deep copy of any type with conversion
        template <typename U>
        bayer_image_s(const bayer_image_s<U> &other)
            : _img(other._img.template convert_to<T>()), _layout(other._layout)
        {
        }

        // Assignment operator deep copies data (with packing)
        bayer_image_s& operator=(const bayer_image_s &other)
        {
            _img = other._img.convert_to();
            _layout = other._layout;
            return *this;
        }

        // Assignment operator deep copies data with conversion
        template <typename U>
        bayer_image_s& operator=(const bayer_image_s<U> &other)
        {
            _img = other._img.template convert_to<T>();
            _layout = other._layout;
            return *this;
        }

        image<T> _img;
        bayer_info_s _layout;
    };

    /// Helper function to get Nth Bayer channel of image as view
    template <typename T>
    image<T> get_sub_channel(uint32_t ch, image<T> &img, bayer_info_s &info)
    {
        uint32_t h_chans = info.get_height();
        uint32_t w_chans = info.get_width();
        if (h_chans == 0 || w_chans == 0)
            throw std::runtime_error("Division by Bayer sensor width = 0");

        return img.subview(h_chans, w_chans, ch / w_chans, ch % w_chans);
    }

    /// Helper function to get first color channel by enum as view
    template <typename T>
    image<T> get_sub_channel(bayer_info_s::color_info_e color, image<T> &img, bayer_info_s &info)
    {
        int idx = info.locate_color(color);
        if (idx < 0)
            throw std::runtime_error("Invalid color channel queried");
        return get_sub_channel(idx, img, info);
    }

    // Returns image dimension in pixels (e.g. 4208 x 3120)
    template <typename T>
    inline roi_point get_dim(bayer_image_s<T> &img) { return { (int)img._img._width, (int)img._img._height }; }

    // returns image dimension of a single color channel of raw image
    template <typename T>
    inline roi_point get_channel_dim(bayer_image_s<T> &img) {
        auto w = img._layout.get_width();
        auto h = img._layout.get_height();
        if (w == 0 || h == 0)
            throw std::runtime_error("Division by Bayer pattern size 0");
        return{ (int)(img._img._width / w), (int)(img._img._height / h) };
    }

    // returns image dimension in pixels (e.g. 4208 x 3120)
    template <typename T>
    inline roi_point get_dim(image<T> &img) { return { (int)img._width, (int)img._height }; }

    // return channel dimensions e.g. 2x2, 4x2, 4x4
    inline roi_point get_dim(bayer_info_s &layout) {
        return { (int)layout.get_width(), (int)layout.get_height() };
    }

}

