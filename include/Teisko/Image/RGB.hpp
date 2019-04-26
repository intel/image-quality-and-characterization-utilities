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
#include <map>
namespace Teisko
{
    /// <summary>
    /// Layouts for RGB images
    /// </summary>
    enum class rgb_layout_e
    {
        paged, ///!< aka Planar, all channels concatenated
        tiled, ///!< all channels concatenated horizontally
        interleaved, ///!< R0G0B0,R1G1B1, etc
        interleaved_stride_4 ///!< Each row padded to multiples of 4 bytes as used in [X]-Windows
    };

    /// <summary>
    /// Color orders of RGB images
    /// </summary>
    enum class rgb_order_e { rgb, rbg, grb, gbr, brg, bgr };

    /// <summary>
    /// Color codes for RGB images
    /// </summary>
    enum class rgb_color_e { red, green, blue };

    /// <summary>
    /// Container class for RGB images
    /// </summary>
    template <typename T>
    class rgb_image_s
    {
    public:
        /// <summary>
        /// Constructs a self contained RGB image
        /// </summary>
        /// <param name="size">Dimensions of image</param>
        /// <param name="layout">Internal layout of the image, defaults to paged/planar</param>
        /// <param name="order">Internal color order of the image, defaults to RGB</param>
        rgb_image_s(roi_point size,
            rgb_layout_e layout = rgb_layout_e::paged,
            rgb_order_e order = rgb_order_e::rgb)
            : _size(size)
            , _layout(layout)
            , _pattern(get_color_pattern(order))
            , _img(layout_to_size(_layout, _size))
        {
            _img.fill(static_cast<T>(0));
        }

        /// <summary>Constructs a view to externally allocated RGB image</summary>
        /// <param name="size">Dimensions of image</param>
        /// <param name="ptr">Pointer to external image raw data</param>
        /// <param name="layout">Internal layout of the image, defaults to paged/planar</param>
        /// <param name="order">Internal color order of the image, defaults to RGB</param>
        rgb_image_s(roi_point size, T *ptr,
            rgb_layout_e layout = rgb_layout_e::paged,
            rgb_order_e order = rgb_order_e::rgb)
            : _size(size)
            , _layout(layout)
            , _pattern(get_color_pattern(order))
            , _img(layout_to_size(_layout, _size), ptr)
        {
        }

        /// <summary>Retrieve a view to internal container</summary>
        /// <returns>image view containing all pixel data</returns>
        image<T> full_image() const { return _img; }

        /// <summary>Returns the rgb color order as vector of color codes</summary>
        /// <param name="order">Enumeration for the color order</param>
        /// <returns>Vector of color codes</returns>
        static std::vector<rgb_color_e> get_color_pattern(rgb_order_e order)
        {
            switch (order)
            {
            case rgb_order_e::rgb: return{ rgb_color_e::red, rgb_color_e::green, rgb_color_e::blue };
            case rgb_order_e::rbg: return{ rgb_color_e::red, rgb_color_e::blue, rgb_color_e::green };
            case rgb_order_e::grb: return{ rgb_color_e::green, rgb_color_e::red, rgb_color_e::blue };
            case rgb_order_e::gbr: return{ rgb_color_e::green, rgb_color_e::blue, rgb_color_e::red };
            case rgb_order_e::brg: return{ rgb_color_e::blue, rgb_color_e::red, rgb_color_e::green };
            case rgb_order_e::bgr: return{ rgb_color_e::blue, rgb_color_e::green, rgb_color_e::red };
            default:
                throw std::runtime_error("nonexistent rgb color order");
            }
        }

        /// <summary>Returns a view to a color channel</summary>
        /// <param name="color">Color code to get</param>
        /// <returns>image view to channel or throws on error</returns>
        image<T> operator [](rgb_color_e color)
        {
            auto it = std::find(_pattern.begin(), _pattern.end(), color);
            if (it == _pattern.end())
                throw std::runtime_error("Asking for unknown color");

            return operator [](static_cast<int>(std::distance(_pattern.begin(), it)));
        }

        /// <summary>Returns a view to a color channel by index</summary>
        /// <param name="i">index to channel</param>
        /// <returns>image view to channel or throws on error</returns>
        image<T> operator [](int i)
        {
            if (i < 0 || i > 2)
                throw std::runtime_error("Max three channels in rgb image");

            switch (_layout)
            {
            case rgb_layout_e::interleaved:
                return _img.subview(1, 3, 0, i);
            case rgb_layout_e::interleaved_stride_4:
                return _img.subview(1, 3, 0, i).region(_size);
            case rgb_layout_e::tiled:
                return _img.region(_size._y, _size._x, 0, _size._x * i);
            case rgb_layout_e::paged:
                return _img.region(_size._y, _size._x, _size._y * i, 0);
            default:
                throw std::runtime_error("Illegal layout");
            }
        }

    private:
        roi_point _size{};        ///<! Dimension of a single channel
        rgb_layout_e _layout;     ///<! Physical layout of image
        std::vector<rgb_color_e>  _pattern; ///<! Color codes of individual channels
        image<T> _img;            ///<! Container holding all channels

        /// <summary>Calculates container size for given layout</summary>
        /// <param name="layout">Internal layout of RGB image</param>
        /// <param name="size">Dimensions of (single color) in RGB image</param>
        /// <returns>Logical container size for given layout</returns>
        static roi_point layout_to_size(rgb_layout_e layout, roi_point size)
        {
            switch (layout)
            {
            case rgb_layout_e::interleaved:
            case rgb_layout_e::tiled:
                return roi_point(size._x * 3, size._y);
            case rgb_layout_e::paged:
                return roi_point(size._x, size._y * 3);
            case rgb_layout_e::interleaved_stride_4:
                return roi_point((size._x * 3 + 3) & ~3, size._y);
            default:
                throw std::runtime_error("Illegal layout");
            }
        }
    };
}
