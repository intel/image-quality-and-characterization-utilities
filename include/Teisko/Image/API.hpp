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
#include "Teisko/Image/Point.hpp"
#include "Teisko/Image/Support.hpp"

#include <cstdint>
#include <memory>
#include <vector>
#include <fstream>
#include <cmath>
#include <algorithm>

namespace Teisko
{
    /// @brief Encodings for border mirroring/replication schemes
    enum mirror_schemes_e
    {
        MIRROR_EVEN = 0x0000,       ///< Mirror   '3 2 1 | 1 2 3 4 | 4 3 2'
        SYMMETRIC = 0x0000,         ///< Matlab naming convention for even mirroring
        MIRROR_ODD = 0x1111,        ///< Mirror   '4 3 2 | 1 2 3 4 | 3 2 1'
        REPLICATE = 0x2222,         ///< Produces '1 1 1 | 1 2 3 4 | 4 4 4'
        ZERO_PAD = 0x3333,          ///< Produces '0 0 0 | 1 2 3 4 | 0 0 0'
        REPLICATE_EVEN = 0x4444     ///< Produces '1 2 3 | 1 2 3 4 | 2 3 4'
    };

    /// @brief Finer control of mirror schemes for each border separately
    ///
    /// @details
    /// e.g. ((LEFT|RIGHT) & MIRROR_ODD) | ((TOP|BOTTOM) & MIRROR_EVEN)
    /// ==     LEFT_RIGHT & MIRROR_ODD | TOP_BOTTOM & MIRROR_EVEN
    /// ==   LEFT & MIRROR_ODD | RIGHT & MIRROR_ODD | TOP & MIRROR_EVEN | BOTTOM & MIRROR_EVEN
    enum border_e
    {
        LEFT = 0x000F,
        RIGHT = 0x00F0,
        LEFT_RIGHT = LEFT | RIGHT,
        TOP = 0x0F00,
        TOP_LEFT = TOP | LEFT,
        BOTTOM = 0xF000,
        BOTTOM_RIGHT = RIGHT | BOTTOM,
        TOP_BOTTOM = TOP | BOTTOM
    };

    /// @brief Container for 2-dimensional image data with built in iterators
    template <typename pixel_type> struct image
    {
        range       _width{};                   ///< Width of the current view in pixels
        range       _height{};                  ///< Height of the current view in pixels
        int32_t     _skip_x{};                  ///< Number of elements between horizontal neighbors
        int32_t     _skip_y{};                  ///< Number of elements between vertical neighbors
        std::shared_ptr <pixel_type> _owned{};  ///< Pointer to owned data
        pixel_type *_begin{};                   ///< Pointer to logical element (0,0)

        /// @brief Default Image with no data
        image() = default;

        /// @brief Creates a new image with self managed uninitialized data
        /// @param height   Height of the image in pixels
        /// @param width    Width of the image in pixels
        image(int32_t height, int32_t width)
            : _width(disallow_negative(width))
            , _height(disallow_negative(height))
            , _skip_x(1)
            , _skip_y(width)
            , _owned(new pixel_type[height * width], [](pixel_type *p) { delete[] p; })
            , _begin((pixel_type*)(_owned.get()))
        { }

        /// @brief Creates a new image with self managed uninitialized data
        /// @param size     Dimensions of the image in pixels
        explicit image(roi_point size)
            : _width(disallow_negative(size._x))
            , _height(disallow_negative(size._y))
            , _skip_x(1)
            , _skip_y(size._x)
            , _owned(new pixel_type[size._x * size._y], [](pixel_type *p) { delete[] p; })
            , _begin((pixel_type*)(_owned.get()))
        { }

        /// @brief Creates a new image view to externally managed contiguous data
        /// @param height   Height of the image in pixels
        /// @param width    Width of the image in pixels
        /// @param begin    Pointer to raw data
        image(int32_t height, int32_t width, pixel_type *begin)
            : _width(disallow_negative(width))
            , _height(disallow_negative(height))
            , _skip_x(1)
            , _skip_y(width)
            , _owned(nullptr)
            , _begin(begin)
        { }

        /// @brief Creates a new image view to externally managed contiguous data
        /// @param size     Dimensions of the image in pixels
        /// @param begin    Pointer to raw data
        image(roi_point size, pixel_type *begin)
            : _width(disallow_negative(size._x))
            , _height(disallow_negative(size._y))
            , _skip_x(1)
            , _skip_y(size._x)
            , _owned(nullptr)
            , _begin(begin)
        { }

        /// @brief Validated input parameter for non-negative or throws
        /// @param value    Integer value to validate
        /// @returns        Input value if legal, throws otherwise
        static int32_t disallow_negative(int32_t value)
        {
            if (value < 0)
                throw std::runtime_error("Negative dimension not allowed");
            return value;
        }

        /// @brief  Returns an xy pixel iterator
        xy_skip_iterator_s<pixel_type> begin() const {
            return{ _begin, { _skip_x, _skip_y }, size() };
        }

        /// @brief  Returns an xy pixel iterator
        xy_skip_iterator_s<pixel_type> end() const {
            auto skip = roi_point{ _skip_x, _skip_y };
            auto sz = size();
            return{ &at(sz._y, 0), skip, sz };
        }

        /// @brief Retrieve image dimensions
        /// @returns Image size as roi_point
        roi_point size() const { return roi_point(_width, _height); }

        /// @brief Reads image from file
        /// @param filename     Full path to file
        /// @returns            Reference to the self
        const image& read(std::string filename) const
        {
            std::ifstream input;
            input.open(filename, std::ios::binary);
            input.seekg(0, std::ios::end);
            std::streampos file_len = input.tellg();
            input.seekg(0, std::ios::beg);
            auto bytecount = _width * _height * sizeof(pixel_type);
            if ((size_t)file_len == (size_t)bytecount)
            {
                if (is_contiguous())
                {
                    input.read(reinterpret_cast<char*>(_begin), file_len);
                }
                else
                {
                    std::vector<pixel_type> data(_width * _height);
                    input.read(reinterpret_cast<char*>(data.data()), file_len);
                    init_from(image<pixel_type>(size(), data.data()));
                }
            }
            return *this;
        }

        /// @brief Writes image to file as contiguous binary representation without header information
        /// @param filename     Full path to file
        /// @returns            Reference to self
        const image& write(std::string filename) const
        {
            std::ofstream output;
            output.open(filename, std::ios::binary);
            auto bytecount = _width * _height * sizeof(pixel_type);
            if (is_contiguous())
            {
                output.write(reinterpret_cast<char*>(_begin), bytecount);
            }
            else
            {
                auto vec = to_vector();
                output.write(reinterpret_cast<char*>(vec.data()), bytecount);
            }
            return *this;
        }

        /// @brief Allows implementation to check if internal image representation is
        /// contiguous to allow optimizations that rely on this fact
        /// @returns        true when image data does not have gaps
        bool is_contiguous() const
        {
            return (_skip_x == 1) && (_skip_y == _width);
        }

        /// @brief Initializes the image with constant value
        /// @param value    Pixel value to fill with
        /// @returns        Reference to the self
        const image& fill(pixel_type value) const
        {
            foreach([value](pixel_type &a) { a = value; });
            return *this;
        }

        /// @brief Initializes the image with constant value -- to be deprecated
        /// @param value    Pixel value to fill with
        /// @returns        Reference to the self
        const image& init(pixel_type value) const { return fill(value); }

        /// @brief Initializes the image by applying a function of type `void (pixel_type &)` for each pixel
        /// @param func     Functor, lambda or function to apply to all pixels
        /// @returns        Reference to the self
        template <typename UnaryOperator>
        const image& generate(UnaryOperator func) const
        {
            foreach(func);
            return *this;
        }

        /// @brief Fills the image from other image having at least the same size as self.
        /// This function allows narrowing conversions e.g from double to float
        /// @param other    Other image to copy the values
        /// @returns        Reference to the self or throws if source image is too small
        template <typename src_type>
        const image& init_from(image<src_type> &other) const
        {
            foreach([](pixel_type &dst, src_type &src) { dst = (pixel_type)src; }, other);
            return *this;
        }

        /// @brief  Fills the image from other image having at least the same size as self.
        /// This function allows narrowing conversions e.g from double to float
        /// @param other    Other image to copy the values
        /// @returns        Reference to the self or throws if source image is too small
        template <typename src_type>
        const image& init_from(image<src_type> &&other) const
        {
            return init_from(other);
        }

        /// @brief Fills the image from raw pointer
        /// This function allows narrowing conversions e.g from double to float
        /// @param ptr      Pointer to raw data to initialize from
        /// @returns        Reference to the self
        template <typename src_type>
        const image& init_from(src_type *ptr) const
        {
            if (ptr == nullptr)
                throw std::runtime_error("Null pointer not allowed");
            return init_from(image<src_type>(size(), ptr));
        }

        /// @brief Populates other image from self allowing narrowing conversions
        /// @param other    Image to copy to
        /// @returns        Reference to the self or throws if target image is too small
        template <typename dst_type>
        const image& copy_to(image<dst_type> &other) const
        {
            foreach([](pixel_type &src, dst_type &dst) { dst = (dst_type)src; }, other);
            return *this;
        }

        /// @brief Populates other image from self allowing narrowing conversions
        /// @param other    Image to copy to
        /// @returns        Reference to the self or throws if target image is too small
        template <typename dst_type>
        const image& copy_to(image<dst_type> &&other) const { return copy_to(other); }


        /// @brief Serializes the image data to raw pointer (allowing narrowing conversions)
        /// @param ptr      Pointer to copy to
        /// @returns        Reference to the self or throws if pointer is null
        template <typename dst_type>
        const image& copy_to(dst_type *ptr) const
        {
            if (ptr == nullptr)
                throw std::runtime_error("Null pointer not allowed");

            return copy_to(image<dst_type>(size(), ptr));
        }

        /// @brief Serializes the image data to vector allowing change of type
        /// @param      container Vector to populate. Resized to exact size
        /// @returns    Reference to the self
        template <typename dst_type>
        const image& to_vector(std::vector<dst_type> &container) const
        {
            container.resize(_width * _height);
            copy_to(container.data());
            return *this;
        }

        /// @brief Serializes the image data to vector (of same type)
        /// @returns    Vector containing serialized data from image
        std::vector<pixel_type>
        to_vector() const
        {
            if (is_contiguous())
                return std::vector<pixel_type>(_begin, _begin + _width * _height);

            std::vector<pixel_type> data(_width * _height);
            copy_to(data.data());
            return data;
        }

        /// @brief Access individual pixels from image without boundary checks
        /// @param row  Row of pixel
        /// @param col  Column of pixel
        /// @returns    Reference to pixel at row and column
        pixel_type &operator()(int row, int col) const
        {
            return _begin[_skip_y * row + _skip_x * col];
        }

        /// @brief Access individual pixels from image without boundary checks
        /// @param row  Row of pixel
        /// @param col  Column of pixel
        /// @returns    Reference to pixel at row and column
        pixel_type &at(int row, int col) const
        {
            return operator()(row, col);
        }

        /// @brief Access individual pixels from image without boundary checks
        /// @param point    Coordinate as roi_point
        /// @returns        Reference to pixel at given coordinate
        pixel_type &operator()(roi_point point) const
        {
            return operator()(point._y, point._x);
        }

        /// @brief Returns an iterator to self with x and y iteration orders exchanged
        /// @returns        Transposed view to self
        image transpose() const
        {
            auto transpose = *this;        // Take copy -- increase shared_ptr ref count
            transpose._width = _height;
            transpose._height = _width;
            transpose._skip_x = _skip_y;
            transpose._skip_y = _skip_x;
            return transpose;
        }

        /// @brief Iterator to self with left/right iteration order swapped
        /// @returns        Horizontally mirrored view to self
        image mirror() const
        {
            auto mirrored = *this;
            mirrored._begin = &at(0, _width - 1);
            mirrored._skip_x = -_skip_x;
            return mirrored;
        }

        /// @brief Iterator to self with top/bottom iteration order swapped
        /// @returns        Vertically mirrored view to self
        image flip() const
        {
            auto flipped = *this;
            flipped._begin = &at(_height - 1, 0);
            flipped._skip_y = -_skip_y;
            return flipped;
        }

        /// @brief Iterator to self skipping every Nth and Mth sample
        /// @param skip_y       Number of rows to advance for each input row
        /// @param skip_x       Number of columns to advance every iteration
        /// @param offset_y     Starting row of iteration, defaults to zero
        /// @param offset_y     Starting column of iteration, defaults to zero
        /// @returns            Sub-sampled / decimated view to self
        image subview(int skip_y, int skip_x, int offset_y = 0, int offset_x = 0) const
        {
            if (skip_y <= 0 || skip_x <= 0)
                throw std::runtime_error("Skip count must be larger than zero");

            auto subsampled = *this;
            subsampled._height = 1 + (_height - offset_y - 1) / skip_y;
            subsampled._width = 1 + (_width - offset_x - 1) / skip_x;
            subsampled._begin = &at(offset_y, offset_x);
            subsampled._skip_y = _skip_y * skip_y;
            subsampled._skip_x = _skip_x * skip_x;
            return subsampled;
        }

        /// @brief Iterator to self skipping every Nth and Mth sample
        /// @param block_size   Skip counts as a tuple
        /// @param offset       Coordinate (in pixels) where the first macro block starts from
        /// @returns            Sub-sampled / decimated view to self
        image subview(roi_point block_size, roi_point offset = { 0, 0 }) const
        {
            return subview(block_size._y, block_size._x, offset._y, offset._x);
        }

        /// @brief Get a view to rectangular area of self
        /// @param  height      Absolute new height (when > 0), relative new height when <= 0
        /// @param  width       Absolute new width (when > 0), relative new width when <= 0
        /// @param  offset_y    Relative position to origin
        /// @param  offset_x    Relative position to origin
        /// @returns            Iterator to rectangular region of self
        image region(int height, int width, int offset_y = 0, int offset_x = 0) const
        {
            auto region = *this;
            region._width = width > 0 ? width : _width + width;
            region._height = height > 0 ? height : _height + height;
            region._begin = &at(offset_y, offset_x);
            return region;
        }

        /// @brief Get a view to rectangular area of self
        /// @param  size        New size as a coordinate pair
        /// @param  offset      Relative position to origin
        /// @returns            Iterator to rectangular region of self
        image region(roi_point size, roi_point offset = roi_point{ 0, 0 }) const
        {
            return region(size._y, size._x, offset._y, offset._x);
        }

        /// @brief Get a view to selected columns from self
        /// @param  width       Number of columns to select
        /// @param  offset_x    Relative position to origin
        /// @returns            Iterator to columns from self
        image<pixel_type> columns(int width, int offset_x = 0) const
        {
            auto columns = *this;
            columns._width = width;
            columns._begin = &at(0, offset_x);
            return columns;
        }

        /// @brief Get a view to selected rows from self
        /// @param  height      Number of rows to select
        /// @param  offset_y    Starting row
        /// @returns            Iterator to rows from self
        image<pixel_type> rows(int height, int offset_y = 0) const
        {
            auto rows = *this;
            rows._height = height;
            rows._begin = &at(offset_y, 0);
            return rows;
        }

        /// @brief Validate that given parameters are not larger than self
        template <typename... ts>
        void validate_sizes(ts... sizes) const
        {
            for (auto &p : std::vector<roi_point>{ sizes... })
            {
                if (_width > p._x || _height > p._y)
                    throw std::runtime_error("Input image dimension too small");
            }
        }

        /// @brief  Applies a function to arbitrary number of rows element wise
        /// @tparam func    Class with prototype of `void(ts... &)`
        /// @tparam ts      Class allowing iteration with `[]`
        /// @param  width   Number of items to iterate
        /// @param  f       function to apply to all iterators as in `f(arg0, arg1, ... argN)`
        /// @param  ptrs    Iterators pointing to beginning of row
        template <class func, typename... ts>
        static void iterate_row(int width, func &&f, ts... ptrs)
        {
            for (int i = 0; i < width; i++)
                f(ptrs[i]...);
        }

        /// @brief Iterator class allowing random access to rows with constant skip.
        template <typename pix_type>
        struct row_iterator
        {
            pix_type *ptr;
            int skipx;
            pix_type& operator[](int column) { return ptr[skipx * column]; }
            row_iterator(pix_type *row, int x) : ptr(row), skipx(x) { }
        };

        /// @brief Iterator maker to allow type deduction from input parameter
        /// @param row  starting row
        /// @returns    row_iterator to given image
        row_iterator<pixel_type> make_row_iterator(int row) const
        {
            return row_iterator<pixel_type>(&at(row, 0), _skip_x);
        }

        /// @brief  Zip-Iterator for arbitrary number of images of any type
        /// @tparam func    Type that can _called_ as `func(pixel0, pixel1, ...)`
        /// @tparam types   Generic types for inputs that bind to everything
        /// @param  f       Function to apply to element wise to each input image
        /// @param  images  Input images - these can be const, ref or refref images of any underlying type
        /// @returns        Input function (functor) or throws if input images are not large enough
        template <class func, typename... types> func&
            foreach(func &&f, types&&... images) const
        {
            // Check that the other input images are at least same size as this is
            validate_sizes(images.size()...);

            auto skipx = std::initializer_list<int>({ _skip_x, images._skip_x... });
            auto all_ones = std::all_of(skipx.begin(), skipx.end(), [](const int skip) { return skip == 1; });
            if (all_ones)
            {
                for (auto y : _height)
                {
                    iterate_row(_width, f, &at(y, 0), &images.at(y, 0)...);
                }
            }
            else
            {
                for (auto y : _height)
                {
                    iterate_row(_width, f, make_row_iterator(y), images.make_row_iterator(y)...);
                }
            }
            return f;
        }

        /// @brief  Iterate all pixels applying a [lambda] function or functor
        /// @tparam UnaryFunction   Class callable as `func(type &pix)`
        /// @param  func            Function to apply to each pixel
        /// @param  skip_rows       Rows to advance -- defaults to 1, i.e. scan every pixel
        /// @param  skip_cols       Columns to advance -- defaults to 1, i.e. scan every pixel
        /// @returns                Reference to given input function (functor)
        template <class UnaryFunction> UnaryFunction&
            for_each(UnaryFunction &&func, const int skip_rows = 1, const int skip_cols = 1) const
        {
            for (int j = 0; j < _height; j += skip_rows)
            {
                pixel_type *ptr = &operator()(j, 0);
                int skip = _skip_x * skip_cols;
                for (int i = 0; i < _width; i += skip_cols, ptr += skip)
                    func(*ptr);
            }
            return func;
        }

        /// @brief Returns position of smallest element
        /// @returns coordinate of first smallest element
        roi_point min_element()
        {
            roi_point result;
            pixel_type minimum = (*this)(result);
            for (auto j: _height)
            {
                pixel_type *src = &at(j, 0);
                const int skip_x = _skip_x;
                for (auto i : _width)
                {
                    pixel_type pix = src[(int)i * skip_x];
                    if (pix < minimum)
                    {
                        minimum = pix;
                        result._x = (int)i;
                        result._y = (int)j;
                    }
                }
            }
            return result;
        }

        /// @brief Returns position of last largest element
        /// @returns coordinate of last largest element
        roi_point max_element()
        {
            roi_point result;
            pixel_type maximum = (*this)(result);
            for (auto j : _height)
            {
                pixel_type *src = &at(j, 0);
                const int skip_x = _skip_x;
                for (auto i : _width)
                {
                    pixel_type pix = src[(int)i * skip_x];
                    if (pix > maximum)
                    {
                        maximum = pix;
                        result._x = (int)i;
                        result._y = (int)j;
                    }
                }
            }
            return result;
        }

        /// @brief Apply a custom transformation to image returning a new image of equal size
        /// Transform functions returns a new range/image from existing image using a transform
        /// with custom function `dst_type UnaryFunction(const src_type &input);`
        template <typename dst_type = pixel_type, class UnaryFunction>
        image<dst_type> transform(UnaryFunction &&func) const
        {
            image<dst_type> buf(size());
            if (is_contiguous())
            {
                auto pixels = range(_height * _width);
                dst_type *dst = buf._begin;
                pixel_type *src = _begin;
                for (auto i : pixels)
                    dst[i] = func(src[i]);
            }
            else
            {
                for (auto i : _height)
                {
                    dst_type *dst = &buf.at(i, 0);
                    pixel_type *src = &at(i, 0);
                    for (auto j : _width)
                        dst[j] = func(src[j * _skip_x]);
                }
            }
            return buf;
        }

        /// @brief Creates a new contiguous image of given type filled with original
        /// @returns A new image of same size with logically the same content
        template <typename dst_type = pixel_type>
        image<dst_type> convert_to() const {
            image<dst_type> buf(size());
            copy_to(buf);
            return buf;
        }

        /// @brief  Fills the borders from self according a replication scheme
        /// @param top          Size of top border
        /// @param left         Size of left border
        /// @param bottom       Size of bottom border
        /// @param right        Size of right border
        /// @param rep_scheme   Replication scheme for each or all borders
        /// @returns            reference to self
        image<pixel_type> &copy_borders(int32_t top, int32_t left, int32_t bottom, int32_t right, uint32_t rep_scheme)
        {
            auto width = _width - left - right;
            auto height = _height - top - bottom;
            // Copy top rows
            if (top > 0)
            {
                auto dst_top = region(top, width, 0, left);
                if ((rep_scheme & TOP) == (TOP & ZERO_PAD))
                {
                    dst_top.fill(0);
                }
                else
                {
                    auto src_top = dst_top.region(0, 0, top, 0);
                    switch (rep_scheme & TOP)
                    {
                    case TOP & REPLICATE:   // Copy closest row to border all over again
                        src_top._skip_y = 0;
                        break;
                    case TOP & MIRROR_ODD:  // Mirror copy a block of rows with offset of 1 row to border
                        src_top._begin += src_top._skip_y;
                        // Fall through
                    case TOP & MIRROR_EVEN: // Mirror copy a block of rows with offset of 0 rows to border
                        src_top = src_top.flip();
                        break;
                    case TOP & REPLICATE_EVEN:  // Block copy a block of rows
                        break;
                    default:
                        throw std::runtime_error("Unknown replication mode");
                    }
                    dst_top.init_from(src_top);
                }
            }
            // Copy bottom rows
            if (bottom > 0)
            {
                auto dst_bottom = region(bottom, width, top + height, left);
                if ((rep_scheme & BOTTOM) == (BOTTOM & ZERO_PAD))
                {
                    dst_bottom.fill(0);
                }
                else
                {
                    auto src_bottom = dst_bottom.region(0, 0, -(int)bottom, 0).flip();
                    switch (rep_scheme & BOTTOM)
                    {
                    case BOTTOM & REPLICATE:
                        src_bottom._skip_y = 0;
                        break;
                    case BOTTOM & MIRROR_ODD:
                        src_bottom._begin += src_bottom._skip_y;
                        // Fall through
                    case BOTTOM & MIRROR_EVEN:
                        break;
                    case BOTTOM & REPLICATE_EVEN:
                        src_bottom = src_bottom.flip();
                        break;
                    default:
                        throw std::runtime_error("Unknown replication mode");
                    }
                    dst_bottom.init_from(src_bottom);
                }
            }
            // Copy Left border columns
            if (left > 0)
            {
                auto dst_left = region(0, left, 0, 0);
                if ((rep_scheme & LEFT) == (LEFT & ZERO_PAD))
                {
                    dst_left.fill(0);
                }
                else
                {
                    auto src_left = dst_left.region(0, 0, 0, left);
                    switch (rep_scheme & LEFT)
                    {
                    case LEFT & REPLICATE:
                        src_left._skip_x = 0;
                        break;
                    case LEFT & MIRROR_ODD:
                        src_left._begin += src_left._skip_x;
                        // Fall through
                    case LEFT & MIRROR_EVEN:
                        src_left = src_left.mirror();
                        break;
                    case LEFT & REPLICATE_EVEN:
                        break;
                    default:
                        throw std::runtime_error("Unknown replication mode");
                    }
                    dst_left.init_from(src_left);
                }
            }
            // Copy right margin
            if (right > 0)
            {
                auto dst_right = region(0, right, 0, left + width);
                if ((rep_scheme & RIGHT) == (RIGHT & ZERO_PAD))
                {
                    dst_right.fill(0);
                }
                else
                {
                    auto src_right = dst_right.region(0, 0, 0, -(int)right).mirror();
                    switch (rep_scheme & RIGHT)
                    {
                    case RIGHT & REPLICATE:
                        src_right._skip_x = 0;
                        break;
                    case RIGHT & MIRROR_ODD:
                        src_right._begin += src_right._skip_x;
                        // Fall through
                    case RIGHT & MIRROR_EVEN:
                        break;
                    case RIGHT & REPLICATE_EVEN:
                        src_right = src_right.mirror();
                        break;
                    default:
                        throw std::runtime_error("Unknown replication mode");
                    }
                    dst_right.init_from(src_right);
                }
            }
            return *this;
        }

        /// @brief Appends borders to self according to a replication scheme
        /// @param top          Size of top border
        /// @param left         Size of left border
        /// @param bottom       Size of bottom border
        /// @param right        Size of right border
        /// @param rep_scheme   Replication scheme for each or all borders
        /// @returns            A new image with borders added
        image<pixel_type> make_borders(int top, int left, int bottom, int right, uint32_t rep_scheme)
        {
            auto result = image<pixel_type>(_height + top + bottom, _width + left + right);
            // Copy center from source
            copy_to(result.region(size(), roi_point{ left, top }));

            result.copy_borders(top, left, bottom, right, rep_scheme);
            return result;
        }

        /// @brief Append borders to self - e.g. for replicating bayer patterns
        /// @param dims         Number of pixels to be added
        /// @param rep_scheme   Replicating scheme defaulting to copying blocks of size `dims`
        /// @returns            A new image with borders added
        image<pixel_type> make_borders(roi_point dims, uint32_t rep_scheme = REPLICATE_EVEN)
        {
            return make_borders(dims._y, dims._x, dims._y, dims._x, rep_scheme);
        }

        /// @brief Filters an image using a fixed support
        /// @details  `this` class works as the target container with arbitrary layout and
        /// can be same as src_image for _logical inplace_ operation
        /// @param  src_image   Image to be filtered
        /// @param  func        Function of prototype `T func(support<W,H,source_type> &)`
        /// @param  rep_scheme  Border replication scheme
        /// @returns            Binary function (to capture side effects)
        template <int HEIGHT, int WIDTH = HEIGHT, typename src_type, class BinaryFunction>
        BinaryFunction& filter(image<src_type> src_image, BinaryFunction &&func, uint32_t rep_scheme = 0)
        {
            static_assert(HEIGHT > 0 && WIDTH > 0, "Kernel dimensions must be positive");
            const auto left = (WIDTH - 1) >> 1;
            const auto top = (HEIGHT - 1) >> 1;

            if (size() == src_image.size())
                src_image = src_image.make_borders(top, left, HEIGHT - top - 1, WIDTH - left - 1, rep_scheme);
            else if (src_image._width < (_width + WIDTH - 1) || src_image._height < (_height + HEIGHT - 1))
                throw std::runtime_error("Src image is not large enough");

            auto matrix = support<HEIGHT, WIDTH, src_type>();
            for (auto i : _height)
            {
                src_type *src = &src_image(i, 0);
                pixel_type *dst = &at(i, 0);

                for (int j = 0; j < WIDTH - 1; j++)
                {
                    matrix.shift(src, src_image._skip_y);
                    src += src_image._skip_x;
                }

                for (auto j : _width)
                {
                    matrix.shift(src, src_image._skip_y);
                    dst[j] = func(matrix);
                    src += src_image._skip_x;
                }
            }
            return func;
        }

        /// @brief Filters an image using a fixed support returning a new image
        /// @details  The output image is not necessarily contiguous -- instead the possible
        /// implementation can allocate the result in excess, use that image _in place_
        /// and return a properly sized image
        /// @param  func        Function of prototype `T func(support<W,H,source_type> &)`
        /// @param  rep_scheme  Border replication scheme
        /// @returns            Filtered image without guarantee of contiguousness
        template <int HEIGHT, int WIDTH = HEIGHT, class BinaryFunction>
            image<pixel_type> filter(BinaryFunction &&func, uint32_t rep_scheme = 0)
        {
            static_assert(HEIGHT > 0 && WIDTH > 0, "Kernel dimensions must be positive");
            const auto left = (WIDTH - 1) >> 1;
            const auto top = (HEIGHT - 1) >> 1;
            // Allocate target (and intermediate) image with excess
            auto result =
                this->make_borders(top, left, HEIGHT - top - 1, WIDTH - left - 1, rep_scheme);

            auto matrix = support<HEIGHT, WIDTH, pixel_type>();

            for (auto i : _height)
            {
                pixel_type *src = &result.at(i, 0);
                pixel_type *dst = src;

                for (int x = 0; x < WIDTH - 1; x++)
                    matrix.shift(src++, result._skip_y);

                for (auto j : _width)
                {
                    matrix.shift(src++, result._skip_y);
                    dst[j] = func(matrix);
                }
            }
            return result.region(size());
        }

        /// @brief Element wise accumulate with other image
        /// @returns    reference to self
        image& operator += (const image &other)
        {
            foreach([](pixel_type &me, pixel_type &src) { me += src; }, other);
            return *this;
        }

        /// @brief Element wise subtract with other image
        /// @returns    reference to self
        image& operator -= (const image &other) /* const */
        {
            foreach([](pixel_type &me, pixel_type &src) { me -= src; }, other);
            return *this;
        }

        /// @brief Element wise multiply with other image
        /// @returns    reference to self
        image& operator *= (const image &other) /* const */
        {
            foreach([](pixel_type &me, pixel_type &src) { me *= src; }, other);
            return *this;
        }

        /// @brief Element wise divide with other image
        /// @returns    reference to self
        image& operator /= (const image &other) /* const */
        {
            foreach([](pixel_type &me, pixel_type &src) { me /= src; }, other);
            return *this;
        }
    };
};
