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
#include <algorithm>
#include <cmath>

#include <type_traits>

namespace Teisko
{
    /// @brief      Ranged for iterators for arithmetic types
    /// range, or range_s<int> is to be preferred / short hand notation
    /// for any member in a class that should work as a terminating count for iterating
    /// @tparam     T   type of value
    template<
        typename T,
        typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type
    >
    struct range_s
    {
        /// @brief Maximum value to iterate (or current value of iterator)
        T value{};

        /// @brief  default constructor
        range_s() = default;

        /// @brief  Alternative constructor
        /// @details omitting the keyword _explicit_ to allow arithmetic operations between class and integers
        range_s(T i) : value(i) { }

        /// @brief returns value when cast to `T`
        operator T() const { return value; }

        /// @brief Returns starting iterator for ranged for
        static range_s begin() { return{}; }

        /// @brief Returns ending iterator for ranged for
        range_s end() const { return *this; }

        /// @brief Dereferencing operator for ranged for
        T& operator *() { return value; }

        /// @brief Inequality operator for ranged for
        bool operator != (const range_s &rhs) const { return value != rhs.value; }

        /// @brief Less than operator for legacy loops
        bool operator < (const range_s &rhs) const { return value < rhs.value; }

        /// @brief Less than operator for legacy loops
        bool operator < (const int rhs) const { return value < rhs; }

        /// @brief Less than or equal operator for legacy loops
        bool operator <= (const range_s &rhs) const { return value <= rhs.value; }

        /// @brief Prefix increment operator for ranged for
        range_s& operator ++() { ++value; return *this; }

        /// @brief Postfix increment operator when used in legacy for loops
        range_s operator ++(int) { auto result = *this; ++value; return result; }

        /// @brief Adds a constant to range -- used in legacy for loops with skip
        range_s operator +=(const int rhs) { value += rhs; return *this; }
    };

    using range = range_s<int>;

    /// @brief iterator class to iterate a 2d range row major order (y first, x second)
    struct xy_iterator_s
    {
        int _x{};
        int _y{};
        /// @brief constructor to setup initial value of iterating at 0,0 to w*h
        xy_iterator_s(int w, int h) : _width(w), _height(h) {}

        /// @brief constructor to setup end value of iterating to w*h
        explicit xy_iterator_s(int rows) : _x(0), _y(rows) { }

        /// @brief prefix increment for ranged for
        xy_iterator_s& operator ++() {
            ++_x;
            if (_x == _width)
            {
                ++_y;
                _x = 0;
            }
            return *this;
        }

        /// @brief dereferencing operator for ranged for
        xy_iterator_s& operator *() { return *this; }

        /// @brief Inequality operator for ranged for
        bool operator != (const xy_iterator_s &rhs) const { return _x != rhs._x || _y != rhs._y; }
    private:
        int _width{};
        int _height{};
    };

    /// Simple container for X/Y --coordinate pair passing
    struct roi_point
    {
        int _x{};
        int _y{};
        roi_point() = default;

        /// @brief  Converting constructor
        /// @details allows expressions of type (1 + roi_point) and (roi_point op int)
        /// to convert implicitly to `roi_point op roi_point` by omitting keyword `explicit`
        roi_point(int k) : _x(k), _y(k) { }
        /// @brief  Constructor for coordinates
        roi_point(int x, int y) : _x(x), _y(y) { }

        /// @brief Constructor taking floating point types rounding to nearest integer
        template <typename T>
        roi_point(T x, T y, typename std::enable_if<std::is_floating_point<T>::value, T>::type = 0)
            : _x(std::lround(x)), _y(std::lround(y)) { }

        /// @brief Converting constructor to allow xy_iterator to implicitly convert to roi_point
        /// for image(roi_point) indexing
        roi_point(const xy_iterator_s &other) : _x(other._x), _y(other._y) { }

        xy_iterator_s begin() { return xy_iterator_s(_x, _y); }
        xy_iterator_s end() { return xy_iterator_s(_y); }
    };

    template <typename pixel_type>
    struct xy_skip_iterator_s
    {
        pixel_type *ptr{};        ///< Current pointer
        pixel_type *eol{};        ///< Marks end of line to allow changing to next line
        int skip_x{};             ///< Spacing between horizontal members
        int skip_y{};             ///< Spacing between vertical members
        int skip_xy{};            ///< Spacing between end of line to start of next line

        /// @brief access the neighbors of current location without boundary checks
        /// @returns    reference to pixel at given offset
        pixel_type& operator()(int y, int x) const {
            return ptr[y * skip_y + x * skip_x];
        }

        /// @brief dereference operator for ranged for
        /// @returns reference to pixel pointed by iterator
        pixel_type& operator *() {
            return *ptr;
        }

        /// @brief  prefix increment operator for ranged for
        /// @returns reference to self
        xy_skip_iterator_s& operator ++()
        {
            if (ptr == eol)
            {
                ptr += skip_xy;
                eol += skip_y;
            }
            else
            {
                ptr += skip_x;
            }
            return *this;
        }

        /// @brief Inequality operator for ranged for
        /// @details the end ptr must be placed at last row, past end of line or `&at(_height - 1, width)`
        /// @returns true if the pointers match
        bool operator != (const xy_skip_iterator_s &rhs) const
        {
            return ptr != rhs.ptr;
        }

        /// @brief Constructs an image iterator from pointer, strides and size
        /// @details The design allows (fast) prefix increment with single comparison
        /// avoiding multiplications in the inner loop
        xy_skip_iterator_s(pixel_type *_ptr, roi_point skip, roi_point size)
            : ptr(_ptr)
            , eol(ptr + (size._x - 1) * skip._x)
            , skip_x(skip._x)
            , skip_y(skip._y)
            , skip_xy(skip._y - (size._x - 1) * skip._x)
        { }
    };

    /// Returns true if points are equal
    inline bool operator == (const roi_point &a, const roi_point &b)
    {
        return a._x == b._x && a._y == b._y;
    }

    /// Returns true if points are not equal
    inline bool operator != (const roi_point &a, const roi_point &b)
    {
        return !(a == b);
    }

    /// Pairwise add coordinates
    inline roi_point operator +(const roi_point &a, const roi_point &b)
    {
        return{ a._x + b._x, a._y + b._y };
    }

    /// Pairwise subtract coordinates
    inline roi_point operator -(const roi_point &a, const roi_point &b)
    {
        return{ a._x - b._x, a._y - b._y };
    }

    /// Pairwise multiply coordinates
    inline roi_point operator *(const roi_point &a, const roi_point &b)
    {
        return{ a._x * b._x, a._y * b._y };
    }

    /// Divide coordinate element wise by another coordinate
    inline roi_point operator /(const roi_point &a, const roi_point &b)
    {
        return{
            b._x == 0 ? 0 : a._x / b._x,
            b._y == 0 ? 0 : a._y / b._y };
    }

    /// Scalar Add value to coordinate
    template <typename T>
    typename std::enable_if<std::is_arithmetic<T>::value, roi_point>::type
        operator +(const roi_point &a, T b)
    {
        return roi_point(a._x + b, a._y + b);
    }

    /// Scalar Subtract value from coordinate
    template <typename T>
    typename std::enable_if<std::is_arithmetic<T>::value, roi_point>::type
        operator -(const roi_point &a, T b)
    {
        return roi_point(a._x - b, a._y - b);
    }

    /// Scalar Multiply coordinate by value
    template <typename T>
    typename std::enable_if<std::is_arithmetic<T>::value, roi_point>::type
        operator *(const roi_point &a, T b)
    {
        return roi_point(a._x * b, a._y * b);
    }

    /// Scalar Divide coordinate by value
    template <typename T>
    typename std::enable_if<std::is_arithmetic<T>::value, roi_point>::type
        operator /(const roi_point &a, T b)
    {
        if (b == (T)0)
            return{ 0, 0 };
        return{ a._x / b, a._y / b };
    }
}