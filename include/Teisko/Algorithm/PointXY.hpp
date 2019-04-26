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

// In alphabetical order
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <stdexcept>
#include <vector>

/// \brief This namespace defines the public interfaces of
/// the \ref libcalc_module module
namespace Teisko
{
    struct point_xy
    {
        double x;
        double y;
        point_xy() : x(0.0), y(0.0) { }
        point_xy(double xy) : x(xy), y(xy) { }
        point_xy(double _x, double _y) : x(_x), y(_y) { }
    };

    inline bool isfinite(const point_xy &a) { return std::isfinite(a.x) && std::isfinite(a.y); }

    inline bool operator ==(const point_xy &a, const point_xy &b) { return a.x == b.x && a.y == b.y; }
    inline point_xy operator -(const point_xy &a) { return point_xy(-a.x, -a.y); }
    inline point_xy operator -(const point_xy &a, const point_xy &b) { return point_xy(a.x - b.x, a.y - b.y); }
    inline point_xy operator +(const point_xy &a, const point_xy &b) { return point_xy(a.x + b.x, a.y + b.y); }
    inline point_xy operator *(const point_xy &a, const point_xy &b) { return point_xy(a.x * b.x, a.y * b.y); }
    inline point_xy operator /(const point_xy &a, const point_xy &b) { return point_xy(a.x / b.x, a.y / b.y); }
    inline point_xy& operator -=(point_xy &a, const point_xy &b) { a.x -= b.x; a.y -= b.y; return a; }
    inline point_xy& operator +=(point_xy &a, const point_xy &b) { a.x += b.x; a.y += b.y; return a; }
    inline point_xy& operator *=(point_xy &a, const point_xy &b) { a.x *= b.x; a.y *= b.y; return a; }
    inline point_xy& operator /=(point_xy &a, const point_xy &b) { a.x /= b.x; a.y /= b.y; return a; }

    inline double dot_prod(const point_xy &a, const point_xy &b) { return a.x*b.x + a.y * b.y; }
    inline double norm2(const point_xy &a) { return dot_prod(a, a); }
    inline double norm(const point_xy &a) { return std::sqrt(norm2(a)); }

    template <typename T>
    typename std::enable_if<std::is_arithmetic<T>::value, point_xy>::type
        operator *(const point_xy &a, T b) { return point_xy(a.x * b, a.y * b); }

    template <typename T>
    typename std::enable_if<std::is_arithmetic<T>::value, point_xy>::type
        operator *(T a, const point_xy &b) { return point_xy(a * b.x, a * b.y); }

    template <typename T>
    typename std::enable_if<std::is_arithmetic<T>::value, point_xy>::type
        operator /(const point_xy &a, T b) { return point_xy(a.x / b, a.y / b); }

    template <typename T>
    typename std::enable_if<std::is_arithmetic<T>::value, point_xy>::type
        operator /(T a, const point_xy &b) { return point_xy(a / b.x, a / b.y); }

    template <typename T>
    typename std::enable_if<std::is_arithmetic<T>::value, point_xy>::type
        operator +(const point_xy &a, T b) { return point_xy(a.x + b, a.y + b); }

    template <typename T>
    typename std::enable_if<std::is_arithmetic<T>::value, point_xy>::type
        operator +(T a, const point_xy &b) { return point_xy(a + b.x, a + b.y); }

    template <typename T>
    typename std::enable_if<std::is_arithmetic<T>::value, point_xy>::type
        operator -(const point_xy &a, T b) { return point_xy(a.x - b, a.y - b); }

    template <typename T>
    typename std::enable_if<std::is_arithmetic<T>::value, point_xy>::type
        operator -(T a, const point_xy &b) { return point_xy(a - b.x, a - b.y); }

    template <typename T>
    typename std::enable_if<std::is_arithmetic<T>::value, point_xy>::type
        operator -=(point_xy &a, T b) { a.x -= b; a.y -= b; return a; }

    template <typename T>
    typename std::enable_if<std::is_arithmetic<T>::value, point_xy>::type
        operator +=(point_xy &a, T b) { a.x += b; a.y += b; return a; }

    template <typename T>
    typename std::enable_if<std::is_arithmetic<T>::value, point_xy>::type
        operator *=(point_xy &a, T b) { a.x *= b; a.y *= b; return a; }

    template <typename T>
    typename std::enable_if<std::is_arithmetic<T>::value, point_xy>::type
        operator /=(point_xy &a, T b) { a.x /= b; a.y /= b; return a; }

    /// \brief Returns a vector of input vector of 2d-points rotated by an angle
    /// \param data     Vector of points to be rotated
    /// \param angle    Angle of rotation in radians
    /// \param center   Center of rotation -- if not zero
    /// \returns        Rotated copy of the input data
    inline std::vector<point_xy> rotate_vector(std::vector<point_xy> data, double angle, point_xy center = point_xy(0, 0))
    {
        auto rot_c = std::cos(angle);
        auto rot_s = std::sin(angle);

        // Rotate each point as (x + iy) * (c + is) = x*c - y*s + i *(sx + yc)
        for (auto &c : data)
        {
            c -= center;
            c = point_xy(c.x * rot_c - c.y * rot_s, rot_s * c.x + rot_c * c.y);
            c += center;
        }

        return data;
    }

    /// \brief Computes the linear regression of data set
    ///
    /// Reference:
    /// mathportal.org/calculators/statistic-calculator/correlation-and-regression-calculator.php
    ///
    /// \param set  Vector of x,y coordinates
    /// \returns    Slope and offset
    inline point_xy linear_regression(std::vector<point_xy> set)
    {
        double x = 0, y = 0, xy = 0, xx = 0;
        for (auto &p : set)
        {
            x += p.x;
            y += p.y;
            xy += p.x * p.y;
            xx += p.x * p.x;
        }
        double n = (double)set.size();
        double d = n * xx - x*x;
        if (n > 1 && d != 0)
            return point_xy((n*xy - x*y) / d, (y*xx - x*xy) / d);
        return point_xy(std::numeric_limits<double>::infinity(), x / n);
    }

    // Calculates weighted squared distance of model y = f(x) to given data set
    // { x0, y0 }, ... { xN, yN }
    // The derived class must provide function to evaluate f(x)
    // according to external parameters `double p[M]`
    struct xy_func_optimizer
    {
        xy_func_optimizer operator =(const xy_func_optimizer &other) = delete;

        std::function<double(double x, double *p)> func;
        std::vector<point_xy> xy;
        std::vector<double> weight;

        xy_func_optimizer(
            std::function<double(double x, double *p)> f,
            std::vector<point_xy> points = {},
            std::vector<double> weight_vec = {})
            : func(f), xy(points), weight(weight_vec) { }

        double operator()(double *p)
        {
            double sum = 0.0;
            size_t n = xy.size();
            size_t w = weight.size();
            for (size_t i = 0; i < n; i++)
            {
                double y = func(xy[i].x, p) - xy[i].y;
                sum += y * y * (i < w ? weight[i] : 1.0);
            }
            return sum;
        }
    };

    // Calculates distance of { x, f(x) } to given point { x0, y0 }
    // according to external parameters `double p[M]`
    struct xy_distance_optimizer
    {
        point_xy target;
        virtual double operator()(double x) = 0;
        double operator()(double *x)
        {
            double y = operator()(x[0]);
            return (x[0] - target.x) * (x[0] - target.x)
                + (y - target.y) * (y - target.y);
        }
    };
}
