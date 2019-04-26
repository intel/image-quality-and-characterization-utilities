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
#include <vector>

namespace Teisko
{
    template <typename T>
    T mean(std::vector<T> &vec)
    {
        average_f<T> avg{};
        for (auto &d : vec)
            avg(d);
        return static_cast<T>(avg);
    }

    template <typename U, typename T>
    std::vector<U> norm(std::vector<T> &data)
    {
        auto result = std::vector<U>(data.size());
        for (size_t i = 0; i < data.size(); i++)
            result[i] = static_cast<U>(sqrt(norm2(data[i])));
        return result;
    }

    // element wise add two vectors
    template <typename T, typename U>
    std::vector<T>& operator +=(std::vector<T> &data, const std::vector<U> &other)
    {
        if (other.size() != data.size())
            throw std::runtime_error("Vector size mismatch");
        for (size_t i = 0; i < data.size(); i++)
            data[i] += other[i];
        return data;
    }

    // scalar add to vector
    template <typename T, typename U>
    std::vector<T> operator +=(std::vector<T> &data, const U other)
    {
        for (auto &d : data)
            d += other;
        return data;
    }

    // generic add vector + (vector or scalar)
    template <typename T, typename U>
    std::vector<T> operator +(const std::vector<T> &data, const U other)
    {
        auto copy = data;
        return copy += other;
    }

    // element wise subtract two vectors
    template <typename T, typename U>
    std::vector<T>& operator -=(std::vector<T> &data, const std::vector<U> &other)
    {
        if (other.size() != data.size())
            throw std::runtime_error("Vector size mismatch");
        for (size_t i = 0; i < data.size(); i++)
            data[i] -= other[i];
        return data;
    }

    // scalar subtract to vector
    template <typename T, typename U>
    std::vector<T> operator -=(std::vector<T> &data, const U other)
    {
        for (auto &d : data)
            d -= other;
        return data;
    }

    // generic add vector + (vector or scalar)
    template <typename T, typename U>
    std::vector<T> operator -(const std::vector<T> &data, const U other)
    {
        auto copy = data;
        return copy -= other;
    }

    // Element wise multiply two vectors
    template <typename T, typename U>
    std::vector<T>& operator *=(std::vector<T> &data, const std::vector<U> &other)
    {
        if (other.size() != data.size())
            throw std::runtime_error("Vector size mismatch");
        for (size_t i = 0; i < data.size(); i++)
            data[i] *= other[i];
        return data;
    }

    // Element wise multiply vector and scalar
    template <typename T, typename U>
    std::vector<T> operator *=(std::vector<T> &data, const U other)
    {
        for (auto &d : data)
            d *= other;
        return data;
    }

    // Generic scale vector * (vector or scalar)
    template <typename T, typename U>
    std::vector<T> operator *(const std::vector<T> &data, const U other)
    {
        auto copy = data;
        return copy *= other;
    }

    // Element wise divide vector
    template <typename T, typename U>
    std::vector<T>& operator /=(std::vector<T> &data, const std::vector<U> &other)
    {
        if (other.size() != data.size())
            throw std::runtime_error("Vector size mismatch");
        for (size_t i = 0; i < data.size(); i++)
            data[i] /= other[i];
        return data;
    }

    // Element wise divide vector by scalar
    template <typename T, typename U>
    std::vector<T> operator /=(std::vector<T> &data, const U other)
    {
        for (auto &d : data)
            d /= other;
        return data;
    }

    // Generic divide vector by (vector or scalar)
    template <typename T, typename U>
    std::vector<T> operator /(const std::vector<T> &data, const U other)
    {
        auto copy = data;
        return copy /= other;
    }

    /// Interpolates/extrapolates linearly from  y(t=0.0) = f0, y(t=1.0) = f1
    template <typename T>
    T linear_mix(const T &f0, const T &f1, double t) { return f0 * (1.0 - t) + f1 * t; }

    /// interpolates the "i"th sample out of n
    /// from a triangular shape defined as f(0) = f(n-1) = edge; f((n-1)/2) = middle;
    inline double triangle_ramp_func(size_t i, size_t n, double edge, double middle)
    {
        double factor = std::abs(2.0 * i - (n - 1)) / (n - 1);
        return linear_mix(middle, edge, factor);
    }

    // 2d-geometry -- returns true if line segment p0->p1 crosses the line segment p2->p3
    inline bool check_line_segment_crossing(point_xy p0, point_xy p1, point_xy p2, point_xy p3)
    {
        auto s1 = p1 - p0;
        auto s2 = p3 - p2;
        auto d = p0 - p2;

        auto denom = s1.x * s2.y - s2.x * s1.y;
        // the line segments are parallel -- should we check if bounding boxes intersect?
        // - we trust that if this happens in the "in_polygon" test, then the next
        //   line segment has same end point and is not parallel
        if (denom == 0.0)
            return false;

        auto s = s1.x * d.y - s1.y * d.x;
        auto t = -s2.y * d.x + s2.x * d.y;
        if (denom < 0.0)
        {
            return (s >= denom) && (s <= 0) && (t > denom) && (t <= 0);
        }
        return (s <= denom) && (s >= 0) && (t <= denom) && (t >= 0);
    }

    /// This polygon test does not assume anything about the points
    /// The polygon can be concave or complex (having self intersections or holes)
    /// The polygon can be closed or open
    /// - each line segment is tested against a ray cast to 'value' from outside
    ///   of the polygons bounding box. If there's an odd amount of intersections
    ///   between that ray and line segments, the 'value' is inside the polygon
    /// - complexity: O(N)
    inline bool is_inside_polygon(const std::vector<point_xy> &vec, point_xy value)
    {
        auto n = vec.size();
        if (n < 2)
            return false;
        auto mins = vec.front();
        auto maxs = vec.front();
        for (auto &p : vec)
        {
            mins.x = std::min(mins.x, p.x);
            mins.y = std::min(mins.y, p.y);
            maxs.x = std::max(maxs.x, p.x);
            maxs.y = std::max(maxs.y, p.y);
        }
        // bounding box
        if (value.x < mins.x || value.y < mins.y ||
            value.x > maxs.x || value.y > maxs.y)
            return false;

        // get a value outside the bounding box -- then check that even number of line segments
        // in vec pass though the line segment from value to point outside the bounding box
        auto outside = point_xy(mins.x - 1.0, mins.y);

        // fetch the last unique point as the "previous"
        // then loop over the n[-1] unique line segments
        auto prev = vec.back();
        if (prev == vec.front())
        {
            prev = vec[n-2];
            --n;
        }

        bool inside = false;
        for (decltype(n) i = 0; i < n; i++)
        {
            if (check_line_segment_crossing(outside, value, prev, vec[i]))
                inside = !inside;
            prev = vec[i];
        }
        return inside;
    }
}
