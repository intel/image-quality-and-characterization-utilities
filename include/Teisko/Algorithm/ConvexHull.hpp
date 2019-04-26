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
    struct convex_hull
    {
    public:
        // Calculates the convex hull from given points using Andrew's Monotone chain algorithm
        // https ://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain
        // The convex hull can have duplicate points and collinear points, which need to be pruned out
        // for minimal convex hull
        // further improvements -- remove the points in output from the input points
        template <typename T>
        static std::vector<T> hull(std::vector<T> &points)
        {
            // Use Andrew's Monotone chain algorithm (O(n log n))
            // - first sort input from left to right (in case of tie, from bottom to top)
            std::sort(points.begin(), points.end(), [](const T &a, const T &b)
            {
                return a.x < b.x || (a.x == b.x && a.y < b.y);
            });

            size_t n = points.size();
            // two points and fewer are either both CW and CCW, collinear or non-unique
            if (n <= 2)
            {
                return std::vector<T>(std::move(points));
            }

            // three points need to be sorted to CW order before returning, let the algorithm do it
            // We allocate just N points -- with just collinear points, the algorithm would add all points
            // on the way right (upper hull) and also on the way left (lower hull)
            // This is limited by adding a guard for hull_pts size
            // other option would be to make a list of indirections, which would only attempt to add
            // points not already added

            auto hull_pts = std::vector<T>(points.size());
            auto hull_ptr = hull_pts.data();        // a raw pointer to vector for sake of efficiency
            auto hull_limit = hull_ptr + 2;         // a lower limit of pruning
            auto hull_limit_end = hull_ptr + hull_pts.size();   // an upper limit of adding points

            // build the upper hull
            for (size_t i = 0; i < n; ++i)
                add_to_hull(hull_ptr, hull_limit, hull_limit_end, points[i]);

            hull_limit = hull_ptr + 1;
            // build lower hull -- starting from the second rightmost point, which was already added
            for (size_t i = n - 2; i > 0; --i)
                add_to_hull(hull_ptr, hull_limit, hull_limit_end, points[i]);

            hull_pts.resize(hull_ptr - hull_pts.data());
            return hull_pts;
        }

        // Generate a convex hull for vector of 'T' (containing elements 'x' and 'y')
        // in clockwise order. The first point is the leftmost point in the original set.
        // The returned polygon is not closed.
        template <typename T>
        static void quick_hull(std::vector<T> &points)
        {
            points = hull(points);
            // Remove successive duplicates
            std::unique(begin(points), end(points), [](const T &a, const T &b) { return (a.x == b.x) && (a.y == b.y); });
        }

        // O(n^2) algo to find smallest feature in a closed simple polygon
        // The polygon must be defined in clockwise order
        template <typename T>
        static void prune_hull(std::vector<T> &points, size_t min_points, size_t max_points, double min_area)
        {
            double area = hull_area(points);
            if (area <= 0)
                return;

            bool is_closed = points.front().x == points.back().x && points.front().y == points.back().y;
            if (!is_closed)
                points.push_back(points.front());

            size_t n = points.size();

            while ((min_points + 1) < (size_t)n)
            {
                double smallest = area;
                size_t index = 0;
                for (size_t i = 1; i + 1 < n; i++)
                {
                    double a = 0;
                    // Shoelace algo to calculate contribution of each triangle formed by consecutive vertices
                    a -= points[i - 1].x * (points[i + 0].y - points[i + 1].y);
                    a -= points[i + 0].x * (points[i + 1].y - points[i - 1].y);
                    a -= points[i + 1].x * (points[i - 1].y - points[i + 0].y);
                    if (a < smallest)
                    {
                        index = i;
                        smallest = a;
                    }
                }
                if (index > 0 && ((n > max_points) || (smallest < area * min_area)))
                {
                    --n;
                    points.erase(points.begin() + index);
                }
                else
                {
                    break;
                }
            }

            // remove the additional point
            if (!is_closed)
                points.resize(n - 1);
        }

        template <typename T>
        static double hull_area(std::vector<T> &points)
        {
            double area = 0;
            size_t n = points.size();
            // Shoelace algo to calculate total area of convex hull
            for (size_t i = 0; i + 1 < n; i++)
                area -= (double)points[i].x * points[i + 1].y - (double)points[i].y * points[i + 1].x;
            return area;
        }

    private:
        // Disallow creating an instance of this object
        convex_hull() {};

        template <typename T>
        static bool is_non_right_turn(T &o, T &a, T &b)
        {
            return (a.x - o.x) * (b.y - o.y) >= (a.y - o.y) * (b.x - o.x);
        }

        template <typename T>
        static void add_to_hull(T*& hull_ptr, T *hull_limit, T *hull_limit_end, T &try_me)
        {
            while (hull_ptr >= hull_limit && is_non_right_turn(hull_ptr[-2], hull_ptr[-1], try_me))
                --hull_ptr;
            if (hull_ptr != hull_limit_end)
                *hull_ptr++ = try_me;
        }
    };
}
