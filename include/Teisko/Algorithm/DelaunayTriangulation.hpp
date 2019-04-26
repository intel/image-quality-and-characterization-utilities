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
    struct triangle_vertex
    {
        int a = 0;
        int b = 0;
        int c = 0;
        triangle_vertex(int v0 = 0, int v1 = 0, int v2 = 0) : a(v0), b(v1), c(v2) { }
    };

    // holds barycentric coordinates for a triangle
    struct interpolation_weights
    {
        double a = 0.0;
        double b = 0.0;
        double c = 0.0;
        triangle_vertex v;
    };

    namespace dt_triangulator
    {
        inline void cyclical_sort_indices(triangle_vertex &vertices)
        {
            if (vertices.a < vertices.b)
            {
                if (vertices.a < vertices.c)
                    return;
            }
            else  // b < a
            {
                if (vertices.b < vertices.c)
                {
                    vertices = triangle_vertex(vertices.b, vertices.c, vertices.a);
                    return;
                }
            }
            vertices = triangle_vertex(vertices.c, vertices.a, vertices.b);
            return;
        }

        // DT specific point containing copy of original data and index to original data
        //  - used to sort input data in x-direction
        struct point_xyi
        {
            point_xy p;
            int index = 0;
            point_xyi() : p(), index(0) { }
            point_xyi(point_xy xy, int i) : p(xy), index(i) { }
            point_xyi(double x, double y, int i) : p(x, y), index(i) { }
            bool operator <(const point_xyi &other) const { return p.x < other.p.x; }
        };

        inline void add_bounding_triangle(std::vector<point_xyi> &v)
        {
            int max_vertex = static_cast<int>(v.size());
            if (max_vertex < 3)
                throw std::runtime_error("There has to be at least three points in the graph");

            // Find a square containing all points
            // we know that v is sorted, so the minimum and maximum x are in beginning and end of 'v'
            // therefore we only need to compare against 'y' in the inner loop
            auto min_xy = v.front().p.x;
            auto max_xy = v.back().p.x;
            for (auto &p : v)
            {
                min_xy = std::min(min_xy, p.p.y);
                max_xy = std::max(max_xy, p.p.y);
            }
            double r = max_xy - min_xy;
            v.emplace_back(-4.5 * r + min_xy, -2.0 * r + min_xy, max_vertex);
            v.emplace_back(5.5 * r + min_xy, -2.0 * r + min_xy, max_vertex + 1);
            v.emplace_back(min_xy, 3.0 * r + min_xy, max_vertex + 2);
        }

        // Sorts the input data by x-coordinate
        // - associates each point with index to original data
        inline std::vector<point_xyi> sort_input(const std::vector<point_xy> &input)
        {
            std::vector<point_xyi> vertices;
            vertices.reserve(input.size());
            int i = 0;
            for (auto &v : input)
                vertices.emplace_back(v, i++);

            std::sort(vertices.begin(), vertices.end());
            return vertices;
        }

        // Triangle as used by DT stores the circumcircle center and radius
        struct triangle
        {
            point_xy center{};      // Center of the triangle circum circle
            double radius2 = 0;     // Radius of the triangle circum circle (squared)
            double max_x = 0;       // Maximum X-coordinate of the triangle circum circle

            // Copy of the coordinate data with index to original location
            triangle_vertex vertices;

            triangle() = default;

            // Generates a triangle from three vertices (in ccw order)
            triangle(int indexa, int indexb, int indexc, std::vector<point_xyi> &points)
                : center(0, 0)
                , radius2(-1.0)
                , max_x(-std::numeric_limits<double>::infinity())
                , vertices(indexa, indexb, indexc)
            {
                auto &A = points[indexa].p;
                auto &B = points[indexb].p;
                auto &C = points[indexc].p;

                auto ba = B - A;
                auto ca = C - A;

                double det = ba.x * ca.y - ba.y * ca.x;

                // Degenerate or CW triangles will have negative or close to zero determinant
                if (det > 1e-15)
                {
                    det = 0.5 / det;
                    auto za = det * norm2(A);
                    auto zb = det * norm2(B) - za;
                    auto zc = det * norm2(C) - za;

                    center = { zb * ca.y - ba.y * zc, ba.x * zc - zb * ca.x };
                    radius2 = norm2(center - A);
                    max_x = center.x + std::sqrt(radius2);
                }
                else
                {
                    throw std::runtime_error("degenerate triangle");
                }
            }

            bool is_non_boundary_triangle(int max) const
            {
                return (vertices.a < max && vertices.b < max && vertices.c < max);
            }

            bool is_inside_circumcircle(point_xy p) const { return norm2(p - center) <= radius2; }
            bool is_unreachable(point_xy p) const { return max_x < p.x; }
        };

        // Given a triangle a,b,c (in CCW order) fills the barycentric
        // interpolation coordinate weights in `w`
        // - returns true if the point p is inside the triangle, false otherwise
        //   in which case `w` is left in unspecified state
        inline bool find_barycentric_coordinates(
            point_xy p, interpolation_weights &w,
            point_xy a, point_xy b, point_xy c)
        {
            w.a = (b.y - c.y) * (p.x - b.x) + (c.x - b.x) * (p.y - b.y);
            w.b = (c.y - a.y) * (p.x - c.x) + (a.x - c.x) * (p.y - c.y);
            w.c = (a.y - b.y) * (p.x - a.x) + (b.x - a.x) * (p.y - a.y);

            double sum = w.a + w.b + w.c;
            if (w.a >= 0.0 && w.b >= 0.0 && w.c >= 0.0 && sum > 0.0)
            {
                sum = 1.0 / sum;
                w.a *= sum;
                w.b *= sum;
                w.c *= sum;
                return true;
            }
            return false;
        }

        inline point_xy center(const triangle &t, const std::vector<point_xyi> &points)
        {
            auto &a = points[t.vertices.a].p;
            auto &b = points[t.vertices.b].p;
            auto &c = points[t.vertices.c].p;
            auto mean = (a + b + c) * (1.0 / 3.0);
            return mean;
        }

        // Organizer of triangle edges -- when a new point conflict more than one
        //    already added triangle in the incremental triangulation
        //    all conflicting triangles must be compared against shared edges
        //    - one triangle contains edge AB and the other contains edge BA
        //      - these edges must be excluded from producing further triangles
        struct edge_buffer
        {
            // Edges a,b are encoded as [smaller(a,b) * 2][larger(a,b) * 2][direction:1]
            bool are_same_edges(uint64_t &a, uint64_t &b)
            {
                return (a ^ b) <= 1;
            }

            std::vector<uint64_t> edges;
            void insert(int vertex_a, int vertex_b)
            {
                uint32_t a = static_cast<uint32_t>(vertex_a)* 2;
                uint32_t b = static_cast<uint32_t>(vertex_b)* 2;
                if (a > b)
                {
                    std::swap(a, b);
                    b |= 1;
                }
                edges.emplace_back(((uint64_t)a << 32) | (uint64_t)b);
            }

            // After all triangle candidates are added, we remove shared edges
            // and add the remaining items to triangles with the common point points[vertex]
            void add_as_triangles(std::vector<triangle> &triangles, std::vector<point_xyi> &points, int vertex)
            {
                std::sort(edges.begin(), edges.end());
                auto len = edges.size();
                for (size_t i = 0; i < len; i++)
                {
                    // skip same/opposite edges
                    if (i + 1 < len && are_same_edges(edges[i], edges[i + 1]))
                        i++;
                    else
                    {
                        int edge_1 = static_cast<int>(edges[i] >> 33);
                        int edge_2 = static_cast<int>(edges[i] >> 1);
                        if (edges[i] & 1)
                            std::swap(edge_1, edge_2);
                        triangles.emplace_back(edge_1, edge_2, vertex, points);
                    }
                }
            }
            void reset() { edges.clear(); }
            edge_buffer() { edges.reserve(100); }
        };

        // returns a list of triangles formed by Delaunay triangulation of given two dimensional set of coordinates
        // Original algorithm from published paper
        // "Efficient Triangulation Algorithm Suitable for Terrain Modeling" by Paul Bourke, 1989
        // precondition: receives data sorted in x dimension and attached with 0 <= index < data.size()
        inline std::vector<triangle> triangulate_internal(std::vector<point_xyi> &sorted_data)
        {
            // Step 0 -- check that there is enough data
            auto triangles = std::vector<triangle>();
            auto data_len = sorted_data.size();
            if (data_len < 3)
                return triangles;

            add_bounding_triangle(sorted_data);

            triangles.reserve(sorted_data.size() * 3);
            triangles.emplace_back((int)data_len, (int)data_len + 1, (int)data_len + 2, sorted_data);

            // Add data points in order of increasing x coordinate
            //  - split all existing triangles that have the new point contained in their circumcircle
            //  - optimize out all triangles from the set of triangles to be inspected
            //    once their maximum x-coordinate (of the circum circle) is smaller than the rest of the data points

            edge_buffer edge_buffer;
            size_t finalized_triangles = 0;
            for (size_t i = 0; i < data_len; i++)
            {
                auto &p = sorted_data[i];
                edge_buffer.reset();
                auto triangles_size = triangles.size();
                for (size_t j = finalized_triangles; j < triangles_size; j++)
                {
                    auto &triangle = triangles[j];
                    if (triangle.is_unreachable(p.p))
                    {
                        if (j != finalized_triangles++)
                            std::swap(triangles[finalized_triangles - 1], triangle);
                    }
                    else if (triangle.is_inside_circumcircle(p.p))
                    {
                        edge_buffer.insert(triangle.vertices.a, triangle.vertices.b);
                        edge_buffer.insert(triangle.vertices.b, triangle.vertices.c);
                        edge_buffer.insert(triangle.vertices.c, triangle.vertices.a);
                        // move last item to current (and decrement j so we can process the same
                        // index again)
                        if (j != --triangles_size)
                            std::swap(triangles[triangles_size], triangles[j--]);
                        triangles.pop_back();
                    }
                }
                edge_buffer.add_as_triangles(triangles, sorted_data, static_cast<int>(i));
            }

            auto max_vertex = static_cast<int>(data_len);
            triangles.erase(remove_if(triangles.begin(), triangles.end(), [max_vertex](const triangle &t)
                {
                    return t.vertices.a >= max_vertex || t.vertices.b >= max_vertex || t.vertices.c >= max_vertex;
                }
            ), triangles.end());

            // Remove also the original bounding triangle
            for (int i = 0; i < 3; i++)
                sorted_data.pop_back();

            return triangles;
        }

        inline std::vector<triangle_vertex> extract_vertices(std::vector<triangle> &triangles)
        {
            auto output = std::vector<triangle_vertex>();
            output.reserve(triangles.size());

            // Convert the indices from sorted order to data order
            for (auto &t : triangles)
            {
                auto &v = t.vertices;
                dt_triangulator::cyclical_sort_indices(v);  // return also in "canonical" form
                output.emplace_back(v);
            }
            return output;
        }

        // Convert vertex index from sorted order back to original order stored in 'mapping'
        inline void remap_indices(std::vector<triangle_vertex> &triangles, std::vector<point_xyi> &mapping)
        {
            for (auto &t : triangles)
            {
                t.a = mapping[t.a].index;
                t.b = mapping[t.b].index;
                t.c = mapping[t.c].index;
            }
        }

        // Convert vertex index from sorted order back to original order stored in 'mapping'
        inline void remap_indices(std::vector<interpolation_weights> &weight, std::vector<point_xyi> &mapping)
        {
            for (auto &t : weight)
            {
                t.v.a = mapping[t.v.a].index;
                t.v.b = mapping[t.v.b].index;
                t.v.c = mapping[t.v.c].index;
            }
        }

        // Finds the closest element in the sorted vector of vertices to a given point `p`
        // Scans left and right from the initial `index` until no further improvements are possible
        inline interpolation_weights find_nearest_neighbor(const std::vector<point_xyi> &vertices, point_xy p, size_t index = 0)
        {
            auto max_index = vertices.size();
            if (index >= max_index)
                throw std::runtime_error("Vertex index out of range");

            double best_dist = norm2(p - vertices[index].p);
            // increments first argument to left / right testing if vertices[idx] either improves
            // current best distance, or is horizontally too far from best distance so that the
            // iteration to left/right can be stopped
            auto trial_advance = [&](size_t &idx, bool increment, size_t &best_index, size_t max)
            {
                auto idx2 = idx;
                if (increment)
                    idx2++;
                else
                    idx2--;

                if (idx2 < max)
                {
                    idx = idx2;
                    auto horizontal_distance = (vertices[idx].p.x - p.x) * (vertices[idx].p.x - p.x);
                    auto vertical_distance = (vertices[idx].p.y - p.y) * (vertices[idx].p.y - p.y);

                    if (horizontal_distance > best_dist)
                        idx = max;
                    else
                    {
                        // Select the closest vertex if improvement found
                        double dist2 = horizontal_distance + vertical_distance;
                        if (dist2 < best_dist)
                        {
                            best_index = idx;
                            best_dist = dist2;
                        }
                    }
                }
                else
                {
                    idx = max;
                }
            };
            auto left_idx = index;
            auto right_idx = left_idx;

            while (left_idx < max_index || right_idx < max_index)
            {
                trial_advance(left_idx, false, index, max_index);
                trial_advance(right_idx, true, index, max_index);
            }

            interpolation_weights w;
            w.a = 1.0;
            w.b = 0.0;
            w.c = 0.0;
            w.v.a = static_cast<int>(index);
            w.v.b = static_cast<int>(index);
            w.v.c = static_cast<int>(index);

            return w;
        }

        // Incrementally advances to next neighbor towards the target
        // until no further advances can be made
        // Casts a line from gravity center of triangle towards target
        // finding out the intersecting neighbor
        //  - the connectivity chart must be built so, that
        //  item.a points to neighbor between vertexes a and b (or -1)
        //  item.b points to neighbor between vertexes b and c (or -1)
        //  item.c points to neighbor between vertexes c and a (or -1)
        inline interpolation_weights find_closest_triangle(point_xy p,
            const std::vector<triangle> &triangles,
            const std::vector<point_xyi> &vertices,
            const std::vector<triangle_vertex> &connectivity,
            int &triangle_hint)
        {
            if (triangle_hint < 0 ||
                static_cast<size_t>(triangle_hint) >= triangles.size() ||
                static_cast<size_t>(triangle_hint) >= connectivity.size())
                throw std::runtime_error("Triangle list must contain at least one triangle");

            // sign of cross product tells if a point pt lies on right/left side of the line
            auto sign_of_cross_product = [](point_xy line, point_xy pt) { return line.x * pt.y >= line.y * pt.x; };
            //
            auto set_neighbor = [](int index, int &dst)
            {
                if (index < 0)
                    return false;
                dst = index;
                return true;
            };

            bool found_improvement = false;
            do
            {
                auto &t = triangles[triangle_hint];
                auto &a = vertices[t.vertices.a].p;
                auto &b = vertices[t.vertices.b].p;
                auto &c = vertices[t.vertices.c].p;

                if (t.is_inside_circumcircle(p))    // A quick rejection
                {
                    interpolation_weights w;
                    if (find_barycentric_coordinates(p, w, a, b, c))
                    {
                        w.v = t.vertices;
                        return w;
                    }
                }
                auto mid_point = (a + b + c) * (1.0 / 3.0);
                auto sa = sign_of_cross_product(a - mid_point, p - mid_point);
                auto sb = sign_of_cross_product(b - mid_point, p - mid_point);
                auto sc = sign_of_cross_product(c - mid_point, p - mid_point);
                if (sa == true && sb == false)
                    found_improvement = set_neighbor(connectivity[triangle_hint].a, triangle_hint);
                else if (sb == true && sc == false)
                    found_improvement = set_neighbor(connectivity[triangle_hint].b, triangle_hint);
                else // (sc == true && sa == false)
                    found_improvement = set_neighbor(connectivity[triangle_hint].c, triangle_hint);
            } while (found_improvement);

            // else -- extrapolate or find NN
            auto closest_idx = static_cast<size_t>(triangles[triangle_hint].vertices.a);
            return dt_triangulator::find_nearest_neighbor(vertices, p, closest_idx);
        }

        // Zig-zag scan the input data in a continuous space filling curve
        // so that successive sampling points are typically found in nearby triangles
        inline std::vector<interpolation_weights> zigzag_scan(
            std::vector<point_xyi> &vertices_sorted,
            std::vector<triangle> &triangles,
            std::vector<triangle_vertex> &connectivity,
            std::vector<double> &x_points,
            std::vector<double> &y_points)
        {
            auto output = std::vector<interpolation_weights>(x_points.size() * y_points.size());

            auto x_size = static_cast<int>(x_points.size());
            int triangle_hint = 0;
            bool left_to_right = true;
            auto *dst = &output[0];
            for (auto y : y_points)
            {
                int start_point = left_to_right ? 0 : x_size - 1;
                int increment_x = left_to_right ? 1 : -1;
                for (int i = 0; i < x_size; i++, start_point += increment_x)
                {
                    point_xy p(x_points[start_point], y);
                    dst[start_point] = triangles.size() == 0
                        ? dt_triangulator::find_nearest_neighbor(vertices_sorted, p)
                        : find_closest_triangle(p, triangles, vertices_sorted, connectivity, triangle_hint);
                }
                dst += x_size;
                left_to_right = !left_to_right;
            }
            return output;
        }
    };

    // Returns connectivity matrix of triangulation
    // input:
    //            A----B-----E
    //             \0/ 1\ 2 /
    //              C-----D
    // input =  [ A B C], [B D C], [D B E]
    // output:  [-1 1 -1], [2 -1 0], [1 -1 -1]
    // where -1 or a value outside [0 input.size()[  means that the triangle is not connected
    // at that edge to another triangle
    inline std::vector<triangle_vertex> calculate_connectivity(std::vector<triangle_vertex> &vertices)
    {
        // associates an edge (in sorted order) with triangle index
        auto add_edge = [](std::vector<triangle_vertex> &edges, int a, int b, int i)
        {
            if (a > b)
                std::swap(a, b);
            edges.emplace_back(a, b, i);
        };

        // inserts an index to a triangle while shifting the values
        auto push_index = [](std::vector<triangle_vertex> &v, int index_a, int index_b)
        {
            switch (index_a & 3)
            {
            case 0: v[index_a >> 2].a = index_b >> 2; break;
            case 1: v[index_a >> 2].b = index_b >> 2; break;
            default:
                v[index_a >> 2].c = index_b >> 2; break;
            }
        };

        int len = static_cast<int>(vertices.size());
        auto result = std::vector<triangle_vertex>(len, triangle_vertex{ -1, -1, -1 });

        // Add all 3*N edges (in sorted order to a buffer)
        // also mark in the triangle index field
        auto edges = std::vector<triangle_vertex>();
        edges.reserve(3 * vertices.size() + 1);
        for (int i = 0; i < len; i++)
        {
            auto &t = vertices[i];
            add_edge(edges, t.a, t.b, i * 4 + 0);
            add_edge(edges, t.b, t.c, i * 4 + 1);
            add_edge(edges, t.c, t.a, i * 4 + 2);
        }

        // Sort the edges to allow easy inspection of duplicates
        std::sort(edges.begin(), edges.end(), [](const triangle_vertex &a, const triangle_vertex &b)
        {
            return (a.a == b.a) ? a.b < b.b : a.a < b.a;
        });

        for (int i = 0; i < len * 3 - 1; i++)
        {
            if (edges[i].a == edges[i + 1].a && edges[i].b == edges[i + 1].b)
            {
                int index_a = edges[i].c;
                int index_b = edges[i+1].c;
                push_index(result, index_a, index_b);
                push_index(result, index_b, index_a);
                i++;
            }
        }
        return result;
    }

    inline std::vector<triangle_vertex> triangulate(const std::vector<point_xy> &data)
    {
        // take a copy of input data to be sorted
        auto vertices_sorted = dt_triangulator::sort_input(data);

        auto triangles = dt_triangulator::triangulate_internal(vertices_sorted);

        auto output = dt_triangulator::extract_vertices(triangles);

        dt_triangulator::remap_indices(output, vertices_sorted);

        return output;
    }

    // Interpolate the point cloud 'data' in orthogonal grid formed from x and y point vector
    // The result has W*H tuples where each record contains
    // three indices to original data and three weights that sum up to 1.0
    // - in case the point is not contained within a triangle, the indexes are for the closest point
    inline std::vector<interpolation_weights> scattered_interpolation(
        const std::vector<point_xy> &data,
        std::vector<double> &x_points,
        std::vector<double> &y_points)
    {
        // take a copy of input data to be sorted -- then triangulate
        auto vertices_sorted = dt_triangulator::sort_input(data);

        // calculate the full triangulation with all the circum circle data as well
        auto triangles = dt_triangulator::triangulate_internal(vertices_sorted);

        // find out the connectivity matrix
        auto vertices_only = dt_triangulator::extract_vertices(triangles);
        auto connectivity = calculate_connectivity(vertices_only);

        // "Walk" in the delaunay tessellation trying to find the next point in nearby cells
        auto output = dt_triangulator::zigzag_scan(vertices_sorted, triangles, connectivity, x_points, y_points);
        dt_triangulator::remap_indices(output, vertices_sorted);

        return output;
    }

    // Given an array of barycentric weights and indexes to `data` to be interpolated
    // return a vector of interpolated data based on the weights
    template <typename T>
    std::vector<T> barycentric_interpolation(std::vector<T> &data, std::vector<interpolation_weights> &weights)
    {
        auto output = std::vector<T>();
        output.reserve(weights.size());
        for (auto &w : weights)
            output.push_back(data[w.v.a] * w.a + data[w.v.b] * w.b + data[w.v.c] * w.c);
        return output;
    }

}
