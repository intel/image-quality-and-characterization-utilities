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

#include "Teisko/Algorithm/DelaunayTriangulation.hpp"
#include "catch.hpp"


/// \page libcalc_specs_delaunay_triangulation Specs: Teisko triangulation
///
/// \snippet this snippet-specs-triangulation

using namespace Teisko;

namespace libcalc_triangulation_tests
{
    SCENARIO("Delaunay triangulation gives correct amount of triangles")
    {
        //    +---+---+
        // N  | \ | / |
        //    +---+---+
        //      M
        GIVEN("A regular grid of N*M items")
        {
            auto const x = 70;
            auto const y = 50;
            std::vector<point_xy> point_cloud;
            for (auto i = 0; i < x; i++)
                for (auto j = 0; j < y; j++)
                    point_cloud.push_back(point_xy(i, j));
            WHEN("The point cloud is triangulated")
            {
                auto t = triangulate(point_cloud);
                THEN("The result is a vector of (N-1)*(M-1) * 2 triangles")
                {
                    CHECK(t.size() == (x - 1)*(y - 1) * 2);
                }
            }
        }
    }

    bool is_same_triangulation(const std::vector<triangle_vertex> &lhs, const std::vector<triangle_vertex> &rhs)
    {
        if (lhs.size() != rhs.size())
            return false;

        for (auto &l: lhs)
        {
            bool found = false;
            for (auto &r : rhs)
            {
                if (l.a == r.a && l.b == r.b && l.c == r.c) { found = true; break; }
                if (l.a == r.b && l.b == r.c && l.c == r.a) { found = true; break; }
                if (l.a == r.c && l.b == r.a && l.c == r.b) { found = true; break; }
            }
            if (!found)
                return found;
        }
        return true;
    }

    SCENARIO("Delaunay triangulation of 10 points gives equal result with matlab")
    {
        GIVEN("A set of 10 points")
        {
            double x[10] = { 0.4868, 0.4359, 0.4468, 0.3063, 0.5085, 0.5108, 0.8176, 0.7948, 0.6443, 0.3786 };
            double y[10] = { 0.8116, 0.5328, 0.3507, 0.9390, 0.8759, 0.5502, 0.6225, 0.5870, 0.2077, 0.3012 };
            std::vector<point_xy> p;
            for (int i = 0; i < 10; i++)
                p.emplace_back(x[i], y[i]);

            WHEN("The set is triangulated")
            {
                auto t = triangulate(p);
                THEN("The set contains all the same triangles as in the matlab reference result")
                {
                    // Any vertical permutation and horizontal rotation is ok
                    auto indices_matlab = std::vector<triangle_vertex>({
                        { 0, 4, 3 },
                        { 9, 2, 1 },
                        { 9, 1, 3 },
                        { 1, 0, 3 },
                        { 5, 2, 8 },
                        { 1, 5, 0 },
                        { 0, 6, 4 },
                        { 0, 5, 7 },
                        { 1, 2, 5 },
                        { 7, 8, 6 },
                        { 2, 9, 8 },
                        { 0, 7, 6 },
                        { 5, 8, 7 }
                    });

                    CHECK(is_same_triangulation(indices_matlab, t) == true);
                }
            }
        }
    }
    SCENARIO("Scattered interpolant locates equal coordinates to matlab")
    {
        auto find_smallest_vertex = [](triangle_vertex &v)
        {
            if (v.a < v.b && v.a < v.c)
                return 0;
            if (v.b < v.a && v.b < v.c)
                return 1;
            return 2;
        };
        // i = x + y*w -> return y + x * h
        auto transpose = [](size_t i, int w, int h)
        {
            return (i / w) + (i % w) * h;
        };

        GIVEN("A set of 10 points and points where the data is interpolated")
        {
            auto p = std::vector<point_xy>({
                { 0.4868, 0.8116 }, { 0.4359, 0.5328 }, { 0.4468, 0.3507 },
                { 0.3063, 0.9390 }, { 0.5085, 0.8759 }, { 0.5108, 0.5502 },
                { 0.8176, 0.6225 }, { 0.7948, 0.5870 }, { 0.6443, 0.2077 },
                { 0.3786, 0.3012 }
            });
            std::vector<double> data({ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 });

            auto x_vec = std::vector<double>({ 0.4, 0.5, 0.6, 0.7 });
            auto y_vec = std::vector<double>({ 0.4, 0.5, 0.6, 0.7 });
            const int w = 4;
            const int h = 4;
            REQUIRE(x_vec.size() == w);
            REQUIRE(y_vec.size() == h);

            WHEN("The input point cloud is interpolated in the chosen x/y grid")
            {
                auto weights = scattered_interpolation(p, x_vec, y_vec);
                THEN("The barycentric coordinates match to matlab generated reference")
                {
                    auto const e = 1e-9; // the required accuracy
                    // x = [0.4868 0.4359 0.4468 0.3063 0.5085 0.5108 0.8176 0.7948 0.6443 0.3786];
                    // y = [0.8116 0.5328 0.3507 0.9390 0.8759 0.5502 0.6225 0.5870 0.2077 0.3012];
                    // DT = DT=delaunayTriangulation(x',y');
                    // [xx,yy]=meshgrid([0.4 0.5 0.6 0.7],[0.4 0.5 0.6 0.7]);
                    // [t,bc] = pointLocation(DT, xx(:), yy(:));
                    // v = DT(t,:)-1; %% convert triangle index to vertices
                    // reference vectors are in column major format
                    auto ref_vertices = std::vector<triangle_vertex>({
                        // x=0.4        x=0.5       x=0.6         x=0.7
                        { 9, 1, 3 }, { 9, 1, 3 }, { 9, 1, 3 }, { 1, 0, 3 },     // y = 0.4
                        { 5, 2, 8 }, { 5, 2, 8 }, { 1, 5, 0 }, { 0, 5, 7 },     // y = 0.5
                        { 5, 8, 7 }, { 5, 8, 7 }, { 0, 5, 7 }, { 0, 5, 7 },     // y = 0.6
                        { 5, 8, 7 }, { 5, 8, 7 }, { 0, 5, 7 }, { 0, 6, 4 }      // y = 0.7
                    });

                    double bc[][3] = {
                        { 0.596605181174473, 0.390165473773809, 0.013229345051718 },
                        { 0.353410412564162, 0.525836629410579, 0.120752958025258 },
                        { 0.110215643953852, 0.661507785047350, 0.228276570998799 },
                        { 0.549255158510958, 0.124745326631467, 0.325999514857575 },
                        { 0.357223254879951, 0.489168078346970, 0.153608666773079 },
                        { 0.763993141550771, 0.214212230901124, 0.021794627548104 },
                        { 0.081410409642507, 0.722657910450149, 0.195931679907343 },
                        { 0.571621175493339, 0.418100978690323, 0.010277845816338 },
                        { 0.447669666519219, 0.449580164176358, 0.102750169304423 },
                        { 0.594954728193003, 0.171646891649084, 0.233398380157913 },
                        { 0.144575670120659, 0.529122160573370, 0.326302169305971 },
                        { 0.522633411784752, 0.119115877360199, 0.358250710855049 },
                        { 0.076472165569939, 0.485594053010878, 0.437933781419182 },
                        { 0.223757227243724, 0.207660780483604, 0.568581992272672 },
                        { 0.095587906412072, 0.230137059243246, 0.674275034344682 },
                        { 0.230338564296950, 0.635711248286133, 0.133950187416916 }
                    };

                    CHECK(weights.size() == x_vec.size() * y_vec.size());
                    REQUIRE(ref_vertices.size() == weights.size());

                    for (size_t i = 0; i < ref_vertices.size(); i++)
                    {
                        auto offset_1 = find_smallest_vertex(weights[i].v);
                        auto offset_2 = find_smallest_vertex(ref_vertices[transpose(i, w, h)]);
                        CHECK(weights[i].a == Approx(bc[transpose(i, w, h)][(3 + offset_2 - offset_1) % 3]).epsilon(e));
                        CHECK(weights[i].b == Approx(bc[transpose(i, w, h)][(4 + offset_2 - offset_1) % 3]).epsilon(e));
                        CHECK(weights[i].c == Approx(bc[transpose(i, w, h)][(5 + offset_2 - offset_1) % 3]).epsilon(e));
                    }

                    AND_WHEN("The `data` is interpolated according to the given weights")
                    {
                        auto output = barycentric_interpolation(data, weights);
                        THEN("The output matches the matlab reference with high accuracy")
                        {
                            // V = 1:10;               %% the data
                            // triVals=V(DT(t,:));     %% Data points of the indices
                            // v=dot(bc', V(DT(t,:))')'; reshape(v,[4 4])'   %% row major form
                            double vals[w * h] = {
                                6.799300139499221, 4.993321765278329, 7.554240831137921, 8.332649721870999,
                                5.068789216563814, 5.422747189940941, 6.981737435263077, 7.760146325996155,
                                3.338278293628409, 4.694699961893255, 5.929725988008647, 6.870610536629003,
                                2.527253703083682, 3.162449814165984, 4.103334362786339, 5.350068239384465
                            };
                            REQUIRE(output.size() == w*h);
                            for (size_t j = 0; j < output.size(); j++)
                            {
                                CHECK(output[j] == Approx(vals[j]).epsilon(e));
                            }
                        }
                    }
                }
            }
        }
    }
    SCENARIO("A Connectivity matrix can be found for a triangulation")
    {
        auto vertex_to_vec = [](triangle_vertex v)
        {
            return std::vector<int>({ v.a, v.b, v.c });
        };

        GIVEN("A triangulation of three triangles, where the middle one is connected to both of the other triangles")
        {
            // input:
            //            A----B-----E
            //             \0/ 1\ 2 /
            //              C-----D
            // input =  [ A B C], [B D C], [D B E]
            // output:  [-1 1 -1], [2 -1 0], [1 -1 -1]
            std::vector<triangle_vertex> triangulation;
            const int A = 889;
            const int B = 13;
            const int C = 15435;
            const int D = -1313;
            const int E = 444;
            triangulation.emplace_back(A, B, C);
            triangulation.emplace_back(B, D, C);
            triangulation.emplace_back(D, B, E);
            WHEN("A connectivity of the triangulation is found")
            {
                auto conn = calculate_connectivity(triangulation);
                THEN("The connectivity matrix equals the predefined one")
                {
                    auto reference = std::vector<triangle_vertex>({
                        { -1, 1, -1 },
                        { 2, -1, 0 },
                        { 1, -1, -1 }
                    });
                    REQUIRE(conn.size() == reference.size());
                    for (size_t i = 0; i < conn.size(); i++)
                    {
                        CHECK(vertex_to_vec(conn[i]) == vertex_to_vec(reference[i]));
                    }
                }
            }
        }
        GIVEN("A triangulation with disconnected triangles")
        {
            std::vector<triangle_vertex> triangulation;
            triangulation.emplace_back(0, 1, 2);
            triangulation.emplace_back(200, 300, 400);
            triangulation.emplace_back(-4, -5, -6);
            triangulation.emplace_back(0, 300, -6);
            WHEN("A connectivity of the triangulation is found")
            {
                auto conn = calculate_connectivity(triangulation);
                THEN("The result is all negative ones")
                {
                    REQUIRE(conn.size() == triangulation.size());
                    for (auto &c : conn)
                    {
                        auto reference = std::vector<int>{-1, -1, -1};
                        CHECK(vertex_to_vec(c) == reference);
                    }
                }
            }
        }
    }

    SCENARIO("Test case for vs2017 oddball bug for LCA")
    {
        GIVEN("5 vertices and the point { 528, 0 }")
        {
            auto vertices_sorted = std::vector<dt_triangulator::point_xyi>({
                { 200.0, 100.0, 0 },
                { 200.0, 500.0, 1 },
                { 500.0, 300.0, 2 },
                { 800.0, 100.0, 3 },
                { 800.0, 500.0, 4 }});
            auto pt = point_xy(528.0, 0.0);
            WHEN("The nearest vertex is located given any index i as a hint for starting value")
            {
                auto resulting_idx = std::vector<interpolation_weights>();
                for (auto i = 0; i < 5; i++)
                    resulting_idx.emplace_back(dt_triangulator::find_nearest_neighbor(vertices_sorted, pt, i));
                THEN("All the best instances found match to the point #3 at { 800, 100 }")
                {
                    CHECK(resulting_idx.size() == 5);
                    for (auto &w : resulting_idx)
                    {
                        CHECK(w.a == 1.0);
                        CHECK(w.b == 0.0);
                        CHECK(w.c == 0.0);
                        CHECK(w.v.a == 3);
                        CHECK(w.v.b == 3);
                        CHECK(w.v.c == 3);
                    }
                }
            }
        }
    }
}