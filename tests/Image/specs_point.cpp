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


#include "Teisko/Image/Point.hpp"
#include "catch.hpp"

#include <vector>
using namespace Teisko;

SCENARIO("range_s allows ranged for syntax starting from 0 and ending at N-1")
{
    GIVEN("A range of 5 integers")
    {
        auto width = range(5);
        WHEN("The range is iterated placing each value to a vector")
        {
            auto result = std::vector<int>();

            for (auto &i : width)
                result.emplace_back(i);

            THEN("The result equals to indices 0,1,2,3,4")
            {
                auto expected_result = std::vector<int>{0, 1, 2, 3, 4};
                CHECK(result == expected_result);
            }
        }
    }

    GIVEN("Two ranges representing a Cartesian space of 4x3")
    {
        auto const y_max = 4;
        auto const x_max = 3;
        auto const rows = range(y_max);
        auto const cols = range(x_max);
        WHEN("The rows and columns are iterated in two nested loops adding each coordinate to a vector")
        {
            auto x_coordinates_ranged_for = std::vector<int>();
            auto y_coordinates_ranged_for = std::vector<int>();

            for (auto y : rows)
            {
                for (auto x : cols)
                {
                    x_coordinates_ranged_for.emplace_back(x);
                    y_coordinates_ranged_for.emplace_back(y);
                }
            }

            THEN("The vectors contain the same values as when looping with the standard notation")
            {
                auto x_coordinates_standard = std::vector<int>();
                auto y_coordinates_standard = std::vector<int>();
                for (auto y = 0; y < y_max; ++y)
                {
                    for (auto x = 0; x < x_max; ++x)
                    {
                        x_coordinates_standard.emplace_back(x);
                        y_coordinates_standard.emplace_back(y);
                    }
                }
                CHECK(x_coordinates_ranged_for == x_coordinates_standard);
                CHECK(y_coordinates_ranged_for == y_coordinates_standard);
            }
        }
    }
}

SCENARIO("Roi_point is an utility class to contain arbitrary integers x and y. "
    "Roi point can be constructed with zero, one or two parameter variants.")
{
    WHEN("A roi point is constructed with default ctor -- i.e. zero parameters")
    {
        roi_point p;
        THEN("The two members of the point p are zeros")
        {
            CHECK(p._x == 0);
            CHECK(p._y == 0);
        }
    }

    WHEN("A roi point is constructed with a single value e.g. -7")
    {
        const int minus_seven = -7;
        roi_point p{ minus_seven };
        THEN("Both coordinates x and y equal to -7")
        {
            CHECK(p._x == minus_seven);
            CHECK(p._y == minus_seven);
        }
    }

    WHEN("A roi point is constructed with pair 5, -3")
    {
        roi_point p{ 5, -3 };
        THEN("The x member of the point is 5 and the y member is -3")
        {
            CHECK(p._x == 5);
            CHECK(p._y == -3);
        }
    }

    WHEN("A roi point is constructed with floating point types")
    {
        roi_point p{ -1.5f, 2.5f };
        THEN("The point is constructed by rounding away from zero")
        {
            CHECK(p._x == -2);
            CHECK(p._y == 3);
        }
    }
}

SCENARIO("Roi_point allows ranged for iteration in row minor order")
{
    GIVEN("A range of 2 rows, 3 columns")
    {
        roi_point range(3, 2);
        WHEN("The range is iterated and each coordinate is added to vector")
        {
            auto coordinates = std::vector<int>();
            for (auto &&p : range)
            {
                coordinates.emplace_back(p._x);
                coordinates.emplace_back(p._y);
            }
            THEN("The result matches the sequence 0,0, 1,0, 2,0, 0,1, 1,1, 2,1")
            {
                auto expected_coordinates = std::vector<int>{0, 0, 1, 0, 2, 0, 0, 1, 1, 1, 2, 1};
                CHECK(coordinates == expected_coordinates);
            }
        }
    }
}

SCENARIO("Roi_point class implements basic mathematical operators with self")
{
    GIVEN("Two roi_points")
    {
        roi_point a(1, 2);
        roi_point b(3, 4);

        WHEN("The points are added")
        {
            auto c = a + b;
            THEN("The sum is the pairwise sum 1+3, 2+4")
            {
                CHECK(c._x == (1 + 3));
                CHECK(c._y == (2 + 4));
            }
        }

        WHEN("The points are subtracted")
        {
            auto d = a - b;
            THEN("The sum is the pairwise difference 1-3, 2-4")
            {
                CHECK(d == roi_point(1 - 3, 2 - 4));
            }
        }

        WHEN("The points are multiplied")
        {
            auto e = a * b;
            THEN("The product is the pairwise product of the members 1*3, 2*4")
            {
                CHECK(e == roi_point(1 * 3, 2 * 4));
            }
        }

        WHEN("The point (3,4) is divided by the point(1,2)")
        {
            auto f = b / a;
            THEN("The results is the pairwise integer ratio of the members 3/1, 4/2")
            {
                CHECK(f == roi_point(3, 2));
            }
        }
    }
}

SCENARIO("Roi_point class implements basic mathematical operators with scalar as the right hand operand")
{
    GIVEN("A roi point of x=6, y=5")
    {
        roi_point p(6, 5);
        WHEN("An integer 5 is added to point")
        {
            auto p2 = p + 5;
            THEN("The new point equals 11,10")
            {
                CHECK(p2._x == (6 + 5));
                CHECK(p2._y == (5 + 5));
            }
        }

        WHEN("A Floating point type is added to point")
        {
            auto p3 = p + 0.55;
            THEN("The resulting coordinate is rounded to next integer (away from zero)")
            {
                CHECK(p3._x == 7);
                CHECK(p3._y == 6);
            }
        }

        WHEN("The coordinate is divided by 3")
        {
            auto p4 = p / 3;
            THEN("The coordinates are divided with truncation")
            {
                CHECK(p4._x == 2);
                CHECK(p4._y == 1);
            }
        }

        WHEN("The coordinate is divided by floating point 3.0")
        {
            auto p5 = p / 3.0;
            THEN("The coordinates are divided with rounding")
            {
                CHECK(p5._x == 2);
                CHECK(p5._y == 2);
            }
        }
    }
}
