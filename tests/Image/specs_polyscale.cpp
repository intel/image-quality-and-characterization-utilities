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

#include "Teisko/Image/Polyscale.hpp"
#include "catch.hpp"
#include <memory>
#include <random>
#include <algorithm>
#include <vector>

using namespace Teisko;

namespace libpolyscale_tests
{
    SCENARIO("Scale 2d-grid using 3th order polynomial scaler")
    {
        GIVEN("A predefined 2d double grid")
        {
            const static double e = 0.01;

            grid_2d<double> A(15, 9);
            A._grid = std::vector<double>(
            { 0.4462, 0.4017, 0.3619, 0.3279, 0.3005, 0.2804, 0.2681, 0.2640, 0.2681, 0.2804, 0.3005, 0.3279, 0.3619, 0.4017, 0.4462
            , 0.4225, 0.3766, 0.3356, 0.3005, 0.2722, 0.2515, 0.2388, 0.2345, 0.2388, 0.2515, 0.2722, 0.3005, 0.3356, 0.3766, 0.4225
            , 0.4052, 0.3582, 0.3162, 0.2804, 0.2515, 0.2302, 0.2173, 0.2129, 0.2173, 0.2302, 0.2515, 0.2804, 0.3162, 0.3582, 0.4052
            , 0.3946, 0.3469, 0.3044, 0.2681, 0.2388, 0.2173, 0.2042, 0.1997, 0.2042, 0.2173, 0.2388, 0.2681, 0.3044, 0.3469, 0.3946
            , 0.3910, 0.3432, 0.3005, 0.2640, 0.2345, 0.2129, 0.1997, 0.1953, 0.1997, 0.2129, 0.2345, 0.2640, 0.3005, 0.3432, 0.3910
            , 0.3946, 0.3469, 0.3044, 0.2681, 0.2388, 0.2173, 0.2042, 0.1997, 0.2042, 0.2173, 0.2388, 0.2681, 0.3044, 0.3469, 0.3946
            , 0.4052, 0.3582, 0.3162, 0.2804, 0.2515, 0.2302, 0.2173, 0.2129, 0.2173, 0.2302, 0.2515, 0.2804, 0.3162, 0.3582, 0.4052
            , 0.4225, 0.3766, 0.3356, 0.3005, 0.2722, 0.2515, 0.2388, 0.2345, 0.2388, 0.2515, 0.2722, 0.3005, 0.3356, 0.3766, 0.4225
            , 0.4462, 0.4017, 0.3619, 0.3279, 0.3005, 0.2804, 0.2681, 0.2640, 0.2681, 0.2804, 0.3005, 0.3279, 0.3619, 0.4017, 0.4462 });

            WHEN("The grid is resized to original size")
            {
                grid_2d<double> B(15, 9);

                poly_3_scaler_2d<double> polynomial;
                polynomial.scale(A, B);

                THEN("Grid values should match to original values")
                {
                    for (size_t i = 0; i < A._grid.size(); i++)
                    {
                        CHECK(A._grid[i] == Approx(B._grid[i]).epsilon(e));
                    }
                }
            }
        }

        GIVEN("A predefined 2d uint32 grid")
        {
            const static double e = 1; // Model fitting differences

            grid_2d<uint32_t> A(15, 9);
            A._grid = std::vector<uint32_t>(
            { 104, 95, 87, 80, 75, 71, 69, 68, 69, 71, 75, 80, 87, 95, 104
            , 99, 90, 82, 75, 70, 66, 63, 62, 63, 66, 70, 75, 82, 90, 99
            , 96, 86, 78, 71, 66, 62, 59, 58, 59, 62, 66, 71, 78, 86, 96
            , 94, 84, 76, 69, 63, 59, 57, 56, 57, 59, 63, 69, 76, 84, 94
            , 93, 84, 75, 68, 62, 58, 56, 55, 56, 58, 62, 68, 75, 84, 93
            , 94, 84, 76, 69, 63, 59, 57, 56, 57, 59, 63, 69, 76, 84, 94
            , 96, 86, 78, 71, 66, 62, 59, 58, 59, 62, 66, 71, 78, 86, 96
            , 99, 90, 82, 75, 70, 66, 63, 62, 63, 66, 70, 75, 82, 90, 99
            , 104, 95, 87, 80, 75, 71, 69, 68, 69, 71, 75, 80, 87, 95, 104 });

            WHEN("The grid is resized to original size")
            {
                grid_2d<uint32_t> B(15, 9);

                poly_3_scaler_2d<uint32_t> polynomial;
                polynomial.scale(A, B);

                THEN("Grid values should match to original values")
                {
                    for (size_t i = 0; i < A._grid.size(); i++)
                    {
                        CHECK(A._grid[i] == Approx(B._grid[i]).epsilon(e));
                    }
                }
            }
        }
    }

    SCENARIO("Scale 2d-grid using 4th order polynomial scaler")
    {
        GIVEN("A predefined 2d double grid")
        {
            const static double e = 0.001;

            grid_2d<double> A(15, 9);
            A._grid = std::vector<double>(
            { 0.4462, 0.4017, 0.3619, 0.3279, 0.3005, 0.2804, 0.2681, 0.2640, 0.2681, 0.2804, 0.3005, 0.3279, 0.3619, 0.4017, 0.4462
            , 0.4225, 0.3766, 0.3356, 0.3005, 0.2722, 0.2515, 0.2388, 0.2345, 0.2388, 0.2515, 0.2722, 0.3005, 0.3356, 0.3766, 0.4225
            , 0.4052, 0.3582, 0.3162, 0.2804, 0.2515, 0.2302, 0.2173, 0.2129, 0.2173, 0.2302, 0.2515, 0.2804, 0.3162, 0.3582, 0.4052
            , 0.3946, 0.3469, 0.3044, 0.2681, 0.2388, 0.2173, 0.2042, 0.1997, 0.2042, 0.2173, 0.2388, 0.2681, 0.3044, 0.3469, 0.3946
            , 0.3910, 0.3432, 0.3005, 0.2640, 0.2345, 0.2129, 0.1997, 0.1953, 0.1997, 0.2129, 0.2345, 0.2640, 0.3005, 0.3432, 0.3910
            , 0.3946, 0.3469, 0.3044, 0.2681, 0.2388, 0.2173, 0.2042, 0.1997, 0.2042, 0.2173, 0.2388, 0.2681, 0.3044, 0.3469, 0.3946
            , 0.4052, 0.3582, 0.3162, 0.2804, 0.2515, 0.2302, 0.2173, 0.2129, 0.2173, 0.2302, 0.2515, 0.2804, 0.3162, 0.3582, 0.4052
            , 0.4225, 0.3766, 0.3356, 0.3005, 0.2722, 0.2515, 0.2388, 0.2345, 0.2388, 0.2515, 0.2722, 0.3005, 0.3356, 0.3766, 0.4225
            , 0.4462, 0.4017, 0.3619, 0.3279, 0.3005, 0.2804, 0.2681, 0.2640, 0.2681, 0.2804, 0.3005, 0.3279, 0.3619, 0.4017, 0.4462 });

            WHEN("The grid is resized to original size")
            {
                grid_2d<double> B(15, 9);

                poly_4_scaler_2d<double> polynomial;
                polynomial.scale(A, B);

                THEN("Grid values should match to original values")
                {
                    for (size_t i = 0; i < A._grid.size(); i++)
                    {
                        CHECK(A._grid[i] == Approx(B._grid[i]).epsilon(e));
                    }
                }
            }
        }

        GIVEN("A predefined 2d uint32 grid")
        {
            const static double e = 1; // Model fitting differences

            grid_2d<uint32_t> A(15, 9);
            A._grid = std::vector<uint32_t>(
            { 104, 95, 87, 80, 75, 71, 69, 68, 69, 71, 75, 80, 87, 95, 104
            , 99, 90, 82, 75, 70, 66, 63, 62, 63, 66, 70, 75, 82, 90, 99
            , 96, 86, 78, 71, 66, 62, 59, 58, 59, 62, 66, 71, 78, 86, 96
            , 94, 84, 76, 69, 63, 59, 57, 56, 57, 59, 63, 69, 76, 84, 94
            , 93, 84, 75, 68, 62, 58, 56, 55, 56, 58, 62, 68, 75, 84, 93
            , 94, 84, 76, 69, 63, 59, 57, 56, 57, 59, 63, 69, 76, 84, 94
            , 96, 86, 78, 71, 66, 62, 59, 58, 59, 62, 66, 71, 78, 86, 96
            , 99, 90, 82, 75, 70, 66, 63, 62, 63, 66, 70, 75, 82, 90, 99
            , 104, 95, 87, 80, 75, 71, 69, 68, 69, 71, 75, 80, 87, 95, 104 });

            WHEN("The grid is resized to original size")
            {
                grid_2d<uint32_t> B(15, 9);

                poly_4_scaler_2d<uint32_t> polynomial;
                polynomial.scale(A, B);

                THEN("Grid values should match to original values")
                {
                    for (size_t i = 0; i < A._grid.size(); i++)
                    {
                        CHECK(A._grid[i] == Approx(B._grid[i]).epsilon(e));
                    }
                }
            }
        }
    }
}
