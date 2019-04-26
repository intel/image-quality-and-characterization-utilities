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


#include "Teisko/Image/Support.hpp"
#include "catch.hpp"

using namespace Teisko;

SCENARIO("Libimage helper function 'support' caches small rectangular areas to linear local array")
{
    GIVEN("A 5x5 array of some type")
    {
        using type = float;
        const int rows = 5;
        const int columns = 5;
        type array[rows][columns] = {
            { 0, 1, 2, 3, 4 },
            { 5, 6, 7, 8, 9 },
            { 10, 11, 12, 13, 14 },
            { 15, 16, 17, 18, 19 },
            { 20, 21, 22, 23, 24 }
        };
        WHEN("We initialize a support of 3x3 pointing to the 5x5 array top left corner")
        {
            auto pointer = &array[0][0];
            auto stride = columns;
            auto sup = support<3, 3, type>(pointer, stride);
            THEN("The support data contains the top left 3x3 matrix from the array")
            {
                type reference[3][3] = { { 0, 1, 2 }, { 5, 6, 7 }, { 10, 11, 12 } };
                CHECK(to_vector(reference) == to_vector(sup.data));
                AND_WHEN("We shift new data to the support from column 3")
                {
                    sup.shift(&array[0][3], columns);
                    THEN("the support contains the next 3x3 sub-matrix starting with {1,2,3}")
                    {
                        type reference2[3][3] = { { 1, 2, 3 }, { 6, 7, 8 }, { 11, 12, 13 } };
                        CHECK(to_vector(reference2) == to_vector(sup.data));
                    }
                }
            }
        }
    }
}

SCENARIO("Libimage 'support' contains methods for generating ordered filters")
{
    GIVEN("Initial values for a 3x3 support")
    {
        auto values = std::vector<int>{ 0, 0, 0, 1, 2, 3, 4, 5, 5};
        auto supp = support<3, 3, int>(values.data(), 3);
        THEN("The support function can be queried for median, minimum and maximum")
        {
            CHECK(supp.maximum() == 5);
            CHECK(supp.minimum() == 0);
            CHECK(supp.median() == 2);
        }
    }
}