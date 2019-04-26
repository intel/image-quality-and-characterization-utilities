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

#include "Teisko/Algorithm/LinearSpace.hpp"
#include "catch.hpp"
#include <memory>
#include <random>
#include <algorithm>
#include <vector>

/// \page libcalc_specs_linear_space Specs: Teisko linear space
///
/// \snippet this snippet-specs-linear-space

/// [snippet-specs-linear-space]
// Helper functions for the scenarios
double any_value()
{
    std::random_device random_device;
    std::mt19937 random_generator(random_device());

    auto lower = std::numeric_limits<int>::min();
    auto upper = std::numeric_limits<int>::max();

    std::uniform_real_distribution<double> distribution(lower, upper);

    return distribution(random_generator);
}

auto reverse = [](std::vector<double> const& x)
{
    std::vector<double> output(x.size());
    std::reverse_copy(x.begin(), x.end(), output.begin());
    return output;
};

std::size_t count_unique_differences(std::vector<double> const &  vec)
{
    auto n_1 = vec.size() - 1;
    std::vector<double> distributed_vector(n_1);
    for (unsigned int i = 0; i < n_1; i++)
    {
        distributed_vector[i] = vec[i + 1] - vec[i];
    }

    //we need this threshold and not using std::numeric_limits<double>::epsilon())
    //since it is too small for std::unique
    double threshold = 0.00001;
    auto last = std::unique(distributed_vector.begin(), distributed_vector.end(), [threshold](double left, double right) {return std::abs(left - right) < threshold; });
    return std::distance(distributed_vector.begin(), last);
}

using namespace Teisko;

SCENARIO("Create linear_space for valid input")
{
    GIVEN("5 and 10 as points and number of cells is 5")
    {
        WHEN("trying to create linear_space")
        {
            std::vector<double> result = linear_space(5, 10, 5);
            THEN("the result should be 5 cells")
            {
                CHECK(result.size() == 5);
                CHECK(result[0] == 5.000);
                CHECK(result[1] == 6.2500);
                CHECK(result[2] == 7.5000);
                CHECK(result[3] == 8.7500);
                CHECK(result[4] == 10.000);
            }
        }
    }
    GIVEN("5 and 5 as points and number of cells is 5")
    {
        WHEN("trying to create linear_space =")
        {
            std::vector<double> result = linear_space(5, 5, 5);
            THEN("the result should be 5 cells")
            {
                CHECK(result.size() == 5);
                CHECK(result[0] == 5.000);
                CHECK(result[1] == 5.000);
                CHECK(result[2] == 5.000);
                CHECK(result[3] == 5.000);
                CHECK(result[4] == 5.000);
            }
        }
    }

    GIVEN("10 and 5 as points and number of cells is 5")
    {
        WHEN("trying to create linear_space")
        {
            std::vector<double> result = linear_space(10, 5, 5);
            THEN("the result should be 5 cells")
            {
                CHECK(result.size() == 5);
                CHECK(result[0] == 10.000);
                CHECK(result[1] == 8.7500);
                CHECK(result[2] == 7.5000);
                CHECK(result[3] == 6.2500);
                CHECK(result[4] == 5.000);
            }
        }
    }
}

SCENARIO("Linear space with zero, one or two points is well defined")
{
    auto number_of_repetitions = 10;

    THEN("any linear space with zero points will be empty")
    {
        while (number_of_repetitions-- > 0)
        {
            auto a = any_value();
            auto b = any_value();

            INFO("a: " << a);
            INFO("b: " << b);
            CHECK(true == linear_space(a, b, 0).empty());
        }
    }

    THEN("any linear space with one point will contain exactly the first input value")
    {
        while (number_of_repetitions-- > 0)
        {
            auto a = any_value();
            auto b = any_value();
            auto result = linear_space(a, b, 1);

            INFO("a: " << a);
            INFO("b: " << b);
            REQUIRE(1 == result.size());
            CHECK(a == result[0]);
        }
    }

    THEN("any linear space with two points will contain exactly the input values")
    {
        while (number_of_repetitions-- > 0)
        {
            auto a = any_value();
            auto b = any_value();
            auto result = linear_space(a, b, 2);

            INFO("a: " << a);
            INFO("b: " << b);
            REQUIRE(2 == result.size());
            CHECK(a == result[0]);
            CHECK(b == result[1]);
        }
    }
}

SCENARIO("Any linear space (a, b) with two or more points will be same as"
    " the reversed linear space of (b, a)")
{
    GIVEN("a range of points to check starting from 2")
    {
        unsigned int starting_point = 2;
        unsigned int ending_point = 10;

        THEN("the linear space (a, b) is same as reversed linear space (b, a)")
        {
            for (unsigned int n = starting_point; n < ending_point; ++n)
            {
                auto a = any_value();
                auto b = any_value();

                auto fwd = linear_space(a, b, n);
                auto rev = reverse(linear_space(b, a, n));
                INFO("a: " << a);
                INFO("b: " << b);
                INFO("n: " << n);
                REQUIRE(fwd.size() == rev.size());
                for (unsigned int i = 0; i < n; i ++)
                {
                    REQUIRE(fwd[i] == Approx(rev[i]));
                }
            }
        }
    }
}

SCENARIO("For any n greater than two the linear space will contain the"
    " given number of equally distributed points")
{
    GIVEN("a range of points to check starting from 2")
    {
        unsigned int starting_point = 2;
        unsigned int ending_point = 10;

        THEN("all the points should have the same space between them i.e equally distributed")
        {
            for (auto n = starting_point; n < ending_point; ++n)
            {
                auto a = any_value();
                auto b = any_value();

                auto x = linear_space(a, b, n);
                auto y = reverse(linear_space(b, a, n));

                auto count_forward = count_unique_differences(x);
                auto count_reverse = count_unique_differences(y);;
                INFO("a: " << a);
                INFO("b: " << b);
                INFO("n: " << n);
                CHECK(count_forward == 1);
                CHECK(count_reverse == 1);
            }
        }
    }
}
/// [snippet-specs-linear-space]
