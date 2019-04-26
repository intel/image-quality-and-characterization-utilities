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

#include "Teisko/Algorithm/VectorMedian.hpp"
#include "catch.hpp"
#include <random>
#include <functional>

/// \page libcalc_specs_vector_median Specs: Teisko vector median
///
/// \snippet this snippet-specs-vector-median

/// [snippet-specs-vector-median]
// Helper functions for the scenarios
int random_int(int lower = 1, int upper = 10)
{
    std::random_device random_device;
    std::mt19937 random_generator(random_device());
    std::uniform_int_distribution<int> distribution(lower, upper);
    return distribution(random_generator);
}

std::vector<int> random_vector(int size, int lower = std::numeric_limits<int>::min(), int upper = std::numeric_limits<int>::max())
{
    //define random generator
    std::random_device random_device;
    std::mt19937 random_generator(random_device());
    std::vector<int> temp_vector(size);

    std::uniform_int_distribution<int> distribution(lower, upper);

    auto gen = std::bind(distribution, random_generator);
    std::generate(temp_vector.begin(), temp_vector.end(), gen);
    return temp_vector;
}

// This function Validates if given number is median of given set
// A median is a pivot point dividing a set two halves, where one half of the points are smaller(or equal)
// and one half of the points are larger(or equal) to.
template <typename T> bool is_valid_median(T result, std::vector<T> &test_vector)
{
    auto size = test_vector.size();
    decltype(size) items = std::count_if(test_vector.begin(), test_vector.end(), [result](T i){return i <= result; });
    if (!(items >= (size + 1) / 2))
    {
        return false;
    }
    items = std::count_if(test_vector.begin(), test_vector.end(), [result](T i){return i >= result; });
    return (items >= (size + 1) / 2);
}

using namespace Teisko;

SCENARIO("Calculating the median value of a vector will not change the input vector "
    "and the returned value separates the higher half of the input data from the lower half")
{
    GIVEN("A vector with same value repeated ANY positive number of times")
    {
        int number_of_times = random_int();
        int value_to_be_repeated = random_int();
        std::vector<int> test_vector(number_of_times, value_to_be_repeated);
        WHEN("The median is calculated for the vector")
        {
            int result = vector_median(test_vector);
            THEN("the result should be the same as the random value")
            {
                CHECK(result == value_to_be_repeated);
                CHECK(true == is_valid_median(result, test_vector));
            }
        }
    }

    GIVEN("a vector with 2 values positive even number of times")
    {
        int rand_size = random_int(3, 17);
        std::vector<int> test_vector(rand_size * 2);// we want an even size
        int a = random_int();
        int b = random_int();
        int result_aaa_bbb;
        WHEN("The median is calculated for the vector aaabbb")
        {
            for (int i = 0; i < rand_size; ++i)
            {
                test_vector[i] = a;
                test_vector[rand_size * 2 - i - 1] = b;
            }
            result_aaa_bbb = vector_median(test_vector);
            AND_WHEN("The median is calculated for the vector ababab")
            {
                for (int i = 0; i < rand_size; ++i)
                {
                    test_vector[i * 2] = a;
                    test_vector[i * 2 + 1] = b;
                }
                int result_ababab = vector_median(test_vector);
                THEN("the result should be the same in both cases")
                {
                    CHECK(result_aaa_bbb == result_ababab);
                }
            }
        }
    }

    GIVEN("a vector of odd number of values using the numbers 1 to 5")
    {
        WHEN("The median is calculated for the vector")
        {
            std::vector<double> test_vector{ 1, 2, 3, 4, 5 };
            double result = vector_median(test_vector);

            THEN("the result should be 3")
            {
                CHECK(result == 3.000);
            }
            AND_THEN("the input vector is not changed")
            {
                std::vector<double> expected_vector{ 1, 2, 3, 4, 5 };
                CHECK(test_vector == expected_vector);
            }
            AND_WHEN("trying to compute for the same vector in reverse")
            {
                std::vector<double> reverse_vector(test_vector.size());
                std::reverse_copy(test_vector.begin(), test_vector.end(), reverse_vector.begin());

                double reverse_result = vector_median(reverse_vector);
                THEN("the result should be the same as before 3")
                {
                    CHECK(reverse_result == 3.000);
                }
                AND_THEN("the input vector is not changed")
                {
                    std::vector<double> expected_vector{ 5, 4, 3, 2, 1 };
                    CHECK(reverse_vector == expected_vector);
                }
            }
        }
    }

    GIVEN("a vector of even number of values using the numbers 1 to 4")
    {
        WHEN("The median is calculated for the vector")
        {
            std::vector<double> test_vector{ 1, 2, 3, 4};
            double result = vector_median(test_vector);

            THEN("the result should be 2.5")
            {
                CHECK(result == 2.500);
            }
            AND_THEN("the input vector is not changed")
            {
                std::vector<double> expected_vector{ 1, 2, 3, 4 };
                CHECK(test_vector == expected_vector);
            }
            AND_WHEN("trying to compute for the same vector in reverse")
            {
                std::vector<double> reverse_vector(test_vector.size());
                std::reverse_copy(test_vector.begin(), test_vector.end(), reverse_vector.begin());

                double reverse_result = vector_median(reverse_vector);
                THEN("the result should be the same as before")
                {
                    CHECK(reverse_result == 2.500);
                }
                AND_THEN("the input vector is not changed")
                {
                    std::vector<double> expected_vector{ 4, 3, 2, 1 };
                    CHECK(reverse_vector == expected_vector);
                }
            }
        }
    }

    GIVEN("a large test vector with random input")
    {
        int size = random_int(50, 150);
        std::vector<int> test_vector = random_vector(size, 1, 60);
        std::vector<int> test_original_vector(test_vector);
        WHEN("the median is calculated")
        {
            int result = vector_median(test_vector);
            THEN("there are at least size/2 elements in the test vector less or equal to reported median")
            {
                CHECK(true == is_valid_median(result, test_vector));
            }
            AND_THEN("the input vector is not changed")
            {
                CHECK(test_vector == test_original_vector);
            }
        }
    }
}

SCENARIO("Calculating median for an empty vector will end in an exception at runtime")
{
    GIVEN("an empty vector")
    {
        std::vector<double> test_vector;
        THEN("calculating median will throw an exception")
        {
            REQUIRE_THROWS_AS(vector_median(test_vector), std::runtime_error);
        }
    }
}
/// [snippet-specs-vector-median]