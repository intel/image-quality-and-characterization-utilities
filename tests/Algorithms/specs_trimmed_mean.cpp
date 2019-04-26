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

#include "Teisko/Algorithm/TrimmedMean.hpp"
#include "catch.hpp"


/// \page libcalc_specs_trimmed_mean Specs: Teisko trimmed_mean
///
/// \snippet this snippet-specs-trimmed_mean

using namespace Teisko;

namespace libcalc_trimmed_mean_tests
{
    SCENARIO("Trimmed mean calculates the mean of samples with the smallest and largest items rejected")
    {
        GIVEN("A vector of ten items")
        {
            auto ten_items = std::vector<int>({ 1, 3, 9, 10, 11, 12, 13, 90, 91, 10001 });
            WHEN("The average of the middle two items, i.e. the middle 20% is calculated")
            {
                auto x = trimmed_mean(ten_items, 0.2);
                THEN("The result is the mean of 11 and 12")
                {
                    auto ref = 0.5 * (11.0 + 12.0);
                    CHECK(x == ref);
                }
            }
        }

        GIVEN("A vector of 100 ones")
        {
            auto one_hundred_ones = std::vector<double>(100, 1.0);
            WHEN("The set is averaged with an increasing percentage of items counted in")
            {
                auto result = std::vector<double>();
                for (auto percentage = 2; percentage <= 100; percentage++)
                {
                    result.push_back(trimmed_mean(one_hundred_ones, percentage / 100.0));
                }
                THEN("The means of all the evaluated subsets match the 1.0")
                {
                    auto ref = std::vector<double>(result.size(), 1.0);
                    CHECK(ref == result);
                }
            }
        }

        GIVEN("A vector of 100 ones with few outliers")
        {
            auto one_hundred_mostly_ones = std::vector<double>(100, 1.0);
            for (int k = 1; k < 100; k += 9)
            {
                one_hundred_mostly_ones[k] = k - 50;
            }
            WHEN("The set is averaged with up to 80 percent of samples included")
            {
                auto result = std::vector<double>();
                for (auto percentage = 2; percentage <= 80; percentage++)
                {
                    result.push_back(trimmed_mean(one_hundred_mostly_ones, percentage / 100.0));
                }
                THEN("The means of all the evaluated subsets match the average of the majority, i.e. 1.0")
                {
                    auto ref = std::vector<double>(result.size(), 1.0);
                    CHECK(ref == result);
                }
            }
        }
    }

    SCENARIO("Trimmed mean can be calculated from a histogram")
    {
        GIVEN("A population of 10 items and it's matching histogram")
        {
            auto pop = std::vector<int>({ 1, 0, 9, 6, 2, 1, 1, 5, 6, 7 });
            auto hist = std::vector<int>({ 1, 3, 1, 0, 0, 1, 2, 1, 0, 1 });

            WHEN("The trimmed mean (of middle 60%) is calculated from the histogram by skipping 2 items and counting 6")
            {
                auto weights = std::vector<double>();
                auto trimmed_mean_hist = trimmed_mean(hist, 2, 6, weights);
                THEN("The result matches the mean calculated from the original set")
                {
                    auto reference_trimmed_mean = trimmed_mean(pop, 0.6);
                    CHECK(trimmed_mean_hist == reference_trimmed_mean);
                }
            }
        }
    }

    SCENARIO("Trimmed mean can be calculated from a histogram with irregular grid centers")
    {
        GIVEN("A histogram of 10 items with it's associated centers")
        {
            auto hist = std::vector<int>({ 1, 3, 1, 0, 0, 1, 2, 1, 0, 1 });
            auto weights = std::vector<int>({ 2, 3, 5, 7, 11, 13, 17, 19, 23, 29 });

            WHEN("The trimmed mean (of middle 60%) is calculated from the histogram by skipping 2 items and counting 6")
            {
                auto trimmed_mean_hist = trimmed_mean(hist, 2, 6, weights);
                THEN("The result matches the mean of samples produces from the weights")
                {
                    auto population = std::vector<int>({
                        2,              // two is represented once
                        3, 3, 3,        // three is there three times in the histogram
                        5,              // five is once
                        13,
                        17, 17,
                        19,
                        29
                    });

                    auto reference_trimmed_mean = trimmed_mean(population, 0.6);
                    CHECK(trimmed_mean_hist == reference_trimmed_mean);
                    CHECK(reference_trimmed_mean == Approx((3 + 3 + 5 + 13 + 17 + 17) / 6.0).epsilon(1e-12));
                }
            }
        }
    }
}