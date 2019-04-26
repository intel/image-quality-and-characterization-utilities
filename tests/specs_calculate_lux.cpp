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

#include "Teisko/CalculateLux.hpp"
#include "catch.hpp"

static const double e = 0.001;

SCENARIO("Verify assumptions about testing for equality")
{
    GIVEN("+zero and -zero")
    {
        WHEN("variables are of type double")
        {
            double zero_plus  = +0.0;
            double zero_minus = -0.0;

            THEN("those should be equal to 0.0")
            {
                CHECK(zero_plus  == 0.0);
                CHECK(zero_plus  == 0.0);
                CHECK(zero_minus == 0.0f);
                CHECK(zero_minus == 0.0f);
                CHECK(zero_plus  == zero_minus);
            }
        }

        WHEN("variables are of type float")
        {
            float zero_plus  = +0.0f;
            float zero_minus = -0.0f;

            THEN("those should be equal to 0.0")
            {
                CHECK(zero_plus  == 0.0);
                CHECK(zero_minus == 0.0);
                CHECK(zero_plus  == 0.0f);
                CHECK(zero_minus == 0.0f);
                CHECK(zero_plus  == zero_minus);
            }
        }
    }
}

SCENARIO("Calculate lux estimates for exposure")
{
    GIVEN("a known set of exposure coefficients")
    {
        double base_iso = 100;
        double low_limit_pct = 0.125;
        double upper_limit_pct = 0.5;
        double low_limit_lux = 10;
        double upper_limit_lux = 100;
        double dark_adaptation_factor = 0.25;

        WHEN("the exposure and total gain are arrays")
        {
            double exposure[7];
            double total_gain[7];
            uint32_t array_length = 7;

            exposure[0] = 0.001;
            exposure[1] = 20.0;
            exposure[2] = 20.0;
            exposure[3] = 33.0;
            exposure[4] = 33.0;
            exposure[5] = 33.333;
            exposure[6] = 33.333;

            total_gain[0] = 1;
            total_gain[1] = 1;
            total_gain[2] = 2;
            total_gain[3] = 2;
            total_gain[4] = 4;
            total_gain[5] = 16;
            total_gain[6] = 64;

            THEN("the lux calculation function should return successfully")
            {
                double lux_average[7];
                double lux_lower_limit[7];
                double lux_upper_limit[7];

                auto res = Teisko::calculate_lux(
                    exposure, total_gain, array_length, base_iso, low_limit_pct,
                    upper_limit_pct, low_limit_lux, upper_limit_lux, dark_adaptation_factor,
                    lux_average, lux_lower_limit, lux_upper_limit);

                CHECK(Teisko::CALCULATELUX_SUCCESS == res);

                AND_THEN("the first lux estimates should be at infinity")
                {
                    CHECK(std::numeric_limits<double>::infinity() == lux_average[0]);
                    CHECK(std::numeric_limits<double>::infinity() == lux_lower_limit[0]);
                    CHECK(std::numeric_limits<double>::infinity() == lux_upper_limit[0]);

                    AND_THEN("the rest of the estimates should match expected values")
                    {
                        CHECK(1212 == lux_average[1]);
                        CHECK(606 == lux_average[2]);
                        CHECK(367 == lux_average[3]);
                        CHECK(184 == lux_average[4]);
                        CHECK(25 == lux_average[5]);
                        CHECK(3 == lux_average[6]);

                        CHECK(485 == lux_lower_limit[1]);
                        CHECK(242 == lux_lower_limit[2]);
                        CHECK(147 == lux_lower_limit[3]);
                        CHECK(74 == lux_lower_limit[4]);
                        CHECK(10 == lux_lower_limit[5]);
                        CHECK(1 == lux_lower_limit[6]);

                        CHECK(1939 == lux_upper_limit[1]);
                        CHECK(970 == lux_upper_limit[2]);
                        CHECK(587 == lux_upper_limit[3]);
                        CHECK(294 == lux_upper_limit[4]);
                        CHECK(40 == lux_upper_limit[5]);
                        CHECK(5 == lux_upper_limit[6]);
                    }
                }
            }
        }

        WHEN("the exposure and total gain values are scalars")
        {
            uint32_t array_length = 1;
            double exposure = 20.0;
            double total_gain = 1;

            THEN("the function call should return successfully")
            {
                double lux_average;
                double lux_lower_limit;
                double lux_upper_limit;

                auto res = Teisko::calculate_lux(
                    &exposure, &total_gain, array_length, base_iso, low_limit_pct,
                    upper_limit_pct, low_limit_lux, upper_limit_lux, dark_adaptation_factor,
                    &lux_average, &lux_lower_limit, &lux_upper_limit);

                CHECK(Teisko::CALCULATELUX_SUCCESS == res);

                AND_THEN("the first lux estimate should not be at infinity")
                {
                    CHECK(std::numeric_limits<double>::infinity() != lux_average);
                    CHECK(std::numeric_limits<double>::infinity() != lux_lower_limit);
                    CHECK(std::numeric_limits<double>::infinity() != lux_upper_limit);

                    AND_THEN("it should match the expected value")
                    {
                        CHECK(1212 == lux_average);
                        CHECK(485 == lux_lower_limit);
                        CHECK(1939 == lux_upper_limit);
                    }
                }
            }
        }
    }
}

SCENARIO("Calculate lux estimates for exposure can handle invalid input values")
{
    GIVEN("a known set of exposure coefficients")
    {
        double base_iso = 100;
        double low_limit_pct = 0.125;
        double upper_limit_pct = 0.5;
        double low_limit_lux = 10;
        double upper_limit_lux = 100;
        double dark_adaptation_factor = 0.25;

        double exposure[2];
        double total_gain[2];
        uint32_t array_length = 2;

        exposure[0] = 0.001;
        exposure[1] = 20.0;

        total_gain[0] = 1;
        total_gain[1] = 1;

        // Place holders for output values
        double lux_average[2];
        double lux_lower_limit[2];
        double lux_upper_limit[2];

        WHEN("the exposure array is not defined")
        {
            double *invalid_exposure = nullptr;

            THEN("the function call should return error value")
            {
                CHECK(Teisko::CALCULATELUX_ERROR_INVALID_PARAMETER == Teisko::calculate_lux(
                    invalid_exposure, total_gain, array_length, base_iso, low_limit_pct,
                    upper_limit_pct, low_limit_lux, upper_limit_lux, dark_adaptation_factor,
                    lux_average, lux_lower_limit, lux_upper_limit));
            }
        }

        WHEN("the total gain array is not defined")
        {
            double *invalid_total_gain = nullptr;

            THEN("the function call should return error value")
            {
                CHECK(Teisko::CALCULATELUX_ERROR_INVALID_PARAMETER == Teisko::calculate_lux(
                    exposure, invalid_total_gain, array_length, base_iso, low_limit_pct,
                    upper_limit_pct, low_limit_lux, upper_limit_lux, dark_adaptation_factor,
                    lux_average, lux_lower_limit, lux_upper_limit));
            }
        }

        WHEN("the exposure and total gain array are not defined")
        {
            double *invalid_exposure = nullptr;
            double *invalid_total_gain = nullptr;

            THEN("the function call should return error value")
            {
                CHECK(Teisko::CALCULATELUX_ERROR_INVALID_PARAMETER == Teisko::calculate_lux(
                    invalid_exposure, invalid_total_gain, array_length, base_iso, low_limit_pct,
                    upper_limit_pct, low_limit_lux, upper_limit_lux, dark_adaptation_factor,
                    lux_average, lux_lower_limit, lux_upper_limit));
            }
        }

        WHEN("the array length value is invalid")
        {
            array_length = 0;

            THEN("the function call should return error value")
            {
                CHECK(Teisko::CALCULATELUX_ERROR_INVALID_PARAMETER == Teisko::calculate_lux(
                    exposure, total_gain, array_length, base_iso, low_limit_pct,
                    upper_limit_pct, low_limit_lux, upper_limit_lux, dark_adaptation_factor,
                    lux_average, lux_lower_limit, lux_upper_limit));
            }
        }

        WHEN("base is not valid")
        {
            base_iso = -1;

            THEN("the function call should return error value")
            {
                CHECK(Teisko::CALCULATELUX_ERROR_INVALID_PARAMETER == Teisko::calculate_lux(
                    exposure, total_gain, array_length, base_iso, low_limit_pct,
                    upper_limit_pct, low_limit_lux, upper_limit_lux, dark_adaptation_factor,
                    lux_average, lux_lower_limit, lux_upper_limit));
            }
        }
    }
}
