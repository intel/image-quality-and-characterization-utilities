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


#include "Teisko/ColorCorrection.hpp"
#include "Teisko/Preprocessing.hpp"
#include "catch.hpp"

using namespace Teisko;

namespace teisko_libacm_tests
{
    static const double e = 0.001;

    ccm_input_params<rgb_cs::sRGB> get_ccm_parameters()
    {
        std::array<rgb<double, rgb_cs::Sensor, wp::None>, 24> input =
        {{{ 0.1339, 0.0861, 0.0692},
          { 0.4474, 0.3033, 0.2538},
          { 0.1537, 0.1998, 0.2646},
          { 0.1141, 0.1173, 0.0774},
          { 0.2602, 0.2592, 0.3556},
          { 0.3059, 0.4710, 0.4397},
          { 0.4656, 0.2092, 0.1086},
          { 0.1111, 0.1499, 0.2796},
          { 0.3987, 0.1570, 0.1425},
          { 0.0959, 0.0682, 0.1086},
          { 0.3743, 0.3921, 0.1954},
          { 0.5737, 0.3332, 0.1534},
          { 0.0548, 0.0861, 0.1859},
          { 0.1385, 0.2217, 0.1344},
          { 0.3150, 0.0990, 0.0746},
          { 0.6756, 0.4733, 0.2022},
          { 0.4093, 0.1909, 0.2674},
          { 0.1339, 0.2734, 0.3515},
          { 0.8369, 0.8150, 0.7926},
          { 0.5554, 0.5504, 0.5469},
          { 0.3515, 0.3515, 0.3515},
          { 0.1948, 0.1923, 0.1941},
          { 0.0913, 0.0919, 0.0950},
          { 0.0335, 0.0330, 0.0339}}
        };

        std::array<rgb<double, rgb_cs::sRGB, wp::D65>, 24> target =
        {{{0.1714, 0.0844, 0.0578},
          {0.5395, 0.3050, 0.2232},
          {0.1221, 0.1946, 0.3372},
          {0.0953, 0.1500, 0.0561},
          {0.2346, 0.2159, 0.4397},
          {0.1356, 0.5089, 0.4020},
          {0.6724, 0.2086, 0.0252},
          {0.0802, 0.1046, 0.3813},
          {0.5333, 0.1022, 0.1248},
          {0.1119, 0.0452, 0.1500},
          {0.3372, 0.5029, 0.0513},
          {0.7454, 0.3663, 0.0273},
          {0.0395, 0.0467, 0.3050},
          {0.0612, 0.2961, 0.0666},
          {0.4287, 0.0369, 0.0452},
          {0.7991, 0.5711, 0.0137},
          {0.4969, 0.0931, 0.3005},
          {0.0024, 0.2346, 0.3564},
          {0.8963, 0.8963, 0.8879},
          {0.5776, 0.5776, 0.5776},
          {0.3515, 0.3515, 0.3515},
          {0.1946, 0.1946, 0.1912},
          {0.0908, 0.0908, 0.0908},
          {0.0343, 0.0343, 0.0343}}
        };

        ccm_input_params<rgb_cs::sRGB> output;
        for (int i = 0; i < 24; ++i)
        {
            output.color_checker_classic[i].input = input[i];
            output.color_checker_classic[i].target = target[i];
            output.color_checker_classic[i].weight = i < 18 ? 1.0 : 0.0;
            output.color_checker_classic[i].achromatic = i < 18 ? false : true;
        }

        return output;
    }

    SCENARIO("It is possible to calculate color correction matrix by "
        "using predefined input to validate internal calculations.", "[CCM]")
    {
        GIVEN("example sensor RGB response given for color checker classic")
        {
            auto parameters = get_ccm_parameters();
            auto optimizer = ccm_optimization<rgb_cs::sRGB>(parameters);
            WHEN("the developer characterizes color matrix")
            {
                auto ccm = optimizer.characterize();
                THEN("results match predefined values")
                {
                    // Do regression testing
                    CHECK(1.6782 == Approx(ccm(0, 0)).epsilon(e));
                    CHECK(-0.6897 == Approx(ccm(0, 1)).epsilon(e));
                    CHECK(0.0115 == Approx(ccm(0, 2)).epsilon(e));
                    CHECK(-0.2772 == Approx(ccm(1, 0)).epsilon(e));
                    CHECK(1.7945 == Approx(ccm(1, 1)).epsilon(e));
                    CHECK(-0.5173 == Approx(ccm(1, 2)).epsilon(e));
                    CHECK(-0.0177 == Approx(ccm(2, 0)).epsilon(e));
                    CHECK(-0.6744 == Approx(ccm(2, 1)).epsilon(e));
                    CHECK(1.6921 == Approx(ccm(2, 2)).epsilon(e));
                }
                AND_THEN("rgb values in target color space can be calculated")
                {
                    auto patch1 = parameters.color_checker_classic.front().input;
                    auto corr = ccm.apply<rgb_cs::sRGB>(parameters.color_checker_classic.front().input);
                    CHECK(patch1.r != corr.r);
                    CHECK(patch1.g != corr.g);
                    CHECK(patch1.b != corr.b);
                }
            }
        }
    }
}