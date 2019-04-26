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


#include "Teisko/BayerInfo.hpp"
#include "catch.hpp"
#include <memory>
#include <random>
#include <algorithm>
#include <vector>

using namespace Teisko;

namespace teisko_libutils_tests
{
    SCENARIO("Bayer pattern info tests")
    {
        GIVEN("The default constructor for bayer info")
        {
            bayer_info_s info;
            THEN("The struct gives properties for RGGB sensor")
            {
                CHECK(info.is_2x2_sensor() == true);
                CHECK(bayer_pattern_e::rggb == info);
                CHECK(0 == (int)info);
            }
        }

        GIVEN("A list of regular sensor patterns")
        {
            auto patterns = bayer_info_s::get_regular_2x2_patterns();

            THEN("These patterns are shown to be non IR, non DP 2x2 sensors")
            {
                for (auto &p : patterns)
                {
                    bayer_info_s info(p);
                    CHECK(info.get_width() == 2);
                    CHECK(info.get_height() == 2);
                    CHECK(info.get_channels() == 4);
                    CHECK(info.is_ir_sensor() == false);
                    CHECK(info.is_dp_4x2_sensor() == false);
                    CHECK(info.is_2x2_sensor() == true);
                }
            }
        }

        GIVEN("A list of regular sensor patterns")
        {
            auto patterns = bayer_info_s::get_ir_2x2_patterns();

            THEN("The these patterns are shown to IR, non DP 2x2 sensors")
            {
                for (auto &p : patterns)
                {
                    bayer_info_s info(p);
                    CHECK(info.get_width() == 2);
                    CHECK(info.get_height() == 2);
                    CHECK(info.get_channels() == 4);
                    CHECK(info.is_ir_sensor() == true);
                    CHECK(info.is_dp_4x2_sensor() == false);
                    CHECK(info.is_2x2_sensor() == false);
                }
            }
        }

        GIVEN("A list of DP sensor patterns")
        {
            auto patterns = bayer_info_s::get_dp_4x2_patterns();

            THEN("The these patterns are shown to DP 4x2 sensors")
            {
                for (auto &p : patterns)
                {
                    bayer_info_s info(p);
                    CHECK(info.get_width() == 4);
                    CHECK(info.get_height() == 2);
                    CHECK(info.get_channels() == 8);
                    CHECK(info.is_ir_sensor() == false);
                    CHECK(info.is_dp_4x2_sensor() == true);
                    CHECK(info.is_2x2_sensor() == false);
                }
            }
        }

        GIVEN("A list of 4x4 RGBIR patterns")
        {
            auto patterns = bayer_info_s::get_ir_4x4_patterns();

            THEN("The these patterns are shown to 4x4 IR sensors")
            {
                for (auto &p : patterns)
                {
                    bayer_info_s info(p);
                    CHECK(info.get_width() == 4);
                    CHECK(info.get_height() == 4);
                    CHECK(info.get_channels() == 16);
                    CHECK(info.is_ir_sensor() == true);
                    CHECK(info.is_dp_4x2_sensor() == false);
                    CHECK(info.is_2x2_sensor() == false);
                    CHECK(info.is_2x2_ir_sensor() == false);
                }
            }
        }

        GIVEN("A list of unknown patterns")
        {
            // This test case is the first thing to fail when new patterns are added
            auto patterns = {
                (int)bayer_pattern_e::ibrg + 1,
                (int)bayer_pattern_e::gbrg_4x2 + 1,
                (int)bayer_pattern_e::igig_grgb_igig_gbgr + 1
            };

            THEN("The these patterns return 0x0 size")
            {
                for (auto &p : patterns)
                {
                    bayer_info_s info((bayer_pattern_e)p);
                    CHECK(info.get_width() == 0);
                    CHECK(info.get_height() == 0);
                    CHECK(info.get_channels() == 0);
                    CHECK(info.is_ir_sensor() == false);
                    CHECK(info.is_dp_4x2_sensor() == false);
                    CHECK(info.is_4x4_ir_sensor() == false);
                    CHECK(info.is_2x2_sensor() == false);
                }
            }
        }
    }
}