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

#include "catch.hpp"
#include "Teisko/SpectralResponse.hpp"

/// \page libmsr_specs Specs: Teisko
///
/// \snippet this snippet-specs-msr

using namespace Teisko;

namespace msr_tests
{
    const std::string standard_illuminant_A = "Standard Illuminant A";
    const std::string standard_illuminant_F12 = "Standard Illuminant F12";
    const std::string standard_illuminant_D50 = "Standard Illuminant D50";
    const std::string standard_illuminant_D55 = "Standard Illuminant D55";
    const std::string standard_illuminant_D65 = "Standard Illuminant D65";
    const std::string standard_illuminant_D75 = "Standard Illuminant D75";

    const std::vector<float> spectral_power_distribution_f12 = { 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.9600f, 0.8960f, 0.8320f, 0.7680f, 0.7040f, 0.6400f, 0.6020f, 0.5640f, 0.5260f, 0.4880f, 0.4500f, 0.4260f, 0.4020f, 0.3780f, 0.3540f, 0.3300f, 0.5020f, 0.6740f, 0.8460f, 1.0180f, 1.1900f, 3.4480f, 5.7060f, 7.9640f, 10.2220f, 12.4800f, 10.2080f, 7.9360f, 5.6640f, 3.3920f, 1.1200f, 1.0840f, 1.0480f, 1.0120f, 0.9760f, 0.9400f, 0.9680f, 0.9960f, 1.0240f, 1.0520f, 1.0800f, 1.1380f, 1.1960f, 1.2540f, 1.3120f, 1.3700f, 1.4520f, 1.5340f, 1.6160f, 1.6980f, 1.7800f, 7.2340f, 12.6880f, 18.1420f, 23.5960f, 29.0500f, 24.8200f, 20.5900f, 16.3600f, 12.1300f, 7.9000f, 6.8500f, 5.8000f, 4.7500f, 3.7000f, 2.6500f, 2.6620f, 2.6740f, 2.6860f, 2.6980f, 2.7100f, 2.6980f, 2.6860f, 2.6740f, 2.6620f, 2.6500f, 2.6180f, 2.5860f, 2.5540f, 2.5220f, 2.4900f, 2.4580f, 2.4260f, 2.3940f, 2.3620f, 2.3300f, 2.2840f, 2.2380f, 2.1920f, 2.1460f, 2.1000f, 2.0620f, 2.0240f, 1.9860f, 1.9480f, 1.9100f, 2.1300f, 2.3500f, 2.5700f, 2.7900f, 3.0100f, 4.5740f, 6.1380f, 7.7020f, 9.2660f, 10.8300f, 11.0400f, 11.2500f, 11.4600f, 11.6700f, 11.8800f, 10.8800f, 9.8800f, 8.8800f, 7.8800f, 6.8800f, 6.1900f, 5.5000f, 4.8100f, 4.1200f, 3.4300f, 3.0420f, 2.6540f, 2.2660f, 1.8780f, 1.4900f, 1.3760f, 1.2620f, 1.1480f, 1.0340f, 0.9200f, 0.8780f, 0.8360f, 0.7940f, 0.7520f, 0.7100f, 0.6880f, 0.6660f, 0.6440f, 0.6220f, 0.6000f, 0.6060f, 0.6120f, 0.6180f, 0.6240f, 0.6300f, 0.7240f, 0.8180f, 0.9120f, 1.0060f, 1.1000f, 1.7920f, 2.4840f, 3.1760f, 3.8680f, 4.5600f, 10.5280f, 16.4960f, 22.4640f, 28.4320f, 34.4000f, 40.6000f, 46.8000f, 53.0000f, 59.2000f, 65.4000f, 58.2160f, 51.0320f, 43.8480f, 36.6640f, 29.4800f, 25.0160f, 20.5520f, 16.0880f, 11.6240f, 7.1600f, 6.3440f, 5.5280f, 4.7120f, 3.8960f, 3.0800f, 2.9580f, 2.8360f, 2.7140f, 2.5920f, 2.4700f, 2.4300f, 2.3900f, 2.3500f, 2.3100f, 2.2700f, 2.8340f, 3.3980f, 3.9620f, 4.5260f, 5.0900f, 6.4640f, 7.8380f, 9.2120f, 10.5860f, 11.9600f, 12.6320f, 13.3040f, 13.9760f, 14.6480f, 15.3200f, 15.1100f, 14.9000f, 14.6900f, 14.4800f, 14.2700f, 13.7880f, 13.3060f, 12.8240f, 12.3420f, 11.8600f, 11.3440f, 10.8280f, 10.3120f, 9.7960f, 9.2800f, 9.8860f, 10.4920f, 11.0980f, 11.7040f, 12.3100f, 23.5540f, 34.7980f, 46.0420f, 57.2860f, 68.5300f, 65.4280f, 62.3260f, 59.2240f, 56.1220f, 53.0200f, 45.3500f, 37.6800f, 30.0100f, 22.3400f, 14.6700f, 14.6120f, 14.5540f, 14.4960f, 14.4380f, 14.3800f, 14.4460f, 14.5120f, 14.5780f, 14.6440f, 14.7100f, 13.0600f, 11.4100f, 9.7600f, 8.1100f, 6.4600f, 5.6820f, 4.9040f, 4.1260f, 3.3480f, 2.5700f, 2.6060f, 2.6420f, 2.6780f, 2.7140f, 2.7500f, 3.0360f, 3.3220f, 3.6080f, 3.8940f, 4.1800f, 4.0320f, 3.8840f, 3.7360f, 3.5880f, 3.4400f, 3.3140f, 3.1880f, 3.0620f, 2.9360f, 2.8100f, 2.7320f, 2.6540f, 2.5760f, 2.4980f, 2.4200f, 2.2640f, 2.1080f, 1.9520f, 1.7960f, 1.6400f, 1.5840f, 1.5280f, 1.4720f, 1.4160f, 1.3600f, 1.3860f, 1.4120f, 1.4380f, 1.4640f, 1.4900f, 1.6200f, 1.7500f, 1.8800f, 2.0100f, 2.1400f, 2.1800f, 2.2200f, 2.2600f, 2.3000f, 2.3400f, 2.1560f, 1.9720f, 1.7880f, 1.6040f, 1.4200f, 1.4580f, 1.4960f, 1.5340f, 1.5720f, 1.6100f, 2.2960f, 2.9820f, 3.6680f, 4.3540f, 5.0400f, 5.4280f, 5.8160f, 6.2040f, 6.5920f, 6.9800f, 6.2220f, 5.4640f, 4.7060f, 3.9480f, 3.1900f, 2.6940f, 2.1980f, 1.7020f, 1.2060f, 0.7100f, 0.6280f, 0.5460f, 0.4640f, 0.3820f, 0.3000f, 0.2920f, 0.2840f, 0.2760f, 0.2680f, 0.2600f, 0.2540f, 0.2480f, 0.2420f, 0.2360f, 0.2300f, 0.2400f, 0.2500f, 0.2600f, 0.2700f, 0.2800f, 0.2800f, 0.2800f, 0.2800f, 0.2800f, 0.2800f, 0.2660f, 0.2520f, 0.2380f, 0.2240f, 0.2100f, 0.2020f, 0.1940f, 0.1860f, 0.1780f, 0.1700f, 0.1780f, 0.1860f, 0.1940f, 0.2020f, 0.2100f, 0.2060f, 0.2020f, 0.1980f, 0.1940f, 0.1900f, 0.1820f, 0.1740f, 0.1660f, 0.1580f, 0.1500f, 0.1400f, 0.1300f, 0.1200f, 0.1100f, 0.1000f, 0.0900f, 0.0800f, 0.0700f, 0.0600f, 0.0500f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f };
    const uint16_t minimum_wavelenght_f12 = 360;
    const uint16_t maximum_wavelenght_f12 = 830;

    const module_spectral_response_data msr_data(
        bayer_pattern_e::rggb,
        std::array<int, 16>{ { 0, 1, -1, -1, 2, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 } },
        std::vector<float>{ 380, 385, 390, 395, 400, 405, 410, 415, 420, 425, 430, 435, 440, 445, 450, 455, 460, 465, 470, 475, 480, 485, 490, 495, 500, 505, 510, 515, 520, 525, 530, 535, 540, 545, 550, 555, 560, 565, 570, 575, 580, 585, 590, 595, 600, 605, 610, 615, 620, 625, 630, 635, 640, 645, 650, 655, 660, 665, 670, 675, 680, 685, 690, 695, 700, 705, 710, 715, 720, 725, 730, 735, 740, 745, 750, 755, 760, 765, 770, 775, 780 },
        std::vector<std::vector<float>>{
            { 0.038216859f, 0.101574881f, 0.043915076f, 0.122957809f, 0.281378238f, 1.046717772f, 2.753346067f, 4.397088966f, 4.056282409f, 3.458097481f, 3.325477778f, 3.002000111f, 2.668543591f, 2.49105402f, 2.47852412f, 2.191135635f, 2.231940657f, 2.588955675f, 3.361769656f, 4.490368405f, 5.489541495f, 6.60048636f, 7.511266772f, 8.227252536f, 9.12842693f, 9.786561853f, 10.9208052f, 12.29122335f, 13.64025794f, 14.93174673f, 14.98290919f, 14.26511617f, 12.98007383f, 12.53712585f, 11.98298521f, 11.82589602f, 12.10591165f, 12.93774877f, 14.64144098f, 22.59904695f, 43.87276944f, 68.70762961f, 82.86342154f, 85.90088736f, 82.82270773f, 78.42303001f, 72.72159875f, 68.84666247f, 64.73342846f, 60.33065496f, 55.96665621f, 51.53745653f, 46.93563105f, 43.02297722f, 41.11368045f, 36.46029687f, 32.90887846f, 31.22024534f, 28.37710555f, 23.19019685f, 16.69878885f, 9.423053721f, 4.37905917f, 1.853533179f, 1.016600106f, 0.505941713f, 0.239427227f, 0.227453901f, 0.019721761f, 0.037292023f, 0.062019683f, 0.023175932f, 0.157789391f, 0.0f, 0.004141147f, 0.067865657f, 0.128281791f, 0.0f, 0.042935523f, 0.0f, 0.118663556f },
            { 0.135351375f, 0.135881629f, 0.193910725f, 0.278537078f, 0.588252033f, 1.458154888f, 3.532427913f, 6.005137088f, 5.737450374f, 5.359044089f, 5.381427007f, 5.813452436f, 6.575417761f, 7.705438609f, 7.916690521f, 7.677146788f, 9.108801829f, 13.80447353f, 24.75270462f, 41.27705918f, 60.39867474f, 76.06474955f, 88.87979307f, 96.26216052f, 99.57673469f, 100.9185254f, 100.3461882f, 99.74130238f, 99.13182504f, 97.91270568f, 94.24873818f, 93.80388811f, 91.40173503f, 91.7974507f, 87.7824115f, 83.02342603f, 78.93014979f, 75.53050128f, 71.07236035f, 66.33084727f, 62.33048663f, 57.65900798f, 53.17922641f, 46.68470127f, 38.98565424f, 31.84232286f, 25.79717595f, 21.2811326f, 18.08210481f, 15.49750967f, 13.51297898f, 12.03521208f, 10.76500507f, 9.648401633f, 9.001764638f, 8.306278847f, 7.876764508f, 7.782906577f, 7.696336862f, 7.182644064f, 5.845316604f, 3.650785487f, 1.758944979f, 0.790324641f, 0.222581918f, 0.47918518f, 0.148833141f, 0.067974729f, 0.0f, 0.0f, 0.08494f, 0.038172123f, 0.0f, 0.151692223f, 0.109050213f, 0.0f, 0.052663051f, 0.05610561f, 0.0f, 0.039383521f, 0.012867133f },
            { 0.086784117f, 0.178933235f, 0.197332679f, 0.393966858f, 0.45428447f, 1.464045507f, 3.618430454f, 6.061960647f, 5.716797204f, 5.318763808f, 5.390928372f, 5.800888846f, 6.581130548f, 7.679843397f, 7.913064002f, 7.629646645f, 9.081202386f, 13.57877476f, 24.60507986f, 41.22216022f, 60.21885388f, 75.87256068f, 88.7018612f, 95.38997103f, 98.97181493f, 100.4655994f, 100.0f, 99.10633032f, 98.62121661f, 97.4631588f, 93.51322301f, 92.99337375f, 91.00956564f, 91.34598641f, 87.04829152f, 82.4378642f, 78.20990485f, 75.07671731f, 70.59057657f, 65.71727304f, 61.78748638f, 57.01478928f, 52.78974423f, 46.20886877f, 38.94513661f, 31.72749718f, 25.85208614f, 21.47171887f, 18.04585094f, 15.67785085f, 13.84362684f, 12.17832015f, 10.88988692f, 9.949719367f, 9.341939294f, 8.413690141f, 8.125928029f, 7.954156868f, 8.002952873f, 7.275370079f, 5.977828115f, 3.722152948f, 1.873340464f, 0.858309557f, 0.556454795f, 0.403780406f, 0.17989397f, 0.0f, 0.022351329f, 0.111876068f, 0.0f, 0.194950484f, 0.0f, 0.0f, 0.098007154f, 0.0f, 0.028357028f, 0.036947597f, 0.036010439f, 0.029537641f, 0.105796423f },
            { 0.115446761f, 0.360557195f, 0.411775128f, 0.630346973f, 1.733697876f, 7.908836362f, 24.10727119f, 48.55778102f, 52.3696776f, 56.42172066f, 61.44934779f, 63.17636296f, 64.99269672f, 71.69218916f, 74.71998728f, 73.2150862f, 72.7353745f, 74.07229529f, 77.67671726f, 74.72442961f, 70.93914314f, 66.68841033f, 63.14004175f, 57.73010798f, 51.40878974f, 44.8220131f, 37.66176516f, 31.13222693f, 26.1618525f, 22.48856012f, 19.28181807f, 17.56674231f, 16.16629743f, 14.94982697f, 13.31452119f, 11.90575832f, 10.80686224f, 10.33615624f, 9.83413223f, 9.68951131f, 9.57033808f, 9.631085486f, 9.397577684f, 8.631380248f, 7.754938943f, 6.68241832f, 5.9614169f, 5.344294019f, 4.971847899f, 4.643173526f, 4.472756229f, 4.32018531f, 4.379411663f, 4.423656036f, 4.439993914f, 4.368837648f, 4.449112348f, 4.472267212f, 4.467931511f, 4.010400155f, 3.065692786f, 1.839411352f, 0.917909373f, 0.468340528f, 0.17121686f, 0.1812147f, 0.055650653f, 0.044445015f, 0.0f, 0.053274318f, 0.0f, 0.0f, 0.078894695f, 0.0f, 0.013803824f, 0.0f, 0.010802677f, 0.150527246f, 0.157891924f, 0.0f, 0.0f }
    });

    template <typename T, typename U>
    bool has_item(T container, U item)
    {
        return std::find(container.begin(), container.end(), item) != container.end();
    }

    template <typename T>
    bool is_between_range(T value, T lower_limit, T upper_limit)
    {
        return value >= lower_limit && value <= upper_limit;
    }

    SCENARIO("Create module spectral response object and query supported illuminants")
    {
        GIVEN("valid module spectral response object")
        {
            module_spectral_response msr;
            WHEN("quering model light source spectrum identifiers")
            {
                auto result = msr.query_light_sources();
                THEN("result should have 16 illuminant spectrums")
                {
                    CHECK(result.size() == 16);
                }
                AND_THEN("result should include predefined standard illuminants")
                {
                    CHECK(has_item(result, standard_illuminant_A));
                    CHECK(has_item(result, standard_illuminant_F12));
                    CHECK(has_item(result, standard_illuminant_D50));
                    CHECK(has_item(result, standard_illuminant_D55));
                    CHECK(has_item(result, standard_illuminant_D65));
                    CHECK(has_item(result, standard_illuminant_D75));
                }
            }
        }
    }

    SCENARIO("Create module spectral response object and get sensor responses for color checker classic ")
    {
        GIVEN("valid module spectral response object")
        {
            module_spectral_response msr;
            WHEN("user calculates sensor response before setting module spectral response data")
            {
                THEN("model should throw expection")
                {
                    CHECK_THROWS(msr.get_color_checker_classic(standard_illuminant_A));
                }
            }

            msr.set_module_spectral_response_data(msr_data);
            WHEN("quering model light source spectrum identifiers")
            {
                auto light_sources = msr.query_light_sources();
                THEN("sensor response for color checker classic can be calculated for each illuminant")
                {
                    for (auto& illuminant : light_sources)
                    {
                        auto sensor_response = msr.get_color_checker_classic(illuminant);
                        CHECK(sensor_response.size() == 24);
                    }
                }
            }
            WHEN("calculating the sensor response for unknown illuminant")
            {
                THEN("model should throw expection")
                {
                    CHECK_THROWS(msr.get_color_checker_classic("test"));
                }
            }
        }
    }

    SCENARIO("Create module spectral response object and get chromaticities")
    {
        GIVEN("valid module spectral response object")
        {
            module_spectral_response msr;
            WHEN("user calculates chromaticity before setting module spectral response data")
            {
                THEN("model should throw expection")
                {
                    CHECK_THROWS(msr.get_chromaticity(standard_illuminant_A));
                }
            }

            msr.set_module_spectral_response_data(msr_data);
            WHEN("quering model light source spectrum identifiers")
            {
                auto light_sources = msr.query_light_sources();
                THEN("calculate chromaticity for each illuminant")
                {
                    std::vector<chromaticity> chromaticities;
                    for (auto& illuminant : light_sources)
                    {
                        chromaticities.push_back(msr.get_chromaticity(illuminant));
                    }
                    THEN("chromaticity for r/g and b/g should be positive number")
                    {
                        for (auto &chromaticity : chromaticities)
                        {
                            CHECK(chromaticity._r_per_g > 0.0);
                            CHECK(chromaticity._b_per_g > 0.0);
                        }
                    }
                    AND_THEN("chromaticity for normal bayer sensor should have zero i/g ratio")
                    {
                        for (auto &chromaticity : chromaticities)
                        {
                            CHECK(chromaticity._i_per_g == 0.0);
                        }
                    }
                }
            }
            WHEN("calculating the sensor response for unknown illuminant")
            {
                THEN("model should throw expection")
                {
                    CHECK_THROWS(msr.get_chromaticity("test"));
                }
            }
        }
    }

    SCENARIO("Create module spectral response object and set new illuminant")
    {
        GIVEN("valid module spectral response object")
        {
            module_spectral_response msr;
            illuminant spectrum_f12 = illuminant(lightsource_type::F12,
                { spectral_power_distribution_f12, minimum_wavelenght_f12, maximum_wavelenght_f12 });
            WHEN("setting new illuminant values with existing id")
            {
                auto xyz_old = msr.get_xyz(standard_illuminant_A);
                msr.set_illuminant(standard_illuminant_A, spectrum_f12);

                THEN("the new values should differ to original")
                {
                    auto is_item_equal = [](const xyz<double> &a, const xyz<double> &b) {
                        return a.x == b.x && a.y == b.y && a.z == b.z;
                    };

                    auto xyz_new = msr.get_xyz(standard_illuminant_A);
                    CHECK_FALSE(is_item_equal(xyz_old, xyz_new));
                }
            }
            WHEN("setting new illuminant id")
            {
                auto light_sources = msr.query_light_sources();
                AND_WHEN("module spectral response does not contain new illuminant")
                {
                    auto illuminant_id = "foo";
                    CHECK_FALSE(has_item(light_sources, illuminant_id));

                    THEN("set operation should succeed")
                    {
                        CHECK_NOTHROW(msr.set_illuminant(illuminant_id, spectrum_f12));
                    }
                }
            }
        }
    }

    SCENARIO("Create module spectral response object and get illuminant XYZ")
    {
        GIVEN("valid module spectral response object")
        {
            module_spectral_response msr;
            WHEN("quering model light source spectrum identifiers")
            {
                auto light_sources = msr.query_light_sources();
                THEN("calculate XYZ for each illuminant")
                {
                    // Not optimal solution. Catch testframework does not allow
                    // dynamic test cases, e.g. for loops and continuing SECTOR
                    // So calculate all values and push to results vector
                    std::vector<xyz<double>> xyz;
                    for (auto& illuminant : light_sources)
                    {
                        xyz.push_back(msr.get_xyz(illuminant));
                    }
                    AND_THEN("all xyz values should be positive numbers")
                    {
                        for (auto &pt : xyz)
                        {
                            CHECK(pt.x > 0.0);
                            CHECK(pt.y > 0.0);
                            CHECK(pt.z > 0.0);
                        }

                        AND_THEN("cie xy calculated from xyz should be in the range of 0.29-0.5")
                        {
                            for (auto &pt : xyz)
                            {
                                CHECK(is_between_range(pt.x / (pt.x + pt.y + pt.z), 0.29, 0.5));
                                CHECK(is_between_range(pt.y / (pt.x + pt.y + pt.z), 0.29, 0.5));
                            }
                        }
                    }
                }
            }
            WHEN("calculating the sensor response for unknown illuminant")
            {
                THEN("model should throw expection")
                {
                    CHECK_THROWS(msr.get_xyz("test"));
                }
            }
        }
    }

    SCENARIO("Create module spectral response object and get illuminant lightsource type")
    {
        GIVEN("valid module spectral response object")
        {
            module_spectral_response msr;
            WHEN("get the light source type for standard illuminants from the model")
            {
                THEN("lightsource type should match to known types")
                {
                    CHECK(msr.get_light_source_type(standard_illuminant_A) == lightsource_type::A);
                    CHECK(msr.get_light_source_type(standard_illuminant_F12) == lightsource_type::F12);
                    CHECK(msr.get_light_source_type(standard_illuminant_D50) == lightsource_type::D50);
                    CHECK(msr.get_light_source_type(standard_illuminant_D55) == lightsource_type::D55);
                    CHECK(msr.get_light_source_type(standard_illuminant_D65) == lightsource_type::D65);
                    CHECK(msr.get_light_source_type(standard_illuminant_D75) == lightsource_type::D75);
                }
            }
        }
    }

    SCENARIO("Create module spectral response object and get illuminant Munsell color patch lab values")
    {
        GIVEN("valid module spectral response object")
        {
            module_spectral_response msr;
            WHEN("get Munsell color patches in Lab")
            {
                auto lab = msr.get_munsell_patches();
                THEN("we should have values for all 1600 Munsell color patches")
                {
                    CHECK(lab.size() == 1600);
                    AND_THEN("lab whitepoint is D65 and "
                        "l values should be in range [0,100] and "
                        "a values should be in range [-100,100] and "
                        "b values should be in range [-110,110]")
                    {
                        for (auto &pt : lab)
                        {
                            //CHECK(pt == wp::D65)
                            CHECK(is_between_range(pt.l, 0.0f, 100.0f));
                            CHECK(is_between_range(pt.a, -100.0f, 100.0f));
                            CHECK(is_between_range(pt.b, -110.0f, 110.0f));
                        }
                    }
                }
            }
        }
    }

    SCENARIO("Create module spectral response object and get sensor responses for Munsell color patches")
    {
        GIVEN("valid module spectral response object")
        {
            module_spectral_response msr;
            WHEN("user calculates sensor response for Munsell color patches "
                "before setting module spectral response data")
            {
                THEN("model should throw expection")
                {
                    CHECK_THROWS(msr.get_munsell_patches(standard_illuminant_A));
                }
            }

            msr.set_module_spectral_response_data(msr_data);
            WHEN("calculating sensor response for Munsell color patches for specific illuminant")
            {
                auto sensor_response = msr.get_munsell_patches(standard_illuminant_A);
                THEN("size of the sensor_response array should match predefined value")
                {
                    CHECK(sensor_response.size() == 1600);

                    AND_THEN("all individual patch values should be within range [0,1]")
                    {
                        for (auto& pt : sensor_response)
                        {
                            CHECK(is_between_range(pt.r, 0.0f, 1.0f));
                            CHECK(is_between_range(pt.g, 0.0f, 1.0f));
                            CHECK(is_between_range(pt.b, 0.0f, 1.0f));
                        }
                    }
                }
            }
            WHEN("calculating the sensor response for unknown illuminant")
            {
                THEN("model should throw expection")
                {
                    CHECK_THROWS(msr.get_munsell_patches("test"));
                }
            }
        }
    }
}