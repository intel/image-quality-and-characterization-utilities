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


#include "Teisko/Algorithm/Bit.hpp"
#include "catch.hpp"
#include <random>

template <typename T>
std::vector<T> make_shuffled(std::vector<T> &src)
{
    auto shuffled = src;
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(shuffled.begin(), shuffled.end(), g);
    return shuffled;
}

/// \page libcalc_specs_libbit Specs: Teisko low level
///
/// \snippet this snippet-specs-low_level

using namespace Teisko;

SCENARIO("Libbit low level library provides fast interfaces to bit counting")
{
    GIVEN("A few 8 bit bit-patterns with one bit set")
    {
        auto one_bits = std::vector<int>{ 1, 2, 4, 8, 128 };

        WHEN("The bit counts are calculated interpreting the values as either int8_t or uint8_t")
        {
            auto vec_a = std::vector<int>();
            auto vec_b = std::vector<int>();
            for (auto &x : one_bits)
            {
                auto a = static_cast<int8_t>(x);
                auto b = static_cast<uint8_t>(x);
                vec_a.emplace_back(popcount(a));
                vec_b.emplace_back(popcount(b));
            }
            THEN("Both vectors produce show only one bit set")
            {
                CHECK(vec_a == std::vector<int>(one_bits.size(), 1));
                CHECK(vec_b == std::vector<int>(one_bits.size(), 1));
            }
        }
    }
}

SCENARIO("Libbit low level library provides a quick implementation of a 'set' using bits")
{
    GIVEN("A container for maximum 32 items and four items to added to set")
    {
        auto my_set = small_bit_set<uint32_t>{};
        auto ref_values = std::vector<int>{ 2, 3, 9, 11};

        WHEN("The set is added with bits in no particular order")
        {
            for (auto &x : make_shuffled(ref_values))
                my_set.set(x);

            THEN("The bitset returns a size of four")
            {
                CHECK(4 == my_set.size());
            }

            THEN("The set queried for items not present in original list gives false")
            {
                for (auto i : { 0, 1, 4, 5, 6, 7, 8, 10, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31 })
                {
                    CHECK(my_set.is_set(i) == false);
                }
                AND_THEN("The set can be iterated for the same bits in numerical order")
                {
                    std::vector<int> indices;
                    for (auto i : my_set)
                    {
                        indices.emplace_back(i);
                    }
                    CHECK(indices == ref_values);
                }
                AND_THEN("The set can be iterated by toggling the rightmost bit out")
                {
                    auto ref_masks = std::vector<uint32_t>();  // 1 << 2, 1 << 3, 1 << 9, 1 << 11
                    for (auto x : ref_values)
                        ref_masks.emplace_back(1 << x);

                    auto lowest_bits = std::vector<uint32_t>{};
                    while (my_set)
                    {
                        auto low_bit_mask = my_set.isolate_lsb();
                        lowest_bits.emplace_back(low_bit_mask);
                        my_set.mask ^= low_bit_mask;    // toggle the lowest bit out
                    }
                    CHECK(ref_masks == lowest_bits);
                }
            }
        }
    }
}