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


#pragma once
#include <cstdint>
#include <type_traits>

#if (defined(WIN32) || defined(_WIN32))
#include "intrin.h"
#else
// Linux
#include "x86intrin.h"
#endif

namespace Teisko
{
#ifdef NO_64_BIT_REGS
    // popcount (i.e. count number of bits set in an integer)
    // defined for all 8,16,32 bit types
    template <typename T>
    typename std::enable_if<std::is_integral<T>::value && sizeof(T) <= 4, T>::type
        popcount(T i) {
        typedef typename std::make_unsigned<T>::type utype;
        auto x = static_cast<utype>(i);
        return static_cast<T>(_mm_popcnt_u32(static_cast<uint32_t>(x)));
    }
    template <typename T>
    typename std::enable_if<std::is_integral<T>::value && sizeof(T) == 8, T>::type
        popcount(T i) {
        typedef typename std::make_unsigned<T>::type utype;
        auto x = static_cast<utype>(i);
        auto lsb_count = _mm_popcnt_u32(static_cast<uint32_t>(x & 0xFFFFFFFF));
        auto msb_count = _mm_popcnt_u32(static_cast<uint32_t>(x >> 32) & 0xffffffff);
        return static_cast<T>(lsb_count + msb_count);
    }
#else
    template <typename T>
    typename std::enable_if<std::is_integral<T>::value && (sizeof(T) <= 8), T>::type
        popcount(T i) {
        typedef typename std::make_unsigned<T>::type utype;
        auto x = static_cast<utype>(i);
        return static_cast<T>(_mm_popcnt_u64(static_cast<uint64_t>(x)));
    }
#endif

    /// container allowing quick erasing, setting, looping and testing of items
    /// The input argument for all set/erase/toggle etc items is never checked for proper range
    template <typename T = uint32_t>
    struct small_bit_set
    {
        T mask;         ///!< Container of bits
        small_bit_set(T initial_mask = 0) : mask(initial_mask) { }

        // Sets individual bit by index
        small_bit_set& set(int i) {
            mask |= static_cast<T>(1) << i;
            return *this;
        }

        // Clears individual bit by index
        small_bit_set& clear(int i) {
            mask &= ~(static_cast<T>(1) << i); return *this;
        }

        // Reverses individual bits by index
        small_bit_set& toggle(int i) {
            mask ^= static_cast<T>(1) << i;
            return *this;
        }

        // Bit hack to isolate rightmost bit -- originally x & (-x)
        // where -x == ~x + 1, because negating unsigned gives warning on MSVC
        T isolate_lsb() { return mask & (~mask + 1); }

        // Tests individual bit by index
        bool is_set(int i) { return (mask >> i) & 1; }

        // return number of bits in set
        T size() { return popcount(mask); }

        // Returns index of the lowest bit set in mask
        int lowest()
        {
            if (mask == 0)
                return -1;
            auto lsb = isolate_lsb();
            return static_cast<int>(popcount(--lsb));
        }

        // Bit hack to remove smallest set bit from mask
        small_bit_set& remove_lowest()
        {
            auto msk_1 = mask;
            mask &= --msk_1;
            return *this;
        }

        operator T() { return mask; }

        // API to allow small_bit_set to be iterated
        small_bit_set begin() { return *this; }
        small_bit_set end() { return{}; }
        small_bit_set& operator++() { return remove_lowest(); }
        int operator *() { return lowest(); }
        bool operator !=(const small_bit_set &rhs) { return mask != rhs.mask; }
    };

    template <int N>
    struct cyclical_index
    {
        int idx;   // index guaranteed to be between 0..N-1
        cyclical_index(int i) : idx(i % N) { }
        cyclical_index& operator ++()
        {
            idx = idx == (N - 1) ? 0 : idx + 1;
            return *this;
        }
        cyclical_index& operator --()
        {
            idx = idx == 0 ? N - 1 : idx - 1;
            return *this;
        }
        operator int() { return idx; }
    };

}