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

// In alphabetical order
#include <algorithm>
#include <limits>
#include <type_traits>

namespace Teisko
{
    // Most internal calculations are done as doubles, but we need to clip the internal result
    // to a final type of T. We use enable_if to allow rounding and clipping for integral types
    // and static casting for floating point types
    // See -- https://stackoverflow.com/a/12073915/1716339
    // For the casual reader -- the return type is defined as the last parameter of enable_if
    //  - without the `, T` part the reduce_to return value would be `void`
    template<typename T, typename FloatType>
    typename std::enable_if<std::is_floating_point<T>::value &&
        std::is_floating_point<FloatType>::value, T>::type reduce_to(FloatType t)
    {
        return (T)t;
    }

    // For integral types T only
    template<typename T, typename FloatType>
    typename std::enable_if<std::is_integral<T>::value &&
        std::is_floating_point<FloatType>::value, T>::type reduce_to(FloatType t)
    {
        return (T)std::max((FloatType)0.0, t + (FloatType)0.5);
    }

    // For integral types T only
    template<typename T>  // T <- T do nothing
    typename std::enable_if<std::is_integral<T>::value, T>::type reduce_to(T t) { return t; }

    /// Most internal calculations are done as doubles, but we need to clip the internal result
    /// to a final type of T. We use enable_if to allow rounding and clipping for integral types
    /// and static casting for floating point types
    /// See -- https://stackoverflow.com/a/12073915/1716339
    /// For the casual reader -- the return type is defined as the last parameter of enable_if
    ///  - without the `, T` part the reduce_to return value would be `void`
    /// For following conversion:
    /// FloatingType -> FloatingType
    template<typename T, typename FloatType>
    typename std::enable_if<std::is_floating_point<T>::value &&
        std::is_floating_point<FloatType>::value, T>::type clamp_to(FloatType t)
    {
        return (T)(t < std::numeric_limits<T>::lowest() ? std::numeric_limits<T>::lowest() : std::numeric_limits<T>::max() < t ? std::numeric_limits<T>::max() : t);
    }

    /**
    * @brief Clamp value between low and high boundary
    *
    * @tparam T0    Type of v, numeric type
    * @tparam T1    Type of lo and high, numeric type
    * @param v      The value to clamp
    * @param lo     Low boundary to clamp
    * @param hi     High boundary to clamp
    * @return       Clamped value
    */
    template <typename T0, typename T1>
    inline T0 clamp(T0 v, T1 lo, T1 hi)
    {
        return v < lo ? (T0)lo : hi < v ? (T0)hi : v;
    }

    /// For following conversions:
    /// IntegralType -> IntegralType
    /// IntegralType -> FloatingType
    template<typename T, typename IntegralType>
    typename std::enable_if<(std::is_integral<T>::value || std::is_floating_point<T>::value) &&
        std::is_integral<IntegralType>::value, T>::type clamp_to(IntegralType t)
    {
        return (T)clamp(t, std::numeric_limits<T>::lowest(), std::numeric_limits<T>::max());
    }

// IQS does not build correctly when the below function is enabled.
// So guard it with preprocessor MACRO that is defined in the project that
// does not work. This can be thought of either as short term solution until
// floating point to integer conversion is needed or long term solution if
// this functionality is not needed in IQS at all.
#ifndef REDUCED_FUNCTIONALITY
    /// For following conversion:
    /// FloatingType -> IntegralType
    template<typename T, typename FloatType>
    typename std::enable_if<std::is_integral<T>::value &&
        std::is_floating_point<FloatType>::value, T>::type clamp_to(FloatType t)
    {
        return (T)clamp(t + (FloatType)0.5, std::numeric_limits<T>::lowest(), std::numeric_limits<T>::max());
    }
#endif

}
