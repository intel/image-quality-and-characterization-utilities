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

#include "Teisko/Algorithm/Pow2.hpp"
#include "Teisko/Algorithm/ReduceTo.hpp"
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <ostream>
#include <stdint.h>

namespace Teisko
{
    /// Color conversion matrices
    static const std::array<double, 9> _rgb2xyz_BT601 = { 0.430619033509700, 0.341541912257496, 0.178309054232804, 0.222037939153439, 0.706638439153439, 0.071323621693122, 0.020185267195767, 0.129550380511464, 0.939094352292769 };
    static const std::array<double, 9> _rgb2xyz_BT709 = { 0.412456439089692, 0.357576077643909, 0.180437483266399, 0.212672851405623, 0.715152155287818, 0.072174993306560, 0.019333895582329, 0.119192025881303, 0.950304078536368 };
    static const std::array<double, 9> _rgb2xyz_BT2100 = { 0.637010191411101, 0.144615027396969, 0.168844781191930, 0.262721717361641, 0.677989275502262, 0.059289007136097, 0.000000000000000, 0.028072328847647, 1.060757671152353 };
    static const std::array<double, 9> _rgb2xyz_DCI_P3 = { 0.486632650000000, 0.265663162500000, 0.198174187500000, 0.229003600000000, 0.691726725000000, 0.079269675000000, -0.000000000000000, 0.045112612500000, 1.043717387500000 };
    static const std::array<double, 9> _xyz2rgb_BT601 = { 3.062897123222697, -1.393179136493678, -0.475751671257954, -0.969266030505187, 1.876010845446694, 0.041556017530350, 0.067877509951752, -0.228854773990332, 1.069348968256285 };
    static const std::array<double, 9> _xyz2rgb_BT709 = { 3.240454162114104, -1.537138512797716, -0.498531409556016, -0.969266030505187, 1.876010845446694, 0.041556017530350, 0.055643430959115, -0.204025913516754, 1.057225188223179 };
    static const std::array<double, 9> _xyz2rgb_BT2100 = { 1.716510669761973, -0.355641669986716, -0.253345541821907, -0.666693001182624, 1.616502208346911, 0.015768750389995, 0.017643638767459, -0.042779781669045, 0.942305072720019 };
    static const std::array<double, 9> _xyz2rgb_DCI_P3 = { 2.493180755328967, -0.931265525497140 ,-0.402659723758882, -0.829503115821079, 1.762694121119792, 0.023625088741740, 0.035853625780072, -0.076188954782652, 0.957092621518022 };

    /// Bradford matrices
    static const std::array<double, 9> _bradford = { 0.8951000,  0.2664000, -0.1614000, -0.7502000,  1.7135000,  0.0367000, 0.0389000, -0.0685000,  1.0296000 };
    static const std::array<double, 9> _bradford_inv = { 0.9869929, -0.1470543, 0.1599627, 0.4323053,  0.5183603, 0.0492912, -0.0085287,  0.0400428, 0.9684867 };

    /// Tristimulus values for each illuminant
    static const std::array<double, 3> _tristimulus_A = { 1.09850, 1.00000, 0.35585 };
    static const std::array<double, 3> _tristimulus_B = { 0.99072, 1.00000, 0.85223 };
    static const std::array<double, 3> _tristimulus_C = { 0.98074, 1.00000, 1.18232 };
    static const std::array<double, 3> _tristimulus_D50 = { 0.96422, 1.00000, 0.82521 };
    static const std::array<double, 3> _tristimulus_D55 = { 0.95682, 1.00000, 0.92149 };
    static const std::array<double, 3> _tristimulus_D65 = { 0.95047, 1.00000, 1.08883 };
    static const std::array<double, 3> _tristimulus_D75 = { 0.94972, 1.00000, 1.22638 };
    static const std::array<double, 3> _tristimulus_E = { 1.00000, 1.00000, 1.00000 };
    static const std::array<double, 3> _tristimulus_F2 = { 0.99186, 1.00000, 0.67393 };
    static const std::array<double, 3> _tristimulus_F7 = { 0.95041, 1.00000, 1.08747 };
    static const std::array<double, 3> _tristimulus_F11 = { 1.00962, 1.00000, 0.64350 };

    /* hcl2rgb inverse matrix precalculated using constants weights of 1/3
     *     |  1    -0.5       -0.5    |
     * M = |  0  sqrt(3)/2 -sqrt(3)/2 |
     *     | 1/3    1/3        1/3    |
     */
    static const std::array<double, 9> _hcl2rgb_inv = { 0.6666666666666667, 0.0000000000000000, 1.0000000000000000, -0.3333333333333333, 0.5773502691896258, 1.0000000000000000, -0.3333333333333333, -0.5773502691896258, 1.0000000000000000 };

    /// Constants for XYZ to Lab conversion
    static const double _epsilon_xyz2lab = 216.0 / 24389.0;
    static const double _kappa_xyz2lab = 24389.0 / 27.0;

    // Constant for pi
    static const double _pi = 3.141592653589793;

    /// RGB color space enum type
    /// Lists supported RGB color spaces
    enum class rgb_cs
    {
        BT601,
        BT709,
        BT2020,
        BT2100,
        sRGB,
        DCI_P3,
        Sensor
    };

    /// Whitepoint enum type
    /// Lists supported whitepoints for rgb,xyz and lab classes
    enum class wp
    {
        A,
        B,
        C,
        D50,
        D55,
        D65,
        D75,
        E,
        F1,
        F2,
        F3,
        F4,
        F5,
        F6,
        F7,
        F8,
        F9,
        F10,
        F11,
        F12,
        None
    };

    /// Gets tristimulus values by whitepoint
    inline std::array<double, 3> get_tristimulus(wp whitepoint)
    {
        switch (whitepoint)
        {
        case wp::A:
            return _tristimulus_A;
        case wp::B:
            return _tristimulus_B;
        case wp::C:
            return _tristimulus_C;
        case wp::D50:
            return _tristimulus_D50;
        case wp::D55:
            return _tristimulus_D55;
        case wp::D65:
            return _tristimulus_D65;
        case wp::D75:
            return _tristimulus_D75;
        case wp::E:
            return _tristimulus_E;
        case wp::F2:
            return _tristimulus_F2;
        case wp::F7:
            return _tristimulus_F7;
        case wp::F11:
            return _tristimulus_F11;
        default:
            throw std::invalid_argument("Unknown tristimulus values");
        }
    }

    // Modulo implementation with always positive output
    template <typename T>
    T mod_p(T m, T n)
    {
        if (n < 0)
            throw std::invalid_argument("denom must be positive");

        return m >= 0 ? (T)std::fmod(m, n) : (T)std::fmod((n - std::abs(std::fmod(m, n))), n);
    }

    /**
    * @brief Matrix multiplication for matrices A x B = C
    *
    * @tparam T0 Type of A, numeric type
    * @tparam T1 Type of B, numeric type
    * @tparam T2 Type of C, numeric type
    * @param A 3x3 matrix
    * @param B 3x1 matrix
    * @param C 3x1 matrix
    */
    template <typename T0, typename T1, typename T2>
    void matrix_mul3x1(const std::array<T0, 9>& A, const std::array<T1, 3>& B, std::array<T2, 3>& C)
    {
        C[0] = clamp_to<T2>(A[0] * B[0] + A[1] * B[1] + A[2] * B[2]);
        C[1] = clamp_to<T2>(A[3] * B[0] + A[4] * B[1] + A[5] * B[2]);
        C[2] = clamp_to<T2>(A[6] * B[0] + A[7] * B[1] + A[8] * B[2]);
    }

    /**
    * @brief Matrix multiplication for matrices A x B = C
    *
    * @tparam T0 Type of A, numeric type
    * @tparam T1 Type of B, numeric type
    * @tparam T2 Type of C, numeric type
    * @param A 3x3 matrix
    * @param B 3x3 matrix
    * @param C 3x3 matrix
    */
    template <typename T0, typename T1, typename T2>
    void matrix_mul3x3(const std::array<T0, 9>& A, const std::array<T1, 9>& B, std::array<T2, 9>& C)
    {
        C[0] = clamp_to<T2>(A[0] * B[0] + A[1] * B[3] + A[2] * B[6]);
        C[1] = clamp_to<T2>(A[0] * B[1] + A[1] * B[4] + A[2] * B[7]);
        C[2] = clamp_to<T2>(A[0] * B[2] + A[1] * B[5] + A[2] * B[8]);
        C[3] = clamp_to<T2>(A[3] * B[0] + A[4] * B[3] + A[5] * B[6]);
        C[4] = clamp_to<T2>(A[3] * B[1] + A[4] * B[4] + A[5] * B[7]);
        C[5] = clamp_to<T2>(A[3] * B[2] + A[4] * B[5] + A[5] * B[8]);
        C[6] = clamp_to<T2>(A[6] * B[0] + A[7] * B[3] + A[8] * B[6]);
        C[7] = clamp_to<T2>(A[6] * B[1] + A[7] * B[4] + A[8] * B[7]);
        C[8] = clamp_to<T2>(A[6] * B[2] + A[7] * B[5] + A[8] * B[8]);
    }

    // Forward declarations
    template <typename T, wp WP>
    struct xyz;
    template <typename T, wp WP>
    struct lab;
    template <typename T>
    struct hcl;

    /**
     * @brief rgb class
     *
     * RGB color class supports decimal and integral value types for
     * red, green and blue components. Class is parametrized by data type,
     * illuminant white point and rgb color space.
     *
     * Class supports conversion routines between different white point and
     * color space conversions. Also conversion routines from xyz, lab and
     * hcl is supported.
     *
     * @tparam T Data type of rgb values
     * @tparam CS Object colorspace
     * @tparam WP Object whitepoint
     */
    template <typename T, rgb_cs CS = rgb_cs::sRGB, wp WP = wp::D65>
    struct rgb
    {
        T r = (T)0;
        T g = (T)0;
        T b = (T)0;

        /**
         * @brief Construct a new rgb object
         */
        rgb() = default;

        /**
         * @brief Construct a new rgb object
         *
         * @param _r Numeric value for red component
         * @param _g Numeric value for green component
         * @param _b Numeric value for blue component
         */
        rgb(T _r, T _g, T _b) : r(_r), g(_g), b(_b)
        {
        }

        /**
         * @brief Construct a new rgb object
         *
         * @tparam N Type to construct rgb object from
         *           Supported types: rgb, xyz and lab
         * @param other Object to copy construct from
         */
        template <typename N>
        rgb(N&& other) : r((T)0), g((T)0), b((T)0)
        {
            convert(std::forward<N>(other));
        }

        /**
         * @brief Copy assignment operator for a new rgb object
         *
         * @tparam T1   Data type of input object
         * @tparam CS1  Colorspace value of other
         * @tparam WP1  Whitepoint value of other
         * @param other Object to copy
         * @return rgb<T, WP, CS>&
         */
        template <typename T1, rgb_cs CS1, wp WP1>
        rgb<T, CS, WP> & operator=(const rgb<T1, CS1, WP1>& other)
        {
            convert(other);
            return *this;
        }

        /**
         * @brief Move assignment operator for a new rgb object
         *
         * @tparam T1   Data type of input object
         * @tparam CS1  Colorspace value of other
         * @tparam WP1  Whitepoint value of other
         * @param other Object to move
         * @return rgb<T, WP, CS>&
         */
        template <typename T1, rgb_cs CS1, wp WP1>
        rgb<T, CS, WP> & operator=(rgb<T1, CS1, WP1>&& other)
        {
            convert(other);
            return *this;
        };

        /**
         * @brief Equality operator for a rgb object
         *
         * @tparam T1   Data type of input object
         * @tparam CS1  Colorspace value of other
         * @tparam WP1  Whitepoint value of other
         * @param other Object to compare
         * @return true when objects rgb values, color space and whitepoint are equal
         */
        template <typename T1, rgb_cs CS1, wp WP1>
        bool operator==(const rgb<T1, CS1, WP1>& other)
        {
            return r == other.r && g == other.g && b == other.b &&
                CS == CS1 && WP == WP1;
        }

        /**
        * @brief Equality operator for a rgb object
        *
        * @tparam T1   Data type of input object
        * @tparam CS1  Colorspace value of other
        * @tparam WP1  Whitepoint value of other
        * @param other Object to compare
        * @return true when objects rgb values, color space and whitepoint are equal
        */
        template <typename T1, rgb_cs CS1, wp WP1>
        bool operator==(rgb<T1, CS1, WP1>&& other)
        {
            return (*this) == other;
        }

        T max() const { return std::max({ r, g, b }); }

        T mean() const { return (r + g + b) / 3.0; }

        /**
         * @brief Operator overload *= for rgb object
         */
        template <typename T1>
        inline rgb<T, CS, WP>& operator *= (T1 v)
        {
            r = clamp_to<T>(r * v);
            g = clamp_to<T>(g * v);
            b = clamp_to<T>(b * v);
            return *this;
        }

        /**
         * @brief Operator overload /= for rgb object
         */
        template <typename T1>
        inline rgb<T, CS, WP>& operator /= (T1 v)
        {
            r = clamp_to<T>(r / v);
            g = clamp_to<T>(g / v);
            b = clamp_to<T>(b / v);
            return *this;
        }

        /**
         * @brief Operator overload * for rgb object
         */
        template <typename T1>
        inline rgb<T, CS, WP> operator * (T1 v)
        {
            return rgb<T, CS, WP> {clamp_to<T>(r * v), clamp_to<T>(g * v), clamp_to<T>(b * v)};
        }

        /**
         * @brief Operator overload / for rgb object
         */
        template <typename T1>
        inline rgb<T, CS, WP> operator / (T1 v)
        {
            return rgb<T, CS, WP> {clamp_to<T>(r / v), clamp_to<T>(g / v), clamp_to<T>(b / v)};
        }

    private:
        /**
         * @brief Conversion routine for rgb to rgb when illuminants
         * and rgb color spaces are same
         *
         * @tparam T1   Data type of input object
         * @tparam CS1  Colorspace value of other
         * @tparam WP1  Whitepoint value of other
         * @param other Object from where rgb is constructed
         */
        template <typename T1, rgb_cs CS1, wp WP1>
        typename std::enable_if<CS == CS1 && WP == WP1, void>::type
            convert(const rgb<T1, CS1, WP1>& other)
        {
            r = clamp_to<T>(other.r);
            g = clamp_to<T>(other.g);
            b = clamp_to<T>(other.b);
        }

        /**
         * @brief Conversion routine for rgb to rgb when illuminants
         * and/or rgb color spaces are different
         *
         * Conversion process to follow:
         *  1. rgb_s to xyz_s
         *  2. xyz_s to xyz_d
         *  3. xyz_d to rgb_d
         *  where *_s is source whitepoint and *_d is destination whitepoint
         *
         * @tparam T1   Data type of input object
         * @tparam CS1  Colorspace value of other
         * @tparam WP1  Whitepoint value of other
         * @param other Object from where rgb is constructed
         */
        template <typename T1, rgb_cs CS1, wp WP1>
        typename std::enable_if<!(CS == CS1 && WP == WP1), void>::type
            convert(const rgb<T1, CS1, WP1>& other)
        {
            // Do chromatic adaptation and rgb color space adaptation
            *this = xyz<double, WP>(other);
        }

        /**
         * @brief Conversion routine for xyz to rgb
         *
         * Conversion process to follow:
         *  1. xyz_s to xyz_d (if not the same illuminant)
         *  2. xyz_d to rgb_d
         *  where *_s is source whitepoint and *_d is destination whitepoint
         *
         * @tparam T1   Data type of input object
         * @tparam WP1  Whitepoint value of other
         * @param other Object from where rgb is constructed
         */
        template <typename T1, wp WP1>
        void convert(const xyz<T1, WP1>& other)
        {
            // Do chromatic adaptation
            xyz<double, WP> xyz_value(other);

            /* RGB to XYZ conversion
            *  |R|       |X|
            *  |G| = M * |Y|
            *  |B|       |Z|
            *
            *  M is the predefined conversion matrix for the current RGB working space
            */

            std::array<double, 9> conversion_matrix;
            switch (CS)
            {
            case rgb_cs::BT601:
                conversion_matrix = _xyz2rgb_BT601;
            case rgb_cs::BT709:
            case rgb_cs::sRGB:
                conversion_matrix = _xyz2rgb_BT709;
                break;
            case rgb_cs::BT2020:
            case rgb_cs::BT2100:
                conversion_matrix = _xyz2rgb_BT2100;
                break;
            case rgb_cs::DCI_P3:
                conversion_matrix = _xyz2rgb_DCI_P3;
                break;
            case rgb_cs::Sensor:
                throw std::invalid_argument("Conversion not supported for rgb_cs::Sensor");
            default:
                throw std::invalid_argument("Unknown rgb color space");
            }

            std::array<decltype(xyz_value.x), 3> input{ xyz_value.x, xyz_value.y, xyz_value.z };
            std::array<T, 3> output{ 0,0,0 };

            matrix_mul3x1(conversion_matrix, input, output);
            r = output[0];
            g = output[1];
            b = output[2];
        }

        /**
         * @brief Conversion routine for lab to rgb
         *
         * Conversion process to follow:
         *  1. lab_s to xyz_s
         *  2. xyz_s to xyz_d (if not the same illuminant)
         *  3. xyz_d to rgb_d
         *  where *_s is source whitepoint and *_d is destination whitepoint
         *
         * @tparam T1   Data type of input object
         * @tparam WP1  Whitepoint value of other
         * @param other Object from where rgb is constructed
         */
        template <typename T1, wp WP1>
        void convert(const lab<T1, WP1>& other)
        {
            *this = xyz<double, WP>(other);
        }

        /**
         * @brief Conversion routine for hcl to rgb
         *
         * @tparam T1   Data type of input object
         * @param other Object from where rgb is constructed
         */
        template <typename T1>
        void convert(const hcl<T1>& other)
        {
            auto alpha = other.c * std::cos(_pi * other.h / 180.0);
            auto beta = other.c * std::sin(_pi * other.h / 180.0);

            r = _hcl2rgb_inv[0] * alpha + _hcl2rgb_inv[1] * beta + _hcl2rgb_inv[2] * other.l;
            g = _hcl2rgb_inv[3] * alpha + _hcl2rgb_inv[4] * beta + _hcl2rgb_inv[5] * other.l;
            b = _hcl2rgb_inv[6] * alpha + _hcl2rgb_inv[7] * beta + _hcl2rgb_inv[8] * other.l;
        }
    };

    /**
     * @brief xyz class
     *
     * XYZ color class supports decimal value types for x, y and z
     * components. Class is parametrized by data type and illuminant
     * white point.
     *
     * Class supports conversion routines between different white points.
     * Also conversion routines from rgb and lab is supported.
     *
     * @tparam T Data type of xyz values
     * @tparam WP Object whitepoint
     */
    template <typename T, wp WP = wp::D65>
    struct xyz
    {
        T x = (T)0;
        T y = (T)0;
        T z = (T)0;

        /**
         * @brief Construct a new xyz object
         *
         */
        xyz() = default;

        /**
         * @brief Construct a new xyz object
         *
         * @param _x Numeric value for x component
         * @param _y Numeric value for y component
         * @param _z Numeric value for z component
         */
        xyz(T _x, T _y, T _z) : x(_x), y(_y), z(_z)
        {
        }

        /**
         * @brief Construct a new xyz object
         *
         * @tparam N Type to construct xyz object from
                     Supported types: rgb, xyz and lab
         * @param other Object to copy construct from
         */
        template <typename N>
        xyz(N&& other) : x((T)0), y((T)0), z((T)0)
        {
            convert(std::forward<N>(other));
        }

        /**
         * @brief Copy assignment operator
         *
         * @tparam T1   Data type of input object
         * @tparam WP1  Whitepoint value of other
         * @param other Object to copy
         * @return xyz<T, WP>&
         */
        template <typename T1, wp WP1>
        xyz<T, WP> & operator=(const xyz<T1, WP1>& other)
        {
            convert(other);
            return *this;
        }

        /**
         * @brief Move assignment operator
         *
         * @tparam T1   Data type of input object
         * @tparam WP1  Whitepoint value of other
         * @param other Object to move
         * @return xyz<T, WP>&
         */
        template <typename T1, wp WP1>
        xyz<T, WP> & operator=(xyz<T1, WP1>&& other)
        {
            convert(other);
            return *this;
        };

        /**
         * @brief Operator overload *= for xyz object
         */
        template <typename T1>
        inline xyz<T, WP>& operator * (T1 value)
        {
            return xyz<T, WP> {clamp_to<T>(x * value), clamp_to<T>(y * value), clamp_to<T>(z * value)};
        }

        /**
         * @brief Operator overload *= for xyz object
         */
        template <typename T1>
        inline xyz<T, WP>& operator *= (T1 value)
        {
            x = (T)(x * value);
            y = (T)(y * value);
            z = (T)(z * value);
            return *this;
        }

    private:
        /**
         * @brief Conversion routine for rgb to xyz
         *
         * Conversion process to follow:
         *  1. rgb_s to xyz_s
         *  2. xyz_s to xyz_d (if not the same illuminant)
         *  where *_s is source whitepoint and *_d is destination whitepoint
         *
         * @tparam T1   Data type of input object
         * @tparam CS1  Colorspace value of other
         * @tparam WP1  Whitepoint value of other
         * @param other Object from where xyz is constructed
         */
        template <typename T1, rgb_cs CS1, wp WP1>
        void convert(const rgb<T1, CS1, WP1>& other)
        {
            /* XYZ to RGB conversion
            *  |X|       |R|
            *  |Y| = M * |G|
            *  |Z|       |B|
            *
            *  M is the predefined conversion matrix for the current RGB working space
            */

            std::array<double, 9> conversion_matrix;
            switch (CS1)
            {
            case rgb_cs::BT601:
                conversion_matrix = _rgb2xyz_BT601;
            case rgb_cs::BT709:
            case rgb_cs::sRGB:
                conversion_matrix = _rgb2xyz_BT709;
                break;
            case rgb_cs::BT2020:
            case rgb_cs::BT2100:
                conversion_matrix = _rgb2xyz_BT2100;
                break;
            case rgb_cs::DCI_P3:
                conversion_matrix = _rgb2xyz_DCI_P3;
                break;
            case rgb_cs::Sensor:
                throw std::invalid_argument("Conversion not supported for rgb_cs::Sensor");
            default:
                throw std::invalid_argument("Unknown rgb color space");
            }

            std::array<T1, 3> input{ other.r, other.g, other.b };
            std::array<T, 3> output{ 0, 0, 0 };
            matrix_mul3x1(conversion_matrix, input, output);

            // Chromatic adaptation
            *this = xyz<T, WP1>(output[0], output[1], output[2]);
        }

        /**
         * @brief Conversion routine for xyz to xyz when illuminants are same
         *
         * @tparam T1   Data type of input object
         * @tparam WP1  Whitepoint value of other
         * @param other Object from where xyz is constructed
         */
        template <typename T1, wp WP1>
        typename std::enable_if<WP == WP1, void>::type
            convert(const xyz<T1, WP1>& other)
        {
            x = clamp_to<T>(other.x);
            y = clamp_to<T>(other.y);
            z = clamp_to<T>(other.z);
        }

        /**
         * @brief Conversion routine for xyz to xyz when illuminants are not same
         *
         * Conversion process to follow:
         *  1. xyz_s to xyz_d
         *  where *_s is source whitepoint and *_d is destination whitepoint
         *
         * @tparam T1   Data type of input object
         * @tparam WP1  Whitepoint value of other
         * @param other Object from where xyz is constructed
         */
        template <typename T1, wp WP1>
        typename std::enable_if<WP != WP1, void>::type
            convert(const xyz<T1, WP1>& other)
        {
            /*
            * Adaptation matrix 3x3 = inverse(M) * WP * M
            *
            *      | Rdw/Rsw  0      0    |
            * WP = |    0  Gdw/Gsw   0    |
            *      |    0     0   Bdw/Bsw |
            *
            * |Rdw|       |Xdw|   |Rsw|       |Xsw|
            * |Gdw| = M * |Ydw| , |Gsw| = M * |Ysw| , sw = source white, dw = destination white
            * |Bdw|       |Zdw|   |Bsw|       |Zsw|
            */

            // Use fixed adaptation matrix (Bradford)
            auto source_wp = get_tristimulus(WP1);
            auto destination_wp = get_tristimulus(WP);

            // Initialize temporary objects
            std::array<double, 3> cone_response_source{ 0,0,0 };
            std::array<double, 3> cone_response_destination{ 0,0,0 };
            std::array<double, 9> cone_response{ 0,0,0,0,0,0,0,0,0 };
            std::array<double, 9> tmp{ 0,0,0,0,0,0,0,0,0 };
            std::array<double, 9> adaptation{ 0,0,0,0,0,0,0,0,0 };

            matrix_mul3x1(_bradford, source_wp, cone_response_source);
            matrix_mul3x1(_bradford, destination_wp, cone_response_destination);

            // Create diagonal matrix (WP)
            cone_response[0] = cone_response_destination[0] / cone_response_source[0];
            cone_response[4] = cone_response_destination[1] / cone_response_source[1];
            cone_response[8] = cone_response_destination[2] / cone_response_source[2];

            // Matrix product from inverse(M) * WP
            matrix_mul3x3(_bradford_inv, cone_response, tmp);

            // Matrix product from inverse(M) * WP * M
            matrix_mul3x3(tmp, _bradford, adaptation);

            std::array<T1, 3> input{ other.x, other.y, other.z };
            std::array<T1, 3> output{ 0,0,0 };

            // Convert
            matrix_mul3x1(adaptation, input, output);

            x = output[0];
            y = output[1];
            z = output[2];
        }

        /**
         * @brief Conversion routine for lab to xyz
         *
         * Conversion process to follow:
         *  1. lab_s to xyz_s
         *  2. xyz_s to xyz_d (if not the same illuminant)
         *  where *_s is source whitepoint and *_d is destination whitepoint
         *
         * @tparam T1   Data type of input object
         * @tparam WP1  Whitepoint value of other
         * @param other Object from where xyz is constructed
         */
        template <typename T1, wp WP1>
        void convert(const lab<T1, WP1>& other)
        {
            /* Lab to XYZ
            * X = xr * Xr                   , Xr = reference white X
            * Y = yr * Yr                   , Yr = reference white Y
            * Z = zr * Zr                   , Zr = reference white Z
            *
            * xr = fx^3                     , fx^3 > epsilon
            *    = (116 * fx - 16) / kappa  , fx^3 <= epsilon
            *
            * yr = ((L + 16)/116)^3         , L > kappa * epsilon
            *    = L / kappa                , L <= kappa * epsilon
            *
            * zr = fz^3                     , fz^3 > epsilon
            *    = (116 * fz - 16) / kappa  , fz^3 <= epsilon
            *
            * fx = a / 500 + fy
            *
            * fy = (L + 16)/116
            *
            * fz = fy - b / 200
            */

            auto fy = (other.l + 16.0) / 116.0;
            auto fx = other.a / 500.0 + fy;
            auto fz = fy - other.b / 200.0;

            auto fx3 = fx * fx * fx;
            auto fz3 = fz * fz * fz;

            auto temp = (other.l + 16.0) / 116.0;

            auto xr = fx3 > _epsilon_xyz2lab ? fx3 : (116.0 * fx - 16.0) / _kappa_xyz2lab;
            auto yr = other.l > _kappa_xyz2lab * _epsilon_xyz2lab ? temp * temp * temp : other.l / _kappa_xyz2lab;
            auto zr = fz3 > _epsilon_xyz2lab ? fz3 : (116.0 * fz - 16.0) / _kappa_xyz2lab;

            auto wp = get_tristimulus(WP1);

            // Conversion and chromatic adaptation
            *this = xyz<T, WP1>(
                clamp_to<T>(xr * wp[0]),
                clamp_to<T>(yr * wp[1]),
                clamp_to<T>(zr * wp[2]));
        }
    };

    /**
     * @brief lab class
     *
     * Lab color class supports decimal value types for l, a and b
     * components. Class is parametrized by data type and illuminant
     * white point.
     *
     * Class supports conversion routines between different white points.
     * Also conversion routines from xyz and lab is supported.
     *
     * @tparam T Data type of lab values
     * @tparam WP Object whitepoint
     */
    template <typename T, wp WP = wp::D65>
    struct lab
    {
        T l = (T)0;
        T a = (T)0;
        T b = (T)0;

        /**
         * @brief Construct a new lab object
         *
         */
        lab() = default;

        /**
         * @brief Construct a new lab object
         *
         * @param _l Numerical value for l component (lightness)
         * @param _a Numerical value for a component (green-red)
         * @param _b Numerical value for b component (blue-yellow)
         */
        lab(T _l, T _a, T _b) : l(_l), a(_a), b(_b)
        {
        }

        /**
         * @brief Construct a new lab object
         *
         * @tparam N Type to construct lab object from
         *           Supported types: rgb, xyz and lab
         * @param other Object to copy construct from
         */
        template <typename N>
        lab(N&& other) : l((T)0), a((T)0), b((T)0)
        {
            convert(std::forward<N>(other));
        }

        /**
         * @brief Copy assignment operator for a new lab object
         *
         * @tparam T1   Data type of input object
         * @tparam WP1  Whitepoint value of other
         * @param other Object to copy
         * @return lab<T, WP>&
         */
        template <typename T1, wp WP1>
        lab<T, WP> & operator=(const lab<T1, WP1>& other)
        {
            convert(other);
            return *this;
        }

        /**
         * @brief Move assignment operator for a new lab object
         *
         * @tparam T1   Data type of input object
         * @tparam WP1  Whitepoint value of other
         * @param other Object to move
         * @return lab<T, WP>&
         */
        template <typename T1, wp WP1>
        lab<T, WP> & operator=(lab<T1, WP1>&& other)
        {
            convert(other);
            return *this;
        };

    private:
        /**
         * @brief Conversion routine for rgb to lab
         * rgb -> xyz -> chromatic adaptation -> lab
         *
         * Conversion process to follow:
         *  1. rgb_s to xyz_s
         *  2. xyz_s to xyz_d (if not the same illuminant)
         *  3. xyz_d to lab_d
         *  where *_s is source whitepoint and *_d is destination whitepoint
         *
         *
         * @tparam T1   Data type of input object
         * @tparam CS1  Colorspace value of other
         * @tparam WP1  Whitepoint value of other
         * @param other Object from where lab is constructed
         */
        template <typename T1, rgb_cs CS1, wp WP1>
        void convert(const rgb<T1, CS1, WP1>& other)
        {
            *this = xyz<T1, WP>(other);
        }

        /**
         * @brief Conversion routine for xyz to lab
         *
         * Conversion process to follow:
         *  1. xyz_s to xyz_d (if not the same illuminant)
         *  2. xyz_d to lab_d
         *  where *_s is source whitepoint and *_d is destination whitepoint
         *
         * @tparam T1   Data type of input object
         * @tparam WP1  Whitepoint value of other
         * @param other Object from where lab is constructed
         */
        template <typename T1, wp WP1>
        void convert(const xyz<T1, WP1>& other)
        {
            /* xyz to lab
            * L = 116 * fy - 16
            * a = 500 * (fx - fy)
            * b = 200 * (fy - fz)
            *
            * fx = xr^(1/3)              , xr > epsilon
            *    = (kappa*xr + 16) / 116 , xr <= epsilon
            *
            * fy = yr^(1/3)              , yr > epsilon
            *    = (kappa*yr + 16) / 116 , yr <= epsilon
            *
            * fz = zr^(1/3)              , zr > epsilon
            *    = (kappa*zr + 16) / 116 , zr <= epsilon
            *
            * xr = X / Xr , Xr = reference white X
            * yr = Y / Yr , Yr = reference white Y
            * zr = Z / Zr , Zr = reference white Z
            */

            auto func = [](double value)
            {
                return value > _epsilon_xyz2lab ? std::cbrt(value) : (_kappa_xyz2lab * value + 16) / 116.0;
            };

            // Chromatic adaptation
            xyz<T1, WP> tmp(other);

            auto ref_wp = get_tristimulus(WP);

            auto xr = tmp.x / ref_wp[0];
            auto yr = tmp.y / ref_wp[1];
            auto zr = tmp.z / ref_wp[2];

            auto fx = func(xr);
            auto fy = func(yr);
            auto fz = func(zr);

            // Convert to Lab
            l = clamp_to<T>(116 * fy - 16);
            a = clamp_to<T>(500 * (fx - fy));
            b = clamp_to<T>(200 * (fy - fz));
        }

        /**
         * @brief Conversion routine for lab when illuminants are same
         *
         * @tparam T1   Data type of input object
         * @tparam WP1  Whitepoint value of other
         * @param other Object from where lab is constructed
         */
        template <typename T1, wp WP1>
        typename std::enable_if<WP == WP1, void>::type
            convert(const lab<T1, WP1>& other)
        {
            l = clamp_to<T>(other.l);
            a = clamp_to<T>(other.a);
            b = clamp_to<T>(other.b);
        }

        /**
         * @brief Conversion routine for lab when illuminants are not same
         *
         * Conversion process to follow:
         *  1. lab_s to xyz_s
         *  2. xyz_s to xyz_d
         *  3. xyz_d to lab_d
         *  where *_s is source whitepoint and *_d is destination whitepoint
         *
         * @tparam T1   Data type of input object
         * @tparam WP1  Whitepoint value of other
         * @param other Object from where lab is constructed
         */
        template <typename T1, wp WP1>
        typename std::enable_if<WP != WP1, void>::type
            convert(const lab<T1, WP1>& other)
        {
            *this = xyz<T1, WP>(other);
        }
    };

    /**
     * @brief hcl class
     *
     * Hcl color class supports decimal value types for h, c and l
     * components. Class is parametrized by data type.
     *
     * Class supports conversion routines between hcl and rgb.
     *
     * @tparam T Data type of lab values
     */
    template <typename T>
    struct hcl
    {
        T h;
        T c;
        T l;

        /**
         * @brief Construct a new hcl object
         *
         */
        hcl() = default;

        /**
         * @brief Construct a new hcl object
         *
         * @param _h Numerical value for h component (hue)
         * @param _c Numerical value for c component (chroma)
         * @param _l Numerical value for l component (luma)
         */
        hcl(T _h, T _c, T _l) : h(_h), c(_c), l(_l)
        {
        }

        /**
         * @brief Construct a new hcl object
         *
         * @tparam N Type to construct hcl object from
         *           Supported types: rgb
         * @param other Object to copy construct from
         */
        template <typename N>
        hcl(N&& other) : h((T)0), c((T)0), l((T)0)
        {
            convert(std::forward<N>(other));
        }

        /**
         * @brief Copy assignment operator for a new hcl object
         *
         * @tparam T1   Data type of input object
         * @param other Object to copy
         * @return hcl<T>&
         */
        template <typename T1>
        hcl<T> & operator=(const hcl<T1>& other)
        {
            convert(other);
            return *this;
        }

        /**
         * @brief Move assignment operator for a new hcl object
         *
         * @tparam T1   Data type of input object
         * @param other Object to move
         * @return hcl<T>&
         */
        template <typename T1>
        hcl<T> & operator=(hcl<T1>&& other)
        {
            convert(other);
            return *this;
        };

    private:
        /**
         * @brief Conversion routine for rgb to hcl
         *
         * @tparam T1   Data type of input object
         * @tparam CS  Colorspace value of other
         * @tparam WP  Whitepoint value of other
         * @param other Object from where lab is constructed
         */
        template <typename T1, rgb_cs CS, wp WP>
        void convert(const rgb<T1, CS, WP>& other)
        {
            auto alpha = other.r - 0.5 * other.g - 0.5 * other.b;
            auto beta = std::sqrt(3.0) / 2 * (other.g - other.b);

            h = clamp_to<T>(mod_p(std::atan2(beta, alpha) * 180 / _pi, 360.0));
            c = clamp_to<T>(std::sqrt(alpha * alpha + beta * beta));
            l = clamp_to<T>(other.mean());
        }

        /**
         * @brief Conversion routine for hcl
         *
         * @tparam T1   Data type of input object
         * @param other Object from where lab is constructed
         */
        template <typename T1>
        void convert(const hcl<T1>& other)
        {
            h = clamp_to<T>(other.h);
            c = clamp_to<T>(other.c);
            l = clamp_to<T>(other.l);
        }
    };

    /**
     * @brief A gamma class.
     *
     * Class uses standard gamma defined for given rgb colorspace to either
     * do forward gamma or reverse gamma.
     *
     * Input values should be floating point precision normalized to range [0,1]
     * Beyond this range the gamma is undefined.
     *
     */
    struct gamma
    {
        /**
         * @brief Forward gamma
         *
         * rgb colorspace standard gamma is applied
         *
         * @tparam T    Data type of input object
         * @tparam CS   Colorspace value of a
         * @tparam WP   Whitepoint value of a
         * @param a     Object to apply gamma
         */
        template <typename T, rgb_cs CS, wp WP>
        static typename std::enable_if<std::is_floating_point<T>::value, rgb<T, CS, WP>>::type
            set(rgb<T, CS, WP>&& a)
        {
            set(a);
            return a;
        }

        /**
         * @brief Forward gamma
         *
         * rgb colorspace standard gamma is applied
         *
         * @tparam T    Data type of input object
         * @tparam CS   Colorspace value of a
         * @tparam WP   Whitepoint value of a
         * @param a     Object to apply gamma
         * @return      rgb object with gamma
         */
        template <typename T, rgb_cs CS, wp WP>
        static typename std::enable_if<std::is_floating_point<T>::value, void>::type
            set(rgb<T, CS, WP>& a)
        {
            std::function<T(T)> eotf;
            switch (CS)
            {
            case rgb_cs::BT601:
            case rgb_cs::BT709:
                eotf = [](T x)
                {
                    if (x < 0.018)
                        return clamp_to<T>(4.5 * x);
                    return clamp_to<T>(1.099 * std::pow(x, 0.45) - 0.099);
                };
                break;
            case rgb_cs::BT2020:
                eotf = [](T x)
                {
                    if (x < 0.0181)
                        return clamp_to<T>(4.5 * x);
                    const auto a = 1.0993;
                    return clamp_to<T>(a * std::pow(x, 0.45) - (a - 1));
                };
                break;
            case rgb_cs::BT2100:
                eotf = [](T x)
                {
                    // 0.083333333333333 == 1.0 / 12
                    if (x <= 0.083333333333333)
                        return clamp_to<T>(std::sqrt(3.0 * x));

                    const auto a = 0.17883277;
                    const auto b = 0.02372241;
                    const auto c = 1.00429347;

                    return clamp_to<T>(a * std::log(x - b) + c);
                };
                break;
            case rgb_cs::sRGB:
                eotf = [](T x)
                {
                    if (x <= 0.0031308)
                        return clamp_to<T>(12.92 * x);
                    const double a = 0.055;
                    // 0.41666666666666669 == 1/2.4
                    return clamp_to<T>((1 + a) *std::pow(x, 0.41666666666666669) - a);
                };
                break;
            case rgb_cs::DCI_P3:
                eotf = [](T x)
                {
                    return clamp_to<T>(std::pow(x, 1 / 2.6));
                };
                break;
            default:
                throw std::invalid_argument("Unknown rgb color space");
            }

            a.r = eotf(a.r);
            a.g = eotf(a.g);
            a.b = eotf(a.b);
        }

        /**
         * @brief Reverse gamma
         *
         * rgb colorspace standard gamma is reversed
         *
         * @tparam T    Data type of input object
         * @tparam CS   Colorspace value of a
         * @tparam WP   Whitepoint value of a
         * @param a     Object to apply reverse gamma
         */
        template <typename T, rgb_cs CS, wp WP>
        static typename std::enable_if<std::is_floating_point<T>::value, rgb<T, CS, WP>>::type
            unset(rgb<T, CS, WP>&& a)
        {
            unset(a);
            return a;
        }

        /**
         * @brief Reverse gamma
         *
         * rgb colorspace standard gamma is reverted
         *
         * @tparam T    Data type of input object
         * @tparam CS   Colorspace value of a
         * @tparam WP   Whitepoint value of a
         * @param a     Object to apply reverse gamma
         * @return      rgb object with gamma reverted
         */
        template <typename T, rgb_cs CS, wp WP>
        static typename std::enable_if<std::is_floating_point<T>::value, void>::type
            unset(rgb<T, CS, WP>& a)
        {
            std::function<T(T)> eotf;
            switch (CS)
            {
            case rgb_cs::BT601:
            case rgb_cs::BT709:
                eotf = [](T x)
                {
                    if (x <= 0.081)
                        return clamp_to<T>(x / 4.5);
                    return clamp_to<T>(std::pow((x + 0.099) / 1.099, 1 / 0.45));
                };
                break;
            case rgb_cs::BT2020:
                eotf = [](T x)
                {
                    if (x <= 0.08145)
                        return clamp_to<T>(x / 4.5);
                    const auto a = 1.0993;
                    return clamp_to<T>(std::pow((x + (1 - a)) / a, 1 / 0.45));
                };
                break;
            case rgb_cs::BT2100:
                eotf = [](T x)
                {
                    const auto threshold = 0.5;
                    const auto a = 0.17883277;
                    const auto b = 0.02372241;
                    const auto c = 1.00429347;

                    return x <= threshold ?
                        clamp_to<T>(x * x / 3.0) :
                        clamp_to<T>(std::exp((x - c) / a) + b);
                };
                break;
            case rgb_cs::sRGB:
                eotf = [](T x)
                {
                    if (x <= 0.04045)
                        return clamp_to<T>(x / 12.92);
                    const auto a = 0.055;
                    return clamp_to<T>(std::pow((x + a) / (1 + a), 2.4));
                };
                break;
            case rgb_cs::DCI_P3:
                eotf = [](T x)
                {
                    return clamp_to<T>(std::pow(x, 1 / 2.6));
                };
                break;
            default:
                throw std::invalid_argument("Unknown rgb color space");
            }

            a.r = eotf(a.r);
            a.g = eotf(a.g);
            a.b = eotf(a.b);
        }

    private:
        /**
         * @brief Block creation of an instance of this object
         */
        gamma() {};
    };

    /**
     * @brief Operator overload for printing class properties in catch testframework
     */
    template <typename T, rgb_cs CS, wp WP>
    std::ostream& operator << (std::ostream& os, rgb<T, CS, WP> const& value)
    {
        // Cast enums to int for simplicity
        os << "r:" << value.r << " g:" << value.g << " b:" << value.b << " wp:"
            << (int)WP << " cs:" << (int)CS;
        return os;
    }

    /**
    * @brief Operator overload for printing class properties in catch testframework
    */
    template <typename T, wp WP>
    std::ostream& operator << (std::ostream& os, xyz<T, WP> const& value)
    {
        // Cast enum to int for simplicity
        os << "x:" << value.x << " y:" << value.y << " z:" << value.z << " wp:"
            << (int)WP;
        return os;
    }

    /**
    * @brief Operator overload for printing class properties in catch testframework
    */
    template <typename T, wp WP>
    std::ostream& operator << (std::ostream& os, lab<T, WP> const& value)
    {
        // Cast enum to int for simplicity
        os << "l:" << value.l << " a:" << value.a << " b:" << value.b << " wp:"
            << (int)WP;
        return os;
    }

    /*
    References:
    http://www.brucelindbloom.com/
    http://www.ece.rochester.edu/~gsharma/ciede2000/ciede2000noteCRNA.pdf
    */
    const double PI_ = 3.14159265358979323846;
    const int K_L = 1;
    const int K_C = 1;
    const int K_H = 1;
    const double X = std::pow(25, 7);

    enum class color_diff_type
    {
        DeltaAB,
        DeltaC2000,
        DeltaE2000
    };

    template <color_diff_type DIFFERENCE>
    struct color_diff
    {
        template <typename T, wp WP>
        static double calculate(const lab<T, WP>& pt1, const lab<T, WP>& pt2)
        {
            switch (DIFFERENCE)
            {
            case Teisko::color_diff_type::DeltaAB:
                return std::hypot(pt2.a - pt1.a, pt2.b - pt1.b);
            case Teisko::color_diff_type::DeltaC2000:
            {
                // Luma component omitted from the final calculation
                return delta2000(pt1, pt2, true);
            }
            case Teisko::color_diff_type::DeltaE2000:
            {
                return delta2000(pt1, pt2);
            }
            default:
                throw std::invalid_argument("Unsupported color difference");
            }
        }

    private:
        template <typename T, wp WP>
        static double inline compute_h1_h2(const lab<T, WP>& pt, double a1)
        {
            auto arctan_h1 = std::atan2(pt.b, a1) * 180 / PI_;
            return arctan_h1 >= 0 ? arctan_h1 : arctan_h1 + 360;
        }

        static double inline compute_H(double h_1, double h_2, double C1, double C2)
        {
            if (C1 * C2 == 0)
                return h_1 + h_2;
            if (std::abs(h_1 - h_2) <= 180)
                return 0.5 * (h_1 + h_2);
            return h_1 + h_2 < 360 ? 0.5 * (h_1 + h_2) + 180 : 0.5 * (h_1 + h_2) - 180;
        }

        static double inline compute_delta_h(double h_1, double h_2, double C1, double C2)
        {
            if (C1 * C2 == 0)
                return 0;
            if (std::abs(h_2 - h_1) <= 180)
                return h_2 - h_1;
            return h_2 <= h_1 ? h_2 - h_1 + 360 : h_2 - h_1 - 360;
        }

        template <typename T, wp WP>
        static double delta2000(const lab<T, WP>& pt1, const lab<T, WP>& pt2, bool omit_luma = false)
        {
            // Computing L', C1, C2, C
            auto L = (pt1.l + pt2.l) / 2;
            auto C1 = std::hypot(pt1.a, pt1.b);
            auto C2 = std::hypot(pt2.a, pt2.b);
            auto C = (C1 + C2) / 2.0;

            // Computing G, a'1, a'2
            auto G = 0.5 * (1.0 - std::sqrt(std::pow(C, 7) / (std::pow(C, 7) + X)));
            auto a1 = pt1.a * (1 + G);
            auto a2 = pt2.a * (1 + G);

            // Computing h'1, h'2 and H'
            auto h_1 = compute_h1_h2(pt1, a1);
            auto h_2 = compute_h1_h2(pt2, a2);

            // Computing C'1, C'2, H, C'
            auto C_1 = std::hypot(a1, pt1.b);
            auto C_2 = std::hypot(a2, pt2.b);
            auto H = compute_H(h_1, h_2, C_1, C_2);
            auto C_ = (C_1 + C_2) * 0.5;

            // Computing T
            auto T_ = 1 - 0.17 * std::cos((H - 30) * PI_ / 180) + 0.24 * std::cos((2 * H) * PI_ / 180) +
                0.32 * std::cos((3 * H + 6) * PI_ / 180) - 0.20 * std::cos((4 * H - 63) *  PI_ / 180);

            // Computing delta_h'
            auto delta_h = compute_delta_h(h_1, h_2, C_1, C_2);

            // Computing delta_L', delta_C', delta_H
            auto delta_L = omit_luma ? 0.0 : pt2.l - pt1.l;
            auto delta_C = C_2 - C_1;
            auto delta_H = 2 * std::sqrt(C_1 * C_2) * std::sin(0.5 * delta_h * PI_ / 180);

            // Computing SL, SC, SH
            auto S_L = 1 + (0.015 * pow2(L - 50)) / (std::sqrt(20 + pow2(L - 50)));
            auto S_C = 1 + 0.045 * C_;
            auto S_H = 1 + 0.015 * C_ * T_;

            // Computing delta_theta, R_C, R_T
            auto delta_theta = 30 * std::exp(-1 * (pow2(((H - 275) / 25))));
            auto R_C = 2 * std::sqrt(std::pow(C_, 7) / (std::pow(C_, 7) + X));
            auto R_T = -1 * R_C * sin(2 * delta_theta * PI_ / 180);

            // Computing the values needed for final error computation
            return std::sqrt(pow2((delta_L / (K_L * S_L))) +
                pow2((delta_C / (K_C * S_C))) +
                pow2((delta_H / (K_H * S_H))) +
                R_T * (delta_C / (K_C * S_C)) * (delta_H / (K_H * S_H)));
        }
    };
};
