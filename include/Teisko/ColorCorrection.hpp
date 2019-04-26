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

#include "Teisko/Algorithm/Bit.hpp"
#include "Teisko/Algorithm/Functors.hpp"
#include "Teisko/Algorithm/NelderMead.hpp"
#include "Teisko/Color.hpp"
#include <Eigen/Dense>
#include <array>
#include <exception>
#include <numeric>
#include <vector>

namespace Teisko
{
    template <int N>
    int operator-(const cyclical_index<N>& a, const cyclical_index<N>& b)
    {
        return a.idx < b.idx ? N - (b.idx - a.idx) : b.idx - a.idx;
    }

    struct color_correction_matrix : Eigen::Matrix3d
    {
        color_correction_matrix()
        {
            setIdentity();
        }

        color_correction_matrix(const Eigen::Matrix3d& v) : Eigen::Matrix3d(v) {}

        // Apply color correction matrix
        template <rgb_cs CS = rgb_cs::sRGB, wp WP = wp::D65>
        rgb<double, CS, WP> apply(const rgb<double, rgb_cs::Sensor, wp::None>& pt) const
        {
            return rgb<double, CS, WP> (
                (*this)(0, 0) * pt.r + (*this)(0, 1) * pt.g + (*this)(0, 2) * pt.b,
                (*this)(1, 0) * pt.r + (*this)(1, 1) * pt.g + (*this)(1, 2) * pt.b,
                (*this)(2, 0) * pt.r + (*this)(2, 1) * pt.g + (*this)(2, 2) * pt.b);
        }
    };

    // Calculate color correction matrix from input and output matrices
    // Either calculate ccm by
    // 1. matrix inverse if matrix is invertible
    // 2. least square method using svd decomposition (high accuracy but slow)
    inline color_correction_matrix calculate_ccm(const Eigen::Matrix3d& in, const Eigen::Matrix3d& out)
    {
        using namespace Eigen;

        const double THRESHOLD = 1e-8;
        if (std::abs(in.determinant()) > THRESHOLD)
        {
            // CCM = M(out) / M(in) = M(out)*inverse(M(in))
            //
            // M(in)  = [RGB(1, in) ^ T, RGB(2, in) ^ T, RGB(grey) ^ T]
            // M(out) = [RGB(1, out) ^ T, RGB(2, out) ^ T, RGB(grey) ^ T]
            return color_correction_matrix(out * in.inverse());
        }
        JacobiSVD<Matrix3d> svd(in.transpose(), ComputeFullU | ComputeFullV);
        return color_correction_matrix(svd.solve(out.transpose()).transpose());
    }

    template <rgb_cs CS, wp WP = wp::D65>
    struct patch
    {
        rgb<double, rgb_cs::Sensor, wp::None> input;
        rgb<double, CS, WP> target;
        double weight;
        bool achromatic;
    };

    template <rgb_cs CS, wp WP = wp::D65>
    struct ccm_input_params
    {
        std::array<patch<CS, WP>, 24> color_checker_classic;
        //std::vector<patch<CS, WP>> munsell; // TODO: Implementation for Munsell

        void white_balance()
        {
            average_f<double> r_per_g;
            average_f<double> b_per_g;

            // Loop over achromatic patches (bottom row), dimmest and brightest omitted
            for (size_t i = 19; i < 23; ++i)
            {
                r_per_g += color_checker_classic[i].input.r / color_checker_classic[i].input.g;
                b_per_g += color_checker_classic[i].input.b / color_checker_classic[i].input.g;
            }

            auto max_gain = std::max({ (double)r_per_g, (double)b_per_g, 1.0 });
            auto gain_r = max_gain / r_per_g;
            auto gain_g = max_gain;
            auto gain_b = max_gain / b_per_g;

            for (auto& patch : color_checker_classic)
            {
                patch.input.r *= gain_r;
                patch.input.g *= gain_g;
                patch.input.b *= gain_b;
            }
        }

        void normalize()
        {
            average_f<double> ratio;
            for (size_t i = 19; i < 23; ++i)
            {
                ratio += color_checker_classic[i].target.mean() / color_checker_classic[i].input.mean();
            }

            for (auto& patch : color_checker_classic)
            {
                patch.input *= ratio;
            }
        }
    };

    template <rgb_cs CS, wp WP = wp::D65>
    struct ccm_optimization_params
    {
        ccm_input_params<CS, WP> input_params;

        ccm_optimization_params() = default;
        ccm_optimization_params(const ccm_input_params<CS, WP>& input) : input_params(input) { }

        void prepare()
        {
            input_params.white_balance();
            input_params.normalize();
        }
    };

    template <rgb_cs CS, wp WP = wp::D65>
    struct ccm_func_optimizer
    {
        ccm_optimization_params<CS, WP> params;

        ccm_func_optimizer(ccm_optimization_params<CS, WP>& p) : params(p) { }

        double operator()(double *p)
        {
            auto sum = 0.0;

            // Create CCM from optimization set. Middle column is constructed from
            // border columns (each row sum up to 1)
            auto ccm = color_correction_matrix();
            ccm <<
                p[0], 1 - p[0] - p[1], p[1],
                p[2], 1 - p[2] - p[3], p[3],
                p[4], 1 - p[4] - p[5], p[5];

            for (const auto& patch : params.input_params.color_checker_classic)
            {
                if (patch.weight > 0)
                {
                    lab<double, WP> pt1(ccm.apply<CS, WP>(patch.input));
                    lab<double, WP> pt2(patch.target); // TODO: Pre calculate

                    sum += color_diff<color_diff_type::DeltaE2000>::calculate(pt1, pt2) * patch.weight;
                }
            }
            return sum;
        }
    };

    template <rgb_cs CS, wp WP = wp::D65>
    class ccm_optimization
    {
    public:
        ccm_optimization(ccm_input_params<CS, WP>& input) : options(input) { }

        color_correction_matrix characterize()
        {
            // Do white balance and adjust input intensity
            options.prepare();

            // Get initial point for nelder mead optimization by using least square CCM
            auto initial_ccm = ls_ccm().data();

            // RinR, BinR, RinG, BinG, RinB, BinB indices
            std::vector<double> initial_pt =
            {
                initial_ccm[0], initial_ccm[2], initial_ccm[3], initial_ccm[5], initial_ccm[6], initial_ccm[8]
            };

            // 1. coarse optimization
            auto solver = ccm_func_optimizer<CS, WP>(options);
            auto result = nelder_mead_simplex(initial_pt, solver, {200, 1e-4, 1e-4});

            // 2. fine optimization
            result = nelder_mead_simplex(result.second, solver, {4000, 1e-12, 1e-4});

            // Create CCM from optimization set. Middle column is constructed from
            // border columns (each row should sum up to 1)
            color_correction_matrix ccm;
            ccm <<
                result.second[0], 1 - result.second[0] - result.second[1], result.second[1],
                result.second[2], 1 - result.second[2] - result.second[3], result.second[3],
                result.second[4], 1 - result.second[4] - result.second[5], result.second[5];
            return ccm;
        }

    private:
        ccm_optimization_params<CS, WP> options;

        // Solve a linear system Ax=b for initial CCM optimization
        // Where A is CCC input matrix and B is CCC target matrix
        // This method does not preserve greys nor does it minimize error in de2000.
        color_correction_matrix ls_ccm()
        {
            Eigen::Matrix<double, 24, 3> A;
            Eigen::Matrix<double, 24, 3> b;
            Eigen::Matrix3d x;

            // Initialize matrices
            size_t n = options.input_params.color_checker_classic.size();
            for (size_t i = 0; i < n; ++i)
            {
                const patch<CS, WP>& patch = options.input_params.color_checker_classic[i];
                A.row(i) << patch.input.r, patch.input.g, patch.input.b;
                b.row(i) << patch.target.r, patch.target.g, patch.target.b;
            }

            // Use HouseholderQR decomposition to solve the linear equation.
            // Prefer speed over accuracy for the decomposition method.
            return color_correction_matrix(A.householderQr().solve(b));
        }
    };
}