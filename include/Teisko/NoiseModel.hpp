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

#include "Eigen/Dense"

// Visual Studio 2013 and 2015 compilers complains about C4714 warning when
// compiling the Eigen/Dense header. However, Visual Studio 2017 compiles
// without warnings.
//
// Eigen provides a utility header file that can be used to disable selected
// warnings and the warning C4714 is one of those identified warnings.
//
// We include the utility header only on Visual Studio 2013 and 2015.
//      VS2013: _MSC_VER 1800
//      VS2015: _MSC_VER 1900
//
#if defined _MSC_VER && (_MSC_VER == 1800 || _MSC_VER == 1900)
#include "Eigen/src/Core/util/DisableStupidWarnings.h"
#endif

#include <algorithm>
#include <cstdint>
#include <vector>

namespace Teisko
{
    // 5 parameter noise model: y = c1*x^2 + AG*c2*x + AG^2*c3 + AG*c4 + c5
    // x is usually mean of image intensities in an image patch
    // y is usually either variance or standard deviation of an image patch
    struct noise_model
    {
        double c1;
        double c2;
        double c3;
        double c4;
        double c5;
    };

    // Helper functions -----------------------------------------------------------
    // Checks if any of the booleans is true
    inline bool any(uint8_t* booleans, int n)
    {
        for (int idx = 0; idx < n; ++idx)
        {
            if (booleans[idx]) return true;
        }
        return false;
    }

    // Calculates n = 1 vector norm for given C matrix
    inline double calculate_C_norm_n1(std::vector<std::vector<double>>& C)
    {
        double norm = 0;
        for (size_t idx = 0; idx < C.size(); ++idx)
        {
            norm = norm + C[idx][0];
        }
        return norm;
    }

    // Estimates a tolerance value from given C matrix
    inline double calculate_tolerance(std::vector<std::vector<double>>& C)
    {
        double epsilon = std::numeric_limits<double>::epsilon();
        return 10 * epsilon * calculate_C_norm_n1(C) * C.size();
    }

    // Checks if any of the values are bigger than tolerance
    // TODO: Why tolerance is not abs from value?
    inline bool check_tolerances(double* w, uint8_t* Z, int n, double tolerance)
    {
        for (int idx = 0; idx < n; ++idx)
        {
            if (Z[idx] && (w[idx] > tolerance))
            {
                return true;
            }
        }
        return false;
    }

    // Check negativity in a positive set
    inline bool check_negativity(double* z, uint8_t* P, int n)
    {
        for (int idx = 0; idx < n; ++idx)
        {
            if (P[idx] && (z[idx] <= 0)) return true;
        }
        return false;
    }

    // Calculates residual = d - C*x
    inline void calculate_residual(std::vector<std::vector<double>>& C, std::vector<double>& d, double* x, std::vector<double>& resid_out, int n)
    {
        for (size_t id_samp = 0; id_samp < resid_out.size(); ++id_samp)
        {
            double mat_mult = 0;
            for (int id_n = 0; id_n < n; ++id_n)
            {
                mat_mult = mat_mult + C[id_samp][id_n] * x[id_n];
            }
            resid_out[id_samp] = d[id_samp] - mat_mult;
        }
    }

    // Calculates lambda = C * residual
    inline void calculate_lambda(std::vector<std::vector<double>>& C, std::vector<double>& residual, double* lambda_out, int n)
    {
        for (int id_n = 0; id_n < n; ++id_n)
        {
            lambda_out[id_n] = 0;
            for (size_t id_samp = 0; id_samp < residual.size(); ++id_samp)
            {
                lambda_out[id_n] = lambda_out[id_n] + C[id_samp][id_n] * residual[id_samp];
            }
        }
    }

    inline void least_squares(int samples, int number_of_coefficients, std::vector<double> &a, std::vector<double> &b)
    {
        // This MKL library call has been replaced with similar call in Eigen
        //
        // LAPACKE_dgels(LAPACK_COL_MAJOR, 'N',
        //      samples, number_of_coefficients, 1,
        //      a.data(), samples,
        //      b.data(), samples);

        Eigen::Map<Eigen::MatrixXd> A(a.data(), samples, number_of_coefficients);
        Eigen::Map<Eigen::VectorXd> B(b.data(), samples);

        Eigen::VectorXd result = A.colPivHouseholderQr().solve(B);

        // The result is reported back to the caller using the vector b.
        // Please note that whatever there was after k (k+1, k+2, ...)
        // still remains.
        for (Eigen::Index k = 0; k < result.rows(); ++k)
        {
            b[k] = result[k];
        }
    }

    // Calculates linear least squares fit for positive set using LAPACKE's dgels routine
    inline void least_squares_fit_for_positive_set(std::vector<std::vector<double>>& C, std::vector<double>& d, uint8_t* P, double* z_out, int n)
    {
        // Number of samples
        auto samples = d.size() > INT_MAX ? INT_MAX : (int)d.size();

        // Calculate max number of coefficients, size of positive set
        int no_coefficients = 0;
        for (int idx = 0; idx < n; ++idx)
        {
            if (P[idx]) no_coefficients++;
        }

        auto A = std::vector<double>(samples * no_coefficients);
        auto B = std::vector<double>(d.begin(), d.begin() + samples);

        // Fill A data from C
        int a_idx = 0;
        for (int id_var = 0; id_var < n; ++id_var)
        {
            if (P[id_var])
            {
                for (size_t id_samp = 0; id_samp < d.size(); ++id_samp)
                {
                    A[a_idx++] = C[id_samp][id_var];
                }
            }
        }

        // Find least squares solution
        // LAPACKE_dgels(LAPACK_COL_MAJOR, 'N', samples, no_coefficients, 1, A.data(), samples, B.data(), samples);
        least_squares(samples, no_coefficients, A, B);

        // Store solution to output
        int var_idx = 0;
        for (int id_var = 0; id_var < n; ++id_var)
        {
            if (P[id_var])
            {
                z_out[id_var] = B[var_idx++];
            }
        }
    }

    // Calculates non-negative linear least squares for given matrix C
    inline void non_neg_least_squares(std::vector<std::vector<double>>& C, std::vector<double>& d, double* x, int n)
    {
        // Calculates noise model from statistics with negativity constrain. Implementation referenced from MATLAB's lsqnonneg.m

        // Initialize set of parameters
        double tolerance = calculate_tolerance(C);
        auto w = std::vector<double>(n, 0.0);
        auto wz = std::vector<double>(n, 0.0);
        auto P = std::vector<uint8_t>(n, 0);
        auto Z = std::vector<uint8_t>(n, 0);
        std::fill(Z.begin(), Z.end(), static_cast<uint8_t>(1));
        auto resid = std::vector<double>(d.size());

        // Calculate initial residual and lambda
        calculate_residual(C, d, &x[0], resid, n);
        calculate_lambda(C, resid, &w[0], n);

        // Set up iteration criterion
        int outeriter = 0;
        int iter = 0;
        int itmax = 3 * n;
        bool continueflag = true;

        // Outer loop to put variables into set to hold positive coefficients
        while (any(&Z[0], n) && check_tolerances(&w[0], &Z[0], n, tolerance) && continueflag)
        {
            outeriter = outeriter + 1;
            // Reset intermediate solution z
            auto z = std::vector<double>(n, 0.0);
            // Create wz, a Lagrange multiplier vector of variables in the zero set.
            // wz must have the same size as w to preserve the correct indices, so
            // set multipliers to - Inf for variables outside of the zero set.
            for (int idx = 0; idx < n; ++idx)
            {
                if (P[idx])
                {
                    wz[idx] = std::numeric_limits<double>::lowest();
                }
            }

            // Fill with active set values
            for (int idx = 0; idx < n; ++idx)
            {
                if (Z[idx])
                {
                    wz[idx] = w[idx];
                }
            }

            // Find variable with largest Lagrange multiplier
            auto t = std::distance(wz.data(), std::max_element(wz.data(), wz.data() + n));

            // Move variable t from zero set to positive set
            P[t] = true;
            Z[t] = false;

            // Compute intermediate solution using only variables in positive set
            least_squares_fit_for_positive_set(C, d, &P[0], &z[0], n);

            // Inner loop to remove elements from the positive set which no longer belong
            while (check_negativity(&z[0], &P[0], n))
            {
                iter = iter + 1;

                // If maximum number of iteration is encountered, return current solution
                if (iter > itmax)
                {
                    continueflag = false;
                    break;
                }

                // Find indices where intermediate solution z is approximately negative
                auto Q = std::vector<uint8_t>(n, 0);
                for (int idx = 0; idx < n; ++idx)
                {
                    if (z[idx] <= 0 && P[idx])
                    {
                        Q[idx] = true;
                    }
                }

                // Choose new x subject to keeping new x nonnegative
                auto alpha = std::vector<double>(n, std::numeric_limits<double>::max());
                for (size_t idx = 0; idx < alpha.size(); ++idx)
                {
                    if (Q[idx])
                    {
                        alpha[idx] = x[idx] / (x[idx] - z[idx]);
                    }
                }
                auto min_alpha = std::min_element(alpha.data(), alpha.data() + alpha.size());
                for (int idx = 0; idx < n; ++idx)
                {
                    x[idx] = x[idx] + *min_alpha * (z[idx] - x[idx]);
                }

                // Reset Z and P given intermediate values of x
                for (int idx = 0; idx < n; ++idx)
                {
                    Z[idx] = ((std::abs(x[idx]) < tolerance) & P[idx]) | Z[idx];
                    P[idx] = !Z[idx];
                }

                // Reset z and resolve
                std::fill(z.begin(), z.end(), 0.0);
                least_squares_fit_for_positive_set(C, d, &P[0], &z[0], n);
            }

            // Replace solution
            for (int idx = 0; idx < n; ++idx) x[idx] = z[idx];
            calculate_residual(C, d, &x[0], resid, n);
            calculate_lambda(C, resid, &w[0], n);
        }
    }

    namespace NoiseModel
    {
        const int32_t SUCCESS = 0;
    }

    inline int32_t calculate_y_noise_model(std::vector<double> const& means, std::vector<double> const& stds, noise_model& noise_model)
    {
        // Create proper data structures for LAPACKE_dgels
        auto samples_size_t = std::min(means.size(), stds.size());
        auto samples = samples_size_t > INT_MAX ? INT_MAX : (int)samples_size_t; // INT_MAX is maximum sample count that LAPACKE_dgels accepts
        const int number_of_coefficients = 3;
        auto A = std::vector<double>(samples * number_of_coefficients);
        auto B = std::vector<double>(stds.begin(), stds.begin() + samples);

        // A = [constant coeff_y coeff_y^2]
        for (int k = 0; k < samples; k++)
        {
            A[k] = 1;
            A[samples + k] = means[k];
            A[samples * 2 + k] = means[k] * means[k];
        }

        // Solve problem and fill noise model
        // LAPACKE_dgels(LAPACK_COL_MAJOR, 'N', samples, number_of_coefficients, 1, A.data(), samples, B.data(), samples);
        least_squares(samples, number_of_coefficients, A, B);
        noise_model.c3 = B[0];
        noise_model.c2 = B[1];
        noise_model.c1 = B[2];

        return NoiseModel::SUCCESS;
    }

    inline int32_t calculate_noise_model_nonneg(std::vector<double> const& means, std::vector<double> const& stds, noise_model& noise_model)
    {
        // Create 2 dimensional vector to store C array
        const int n = 3;
        auto C = std::vector<std::vector<double>>(means.size(), std::vector<double>(n));
        auto d = std::vector<double>(stds.size());
        double x[n] = {};

        // Fill data structures
        for (size_t idx = 0; idx < stds.size(); ++idx)
        {
            C[idx][0] = 1;
            C[idx][1] = means[idx];
            C[idx][2] = means[idx] * means[idx];
            d[idx] = stds[idx] * stds[idx]; // Use variance space for noise model (better would be to store only variance in the first place!!!)
        }

        // Find minimum variance value and remove it from the data (noise floor), add it back after solution is found
        auto min_var = *(std::min_element(d.data(), d.data() + d.size()));

        for (size_t idx = 0; idx < d.size(); ++idx)
        {
            d[idx] = d[idx] - min_var;
        }

        // C and d needs to be refactored elsewhere, also x since it is input for the result
        non_neg_least_squares(C, d, &x[0], n);

        // Fill solution
        noise_model.c3 = x[0] + min_var; // Add minimum variance back to the solution
        noise_model.c2 = x[1];
        noise_model.c1 = x[2];

        return NoiseModel::SUCCESS;
    }

    inline int32_t fit_surface(std::vector<double> const& values, std::vector<double> const& xcoords, std::vector<double> const& ycoords, std::vector<double>& surface)
    {
        // Fits a 2D surface to given values and coordinates
        // Model: z = p1 + p2*x + p3*y + p4*x^2 + p5*y^2 + p6*x*y;
        auto samples_size_t = values.size();
        auto samples = samples_size_t > INT_MAX ? INT_MAX : (int)samples_size_t; // INT_MAX is maximum sample count that LAPACKE_dgels accepts
        const int number_of_coefficients = 6;

        auto A = std::vector<double>(samples * number_of_coefficients);
        auto B = std::vector<double>(values.begin(), values.begin() + samples);

        // A = [constant x y x^2 y^2 x*y]
        for (int k = 0; k < samples; k++)
        {
            A[k] = 1;
            A[samples + k] = xcoords[k];
            A[samples * 2 + k] = ycoords[k];
            A[samples * 3 + k] = xcoords[k] * xcoords[k];
            A[samples * 4 + k] = ycoords[k] * ycoords[k];
            A[samples * 5 + k] = xcoords[k] * ycoords[k];
        }

        // Solve surface polynomial
        // LAPACKE_dgels(LAPACK_COL_MAJOR, 'N', samples, number_of_coefficients, 1, A.data(), samples, B.data(), samples);
        least_squares(samples, number_of_coefficients, A, B);

        // Fill surface
        for (int idx = 0; idx < samples; ++idx)
        {
            surface[idx] = B[0] + B[1] * xcoords[idx] + B[2] * ycoords[idx] + B[3] * xcoords[idx] * xcoords[idx] + B[4] * ycoords[idx] * ycoords[idx] + B[5] * xcoords[idx] * ycoords[idx];
        }

        return NoiseModel::SUCCESS;
    }

    inline int32_t calculate_unified_non_neg_noise_model(std::vector<double> const& means, std::vector<double> const& vars, std::vector<double> const& weights, std::vector<double> const& gains, noise_model& noise_model)
    {
        // Create 2 dimensional vector to store C array
        const int n = 5; // y = c1*x ^ 2 + AG*c2*x + AG ^ 2 * c3 + AG*c4 + c5
        auto C = std::vector<std::vector<double>>(means.size(), std::vector<double>(n));
        auto d = std::vector<double>(vars.size());
        double x[n] = {};

        // Fill data structures
        for (size_t idx = 0; idx < vars.size(); ++idx)
        {
            C[idx][0] = means[idx] * means[idx] * weights[idx];
            C[idx][1] = means[idx] * weights[idx] * gains[idx];
            C[idx][2] = weights[idx] * gains[idx] * gains[idx];
            C[idx][3] = weights[idx] * gains[idx];
            C[idx][4] = weights[idx];
            d[idx] = vars[idx] * weights[idx];
        }

        // Use non-negative least squares solution
        non_neg_least_squares(C, d, &x[0], n);

        // Fill solution
        noise_model.c1 = x[0];
        noise_model.c2 = x[1];
        noise_model.c3 = x[2];
        noise_model.c4 = x[3];
        noise_model.c5 = x[4];

        return NoiseModel::SUCCESS;
    }

    inline int32_t calculate_2param_non_neg_noise_model(std::vector<double> const& means, std::vector<double> const& vars, std::vector<double> const& gains, noise_model& noise_model)
    {
        // Create 2 dimensional vector to store C array
        const int n = 2; // y = AG * c2 * x + AG ^ 2 * c3
        auto C = std::vector<std::vector<double>>(means.size(), std::vector<double>(n));
        auto d = std::vector<double>(vars.size());
        double x[n] = {};

        // Fill data structures
        for (size_t idx = 0; idx < vars.size(); ++idx)
        {
            C[idx][0] = means[idx] * gains[idx];
            C[idx][1] = gains[idx] * gains[idx];
            d[idx] = vars[idx];
        }

        // Use non-negative least squares solution
        non_neg_least_squares(C, d, &x[0], n);

        // Fill solution
        noise_model.c2 = x[0];
        noise_model.c3 = x[1];

        return NoiseModel::SUCCESS;
    }
}
