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

#include "Teisko/Algorithm/PointXY.hpp"
// In alphabetical order
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <stdexcept>
#include <vector>

/// \brief This namespace defines the public interfaces of
/// the \ref libcalc_module module
namespace Teisko
{
    /// Constraints for nelder mead simplex calculation
    struct nelder_mead_properties
    {
        int max_iterations;     // max suggested number of iterations (200*N) if missing
        double tol_x;           // terminating condition for all parameters are within this distance
        double tol_fun;         // terminating condition for parameters evaluated within this value
        int max_func_evals;     // max number of function evaluations

        nelder_mead_properties(int iters = -1, double tx = 1.0e-6, double tf = 1.0e-6, int funcs = -1)
            : max_iterations(iters), tol_x(tx), tol_fun(tf), max_func_evals(funcs) { }
    };

    // Intermediate and return value of Nelder Mead Simplex
    // first: evaluated value; second: parameters for the evaluation
    struct nelder_mead_param
    {
        double first;
        std::vector<double> second;

        nelder_mead_param(size_t n = 0) : first(0.0), second(n) { }
        nelder_mead_param(double value, std::vector<double> &vec) : first(value), second(vec) { }
        bool operator <(const nelder_mead_param &rhs) const { return first < rhs.first; }
        bool operator >=(const nelder_mead_param &rhs) const { return first >= rhs.first; }
    };

    struct nelder_mead_param_generator
    {
        // Remove assignment operator to allow compilation
        // One should never need to make copies anyway
        nelder_mead_param_generator operator =(const nelder_mead_param_generator &other) = delete;

        std::function<double(double *)> &&func;     // Function to evaluate the cost
        size_t n;            // Number of items per vector
        size_t vs;           // Index to smallest evaluated parameter set (best)
        size_t vh;           // Index to second highest evaluated parameter set
        size_t vg;           // Index to largest evaluated parameter set (worst)

        nelder_mead_properties props;

        std::vector<double> avg;                    // Average of N best parameters
        std::vector<nelder_mead_param> full_set;    // Set of parameters and the evaluated costs

        // Replaces the worst parameter set with the argument
        // - equivalent to  full_set[vg] = a -- but without the constructor
        void assign_highest(nelder_mead_param &a)
        {
            auto &highest = full_set[vg];
            highest.first = a.first;

            double *s = a.second.data();
            double *d = highest.second.data();

            for (size_t i = 0; i < n; i++)
                d[i] = s[i];
        }

        // Locates the smallest, highest and second highest value from the full_set
        void sort()
        {
            vs = 0;
            vg = 0;
            vh = 1;

            for (size_t idx = 1; idx <= n; idx++)
            {
                if (full_set[idx] < full_set[vs])
                    vs = idx;
                // Maximum and second largest
                if (full_set[idx] >= full_set[vg])
                {
                    vh = vg;
                    vg = idx;
                }
                else if (full_set[idx] >= full_set[vh])
                    vh = idx;
            }
        }

        // calculates the average value of parameters except the set evaluated highest (worst)
        // Theoretically we could improve this by evaluating it incrementally
        // This would however suffer from random drift, as (a + b - a - b) is not necessarily 0
        void update_avg()
        {
            avg.clear();
            avg.resize(n);
            double *a = avg.data();
            for (size_t i = 0; i <= n; i++)
            {
                if (i == vg) continue;
                double *s = full_set[i].second.data();
                for (size_t j = 0; j < n; j++)
                    a[j] += s[j];
            }
            for (size_t i = 0; i < n; i++)
                a[i] *= (1.0 / static_cast<double>(n));
        }

        // Initializes the full_set with variation in the Nth element of parameter space
        // After initialization there are N+1 vectors with one original vector and N altered vectors
        nelder_mead_param_generator(
            std::vector<double> &init,
            std::function<double(double *)> &&f,
            nelder_mead_properties &options)
            : func(std::move(f)), n(init.size()), vs(0), vh(0), vg(0), props(options), avg(n)
            // copies the init vector and its evaluated value to n+1 positions
            , full_set(std::vector<nelder_mead_param>(n + 1, nelder_mead_param(func(init.data()), init)))
        {
            const double zero_term_delta = 0.00025; // absolute increase in term if it's initially zero
            const double delta = 0.05; // initial five percent difference for each parameter
            // tweak the first N vectors in `full_set` -- each vector tweaked at different element
            for (size_t i = 0; i < n; i++)
            {
                auto &p = full_set[i];
                double *d = p.second.data();
                d[i] = std::abs(d[i]) == 0 ? zero_term_delta : (1 + delta) * d[i];
                p.first = func(d);
            }

            if (props.max_iterations < 0)
                props.max_iterations = 200 * static_cast<int>(n);
            if (props.max_func_evals < 0)
                props.max_func_evals = 400 * static_cast<int>(n);

            props.max_func_evals -= static_cast<int>(n + 1);
        }

        // Linear combination of `avg` and the worst so far
        void mix(nelder_mead_param &result, double factor)
        {
            double *a = avg.data();
            double *s = full_set[vg].second.data();
            double *d = result.second.data();
            for (size_t i = 0; i < n; i++)
                d[i] = a[i] * (1.0 - factor) + s[i] * factor;
            result.first = func(d);
            --props.max_func_evals;
        }

        void shrink()
        {
            auto best = full_set[vs].second.data();
            for (size_t i = 0; i <= n; i++)
            {
                if (i == vs)
                    continue;
                auto dst = full_set[i].second.data();
                const double sigma = 0.5;
                for (size_t j = 0; j < n; j++)
                {
                    dst[j] = (dst[j] - best[j]) * sigma + best[j];
                }
            }
        }

        // Called once and only once per iteration
        bool should_terminate()
        {
            if (props.max_iterations-- <= 0 || props.max_func_evals <= 0)
                return true;

            sort();

            // check difference between largest and smallest keys (evaluation of objective function)
            if (full_set[vg].first - full_set[vs].first > props.tol_fun)
                return false;

            double *s = full_set[vs].second.data();
            // check maximum difference between function 'x' values to currently best vector
            for (size_t i = 0; i <= n; i++)
            {
                if (i == vs)
                    continue;
                double *d = full_set[i].second.data();
                for (size_t j = 0; j < n; j++)
                {
                    if (std::abs(s[j] - d[j]) > props.tol_x)
                        return false;
                }
            }
            return true;
        }
    };

    // Function minimizer using Nelder Mead method
    // initial -- starting point of evaluation in N dimensions
    // func -- functor, lambda or function taking the N parameters in a pointer
    // props -- options to manage the terminating conditions and iteration count
    // returns tuple of minimum value as evaluated by the `func`, vector of parameters
    inline nelder_mead_param nelder_mead_simplex(std::vector<double> initial,
        std::function<double(double *p)> &&func, nelder_mead_properties props = {})
    {
        size_t n = initial.size();
        if (n == 0)
            return{};

        nelder_mead_param_generator params(initial, std::move(func), props);

        // pre-allocate these parameters for efficiency -- and re-use in the loop
        auto ref = nelder_mead_param(n);
        auto exp = nelder_mead_param(n);
        auto con = nelder_mead_param(n);

        while (!params.should_terminate())
        {
            params.update_avg();

            auto &worst = params.full_set[params.vg];
            auto &best = params.full_set[params.vs];

            const double rho = 1.0;
            const double psi = 0.5;
            const double chi = 2.0;

            // reflection
            params.mix(ref, -rho);

            if (ref < best)
            {
                // expansion
                params.mix(exp, -rho*chi);
                params.assign_highest(exp < ref ? exp : ref);
            }
            else
            {
                if (ref < params.full_set[params.vh])
                    params.assign_highest(ref);
                else
                {
                    if (ref < worst)
                    {
                        // outside contraction
                        params.mix(con, -rho*psi);
                        if (con < ref)
                            params.assign_highest(con);
                        else
                            params.shrink();
                    }
                    else
                    {
                        // inside contraction
                        params.mix(con, psi);
                        if (con < worst)
                            params.assign_highest(con);
                        else
                            params.shrink();
                    }
                }
            }
        }
        params.sort();
        return params.full_set[params.vs];
    }
}
