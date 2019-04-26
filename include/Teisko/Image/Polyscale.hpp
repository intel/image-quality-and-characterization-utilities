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
#include "Teisko/Algorithm/ReduceTo.hpp"
#include "Teisko/Image/API.hpp"
#include <cstdint>
#include <vector>

namespace Teisko
{
    template <typename T>
    struct grid_2d
    {
        size_t _cols;
        size_t _rows;
        std::vector<T> _grid;

        explicit grid_2d(image<T> &src_image)
            : _cols(src_image._width)
            , _rows(src_image._height)
            , _grid(src_image.to_vector()) { }

        grid_2d() : _cols(0), _rows(0) { };

        grid_2d(size_t c, size_t r)
            : _cols(c), _rows(r), _grid(c *r) { }

        grid_2d(size_t c, size_t r, const std::vector<T> &i)
            : _cols(c), _rows(r), _grid(i) { }

        grid_2d &resize(size_t c, size_t r)
        {
            _cols = c;
            _rows = r;
            _grid.resize(c*r);
            return *this;
        }

        grid_2d transpose()
        {
            grid_2d<T> trans(_rows, _cols);
            for (size_t j = 0; j < _rows; j++)
            {
                for (size_t i = 0; i < _cols; i++)
                {
                    trans._grid[i * trans._cols + j] = _grid[j * _cols + i];
                }
            }
            return trans;
        }

        image<T> to_view()
        {
            return image<T>((uint32_t)_rows, (uint32_t)_cols, _grid.data());
        }

        /// Multiply this * B
        template <typename V>
        grid_2d mul_mm(grid_2d<V> &B) const
        {
            auto &A = *this;
            grid_2d output(B._cols, A._rows);
            auto Bt = B.transpose();

            auto t = 64; // Cache line
            for (size_t j = 0; j < B._cols; j += t)
            {
                for (size_t k = 0; k < A._cols; k += t)
                {
                    for (size_t i = 0; i < A._rows; i += t)
                    {
                        for (size_t ii = i; ii < std::min((size_t)i + t, A._rows); ii++)
                        {
                            for (size_t jj = j; jj < std::min((size_t)j + t, B._cols); jj++)
                            {
                                auto tmp = output._grid[ii * output._cols + jj];
                                for (size_t kk = k; kk < std::min((size_t)k + t, A._cols); kk++)
                                {
                                    auto product = A._grid[ii *A._cols + kk] * Bt._grid[jj * Bt._cols + kk];
                                    if (!std::isnan(product))
                                        tmp += product;
                                }
                                output._grid[ii *output._cols + jj] = tmp;
                            }
                        }
                    }
                }
            }
            return output;
        }

        // Matrix multiplication by its transpose, A * At
        grid_2d mul_mm() const
        {
            auto &A = *this;
            grid_2d<double> output(A._rows, A._rows);

            auto t = 64; // Cache line
            for (size_t j = 0; j < A._rows; j += t)
            {
                for (size_t k = 0; k < A._cols; k += t)
                {
                    for (size_t i = 0; i < A._rows; i += t)
                    {
                        for (size_t ii = i; ii < std::min((size_t)i + t, A._rows); ii++)
                        {
                            for (size_t jj = j; jj < std::min((size_t)j + t, A._rows); jj++)
                            {
                                auto tmp = output._grid[ii * output._cols + jj];
                                for (size_t kk = k; kk < std::min((size_t)k + t, A._cols); kk++)
                                {
                                    auto product = A._grid[ii *A._cols + kk] * A._grid[jj * A._cols + kk];
                                    if (!std::isnan(product))
                                        tmp += product;
                                }
                                output._grid[ii *output._cols + jj] = tmp;
                            }
                        }
                    }
                }
            }
            return output;
        }
    };

    struct poly_coeffs
    {
        double x, y, x2, y2, x3, y3;
    };

    template<typename T>
    class scaler_2d
    {
    public:
        virtual void scale(image<T> &input, image<T> &output) = 0;
        virtual ~scaler_2d() {}
    };


    // Polynomial scaler base class
    template <typename T>
    class poly_scaler_2d : public scaler_2d<T>
    {
    public:
        virtual void scale(image<T> &input, image<T> &output);
        void scale(grid_2d<T> &input, grid_2d<T> &output)
        {
            auto in_grid = input.to_view();
            auto out_grid = output.to_view();
            scale(in_grid, out_grid);
        }

        ~poly_scaler_2d() {}

        poly_scaler_2d(int coeffs) : num_coeff(coeffs)
        {
            LU = grid_2d<double>(num_coeff, num_coeff);
            x = std::vector<double>(num_coeff);
        }

        void poly_fit_2d(image<T> &input);
        void interp(image<T> &output);

    protected:
        // Polynomial scale protected members
        int num_coeff;

        void lu_decomposition(grid_2d<double> &A);
        void lu_solve();
        virtual void get_basis_function(double* ptr, const poly_coeffs &coeffs) = 0;
        virtual void set_polynomial_value(T &ptr, std::vector<double> &arr, const poly_coeffs &coeffs) = 0;

        grid_2d<double> LU;
        grid_2d<double> b;
        std::vector<double> x;
    };

    // Computes 3th order polynomial model for 2d grid
    // Interpolates the model to the output grid
    template <typename T>
    class poly_3_scaler_2d : public poly_scaler_2d<T>
    {
    public:
        ~poly_3_scaler_2d() {}

        poly_3_scaler_2d() : poly_scaler_2d<T>(10)
        {
        }

    protected:
        virtual void get_basis_function(double* ptr, const poly_coeffs &coeffs)
        {
            ptr[0] = 1.0;
            ptr[1] = coeffs.x;
            ptr[2] = coeffs.y;
            ptr[3] = coeffs.x2;
            ptr[4] = coeffs.x * coeffs.y;
            ptr[5] = coeffs.y2;
            ptr[6] = coeffs.x3;
            ptr[7] = coeffs.x2 * coeffs.y;
            ptr[8] = coeffs.x * coeffs.y2;
            ptr[9] = coeffs.y3;
        }
        virtual void set_polynomial_value(T &ptr, std::vector<double> &arr, const poly_coeffs &coeffs)
        {
            ptr = reduce_to<T>(get_cubic_polynomial_value(arr, coeffs));
        }

        static double get_cubic_polynomial_value(std::vector<double> &arr, const poly_coeffs &coeffs)
        {
            return arr[0] +
                arr[1] * coeffs.x +
                arr[2] * coeffs.y +
                arr[3] * coeffs.x2 +
                arr[4] * coeffs.x * coeffs.y +
                arr[5] * coeffs.y2 +
                arr[6] * coeffs.x3 +
                arr[7] * coeffs.x2 * coeffs.y +
                arr[8] * coeffs.x * coeffs.y2 +
                arr[9] * coeffs.y3;
        }
    };

    // Computes 4th order polynomial model for 2d grid
    // Interpolates the model to the output grid
    template <typename T>
    class poly_4_scaler_2d : public poly_scaler_2d<T>
    {
    public:
        ~poly_4_scaler_2d() {}

        poly_4_scaler_2d() : poly_scaler_2d<T>(15)
        {
        }

    protected:
        virtual void get_basis_function(double* ptr, const poly_coeffs &coeffs)
        {
            ptr[0] = 1.0;
            ptr[1] = coeffs.x;
            ptr[2] = coeffs.y;
            ptr[3] = coeffs.x2;
            ptr[4] = (coeffs.x * coeffs.y);
            ptr[5] = coeffs.y2;
            ptr[6] = coeffs.x3;
            ptr[7] = coeffs.x2 * coeffs.y;
            ptr[8] = coeffs.x * coeffs.y2;
            ptr[9] = coeffs.y3;
            ptr[10] = coeffs.x3 * coeffs.x;
            ptr[11] = coeffs.x3 * coeffs.y;
            ptr[12] = coeffs.x2 * coeffs.y2;
            ptr[13] = coeffs.x * coeffs.y3;
            ptr[14] = coeffs.y3 * coeffs.y;
        }
        virtual void set_polynomial_value(T &ptr, std::vector<double> &arr, const poly_coeffs &coeffs)
        {
            ptr = reduce_to<T>(get_quadric_polynomial_value(arr, coeffs));
        }
        static double get_quadric_polynomial_value(std::vector<double> &arr, const poly_coeffs &coeffs)
        {
            return arr[0] +
                arr[1] * coeffs.x +
                arr[2] * coeffs.y +
                arr[3] * coeffs.x2 +
                arr[4] * coeffs.x * coeffs.y +
                arr[5] * coeffs.y2 +
                arr[6] * coeffs.x3 +
                arr[7] * coeffs.x2 * coeffs.y +
                arr[8] * coeffs.x * coeffs.y2 +
                arr[9] * coeffs.y3 +
                arr[10] * coeffs.x3 * coeffs.x +
                arr[11] * coeffs.x3 * coeffs.y +
                arr[12] * coeffs.x2 * coeffs.y2 +
                arr[13] * coeffs.x * coeffs.y3 +
                arr[14] * coeffs.y3 * coeffs.y;
        }
    };

    template <typename T>
    void poly_scaler_2d<T>::poly_fit_2d(image<T> &input)
    {
        size_t cols = input._width;
        size_t rows = input._height;
        poly_coeffs coeffs;
        // Basis function
        grid_2d<double> A(num_coeff, cols * rows);

        // Fit model as center symmetric
        auto h2 = (rows - 1) * 0.5;
        auto w2 = (cols - 1) * 0.5;
        auto inv_h2 = h2 == 0.0 ? 0.0 : 1.0 / h2;
        auto inv_w2 = w2 == 0.0 ? 0.0 : 1.0 / w2;

        for (size_t j = 0; j < rows; j++)
        {
            coeffs.y = (2 * j - h2) * inv_h2;
            coeffs.y2 = coeffs.y * coeffs.y;
            coeffs.y3 = coeffs.y2 * coeffs.y;
            for (size_t i = 0; i < cols; i++)
            {
                coeffs.x = (2 * i - w2) * inv_w2;
                coeffs.x2 = coeffs.x * coeffs.x;
                coeffs.x3 = coeffs.x2 * coeffs.x;
                get_basis_function(&A._grid[j * cols * num_coeff + i * num_coeff], coeffs);
            }
        }

        auto A_trans = A.transpose();
        auto product = A_trans.mul_mm();    // A_trans * A
        lu_decomposition(product);

        // Treat input as vector for matrix multiplication
        grid_2d<T> InputAsVector(input);
        InputAsVector._rows *= InputAsVector._cols;
        InputAsVector._cols = 1;

        // b = A' * grid(:)
        b = A_trans.mul_mm(InputAsVector);

        lu_solve();
    }

    template <typename T>
    void poly_scaler_2d<T>::lu_decomposition(grid_2d<double> &A)
    {
        auto w = A._cols;

        for (size_t k = 0; k < w; k++)
        {
            // U(k,k) = A(k,k);
            LU._grid[w*k + k] = A._grid[w*k + k];
            for (size_t j = k + 1; j < w; j++)
            {
                // L(k,j) = A(k,j) / U(k,k);
                LU._grid[w*j + k] = A._grid[w*j + k] / (double)LU._grid[w*k + k];

                // U(j,k) = A(j,k);
                LU._grid[w*k + j] = A._grid[w*k + j];
            }
            for (size_t j = k + 1; j < w; j++)
            {
                for (size_t i = k + 1; i < w; i++)
                {
                    // A(i,j) = A(i, j) - L(k, j) * U(i, k);
                    A._grid[w*j + i] -= LU._grid[w*j + k] * LU._grid[w*k + i];
                }
            }
        }
    }

    template <typename T>
    void poly_scaler_2d<T>::lu_solve()
    {
        std::vector<double> y(num_coeff);

        for (auto j = 0; j < num_coeff; j++)
        {
            double sigma = 0.0;
            for (auto i = 0; i <= j - 1; i++)
            {
                sigma += LU._grid[num_coeff * j + i] * y[i];
            }
            y[j] = b._grid[j] - sigma;
        }

        for (int j = num_coeff - 1; j >= 0; j--)
        {
            double sigma = 0.0;
            for (auto i = j + 1; i < num_coeff; i++)
            {
                sigma += LU._grid[j * num_coeff + i] * x[i];
            }
            x[j] = (y[j] - sigma) / LU._grid[num_coeff * j + j];
        }
    }

    template <typename T>
    void poly_scaler_2d<T>::interp(image<T> &output)
    {
        poly_coeffs coeffs;

        auto h2 = (output._height - 1) * 0.5;
        auto w2 = (output._width - 1) * 0.5;
        auto rows = output._height;
        auto cols = output._width;

        for (decltype(rows) i = 0; i < rows; i++)
        {
            // Code working, Aki to optimize the runtime speed
            coeffs.y = (2 * i - h2) / h2;
            coeffs.y2 = coeffs.y * coeffs.y;
            coeffs.y3 = coeffs.y2 * coeffs.y;

            T *ptr = &output.at(i, 0);
            int skipx = output._skip_x;

            for (decltype(cols) j = 0; j < cols; j++)
            {
                coeffs.x = (2 * j - w2) / w2;
                coeffs.x2 = coeffs.x * coeffs.x;
                coeffs.x3 = coeffs.x2 * coeffs.x;

                set_polynomial_value(*ptr, x, coeffs);
                ptr += skipx;
            }
        }
    }

    template <typename T>
    void poly_scaler_2d<T>::scale(image<T> &input, image<T> &output)
    {
        // LU decomposition
        poly_fit_2d(input);

        // Do resizing
        interp(output);
    }
}