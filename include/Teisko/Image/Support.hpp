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
#include <algorithm>
#include <vector>

namespace Teisko
{
    /// Auxiliary function to manage supports for filtering kernels
    /// the matrix is first shifted, then W-1 th column is read from given address
    template <int ROWS, int COLUMNS, typename T>
    struct support {
        T data[ROWS][COLUMNS];

        // Default constructor initializing the support to zeros
        explicit support(T value = 0)
        {
            for (int i = 0; i < ROWS; ++i)
                for (int j = 0; j < COLUMNS; ++j)
                    data[i][j] = value;
        }

        // Initialize from 2d-array
        support(T *begin, int skip_y, int skip_x = 1)
        {
            for (int i = 0; i < ROWS; ++i)
                for (int j = 0; j < COLUMNS; ++j)
                    data[i][j] = begin[i*skip_y + j * skip_x];
        }

        // fill from 2d-array (same functionality as in the constructor)
        void fill(T *begin, int skip_y, int skip_x = 1)
        {
            for (int i = 0; i < ROWS; ++i)
                for (int j = 0; j < COLUMNS; ++j)
                    data[i][j] = begin[i*skip_y + j * skip_x];
        }

        // Returns pointer to row, allowing 'support[row][column]' access
        // instead of my_support.data[row][column]
        T* operator[](int row) { return &data[row][0]; }

        // Shift in new data to rightmost column data[j][COL-1] from begin[j * skip]
        template <typename U>
        void shift(U *begin, int skip_y)
        {
            for (int i = 0; i < ROWS; ++i)
            {
                for (int j = 0; j < COLUMNS - 1; ++j)
                {
                    data[i][j] = data[i][j + 1];
                }
                data[i][COLUMNS - 1] = (T)begin[skip_y * i];
            }
        }

        // quick hack -- median is properly defined only when ROWS/COLS are both odd
        // when ROWS, COLS == 7, please use the method in libprepro/libimage_algorithms.hpp
        T median() const
        {
            auto copy = *this;
            auto begin = reinterpret_cast<T*>(const_cast<T*>(&copy[0][0]));
            std::sort(begin, begin + ROWS*COLUMNS);
            return copy[ROWS / 2][COLUMNS / 2];
        }

        // returns minimum value within the support
        T minimum() const
        {
            auto begin = reinterpret_cast<T*>(const_cast<T*>(&data[0][0]));
            return *std::min_element(begin, begin + ROWS*COLUMNS);
        }

        // returns maximum value within the support
        T maximum() const
        {
            auto begin = reinterpret_cast<T*>(const_cast<T*>(&data[0][0]));
            return *std::max_element(begin, begin + ROWS*COLUMNS);
        }
    };

    /// Helper function to serialize 2-d array
    /// - used primarily in specs
    template <int ROWS, int COLS, typename T>
    std::vector<T> to_vector(T(&arr)[ROWS][COLS])
    {
        return std::vector<T>(&arr[0][0], &arr[0][0] + ROWS*COLS);
    }
}