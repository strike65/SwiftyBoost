//
//  Created by Volker Thieme 2025.
//  Copyright © 2025 Volker Thieme.
//
//  Permission is hereby granted, free of charge, to any person obtaining a copy
//  of this software and associated documentation files (the "Software"), to deal
//  in the Software without restriction, including without limitation the rights
//  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//  copies of the Software, and to permit persons to whom the Software is
//  furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included in
//  all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
//  THE SOFTWARE.
//
#pragma once

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

// Chebyshev
double bs_chebyshev_T(unsigned int n, double x);
double bs_chebyshev_U(unsigned int n, double x);
float bs_chebyshev_T_f(unsigned int n, float x);
float bs_chebyshev_U_f(unsigned int n, float x);
long double bs_chebyshev_T_l(unsigned int n, long double x);
long double bs_chebyshev_U_l(unsigned int n, long double x);

// Chebyshev series evaluation via Clenshaw (first kind)
double bs_chebyshev_clenshaw(const double *c, size_t count, double x);
float bs_chebyshev_clenshaw_f(const float *c, size_t count, float x);
long double bs_chebyshev_clenshaw_l(const long double *c, size_t count, long double x);

#ifdef __cplusplus
}
#endif

