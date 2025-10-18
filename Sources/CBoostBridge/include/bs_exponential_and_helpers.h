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

// Exponential integrals and related helpers
double bs_expint_Ei(double x);
double bs_expint_En(int n, double x);
double bs_expm1(double x);
double bs_log1p(double x);
double bs_log1pmx(double x);
double bs_powm1(double x, double y);
double bs_cbrt(double x);
double bs_sqrt1pm1(double x);
double bs_hypot(double x, double y);
double bs_rsqrt(const double x);

float bs_expint_Ei_f(float x);
float bs_expint_En_f(int n, float x);
float bs_expm1_f(float x);
float bs_log1p_f(float x);
float bs_log1pmx_f(float x);
float bs_powm1_f(float x, float y);
float bs_cbrt_f(float x);
float bs_sqrt1pm1_f(float x);
float bs_hypot_f(float x, float y);
float bs_rsqrt_f(const float x);

long double bs_expint_Ei_l(long double x);
long double bs_expint_En_l(int n, long double x);
long double bs_expm1_l(long double x);
long double bs_log1p_l(long double x);
long double bs_log1pmx_l(long double x);
long double bs_powm1_l(long double x, long double y);
long double bs_cbrt_l(long double x);
long double bs_sqrt1pm1_l(long double x);
long double bs_hypot_l(long double x, long double y);
long double bs_rsqrt_l(const long double x);

#ifdef __cplusplus
}
#endif

