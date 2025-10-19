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

// Hypergeometric functions (Gauss/confluent/general)
double bs_hypergeometric_1F0(double a, double z);
float bs_hypergeometric_1F0_f(float a, float z);
long double bs_hypergeometric_1F0_l(long double a, long double z);

double bs_hypergeometric_0F1(double b, double z);
float bs_hypergeometric_0F1_f(float b, float z);
long double bs_hypergeometric_0F1_l(long double b, long double z);

double bs_hypergeometric_2F0(double a, double b, double z);
float bs_hypergeometric_2F0_f(float a, float b, float z);
long double bs_hypergeometric_2F0_l(long double a, long double b, long double z);

double bs_hypergeometric_1F1(double a, double b, double z);
float bs_hypergeometric_1F1_f(float a, float b, float z);
long double bs_hypergeometric_1F1_l(long double a, long double b, long double z);

// General pFq with arrays of parameters
double bs_hypergeometric_pFq(const double *a, size_t p, const double *b, size_t q, double z);
float bs_hypergeometric_pFq_f(const float *a, size_t p, const float *b, size_t q, float z);
long double bs_hypergeometric_pFq_l(const long double *a, size_t p, const long double *b, size_t q, long double z);

#ifdef __cplusplus
}
#endif

