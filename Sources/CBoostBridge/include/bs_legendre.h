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

#ifdef __cplusplus
extern "C" {
#endif

// Legendre
double bs_legendre_p(int n, double x);
double bs_assoc_legendre_p(int n, int m, double x);
double bs_legendre_p_prime(int n, double x);
void bs_legendre_p_zeros(int l, double *out);

float bs_legendre_p_f(int n, float x);
float bs_assoc_legendre_p_f(int n, int m, float x);
float bs_legendre_p_prime_f(int n, float x);
void bs_legendre_p_zeros_f(int l, float *out);

long double bs_legendre_p_l(int n, long double x);
long double bs_assoc_legendre_p_l(int n, int m, long double x);
long double bs_legendre_p_prime_l(int n, long double x);
void bs_legendre_p_zeros_l(int l, long double *out);

#ifdef __cplusplus
}
#endif

