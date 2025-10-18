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

// Gegenbauer polynomials
double bs_gegenbauer(unsigned int n, double lambda, double x);
double bs_gegenbauer_prime(unsigned int n, double lambda, double x);
double bs_gegenbauer_derivative(unsigned int n, double lambda, double x, unsigned int k);
float bs_gegenbauer_f(unsigned int n, float lambda, float x);
float bs_gegenbauer_prime_f(unsigned int n, float lambda, float x);
float bs_gegenbauer_derivative_f(unsigned int n, float lambda, float x, unsigned int k);
long double bs_gegenbauer_l(unsigned int n, long double lambda, long double x);
long double bs_gegenbauer_prime_l(unsigned int n, long double lambda, long double x);
long double bs_gegenbauer_derivative_l(unsigned int n, long double lambda, long double x, unsigned int k);

#ifdef __cplusplus
}
#endif

