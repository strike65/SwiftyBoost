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

// Factorials
double bs_factorial(unsigned int i);
float bs_factorial_f(unsigned int i);
long double bs_factorial_l(unsigned int i);

// Pochhammer (rising factorial)
double bs_rising_factorial(double x, unsigned int i);
float bs_rising_factorial_f(float x, unsigned int i);
long double bs_rising_factorial_l(long double x, unsigned int i);

// Binomial coefficients
double bs_binomial_coefficient(unsigned int n, unsigned int k);
float bs_binomial_coefficient_f(unsigned int n, unsigned int k);
long double bs_binomial_coefficient_l(unsigned int n, unsigned int k);

// Double factorial
double bs_double_factorial(unsigned int i);
float bs_double_factorial_f(unsigned int i);
long double bs_double_factorial_l(unsigned int i);

#ifdef __cplusplus
}
#endif

