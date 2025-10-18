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

#include "bs_complex.h"

double bs_sinc_pi(double x);
bs_complex_d bs_sincc_pi(bs_complex_d x);
float bs_sinc_pi_f(float x);
bs_complex_f bs_sincc_pi_f(bs_complex_f x);
long double bs_sinc_pi_l(long double x);
bs_complex_l bs_sincc_pi_l(bs_complex_l x);

double bs_sinhc_pi(double x);
bs_complex_d bs_sinhcc_pi(bs_complex_d x);
float bs_sinhc_pi_f(float x);
bs_complex_f bs_sinhcc_pi_f(bs_complex_f x);
long double bs_sinhc_pi_l(long double x);
bs_complex_l bs_sinhcc_pi_l(bs_complex_l x);

#ifdef __cplusplus
}
#endif

