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

// POD layout compatible with Swift's import.
typedef struct {
    float re, im;
} complex_f;
typedef struct {
    double re, im;
} complex_d;
typedef struct {
    long double re, im;
} complex_l;

// Elementary complex arithmetic
complex_d bs_cadd_d(complex_d a, complex_d b);
complex_d bs_csub_d(complex_d a, complex_d b);
complex_d bs_cmul_d(complex_d a, complex_d b);
complex_d bs_cdiv_d(complex_d a, complex_d b);

complex_f bs_cadd_f(complex_f a, complex_f b);
complex_f bs_csub_f(complex_f a, complex_f b);
complex_f bs_cmul_f(complex_f a, complex_f b);
complex_f bs_cdiv_f(complex_f a, complex_f b);

complex_l bs_cadd_l(complex_l a, complex_l b);
complex_l bs_csub_l(complex_l a, complex_l b);
complex_l bs_cmul_l(complex_l a, complex_l b);
complex_l bs_cdiv_l(complex_l a, complex_l b);

// Elementary complex functions
complex_d bs_cexp_d(complex_d z);
complex_d bs_clog_d(complex_d z);
complex_d bs_csqrt_d(complex_d z);
complex_d bs_csin_d(complex_d z);
complex_d bs_ccos_d(complex_d z);
complex_d bs_ctan_d(complex_d z);
complex_d bs_csinh_d(complex_d z);
complex_d bs_ccosh_d(complex_d z);
complex_d bs_ctanh_d(complex_d z);
complex_d bs_catan_d(complex_d z);

complex_f bs_cexp_f(complex_f z);
complex_f bs_clog_f(complex_f z);
complex_f bs_csqrt_f(complex_f z);
complex_f bs_csin_f(complex_f z);
complex_f bs_ccos_f(complex_f z);
complex_f bs_ctan_f(complex_f z);
complex_f bs_csinh_f(complex_f z);
complex_f bs_ccosh_f(complex_f z);
complex_f bs_ctanh_f(complex_f z);
complex_f bs_catan_f(complex_f z);

complex_l bs_cexp_l(complex_l z);
complex_l bs_clog_l(complex_l z);
complex_l bs_csqrt_l(complex_l z);
complex_l bs_csin_l(complex_l z);
complex_l bs_ccos_l(complex_l z);
complex_l bs_ctan_l(complex_l z);
complex_l bs_csinh_l(complex_l z);
complex_l bs_ccosh_l(complex_l z);
complex_l bs_ctanh_l(complex_l z);
complex_l bs_catan_l(complex_l z);

#ifdef __cplusplus
}
#endif

