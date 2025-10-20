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
} bs_complex_f;
typedef struct {
    double re, im;
} bs_complex_d;
typedef struct {
    long double re, im;
} bs_complex_l;

// Elementary complex arithmetic
bs_complex_d bs_cadd_d(bs_complex_d a, bs_complex_d b);
bs_complex_d bs_csub_d(bs_complex_d a, bs_complex_d b);
bs_complex_d bs_cmul_d(bs_complex_d a, bs_complex_d b);
bs_complex_d bs_cdiv_d(bs_complex_d a, bs_complex_d b);

bs_complex_f bs_cadd_f(bs_complex_f a, bs_complex_f b);
bs_complex_f bs_csub_f(bs_complex_f a, bs_complex_f b);
bs_complex_f bs_cmul_f(bs_complex_f a, bs_complex_f b);
bs_complex_f bs_cdiv_f(bs_complex_f a, bs_complex_f b);

bs_complex_l bs_cadd_l(bs_complex_l a, bs_complex_l b);
bs_complex_l bs_csub_l(bs_complex_l a, bs_complex_l b);
bs_complex_l bs_cmul_l(bs_complex_l a, bs_complex_l b);
bs_complex_l bs_cdiv_l(bs_complex_l a, bs_complex_l b);

// Elementary complex functions
bs_complex_d bs_cexp_d(bs_complex_d z);
bs_complex_d bs_clog_d(bs_complex_d z);
bs_complex_d bs_csqrt_d(bs_complex_d z);
bs_complex_d bs_csin_d(bs_complex_d z);
bs_complex_d bs_ccos_d(bs_complex_d z);
bs_complex_d bs_ctan_d(bs_complex_d z);
bs_complex_d bs_csinh_d(bs_complex_d z);
bs_complex_d bs_ccosh_d(bs_complex_d z);
bs_complex_d bs_ctanh_d(bs_complex_d z);
bs_complex_d bs_catan_d(bs_complex_d z);

bs_complex_f bs_cexp_f(bs_complex_f z);
bs_complex_f bs_clog_f(bs_complex_f z);
bs_complex_f bs_csqrt_f(bs_complex_f z);
bs_complex_f bs_csin_f(bs_complex_f z);
bs_complex_f bs_ccos_f(bs_complex_f z);
bs_complex_f bs_ctan_f(bs_complex_f z);
bs_complex_f bs_csinh_f(bs_complex_f z);
bs_complex_f bs_ccosh_f(bs_complex_f z);
bs_complex_f bs_ctanh_f(bs_complex_f z);
bs_complex_f bs_catan_f(bs_complex_f z);

bs_complex_l bs_cexp_l(bs_complex_l z);
bs_complex_l bs_clog_l(bs_complex_l z);
bs_complex_l bs_csqrt_l(bs_complex_l z);
bs_complex_l bs_csin_l(bs_complex_l z);
bs_complex_l bs_ccos_l(bs_complex_l z);
bs_complex_l bs_ctan_l(bs_complex_l z);
bs_complex_l bs_csinh_l(bs_complex_l z);
bs_complex_l bs_ccosh_l(bs_complex_l z);
bs_complex_l bs_ctanh_l(bs_complex_l z);
bs_complex_l bs_catan_l(bs_complex_l z);

#ifdef __cplusplus
}
#endif

