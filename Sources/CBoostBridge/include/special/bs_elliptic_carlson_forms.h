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

// Elliptic integrals (Carlson symmetric forms)
double bs_ellint_rc(double x, double y);
double bs_ellint_rf(double x, double y, double z);
double bs_ellint_rd(double x, double y, double z);
double bs_ellint_rj(double x, double y, double z, double p);
double bs_ellint_rg(double x, double y, double z);

float bs_ellint_rc_f(float x, float y);
float bs_ellint_rf_f(float x, float y, float z);
float bs_ellint_rd_f(float x, float y, float z);
float bs_ellint_rj_f(float x, float y, float z, float p);
float bs_ellint_rg_f(float x, float y, float z);

long double bs_ellint_rc_l(long double x, long double y);
long double bs_ellint_rf_l(long double x, long double y, long double z);
long double bs_ellint_rd_l(long double x, long double y, long double z);
long double bs_ellint_rj_l(long double x, long double y, long double z, long double p);
long double bs_ellint_rg_l(long double x, long double y, long double z);

#ifdef __cplusplus
}
#endif

