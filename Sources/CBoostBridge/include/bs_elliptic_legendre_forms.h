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

// Elliptic integrals (Legendre forms)
double bs_ellint_1_complete(double k);
double bs_ellint_1(double k, double phi);
double bs_ellint_2_complete(double k);
double bs_ellint_2(double k, double phi);
double bs_ellint_3(double k, double nu, double phi);
double bs_ellint_3_complete(double k, double nu);

float bs_ellint_1_complete_f(float k);
float bs_ellint_1_f(float k, float phi);
float bs_ellint_2_complete_f(float k);
float bs_ellint_2_f(float k, float phi);
float bs_ellint_3_f(float k, float nu, float phi);
float bs_ellint_3_complete_f(float k, float nu);

long double bs_ellint_1_complete_l(long double k);
long double bs_ellint_1_l(long double k, long double phi);
long double bs_ellint_2_complete_l(long double k);
long double bs_ellint_2_l(long double k, long double phi);
long double bs_ellint_3_l(long double k, long double nu, long double phi);
long double bs_ellint_3_complete_l(long double k, long double nu);

#ifdef __cplusplus
}
#endif

