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

// Beta family
double bs_beta(double a, double b);
double bs_fullBeta(double a, double b, double x);
double bs_ibeta(double a, double b, double x);
double bs_ibetac(double a, double b, double x);
double bs_ibeta_inv(double a, double b, double p);
double bs_ibetac_inv(double a, double b, double p);
double bs_ibeta_inva(double b, double x, double p);
double bs_ibeta_invb(double a, double x, double p);
// Derivative of regularized incomplete Beta I_x(a, b) w.r.t. x
double bs_ibeta_derivative(double a, double b, double x);

float bs_beta_f(float a, float b);
float bs_fullBeta_f(float a, float b, float x);
float bs_ibeta_f(float a, float b, float x);
float bs_ibetac_f(float a, float b, float x);
float bs_ibeta_inv_f(float a, float b, float p);
float bs_ibetac_inv_f(float a, float b, float p);
float bs_ibeta_inva_f(float b, float x, float p);
float bs_ibeta_invb_f(float a, float x, float p);
// Derivative (float)
float bs_ibeta_derivative_f(float a, float b, float x);

long double bs_beta_l(long double a, long double b);
long double bs_fullBeta_l(long double a, long double b, long double x);
long double bs_ibeta_l(long double a, long double b, long double x);
long double bs_ibetac_l(long double a, long double b, long double x);
long double bs_ibeta_inv_l(long double a, long double b, long double p);
long double bs_ibetac_inv_l(long double a, long double b, long double p);
long double bs_ibeta_inva_l(long double b, long double x, long double p);
long double bs_ibeta_invb_l(long double a, long double x, long double p);
// Derivative (long double)
long double bs_ibeta_derivative_l(long double a, long double b, long double x);

#ifdef __cplusplus
}
#endif

