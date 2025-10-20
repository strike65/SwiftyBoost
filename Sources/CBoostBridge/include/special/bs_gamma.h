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

// Gamma and log-gamma
double bs_tgamma_d(double x);
double bs_lgamma_d(double x);
double bs_tgamma_ratio_d(double a, double b);
double bs_tgamma_delta_ratio_d(double a, double delta);

float bs_tgamma_f(float x);
float bs_lgamma_f(float x);
float bs_tgamma_ratio_f(float a, float b);
float bs_tgamma_delta_ratio_f(float a, float delta);

long double bs_tgamma_l(long double x);
long double bs_lgamma_l(long double x);
long double bs_tgamma_ratio_l(long double a, long double b);
long double bs_tgamma_delta_ratio_l(long double a, long double delta);

// Incomplete gamma (lower/upper, regularized, and inverses)
double bs_tgamma_lower_d(double a, double x);
double bs_tgamma_upper_d(double a, double x);
double bs_gamma_p_d(double a, double x);
double bs_gamma_q_d(double a, double x);
double bs_gamma_p_inv_d(double a, double p);
double bs_gamma_q_inv_d(double a, double q);

// Derivatives of regularized incomplete gamma (w.r.t. x)
double bs_gamma_p_derivative_d(double a, double x);

float bs_tgamma_lower_f(float a, float x);
float bs_tgamma_upper_f(float a, float x);
float bs_gamma_p_f(float a, float x);
float bs_gamma_q_f(float a, float x);
float bs_gamma_p_inv_f(float a, float p);
float bs_gamma_q_inv_f(float a, float q);

// Derivatives (float)
float bs_gamma_p_derivative_f(float a, float x);

long double bs_tgamma_lower_l(long double a, long double x);
long double bs_tgamma_upper_l(long double a, long double x);
long double bs_gamma_p_l(long double a, long double x);
long double bs_gamma_q_l(long double a, long double x);
long double bs_gamma_p_inv_l(long double a, long double p);
long double bs_gamma_q_inv_l(long double a, long double q);

// Derivatives (long double)
long double bs_gamma_p_derivative_l(long double a, long double x);

#ifdef __cplusplus
}
#endif

