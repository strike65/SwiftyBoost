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
// Gamma, Error functions, and incomplete gamma helpers
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/erf.hpp>
#include "../internal/bs_internal.hpp"

extern "C" {

// Gamma / Error
double bs_tgamma_d(double x)              { return bs_wrap<double>([&] { return boost::math::tgamma(x); }); }
double bs_lgamma_d(double x)              { return bs_wrap<double>([&] { return boost::math::lgamma(x); }); }
double bs_erf_d(double x)                 { return bs_wrap<double>([&] { return boost::math::erf(x); }); }
double bs_erfc_d(double x)                { return bs_wrap<double>([&] { return boost::math::erfc(x); }); }
double bs_erf_inv_d(double p)             { return bs_wrap<double>([&] { return boost::math::erf_inv(p); }); }
double bs_erfc_inv_d(double p)            { return bs_wrap<double>([&] { return boost::math::erfc_inv(p); }); }

float bs_tgamma_f(float x)              { return bs_wrap<float>([&] { return boost::math::tgamma(x); }); }
float bs_lgamma_f(float x)              { return bs_wrap<float>([&] { return boost::math::lgamma(x); }); }
float bs_erf_f(float x)                 { return bs_wrap<float>([&] { return boost::math::erf(x); }); }
float bs_erfc_f(float x)                { return bs_wrap<float>([&] { return boost::math::erfc(x); }); }
float bs_erf_inv_f(float p)             { return bs_wrap<float>([&] { return boost::math::erf_inv(p); }); }
float bs_erfc_inv_f(float p)            { return bs_wrap<float>([&] { return boost::math::erfc_inv(p); }); }

long double bs_tgamma_l(long double x)  { return bs_wrap<long double>([&] { return boost::math::tgamma(x); }); }
long double bs_lgamma_l(long double x)  { return bs_wrap<long double>([&] { return boost::math::lgamma(x); }); }
long double bs_erf_l(long double x)     { return bs_wrap<long double>([&] { return boost::math::erf(x); }); }
long double bs_erfc_l(long double x)    { return bs_wrap<long double>([&] { return boost::math::erfc(x); }); }
long double bs_erf_inv_l(long double p) { return bs_wrap<long double>([&] { return boost::math::erf_inv(p); }); }
long double bs_erfc_inv_l(long double p){ return bs_wrap<long double>([&] { return boost::math::erfc_inv(p); }); }

// Incomplete gamma (lower/upper, regularized, and inverses)
double bs_tgamma_lower_d(double a, double x)      { return bs_wrap<double>([&] { return boost::math::tgamma_lower(a, x); }); }
double bs_tgamma_upper_d(double a, double x)      { return bs_wrap<double>([&] { return boost::math::tgamma(a, x); }); }
double bs_gamma_p_d(double a, double x)           { return bs_wrap<double>([&] { return boost::math::gamma_p(a, x); }); }
double bs_gamma_q_d(double a, double x)           { return bs_wrap<double>([&] { return boost::math::gamma_q(a, x); }); }
double bs_gamma_p_inv_d(double a, double p)       { return bs_wrap<double>([&] { return boost::math::gamma_p_inv(a, p); }); }
double bs_gamma_q_inv_d(double a, double q)       { return bs_wrap<double>([&] { return boost::math::gamma_q_inv(a, q); }); }

// Derivatives of regularized incomplete gamma functions
double bs_gamma_p_derivative_d(double a, double x) { return bs_wrap<double>([&] { return boost::math::gamma_p_derivative(a, x); }); }

float bs_tgamma_lower_f(float a, float x)       { return bs_wrap<float>([&] { return boost::math::tgamma_lower(a, x); }); }
float bs_tgamma_upper_f(float a, float x)       { return bs_wrap<float>([&] { return boost::math::tgamma(a, x); }); }
float bs_gamma_p_f(float a, float x)            { return bs_wrap<float>([&] { return boost::math::gamma_p(a, x); }); }
float bs_gamma_q_f(float a, float x)            { return bs_wrap<float>([&] { return boost::math::gamma_q(a, x); }); }
float bs_gamma_p_inv_f(float a, float p)        { return bs_wrap<float>([&] { return boost::math::gamma_p_inv(a, p); }); }
float bs_gamma_q_inv_f(float a, float q)        { return bs_wrap<float>([&] { return boost::math::gamma_q_inv(a, q); }); }

// Derivatives (float)
float bs_gamma_p_derivative_f(float a, float x) { return bs_wrap<float>([&] { return boost::math::gamma_p_derivative(a, x); }); }

long double bs_tgamma_lower_l(long double a, long double x) { return bs_wrap<long double>([&] { return boost::math::tgamma_lower(a, x); }); }
long double bs_tgamma_upper_l(long double a, long double x) { return bs_wrap<long double>([&] { return boost::math::tgamma(a, x); }); }
long double bs_gamma_p_l(long double a, long double x)      { return bs_wrap<long double>([&] { return boost::math::gamma_p(a, x); }); }
long double bs_gamma_q_l(long double a, long double x)      { return bs_wrap<long double>([&] { return boost::math::gamma_q(a, x); }); }
long double bs_gamma_p_inv_l(long double a, long double p)  { return bs_wrap<long double>([&] { return boost::math::gamma_p_inv(a, p); }); }
long double bs_gamma_q_inv_l(long double a, long double q)  { return bs_wrap<long double>([&] { return boost::math::gamma_q_inv(a, q); }); }

// Derivatives (long double)
long double bs_gamma_p_derivative_l(long double a, long double x) { return bs_wrap<long double>([&] { return boost::math::gamma_p_derivative(a, x); }); }

// Γ(a) / Γ(b) and Γ(a) / Γ(a+delta)
double bs_tgamma_ratio_d(double a, double b) { return bs_wrap<double>([&] { return boost::math::tgamma_ratio(a, b); }); }
float  bs_tgamma_ratio_f(float a, float b) { return bs_wrap<float>([&] { return boost::math::tgamma_ratio(a, b); }); }
long double bs_tgamma_ratio_l(long double a, long double b) { return bs_wrap<long double>([&] { return boost::math::tgamma_ratio(a, b); }); }
double bs_tgamma_delta_ratio_d(double a, double delta)  { return bs_wrap<double>([&] { return boost::math::tgamma_delta_ratio(a, delta); }); }
float  bs_tgamma_delta_ratio_f(float a, float delta)  { return bs_wrap<float>([&] { return boost::math::tgamma_delta_ratio(a, delta); }); }
long double bs_tgamma_delta_ratio_l(long double a, long double delta)  { return bs_wrap<long double>([&] { return boost::math::tgamma_delta_ratio(a, delta); }); }

} // extern "C"

