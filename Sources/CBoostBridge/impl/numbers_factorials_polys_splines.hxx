// Numbers (Bernoulli, Tangent, Fibonacci, Prime), Factorials/Combinatorics,
// Laguerre, Chebyshev, Spherical Harmonics, Cardinal B-splines, Gegenbauer
#include <boost/math/special_functions/bernoulli.hpp>
#include <boost/math/special_functions/fibonacci.hpp>
#include <boost/math/special_functions/prime.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/chebyshev.hpp>
#include <boost/math/special_functions/laguerre.hpp>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <boost/math/special_functions/cardinal_b_spline.hpp>
#include <boost/math/special_functions/gegenbauer.hpp>
#include <algorithm>
#include <vector>
#include "../internal/bs_internal.hpp"
#include "../include/bs_complex.h"

// Internal dispatchers for cardinal B-spline (mirroring original)
template <typename Real>
static inline Real dispatch_cardinal_b_spline(unsigned int n, Real x) {
    switch (n) {
        case 0:  return boost::math::cardinal_b_spline<0,  Real>(x);
        case 1:  return boost::math::cardinal_b_spline<1,  Real>(x);
        case 2:  return boost::math::cardinal_b_spline<2,  Real>(x);
        case 3:  return boost::math::cardinal_b_spline<3,  Real>(x);
        case 4:  return boost::math::cardinal_b_spline<4,  Real>(x);
        case 5:  return boost::math::cardinal_b_spline<5,  Real>(x);
        case 6:  return boost::math::cardinal_b_spline<6,  Real>(x);
        case 7:  return boost::math::cardinal_b_spline<7,  Real>(x);
        case 8:  return boost::math::cardinal_b_spline<8,  Real>(x);
        case 9:  return boost::math::cardinal_b_spline<9,  Real>(x);
        case 10: return boost::math::cardinal_b_spline<10, Real>(x);
        case 11: return boost::math::cardinal_b_spline<11, Real>(x);
        case 12: return boost::math::cardinal_b_spline<12, Real>(x);
        case 13: return boost::math::cardinal_b_spline<12, Real>(x);
        case 14: return boost::math::cardinal_b_spline<12, Real>(x);
        case 15: return boost::math::cardinal_b_spline<12, Real>(x);
        case 16: return boost::math::cardinal_b_spline<12, Real>(x);
        case 17: return boost::math::cardinal_b_spline<12, Real>(x);
        case 18: return boost::math::cardinal_b_spline<12, Real>(x);
        case 19: return boost::math::cardinal_b_spline<12, Real>(x);
        case 20: return boost::math::cardinal_b_spline<12, Real>(x);
        default: return std::numeric_limits<Real>::quiet_NaN();
    }
}

template <typename Real>
static inline Real dispatch_cardinal_b_spline_prime(unsigned int n, Real x) {
    switch (n) {
        case 3:  return boost::math::cardinal_b_spline_prime<3,  Real>(x);
        case 4:  return boost::math::cardinal_b_spline_prime<4,  Real>(x);
        case 5:  return boost::math::cardinal_b_spline_prime<5,  Real>(x);
        case 6:  return boost::math::cardinal_b_spline_prime<6,  Real>(x);
        case 7:  return boost::math::cardinal_b_spline_prime<7,  Real>(x);
        case 8:  return boost::math::cardinal_b_spline_prime<8,  Real>(x);
        case 9:  return boost::math::cardinal_b_spline_prime<9,  Real>(x);
        case 10: return boost::math::cardinal_b_spline_prime<10, Real>(x);
        case 11: return boost::math::cardinal_b_spline_prime<11, Real>(x);
        case 12: return boost::math::cardinal_b_spline_prime<12, Real>(x);
        case 13: return boost::math::cardinal_b_spline_prime<12, Real>(x);
        case 14: return boost::math::cardinal_b_spline_prime<12, Real>(x);
        case 15: return boost::math::cardinal_b_spline_prime<12, Real>(x);
        case 16: return boost::math::cardinal_b_spline_prime<12, Real>(x);
        case 17: return boost::math::cardinal_b_spline_prime<12, Real>(x);
        case 18: return boost::math::cardinal_b_spline_prime<12, Real>(x);
        case 19: return boost::math::cardinal_b_spline_prime<12, Real>(x);
        case 20: return boost::math::cardinal_b_spline_prime<12, Real>(x);
        default: return std::numeric_limits<Real>::quiet_NaN();
    }
}

template <typename Real>
static inline Real dispatch_cardinal_b_spline_double_prime(unsigned int n, Real x) {
    switch (n) {
        case 3:  return boost::math::cardinal_b_spline_double_prime<3,  Real>(x);
        case 4:  return boost::math::cardinal_b_spline_double_prime<4,  Real>(x);
        case 5:  return boost::math::cardinal_b_spline_double_prime<5,  Real>(x);
        case 6:  return boost::math::cardinal_b_spline_double_prime<6,  Real>(x);
        case 7:  return boost::math::cardinal_b_spline_double_prime<7,  Real>(x);
        case 8:  return boost::math::cardinal_b_spline_double_prime<8,  Real>(x);
        case 9:  return boost::math::cardinal_b_spline_double_prime<9,  Real>(x);
        case 10: return boost::math::cardinal_b_spline_double_prime<10, Real>(x);
        case 11: return boost::math::cardinal_b_spline_double_prime<11, Real>(x);
        case 12: return boost::math::cardinal_b_spline_double_prime<12, Real>(x);
        case 13: return boost::math::cardinal_b_spline_double_prime<12, Real>(x);
        case 14: return boost::math::cardinal_b_spline_double_prime<12, Real>(x);
        case 15: return boost::math::cardinal_b_spline_double_prime<12, Real>(x);
        case 16: return boost::math::cardinal_b_spline_double_prime<12, Real>(x);
        case 17: return boost::math::cardinal_b_spline_double_prime<12, Real>(x);
        case 18: return boost::math::cardinal_b_spline_double_prime<12, Real>(x);
        case 19: return boost::math::cardinal_b_spline_double_prime<12, Real>(x);
        case 20: return boost::math::cardinal_b_spline_double_prime<12, Real>(x);
        default: return std::numeric_limits<Real>::quiet_NaN();
    }
}

template <typename Real>
static inline Real dispatch_forward_cardinal_b_spline(unsigned int n, Real x) {
    switch (n) {
        case 0:  return boost::math::forward_cardinal_b_spline<0,  Real>(x);
        case 1:  return boost::math::forward_cardinal_b_spline<1,  Real>(x);
        case 2:  return boost::math::forward_cardinal_b_spline<2,  Real>(x);
        case 3:  return boost::math::forward_cardinal_b_spline<3,  Real>(x);
        case 4:  return boost::math::forward_cardinal_b_spline<4,  Real>(x);
        case 5:  return boost::math::forward_cardinal_b_spline<5,  Real>(x);
        case 6:  return boost::math::forward_cardinal_b_spline<6,  Real>(x);
        case 7:  return boost::math::forward_cardinal_b_spline<7,  Real>(x);
        case 8:  return boost::math::forward_cardinal_b_spline<8,  Real>(x);
        case 9:  return boost::math::forward_cardinal_b_spline<9,  Real>(x);
        case 10: return boost::math::forward_cardinal_b_spline<10, Real>(x);
        case 11: return boost::math::forward_cardinal_b_spline<11, Real>(x);
        case 12: return boost::math::forward_cardinal_b_spline<12, Real>(x);
        case 13: return boost::math::forward_cardinal_b_spline<12, Real>(x);
        case 14: return boost::math::forward_cardinal_b_spline<12, Real>(x);
        case 15: return boost::math::forward_cardinal_b_spline<12, Real>(x);
        case 16: return boost::math::forward_cardinal_b_spline<12, Real>(x);
        case 17: return boost::math::forward_cardinal_b_spline<12, Real>(x);
        case 18: return boost::math::forward_cardinal_b_spline<12, Real>(x);
        case 19: return boost::math::forward_cardinal_b_spline<12, Real>(x);
        case 20: return boost::math::forward_cardinal_b_spline<12, Real>(x);
        default: return std::numeric_limits<Real>::quiet_NaN();
    }
}

extern "C" {

    //
    //  Created by VT on 18.10.25.
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
// Bernoulli Numbers
double bs_bernoulli_b2n(const int n) { return bs_wrap<double>([&] { return boost::math::bernoulli_b2n<double>(n); }); }
float bs_bernoulli_b2n_f(const int n) { return bs_wrap<float>([&] { return boost::math::bernoulli_b2n<float>(n); }); }
long double bs_bernoulli_b2n_l(const int n) { return bs_wrap<long double>([&] { return boost::math::bernoulli_b2n<long double>(n); }); }

// Tangent Numbers (scalar)
double bs_tangent_t2n(const int n) { return bs_wrap<double>([&] { return boost::math::tangent_t2n<double>(n); }); }
float bs_tangent_t2n_f(const int n) { return bs_wrap<float>([&] { return boost::math::tangent_t2n<float>(n); }); }
long double bs_tangent_t2n_l(const int n) { return bs_wrap<long double>([&] { return boost::math::tangent_t2n<long double>(n); }); }

// Tangent Numbers (bulk sequence)
void bs_tangent_t2n_seq(int start_index, unsigned int count, double* out) {
    if (!out || count == 0) return;
    try { boost::math::tangent_t2n<double>(start_index, count, out); } catch (...) {}
}
void bs_tangent_t2n_seq_f(int start_index, unsigned int count, float* out) {
    if (!out || count == 0) return;
    try { boost::math::tangent_t2n<float>(start_index, count, out); } catch (...) {}
}
void bs_tangent_t2n_seq_l(int start_index, unsigned int count, long double* out) {
    if (!out || count == 0) return;
    try { boost::math::tangent_t2n<long double>(start_index, count, out); } catch (...) {}
}

// Fibonacci numbers
unsigned long long bs_fibonacci_ull(unsigned long long n) { return bs_wrap<unsigned long long>([&] { return boost::math::fibonacci<unsigned long long>(n); }); }

// prime numbers up to 10000th
unsigned int bs_prime(unsigned int n) { return bs_wrap<unsigned int>([&] { return boost::math::prime(n); }); }

// Factorials
double bs_factorial_d(unsigned int i) { return bs_wrap<double>([&] { return boost::math::factorial<double>(i); }); }
float bs_factorial_f(unsigned int i) { return bs_wrap<float>([&] { return boost::math::factorial<float>(i); }); }
long double bs_factorial_l(unsigned int i) { return bs_wrap<long double>([&] { return boost::math::factorial<long double>(i); }); }

// Pochhammer (rising factorial)
double bs_rising_factorial_d(double x, unsigned int i) { return bs_wrap<double>([&] { return boost::math::rising_factorial<double>(x, i); }); }
float bs_rising_factorial_f(float x, unsigned int i) { return bs_wrap<float>([&] { return boost::math::rising_factorial<float>(x, i); }); }
long double bs_rising_factorial_l(long double x, unsigned int i) { return bs_wrap<long double>([&] { return boost::math::rising_factorial<long double>(x, i); }); }

// Binomial coefficients
double bs_binomial_coefficient_d(unsigned int n, unsigned int k) { return bs_wrap<double>([&] { return boost::math::binomial_coefficient<double>(n, k); }); }
float bs_binomial_coefficient_f(unsigned int n, unsigned int k) { return bs_wrap<float>([&] { return boost::math::binomial_coefficient<float>(n, k); }); }
long double bs_binomial_coefficient_l(unsigned int n, unsigned int k) { return bs_wrap<long double>([&] { return boost::math::binomial_coefficient<long double>(n, k); }); }

// double factorial
double bs_double_factorial_d(unsigned int i) { return bs_wrap<double>([&] { return boost::math::double_factorial<double>(i); }); }
float bs_double_factorial_f(unsigned int i) { return bs_wrap<float>([&] { return boost::math::double_factorial<float>(i); }); }
long double bs_double_factorial_l(unsigned int i) { return bs_wrap<long double>([&] { return boost::math::double_factorial<long double>(i); }); }

// Laguerre
double bs_laguerre_d(unsigned int n, double x) { return bs_wrap<double>([&] { return boost::math::laguerre(n, x); }); }
float bs_laguerre_f(unsigned int n, float x) { return bs_wrap<float>([&] { return boost::math::laguerre(n, x); }); }
long double bs_laguerre_l(unsigned int n, long double x) { return bs_wrap<long double>([&] { return boost::math::laguerre(n, x); }); }
double bs_assoc_laguerre_d(unsigned int n, unsigned int m, double x) { return bs_wrap<double>([&] { return boost::math::laguerre(n, m, x); }); }
float bs_assoc_laguerre_f(unsigned int n, unsigned int m, float x) { return bs_wrap<float>([&] { return boost::math::laguerre(n, m, x); }); }
long double bs_assoc_laguerre_l(unsigned int n, unsigned int m, long double x) { return bs_wrap<long double>([&] { return boost::math::laguerre(n, m, x); }); }

// Chebyshev
double bs_chebyshev_T_d(unsigned int n, double x) { return bs_wrap<double>([&] { return boost::math::chebyshev_t(n, x); }); }
double bs_chebyshev_U_d(unsigned int n, double x) { return bs_wrap<double>([&] { return boost::math::chebyshev_u(n, x); }); }
float bs_chebyshev_T_f(unsigned int n, float x) { return bs_wrap<float>([&] { return boost::math::chebyshev_t(n, x); }); }
float bs_chebyshev_U_f(unsigned int n, float x) { return bs_wrap<float>([&] { return boost::math::chebyshev_u(n, x); }); }
long double bs_chebyshev_T_l(unsigned int n, long double x) { return bs_wrap<long double>([&] { return boost::math::chebyshev_t(n, x); }); }
long double bs_chebyshev_U_l(unsigned int n, long double x) { return bs_wrap<long double>([&] { return boost::math::chebyshev_u(n, x); }); }

// Chebyshev series evaluation via Clenshaw (first kind)
double bs_chebyshev_clenshaw_d(const double* c, size_t count, double x) { return bs_wrap<double>([&] { return boost::math::chebyshev_clenshaw_recurrence(c, count, x);} ); }
float bs_chebyshev_clenshaw_f(const float* c, size_t count, float x) { return bs_wrap<float>([&] { return boost::math::chebyshev_clenshaw_recurrence(c, count, x);} ); }
long double bs_chebyshev_clenshaw_l(const long double* c, size_t count, long double x) { return bs_wrap<long double>([&] { return boost::math::chebyshev_clenshaw_recurrence(c, count, x);} ); }

// Spherical harmonics
bs_complex_d bs_spherical_harmonic_d(unsigned int n, int m, double theta, double phi) {
    auto z = boost::math::spherical_harmonic(n, m, theta, phi);
    bs_complex_d r{ z.real(), z.imag() };
    return r;
}
bs_complex_f bs_spherical_harmonic_f(unsigned int n, int m, float theta, float phi) {
    auto z = boost::math::spherical_harmonic(n, m, theta, phi);
    bs_complex_f r{ static_cast<float>(z.real()), static_cast<float>(z.imag()) };
    return r;
}
bs_complex_l bs_spherical_harmonic_l(unsigned int n, int m, long double theta, long double phi) {
    auto z = boost::math::spherical_harmonic(n, m, theta, phi);
    bs_complex_l r{ static_cast<long double>(z.real()), static_cast<long double>(z.imag()) };
    return r;
}

// Cardinal B-spline wrappers
double bs_cardinal_b_spline_d(unsigned int n, double x) { return bs_wrap<double>([&]{ return dispatch_cardinal_b_spline<double>(n, x); }); }
double bs_cardinal_b_spline_prime_d(unsigned int n, double x) { return bs_wrap<double>([&]{ return dispatch_cardinal_b_spline_prime<double>(n, x); }); }
double bs_cardinal_b_spline_double_prime_d(unsigned int n, double x) { return bs_wrap<double>([&]{ return dispatch_cardinal_b_spline_double_prime<double>(n, x); }); }
double bs_forward_cardinal_b_spline_d(unsigned int n, double x) { return bs_wrap<double>([&]{ return dispatch_forward_cardinal_b_spline<double>(n, x); }); }

float bs_cardinal_b_spline_f(unsigned int n, float x) { return bs_wrap<float>([&]{ return dispatch_cardinal_b_spline<float>(n, x); }); }
float bs_cardinal_b_spline_prime_f(unsigned int n, float x) { return bs_wrap<float>([&]{ return dispatch_cardinal_b_spline_prime<float>(n, x); }); }
float bs_cardinal_b_spline_double_prime_f(unsigned int n, float x) { return bs_wrap<float>([&]{ return dispatch_cardinal_b_spline_double_prime<float>(n, x); }); }
float bs_forward_cardinal_b_spline_f(unsigned int n, float x) { return bs_wrap<float>([&]{ return dispatch_forward_cardinal_b_spline<float>(n, x); }); }

long double bs_cardinal_b_spline_l(unsigned int n, long double x) { return bs_wrap<long double>([&]{ return dispatch_cardinal_b_spline<long double>(n, x); }); }
long double bs_cardinal_b_spline_prime_l(unsigned int n, long double x) { return bs_wrap<long double>([&]{ return dispatch_cardinal_b_spline_prime<long double>(n, x); }); }
long double bs_cardinal_b_spline_double_prime_l(unsigned int n, long double x) { return bs_wrap<long double>([&]{ return dispatch_cardinal_b_spline_double_prime<long double>(n, x); }); }
long double bs_forward_cardinal_b_spline_l(unsigned int n, long double x) { return bs_wrap<long double>([&]{ return dispatch_forward_cardinal_b_spline<long double>(n, x); }); }

// Gegenbauer polynomials
double bs_gegenbauer_d(unsigned int n, double lambda, double x) { return bs_wrap<double>([&]{ return boost::math::gegenbauer(n, lambda, x); }); }
double bs_gegenbauer_prime_d(unsigned int n, double lambda, double x) { return bs_wrap<double>([&]{ return boost::math::gegenbauer_prime(n, lambda, x); }); }
double bs_gegenbauer_derivative_d(unsigned int n, double lambda, double x, unsigned int k) { return bs_wrap<double>([&]{ return boost::math::gegenbauer_derivative(n, lambda, x, k); }); }

float bs_gegenbauer_f(unsigned int n, float lambda, float x) { return bs_wrap<float>([&]{ return boost::math::gegenbauer(n, lambda, x); }); }
float bs_gegenbauer_prime_f(unsigned int n, float lambda, float x) { return bs_wrap<float>([&]{ return boost::math::gegenbauer_prime(n, lambda, x); }); }
float bs_gegenbauer_derivative_f(unsigned int n, float lambda, float x, unsigned int k) { return bs_wrap<float>([&]{ return boost::math::gegenbauer_derivative(n, lambda, x, k); }); }

long double bs_gegenbauer_l(unsigned int n, long double lambda, long double x) { return bs_wrap<long double>([&]{ return boost::math::gegenbauer(n, lambda, x); }); }
long double bs_gegenbauer_prime_l(unsigned int n, long double lambda, long double x) { return bs_wrap<long double>([&]{ return boost::math::gegenbauer_prime(n, lambda, x); }); }
long double bs_gegenbauer_derivative_l(unsigned int n, long double lambda, long double x, unsigned int k) { return bs_wrap<long double>([&]{ return boost::math::gegenbauer_derivative(n, lambda, x, k); }); }

} // extern "C"

