//
//  Created by VT on 11.10.25.
//  Copyright © 2025 Volker Thieme.
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

#include "CBoostBridge.h"

#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/trigamma.hpp>
#include <boost/math/special_functions/polygamma.hpp>
#include <boost/math/special_functions/zeta.hpp>
#include <boost/math/special_functions/owens_t.hpp>
#include <boost/math/special_functions/expint.hpp>
#include <boost/math/special_functions/expm1.hpp>
#include <boost/math/special_functions/log1p.hpp>
#include <boost/math/special_functions/powm1.hpp>
#include <boost/math/special_functions/cbrt.hpp>
#include <boost/math/special_functions/sinc.hpp>
#include <boost/math/special_functions/sinhc.hpp>
#include <boost/math/special_functions/sin_pi.hpp>
#include <boost/math/special_functions/cos_pi.hpp>
#include <boost/math/special_functions/airy.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/legendre.hpp>
#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>
#include <boost/math/special_functions/ellint_3.hpp>
#include <boost/math/special_functions/ellint_rc.hpp>
#include <boost/math/special_functions/ellint_rf.hpp>
#include <boost/math/special_functions/ellint_rd.hpp>
#include <boost/math/special_functions/ellint_rj.hpp>
#include <boost/math/special_functions/ellint_rg.hpp>
#include <boost/math/special_functions/lambert_w.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/hypergeometric_1F0.hpp>
#include <boost/math/special_functions/hypergeometric_0F1.hpp>
#include <boost/math/special_functions/hypergeometric_2F0.hpp>
#include <boost/math/special_functions/hypergeometric_1F1.hpp>
#include <boost/math/special_functions/hypergeometric_pFq.hpp>
#include <boost/math/special_functions/bernoulli.hpp>
#include <boost/math/special_functions/fibonacci.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/prime.hpp>
#include <boost/math/special_functions/chebyshev.hpp>
#include <boost/math/special_functions/cardinal_b_spline.hpp>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <boost/math/special_functions/gegenbauer.hpp>
#include <vector>
#include <limits>
#include <stdexcept>
#include <algorithm> // for std::min, std::copy
#include <complex>   // for std::complex
#include <boost/math/complex/atan.hpp>
#include <boost/math/complex/acos.hpp>
#include <boost/math/complex/acosh.hpp>
#include <boost/math/complex/asin.hpp>
#include <boost/math/complex/asinh.hpp>
#include <boost/math/complex/atanh.hpp>
#include <boost/math/complex/fabs.hpp>
#include <boost/math/special_functions/sqrt1pm1.hpp>
#include <boost/math/special_functions/powm1.hpp>
#include <boost/math/special_functions/hypot.hpp>
#include <boost/math/special_functions/bessel_prime.hpp>



// Exception-safe wrapper: translate Boost throws to numeric sentinels.
// - Overflow -> +infinity
// - Domain/other errors -> quiet NaN
template <class T, class F>
inline T bs_wrap(F&& f) noexcept {
    try {
        return static_cast<T>(f());
    } catch (const std::overflow_error&) {
        return std::numeric_limits<T>::infinity();
    } catch (const std::domain_error&) {
        return std::numeric_limits<T>::quiet_NaN();
    } catch (...) {
        return std::numeric_limits<T>::quiet_NaN();
    }
}

// Complex wrapper with same policy: on overflow => {+Inf,+Inf}, on other => {NaN,NaN}.
template <class T, class F>
static inline std::complex<T> bs_wrap_complex(F&& f) noexcept {
    try {
        return f();
    } catch (const std::overflow_error&) {
        const T inf = std::numeric_limits<T>::infinity();
        return std::complex<T>(inf, inf);
    } catch (const std::domain_error&) {
        const T nan = std::numeric_limits<T>::quiet_NaN();
        return std::complex<T>(nan, nan);
    } catch (...) {
        const T nan = std::numeric_limits<T>::quiet_NaN();
        return std::complex<T>(nan, nan);
    }
}

extern "C" {

// Boost constants (double)
double bs_const_e(void)                 { return bs_wrap<double>([] { return boost::math::constants::e<double>(); }); }
double bs_const_pi(void)                { return bs_wrap<double>([] { return boost::math::constants::pi<double>(); }); }
double bs_const_two_pi(void)            { return bs_wrap<double>([] { return boost::math::constants::two_pi<double>(); }); }
double bs_const_half_pi(void)           { return bs_wrap<double>([] { return boost::math::constants::half_pi<double>(); }); }
double bs_const_third_pi(void)          { return bs_wrap<double>([] { return boost::math::constants::third_pi<double>(); }); }
double bs_const_two_thirds_pi(void)     { return bs_wrap<double>([] { return boost::math::constants::two_thirds_pi<double>(); }); }
double bs_const_quarter_pi(void)        { return bs_wrap<double>([] { return boost::math::constants::quarter_pi<double>(); }); }
double bs_const_three_quarters_pi(void) { return bs_wrap<double>([] { return boost::math::constants::three_quarters_pi<double>(); }); }
double bs_const_sixth_pi(void)          { return bs_wrap<double>([] { return boost::math::constants::sixth_pi<double>(); }); }
double bs_const_pi_sqr(void)            { return bs_wrap<double>([] { return boost::math::constants::pi_sqr<double>(); }); }

double bs_const_root_two(void)          { return bs_wrap<double>([] { return boost::math::constants::root_two<double>(); }); }
double bs_const_root_three(void)        { return bs_wrap<double>([] { return boost::math::constants::root_three<double>(); }); }
double bs_const_root_pi(void)           { return bs_wrap<double>([] { return boost::math::constants::root_pi<double>(); }); }

double bs_const_one_div_pi(void)        { return bs_wrap<double>([] { return boost::math::constants::one_div_pi<double>(); }); }
double bs_const_one_div_two_pi(void)    { return bs_wrap<double>([] { return boost::math::constants::one_div_two_pi<double>(); }); }
double bs_const_one_div_root_pi(void)   { return bs_wrap<double>([] { return boost::math::constants::one_div_root_pi<double>(); }); }
double bs_const_two_div_pi(void)        { return bs_wrap<double>([] { return boost::math::constants::two_div_pi<double>(); }); }
double bs_const_two_div_root_pi(void)   { return bs_wrap<double>([] { return boost::math::constants::two_div_root_pi<double>(); }); }

double bs_const_ln_two(void)            { return bs_wrap<double>([] { return boost::math::constants::ln_two<double>(); }); }
double bs_const_ln_ten(void)            { return bs_wrap<double>([] { return boost::math::constants::ln_ten<double>(); }); }
double bs_const_ln_ln_two(void)         { return bs_wrap<double>([] { return boost::math::constants::ln_ln_two<double>(); }); }

double bs_const_euler(void)             { return bs_wrap<double>([] { return boost::math::constants::euler<double>(); }); }
double bs_const_catalan(void)           { return bs_wrap<double>([] { return boost::math::constants::catalan<double>(); }); }
double bs_const_zeta_three(void)        { return bs_wrap<double>([] { return boost::math::constants::zeta_three<double>(); }); }
double bs_const_phi(void)               { return bs_wrap<double>([] { return boost::math::constants::phi<double>(); }); }

// Boost constants (float)
float bs_const_e_f(void)                 { return bs_wrap<float>([] { return boost::math::constants::e<float>(); }); }
float bs_const_pi_f(void)                { return bs_wrap<float>([] { return boost::math::constants::pi<float>(); }); }
float bs_const_two_pi_f(void)            { return bs_wrap<float>([] { return boost::math::constants::two_pi<float>(); }); }
float bs_const_half_pi_f(void)           { return bs_wrap<float>([] { return boost::math::constants::half_pi<float>(); }); }
float bs_const_third_pi_f(void)          { return bs_wrap<float>([] { return boost::math::constants::third_pi<float>(); }); }
float bs_const_two_thirds_pi_f(void)     { return bs_wrap<float>([] { return boost::math::constants::two_thirds_pi<float>(); }); }
float bs_const_quarter_pi_f(void)        { return bs_wrap<float>([] { return boost::math::constants::quarter_pi<float>(); }); }
float bs_const_three_quarters_pi_f(void) { return bs_wrap<float>([] { return boost::math::constants::three_quarters_pi<float>(); }); }
float bs_const_sixth_pi_f(void)          { return bs_wrap<float>([] { return boost::math::constants::sixth_pi<float>(); }); }
float bs_const_pi_sqr_f(void)            { return bs_wrap<float>([] { return boost::math::constants::pi_sqr<float>(); }); }

float bs_const_root_two_f(void)          { return bs_wrap<float>([] { return boost::math::constants::root_two<float>(); }); }
float bs_const_root_three_f(void)        { return bs_wrap<float>([] { return boost::math::constants::root_three<float>(); }); }
float bs_const_root_pi_f(void)           { return bs_wrap<float>([] { return boost::math::constants::root_pi<float>(); }); }

float bs_const_one_div_pi_f(void)        { return bs_wrap<float>([] { return boost::math::constants::one_div_pi<float>(); }); }
float bs_const_one_div_two_pi_f(void)    { return bs_wrap<float>([] { return boost::math::constants::one_div_two_pi<float>(); }); }
float bs_const_one_div_root_pi_f(void)   { return bs_wrap<float>([] { return boost::math::constants::one_div_root_pi<float>(); }); }
float bs_const_two_div_pi_f(void)        { return bs_wrap<float>([] { return boost::math::constants::two_div_pi<float>(); }); }
float bs_const_two_div_root_pi_f(void)   { return bs_wrap<float>([] { return boost::math::constants::two_div_root_pi<float>(); }); }

float bs_const_ln_two_f(void)            { return bs_wrap<float>([] { return boost::math::constants::ln_two<float>(); }); }
float bs_const_ln_ten_f(void)            { return bs_wrap<float>([] { return boost::math::constants::ln_ten<float>(); }); }
float bs_const_ln_ln_two_f(void)         { return bs_wrap<float>([] { return boost::math::constants::ln_ln_two<float>(); }); }

float bs_const_euler_f(void)             { return bs_wrap<float>([] { return boost::math::constants::euler<float>(); }); }
float bs_const_catalan_f(void)           { return bs_wrap<float>([] { return boost::math::constants::catalan<float>(); }); }
float bs_const_zeta_three_f(void)        { return bs_wrap<float>([] { return boost::math::constants::zeta_three<float>(); }); }
float bs_const_phi_f(void)               { return bs_wrap<float>([] { return boost::math::constants::phi<float>(); }); }

// Boost constants (long double)
long double bs_const_e_l(void)                 { return bs_wrap<long double>([] { return boost::math::constants::e<long double>(); }); }
long double bs_const_pi_l(void)                { return bs_wrap<long double>([] { return boost::math::constants::pi<long double>(); }); }
long double bs_const_two_pi_l(void)            { return bs_wrap<long double>([] { return boost::math::constants::two_pi<long double>(); }); }
long double bs_const_half_pi_l(void)           { return bs_wrap<long double>([] { return boost::math::constants::half_pi<long double>(); }); }
long double bs_const_third_pi_l(void)          { return bs_wrap<long double>([] { return boost::math::constants::third_pi<long double>(); }); }
long double bs_const_two_thirds_pi_l(void)     { return bs_wrap<long double>([] { return boost::math::constants::two_thirds_pi<long double>(); }); }
long double bs_const_quarter_pi_l(void)        { return bs_wrap<long double>([] { return boost::math::constants::quarter_pi<long double>(); }); }
long double bs_const_three_quarters_pi_l(void) { return bs_wrap<long double>([] { return boost::math::constants::three_quarters_pi<long double>(); }); }
long double bs_const_sixth_pi_l(void)          { return bs_wrap<long double>([] { return boost::math::constants::sixth_pi<long double>(); }); }
long double bs_const_pi_sqr_l(void)            { return bs_wrap<long double>([] { return boost::math::constants::pi_sqr<long double>(); }); }

long double bs_const_root_two_l(void)          { return bs_wrap<long double>([] { return boost::math::constants::root_two<long double>(); }); }
long double bs_const_root_three_l(void)        { return bs_wrap<long double>([] { return boost::math::constants::root_three<long double>(); }); }
long double bs_const_root_pi_l(void)           { return bs_wrap<long double>([] { return boost::math::constants::root_pi<long double>(); }); }

long double bs_const_one_div_pi_l(void)        { return bs_wrap<long double>([] { return boost::math::constants::one_div_pi<long double>(); }); }
long double bs_const_one_div_two_pi_l(void)    { return bs_wrap<long double>([] { return boost::math::constants::one_div_two_pi<long double>(); }); }
long double bs_const_one_div_root_pi_l(void)   { return bs_wrap<long double>([] { return boost::math::constants::one_div_root_pi<long double>(); }); }
long double bs_const_two_div_pi_l(void)        { return bs_wrap<long double>([] { return boost::math::constants::two_div_pi<long double>(); }); }
long double bs_const_two_div_root_pi_l(void)   { return bs_wrap<long double>([] { return boost::math::constants::two_div_root_pi<long double>(); }); }

long double bs_const_ln_two_l(void)            { return bs_wrap<long double>([] { return boost::math::constants::ln_two<long double>(); }); }
long double bs_const_ln_ten_l(void)            { return bs_wrap<long double>([] { return boost::math::constants::ln_ten<long double>(); }); }
long double bs_const_ln_ln_two_l(void)         { return bs_wrap<long double>([] { return boost::math::constants::ln_ln_two<long double>(); }); }

long double bs_const_euler_l(void)             { return bs_wrap<long double>([] { return boost::math::constants::euler<long double>(); }); }
long double bs_const_catalan_l(void)           { return bs_wrap<long double>([] { return boost::math::constants::catalan<long double>(); }); }
long double bs_const_zeta_three_l(void)        { return bs_wrap<long double>([] { return boost::math::constants::zeta_three<long double>(); }); }
long double bs_const_phi_l(void)               { return bs_wrap<long double>([] { return boost::math::constants::phi<long double>(); }); }


// MARK: Complex numbers

// Converters for C POD <-> std::complex
static inline std::complex<double> to_std(bs_complex_d z) noexcept { return { z.re, z.im }; }
static inline bs_complex_d from_std(std::complex<double> z) noexcept { return { static_cast<double>(z.real()), static_cast<double>(z.imag()) }; }

static inline std::complex<float> to_std(bs_complex_f z) noexcept { return { z.re, z.im }; }
static inline bs_complex_f from_std(std::complex<float> z) noexcept { return { static_cast<float>(z.real()), static_cast<float>(z.imag()) }; }

static inline std::complex<long double> to_std(bs_complex_l z) noexcept { return { z.re, z.im }; }
static inline bs_complex_l from_std(std::complex<long double> z) noexcept { return { static_cast<long double>(z.real()), static_cast<long double>(z.imag()) }; }

// Elementary arithmetic
bs_complex_d bs_cadd(bs_complex_d a, bs_complex_d b) {
    auto r = bs_wrap_complex<double>([&]{ return to_std(a) + to_std(b); });
    return from_std(r);
}
bs_complex_d bs_csub(bs_complex_d a, bs_complex_d b) {
    auto r = bs_wrap_complex<double>([&]{ return to_std(a) - to_std(b); });
    return from_std(r);
}
bs_complex_d bs_cmul(bs_complex_d a, bs_complex_d b) {
    auto r = bs_wrap_complex<double>([&]{ return to_std(a) * to_std(b); });
    return from_std(r);
}
bs_complex_d bs_cdiv(bs_complex_d a, bs_complex_d b) {
    auto r = bs_wrap_complex<double>([&]{ return to_std(a) / to_std(b); });
    return from_std(r);
}

bs_complex_f bs_cadd_f(bs_complex_f a, bs_complex_f b) {
    auto r = bs_wrap_complex<float>([&]{ return to_std(a) + to_std(b); });
    return from_std(r);
}
bs_complex_f bs_csub_f(bs_complex_f a, bs_complex_f b) {
    auto r = bs_wrap_complex<float>([&]{ return to_std(a) - to_std(b); });
    return from_std(r);
}
bs_complex_f bs_cmul_f(bs_complex_f a, bs_complex_f b) {
    auto r = bs_wrap_complex<float>([&]{ return to_std(a) * to_std(b); });
    return from_std(r);
}
bs_complex_f bs_cdiv_f(bs_complex_f a, bs_complex_f b) {
    auto r = bs_wrap_complex<float>([&]{ return to_std(a) / to_std(b); });
    return from_std(r);
}

bs_complex_l bs_cadd_l(bs_complex_l a, bs_complex_l b) {
    auto r = bs_wrap_complex<long double>([&]{ return to_std(a) + to_std(b); });
    return from_std(r);
}
bs_complex_l bs_csub_l(bs_complex_l a, bs_complex_l b) {
    auto r = bs_wrap_complex<long double>([&]{ return to_std(a) - to_std(b); });
    return from_std(r);
}
bs_complex_l bs_cmul_l(bs_complex_l a, bs_complex_l b) {
    auto r = bs_wrap_complex<long double>([&]{ return to_std(a) * to_std(b); });
    return from_std(r);
}
bs_complex_l bs_cdiv_l(bs_complex_l a, bs_complex_l b) {
    auto r = bs_wrap_complex<long double>([&]{ return to_std(a) / to_std(b); });
    return from_std(r);
}

// Elementary functions: prefer std::complex; Boost.Math provides some complex overloads too.
bs_complex_d bs_cexp(bs_complex_d z) {
    auto r = bs_wrap_complex<double>([&]{ return std::exp(to_std(z)); });
    return from_std(r);
}
bs_complex_d bs_clog(bs_complex_d z) {
    auto r = bs_wrap_complex<double>([&]{ return std::log(to_std(z)); });
    return from_std(r);
}
bs_complex_d bs_csqrt(bs_complex_d z) {
    auto r = bs_wrap_complex<double>([&]{ return std::sqrt(to_std(z)); });
    return from_std(r);
}
bs_complex_d bs_csin(bs_complex_d z) {
    auto r = bs_wrap_complex<double>([&]{ return std::sin(to_std(z)); });
    return from_std(r);
}
bs_complex_d bs_ccos(bs_complex_d z) {
    auto r = bs_wrap_complex<double>([&]{ return std::cos(to_std(z)); });
    return from_std(r);
}

bs_complex_d bs_ctan(bs_complex_d z) {
    auto r = bs_wrap_complex<double>([&]{ return std::tan(to_std(z)); });
    return from_std(r);
}

bs_complex_d bs_csinh(bs_complex_d z) {
    auto r = bs_wrap_complex<double>([&]{ return std::sinh(to_std(z)); });
    return from_std(r);
}
bs_complex_d bs_ccosh(bs_complex_d z) {
    auto r = bs_wrap_complex<double>([&]{ return std::cosh(to_std(z)); });
    return from_std(r);
}
bs_complex_d bs_ctanh(bs_complex_d z) {
    auto r = bs_wrap_complex<double>([&]{ return std::tanh(to_std(z)); });
    return from_std(r);
}

bs_complex_d bs_catan(bs_complex_d z) {
    auto r = bs_wrap_complex<double>([&]{ return std::atan(to_std(z)); });
    return from_std(r);
}


bs_complex_f bs_cexp_f(bs_complex_f z) {
    auto r = bs_wrap_complex<float>([&]{ return std::exp(to_std(z)); });
    return from_std(r);
}
bs_complex_f bs_clog_f(bs_complex_f z) {
    auto r = bs_wrap_complex<float>([&]{ return std::log(to_std(z)); });
    return from_std(r);
}
bs_complex_f bs_csqrt_f(bs_complex_f z) {
    auto r = bs_wrap_complex<float>([&]{ return std::sqrt(to_std(z)); });
    return from_std(r);
}
bs_complex_f bs_csin_f(bs_complex_f z) {
    auto r = bs_wrap_complex<float>([&]{ return std::sin(to_std(z)); });
    return from_std(r);
}
bs_complex_f bs_ccos_f(bs_complex_f z) {
    auto r = bs_wrap_complex<float>([&]{ return std::cos(to_std(z)); });
    return from_std(r);
}

bs_complex_f bs_ctan_f(bs_complex_f z) {
    auto r = bs_wrap_complex<float>([&]{ return std::tan(to_std(z)); });
    return from_std(r);
}
bs_complex_f bs_csinh_f(bs_complex_f z)  {
    auto r = bs_wrap_complex<float>([&]{ return std::sinh(to_std(z)); });
    return from_std(r);
}

bs_complex_f bs_ccosh_f(bs_complex_f z)  {
    auto r = bs_wrap_complex<float>([&]{ return std::cosh(to_std(z)); });
    return from_std(r);
}

bs_complex_f bs_ctanh_f(bs_complex_f z)  {
    auto r = bs_wrap_complex<float>([&]{ return std::tanh(to_std(z)); });
    return from_std(r);
}
bs_complex_f bs_catan_f(bs_complex_f z)  {
    auto r = bs_wrap_complex<float>([&]{ return std::atan(to_std(z)); });
    return from_std(r);
}

bs_complex_l bs_cexp_l(bs_complex_l z) {
    auto r = bs_wrap_complex<long double>([&]{ return std::exp(to_std(z)); });
    return from_std(r);
}
bs_complex_l bs_clog_l(bs_complex_l z) {
    auto r = bs_wrap_complex<long double>([&]{ return std::log(to_std(z)); });
    return from_std(r);
}
bs_complex_l bs_csqrt_l(bs_complex_l z) {
    auto r = bs_wrap_complex<long double>([&]{ return std::sqrt(to_std(z)); });
    return from_std(r);
}
bs_complex_l bs_csin_l(bs_complex_l z) {
    auto r = bs_wrap_complex<long double>([&]{ return std::sin(to_std(z)); });
    return from_std(r);
}
bs_complex_l bs_ccos_l(bs_complex_l z) {
    auto r = bs_wrap_complex<long double>([&]{ return std::cos(to_std(z)); });
    return from_std(r);
}

bs_complex_l bs_ctan_l(bs_complex_l z) {
    auto r = bs_wrap_complex<long double>([&]{ return std::tan(to_std(z)); });
    return from_std(r);
}

bs_complex_l bs_csinh_l(bs_complex_l z) {
    auto r = bs_wrap_complex<long double>([&]{ return std::sinh(to_std(z)); });
    return from_std(r);
}
bs_complex_l bs_ccosh_l(bs_complex_l z) {
    auto r = bs_wrap_complex<long double>([&]{ return std::cosh(to_std(z)); });
    return from_std(r);
}
bs_complex_l bs_ctanh_l(bs_complex_l z) {
    auto r = bs_wrap_complex<long double>([&]{ return std::tanh(to_std(z)); });
    return from_std(r);
}

bs_complex_l bs_catan_l(bs_complex_l z) {
    auto r = bs_wrap_complex<long double>([&]{ return std::atan(to_std(z)); });
    return from_std(r);
}

// Gamma / Error
double bs_tgamma(double x)              { return bs_wrap<double>([&] { return boost::math::tgamma(x); }); }
double bs_lgamma(double x)              { return bs_wrap<double>([&] { return boost::math::lgamma(x); }); }
double bs_erf(double x)                 { return bs_wrap<double>([&] { return boost::math::erf(x); }); }
double bs_erfc(double x)                { return bs_wrap<double>([&] { return boost::math::erfc(x); }); }
double bs_erf_inv(double p)             { return bs_wrap<double>([&] { return boost::math::erf_inv(p); }); }
double bs_erfc_inv(double p)            { return bs_wrap<double>([&] { return boost::math::erfc_inv(p); }); }

float bs_tgamma_f(float x)              { return bs_wrap<float>([&] { return boost::math::tgamma(x); }); }
float bs_lgamma_f(float x)              { return bs_wrap<float>([&] { return boost::math::lgamma(x); }); }
float bs_erf_f(float x)                 { return bs_wrap<float>([&] { return boost::math::erf(x); }); }
float bs_erfc_f(float x)                { return bs_wrap<float>([&] { return boost::math::erfc(x); }); }
float bs_erf_inv_f(float p)            { return bs_wrap<float>([&] { return boost::math::erf_inv(p); }); }
float bs_erfc_inv_f(float p)           { return bs_wrap<float>([&] { return boost::math::erfc_inv(p); }); }

long double bs_tgamma_l(long double x)  { return bs_wrap<long double>([&] { return boost::math::tgamma(x); }); }
long double bs_lgamma_l(long double x)  { return bs_wrap<long double>([&] { return boost::math::lgamma(x); }); }
long double bs_erf_l(long double x)     { return bs_wrap<long double>([&] { return boost::math::erf(x); }); }
long double bs_erfc_l(long double x)    { return bs_wrap<long double>([&] { return boost::math::erfc(x); }); }
long double bs_erf_inv_l(long double p) { return bs_wrap<long double>([&] { return boost::math::erf_inv(p); }); }
long double bs_erfc_inv_l(long double p){ return bs_wrap<long double>([&] { return boost::math::erfc_inv(p); }); }



// Incomplete gamma (lower/upper, regularized, and inverses)
double bs_tgamma_lower(double a, double x)      { return bs_wrap<double>([&] { return boost::math::tgamma_lower(a, x); }); }
double bs_tgamma_upper(double a, double x)      { return bs_wrap<double>([&] { return boost::math::tgamma(a, x); }); }
double bs_gamma_p(double a, double x)           { return bs_wrap<double>([&] { return boost::math::gamma_p(a, x); }); }
double bs_gamma_q(double a, double x)           { return bs_wrap<double>([&] { return boost::math::gamma_q(a, x); }); }
double bs_gamma_p_inv(double a, double p)       { return bs_wrap<double>([&] { return boost::math::gamma_p_inv(a, p); }); }
double bs_gamma_q_inv(double a, double q)       { return bs_wrap<double>([&] { return boost::math::gamma_q_inv(a, q); }); }

// Derivatives of regularized incomplete gamma functions
double bs_gamma_p_derivative(double a, double x) { return bs_wrap<double>([&] { return boost::math::gamma_p_derivative(a, x); }); }

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

// Beta family
double bs_beta(double a, double b)                     { return bs_wrap<double>([&] { return boost::math::beta(a, b); }); }
double bs_fullBeta(double a, double b, double x)       { return bs_wrap<double>([&] { return boost::math::beta(a, b, x); }); }
double bs_ibeta(double a, double b, double x)          { return bs_wrap<double>([&] { return boost::math::ibeta(a, b, x); }); }
double bs_ibetac(double a, double b, double x)         { return bs_wrap<double>([&] { return boost::math::ibetac(a, b, x); }); }
double bs_ibeta_inv(double a, double b, double p)      { return bs_wrap<double>([&] { return boost::math::ibeta_inv(a, b, p); }); }
double bs_ibetac_inv(double a, double b, double p)     { return bs_wrap<double>([&] { return boost::math::ibetac_inv(a, b, p); }); }
double bs_ibeta_inva(double b, double x, double p)     { return bs_wrap<double>([&] { return boost::math::ibeta_inva(b, x, p); }); }
double bs_ibeta_invb(double a, double x, double p)     { return bs_wrap<double>([&] { return boost::math::ibeta_invb(a, x, p); }); }
double bs_ibeta_derivative(double a, double b, double x) { return bs_wrap<double>([&] { return boost::math::ibeta_derivative(a, b, x); }); }

float bs_beta_f(float a, float b)                      { return bs_wrap<float>([&] { return boost::math::beta(a, b); }); }
float bs_fullBeta_f(float a, float b, float x)         { return bs_wrap<float>([&] { return boost::math::beta(a, b, x); }); }
float bs_ibeta_f(float a, float b, float x)            { return bs_wrap<float>([&] { return boost::math::ibeta(a, b, x); }); }
float bs_ibetac_f(float a, float b, float x)           { return bs_wrap<float>([&] { return boost::math::ibetac(a, b, x); }); }
float bs_ibeta_inv_f(float a, float b, float p)        { return bs_wrap<float>([&] { return boost::math::ibeta_inv(a, b, p); }); }
float bs_ibetac_inv_f(float a, float b, float p)       { return bs_wrap<float>([&] { return boost::math::ibetac_inv(a, b, p); }); }
float bs_ibeta_inva_f(float b, float x, float p)       { return bs_wrap<float>([&] { return boost::math::ibeta_inva(b, x, p); }); }
float bs_ibeta_invb_f(float a, float x, float p)       { return bs_wrap<float>([&] { return boost::math::ibeta_invb(a, x, p); }); }
float bs_ibeta_derivative_f(float a, float b, float x) { return bs_wrap<float>([&] { return boost::math::ibeta_derivative(a, b, x); }); }

long double bs_beta_l(long double a, long double b)                      { return bs_wrap<long double>([&] { return boost::math::beta(a, b); }); }
long double bs_fullBeta_l(long double a, long double b, long double x)   { return bs_wrap<long double>([&] { return boost::math::beta(a, b, x); }); }
long double bs_ibeta_l(long double a, long double b, long double x)      { return bs_wrap<long double>([&] { return boost::math::ibeta(a, b, x); }); }
long double bs_ibetac_l(long double a, long double b, long double x)     { return bs_wrap<long double>([&] { return boost::math::ibetac(a, b, x); }); }
long double bs_ibeta_inv_l(long double a, long double b, long double p)  { return bs_wrap<long double>([&] { return boost::math::ibeta_inv(a, b, p); }); }
long double bs_ibetac_inv_l(long double a, long double b, long double p) { return bs_wrap<long double>([&] { return boost::math::ibetac_inv(a, b, p); }); }
long double bs_ibeta_inva_l(long double b, long double x, long double p) { return bs_wrap<long double>([&] { return boost::math::ibeta_inva(b, x, p); }); }
long double bs_ibeta_invb_l(long double a, long double x, long double p) { return bs_wrap<long double>([&] { return boost::math::ibeta_invb(a, x, p); }); }
long double bs_ibeta_derivative_l(long double a, long double b, long double x) { return bs_wrap<long double>([&] { return boost::math::ibeta_derivative(a, b, x); }); }

// Digamma / Polygamma / Zeta
double bs_digamma(double x)            { return bs_wrap<double>([&] { return boost::math::digamma(x); }); }
double bs_trigamma(double x)           { return bs_wrap<double>([&] { return boost::math::trigamma(x); }); }
double bs_polygamma(int n, double x)   { return bs_wrap<double>([&] { return boost::math::polygamma(n, x); }); }
double bs_riemann_zeta(double x)       { return bs_wrap<double>([&] { return boost::math::zeta(x); }); }

float bs_digamma_f(float x)            { return bs_wrap<float>([&] { return boost::math::digamma(x); }); }
float bs_trigamma_f(float x)           { return bs_wrap<float>([&] { return boost::math::trigamma(x); }); }
float bs_polygamma_f(int n, float x)   { return bs_wrap<float>([&] { return boost::math::polygamma(n, x); }); }
float bs_riemann_zeta_f(float x)       { return bs_wrap<float>([&] { return boost::math::zeta(x); }); }

long double bs_digamma_l(long double x)          { return bs_wrap<long double>([&] { return boost::math::digamma(x); }); }
long double bs_trigamma_l(long double x)         { return bs_wrap<long double>([&] { return boost::math::trigamma(x); }); }
long double bs_polygamma_l(int n, long double x) { return bs_wrap<long double>([&] { return boost::math::polygamma(n, x); }); }
long double bs_riemann_zeta_l(long double x)     { return bs_wrap<long double>([&] { return boost::math::zeta(x); }); }

// Owen's T
double bs_owens_t(double h, double a)           { return bs_wrap<double>([&] { return boost::math::owens_t(h, a); }); }
float bs_owens_t_f(float h, float a)            { return bs_wrap<float>([&] { return boost::math::owens_t(h, a); }); }
long double bs_owens_t_l(long double h, long double a) { return bs_wrap<long double>([&] { return boost::math::owens_t(h, a); }); }

// Exponential integrals and related
double bs_expint_Ei(double x)        { return bs_wrap<double>([&] { return boost::math::expint(x); }); }
double bs_expint_En(int n, double x) { return bs_wrap<double>([&] { return boost::math::expint(n, x); }); }
double bs_expm1(double x)            { return bs_wrap<double>([&] { return boost::math::expm1(x); }); }
double bs_log1p(double x)            { return bs_wrap<double>([&] { return boost::math::log1p(x); }); }
double bs_log1pmx(double x)          { return bs_wrap<double>([&] { return boost::math::log1pmx(x); }); }
double bs_powm1(double x, double y)  { return bs_wrap<double>([&] { return boost::math::powm1(x, y); }); }
double bs_cbrt(double x)             { return bs_wrap<double>([&] { return boost::math::cbrt(x); }); }
double bs_sqrt1pm1(double x)         { return bs_wrap<double>([&] { return boost::math::sqrt1pm1(x); }); }
double bs_hypot(double x, double y)  { return bs_wrap<double>([&] { return boost::math::hypot(x, y); }); }

float bs_expint_Ei_f(float x)        { return bs_wrap<float>([&] { return boost::math::expint(x); }); }
float bs_expint_En_f(int n, float x) { return bs_wrap<float>([&] { return boost::math::expint(n, x); }); }
float bs_expm1_f(float x)            { return bs_wrap<float>([&] { return boost::math::expm1(x); }); }
float bs_log1p_f(float x)            { return bs_wrap<float>([&] { return boost::math::log1p(x); }); }
float bs_log1pmx_f(float x)          { return bs_wrap<float>([&] { return boost::math::log1pmx(x); }); }
float bs_powm1_f(float x, float y)   { return bs_wrap<float>([&] { return boost::math::powm1(x, y); }); }
float bs_cbrt_f(float x)             { return bs_wrap<float>([&] { return boost::math::cbrt(x); }); }
float bs_sqrt1pm1_f(float x)         { return bs_wrap<float>([&] { return boost::math::sqrt1pm1(x); }); }
float bs_hypot_f(float x, float y)  { return bs_wrap<float>([&] { return boost::math::hypot(x, y); }); }

long double bs_expint_Ei_l(long double x)          { return bs_wrap<long double>([&] { return boost::math::expint(x); }); }
long double bs_expint_En_l(int n, long double x)   { return bs_wrap<long double>([&] { return boost::math::expint(n, x); }); }
long double bs_expm1_l(long double x)              { return bs_wrap<long double>([&] { return boost::math::expm1(x); }); }
long double bs_log1p_l(long double x)              { return bs_wrap<long double>([&] { return boost::math::log1p(x); }); }
long double bs_log1pmx_l(long double x)            { return bs_wrap<long double>([&] { return boost::math::log1pmx(x); }); }
long double bs_powm1_l(long double x, long double y) { return bs_wrap<long double>([&] { return boost::math::powm1(x, y); }); }
long double bs_cbrt_l(long double x)               { return bs_wrap<long double>([&] { return boost::math::cbrt(x); }); }
long double bs_sqrt1pm1_l(long double x)         { return bs_wrap<long double>([&] { return boost::math::sqrt1pm1(x); }); }
long double bs_hypot_l(long double x, long double y)  { return bs_wrap<long double>([&] { return boost::math::hypot(x, y); }); }

// Trig helpers
double bs_sin_pi(double x)             { return bs_wrap<double>([&] { return boost::math::sin_pi(x); }); }
double bs_cos_pi(double x)             { return bs_wrap<double>([&] { return boost::math::cos_pi(x); }); }

float bs_sin_pi_f(float x)             { return bs_wrap<float>([&] { return boost::math::sin_pi(x); }); }
float bs_cos_pi_f(float x)             { return bs_wrap<float>([&] { return boost::math::cos_pi(x); }); }

long double bs_sin_pi_l(long double x) { return bs_wrap<long double>([&] { return boost::math::sin_pi(x); }); }
long double bs_cos_pi_l(long double x) { return bs_wrap<long double>([&] { return boost::math::cos_pi(x); }); }

// Airy
double bs_airy_ai(double x)            { return bs_wrap<double>([&] { return boost::math::airy_ai(x); }); }
double bs_airy_bi(double x)            { return bs_wrap<double>([&] { return boost::math::airy_bi(x); }); }
double bs_airy_ai_prime(double x)      { return bs_wrap<double>([&] { return boost::math::airy_ai_prime(x); }); }
double bs_airy_bi_prime(double x)      { return bs_wrap<double>([&] { return boost::math::airy_bi_prime(x); }); }

float bs_airy_ai_f(float x)            { return bs_wrap<float>([&] { return boost::math::airy_ai(x); }); }
float bs_airy_bi_f(float x)            { return bs_wrap<float>([&] { return boost::math::airy_bi(x); }); }
float bs_airy_ai_prime_f(float x)      { return bs_wrap<float>([&] { return boost::math::airy_ai_prime(x); }); }
float bs_airy_bi_prime_f(float x)      { return bs_wrap<float>([&] { return boost::math::airy_bi_prime(x); }); }

long double bs_airy_ai_l(long double x)       { return bs_wrap<long double>([&] { return boost::math::airy_ai(x); }); }
long double bs_airy_bi_l(long double x)       { return bs_wrap<long double>([&] { return boost::math::airy_bi(x); }); }
long double bs_airy_ai_prime_l(long double x) { return bs_wrap<long double>([&] { return boost::math::airy_ai_prime(x); }); }
long double bs_airy_bi_prime_l(long double x) { return bs_wrap<long double>([&] { return boost::math::airy_bi_prime(x); }); }

// Bessel (cylindrical, real)
double bs_cyl_bessel_j(double v, double x) { return bs_wrap<double>([&] { return boost::math::cyl_bessel_j(v, x); }); }
double bs_cyl_neumann(double v, double x)  { return bs_wrap<double>([&] { return boost::math::cyl_neumann(v, x); }); }
double bs_cyl_bessel_i(double v, double x) { return bs_wrap<double>([&] { return boost::math::cyl_bessel_i(v, x); }); }
double bs_cyl_bessel_k(double v, double x) { return bs_wrap<double>([&] { return boost::math::cyl_bessel_k(v, x); }); }
double bs_cyl_bessel_j_zero ( double v, int m) { return bs_wrap<double>([&] { return boost::math::cyl_bessel_j_zero<double>(v, m); }); }
void bs_cyl_bessel_j_zeros(double v, int start_index, unsigned int number_of_zeros, double* out) {
    if (!out || number_of_zeros == 0) return;
        // Fill the caller-provided buffer; the function returns the output iterator, which we can ignore.
    boost::math::cyl_bessel_j_zero(v, start_index, number_of_zeros, out);
}

float bs_cyl_bessel_j_f(float v, float x) { return bs_wrap<float>([&] { return boost::math::cyl_bessel_j(v, x); }); }
float bs_cyl_neumann_f(float v, float x)  { return bs_wrap<float>([&] { return boost::math::cyl_neumann(v, x); }); }
float bs_cyl_bessel_i_f(float v, float x) { return bs_wrap<float>([&] { return boost::math::cyl_bessel_i(v, x); }); }
float bs_cyl_bessel_k_f(float v, float x) { return bs_wrap<float>([&] { return boost::math::cyl_bessel_k(v, x); }); }
float bs_cyl_bessel_j_zero_f ( float v, int m) { return bs_wrap<float>([&] { return boost::math::cyl_bessel_j_zero<float>(v, m); }); }
void bs_cyl_bessel_j_zeros_f(float v, int start_index, unsigned int number_of_zeros, float* out) {
    if (!out || number_of_zeros == 0) return;
        // Fill the caller-provided buffer; the function returns the output iterator, which we can ignore.
    boost::math::cyl_bessel_j_zero(v, start_index, number_of_zeros, out);
}

long double bs_cyl_bessel_j_l(long double v, long double x) { return bs_wrap<long double>([&] { return boost::math::cyl_bessel_j(v, x); }); }
long double bs_cyl_neumann_l(long double v, long double x)  { return bs_wrap<long double>([&] { return boost::math::cyl_neumann(v, x); }); }
long double bs_cyl_bessel_i_l(long double v, long double x) { return bs_wrap<long double>([&] { return boost::math::cyl_bessel_i(v, x); }); }
long double bs_cyl_bessel_k_l(long double v, long double x) { return bs_wrap<long double>([&] { return boost::math::cyl_bessel_k(v, x); }); }
long double bs_cyl_bessel_j_zero_l( long double v, int m) { return bs_wrap<long double>([&] { return boost::math::cyl_bessel_j_zero<long double>(v, m); }); }
void bs_cyl_bessel_j_zeros_l(long double v, int start_index, unsigned int number_of_zeros, long double* out) {
    if (!out || number_of_zeros == 0) return;
        // Fill the caller-provided buffer; the function returns the output iterator, which we can ignore.
    boost::math::cyl_bessel_j_zero(v, start_index, number_of_zeros, out);
}
// Spherical Bessel/Neumann bridges (no helpers)
double bs_sph_bessel(unsigned int n, double x) {
    return boost::math::sph_bessel(n, x);
}
float bs_sph_bessel_f(unsigned int n, float x) {
    return boost::math::sph_bessel(n, x);
}
long double bs_sph_bessel_l(unsigned int n, long double x) {
    return boost::math::sph_bessel(n, x);
}

double bs_sph_neumann(unsigned int n, double x) {
    return boost::math::sph_neumann(n, x);
}
float bs_sph_neumann_f(unsigned int n, float x) {
    return boost::math::sph_neumann(n, x);
}
long double bs_sph_neumann_l(unsigned int n, long double x) {
    return boost::math::sph_neumann(n, x);
}

double bs_cyl_bessel_j_prime(double v, double x) {
    return boost::math::cyl_bessel_j_prime(v, x);
}

double bs_cyl_bessel_i_prime(double v, double x) {
    return boost::math::cyl_bessel_i_prime(v, x);
}
double bs_cyl_bessel_k_prime(double v, double x) {
    return boost::math::cyl_bessel_k_prime(v, x);
}
double bs_sph_bessel_prime(unsigned int n, double x) {
    return boost::math::sph_bessel_prime(n, x);
}
double bs_sph_neumann_prime(unsigned int n, double x) {
    return boost::math::sph_neumann_prime(n, x);
}
float bs_cyl_bessel_j_prime_f(float v, float x) {
    return boost::math::cyl_bessel_j_prime(v, x);
}
float bs_cyl_bessel_i_prime_f(float v, float x) {
    return boost::math::cyl_bessel_i_prime(v, x);
}
float bs_cyl_bessel_k_prime_f(float v, float x) {
    return boost::math::cyl_bessel_k_prime(v, x);
}
float bs_sph_bessel_prime_f(unsigned int n, float x) {
    return boost::math::sph_bessel_prime(n, x);
}
float bs_sph_neumann_prime_f(unsigned int n, float x) {
    return boost::math::sph_neumann_prime(n, x);
}
long double bs_cyl_bessel_j_prime_l(long double v, long double x) {
    return boost::math::cyl_bessel_j_prime(v, x);
}

long double bs_cyl_bessel_i_prime_l(long double v, long double x) {
    return boost::math::cyl_bessel_i_prime(v, x);
}
long double bs_cyl_bessel_k_prime_l(long double v, long double x) {
    return boost::math::cyl_bessel_k_prime(v, x);
}
long double bs_sph_bessel_prime_l(unsigned int n, long double x) {
    return boost::math::sph_bessel_prime(n, x);
}
long double bs_sph_neumann_prime_l(unsigned int n, long double x) {
    return boost::math::sph_neumann_prime(n, x);
}

                                           
// Legendre
double bs_legendre_p(int n, double x)                 { return bs_wrap<double>([&] { return boost::math::legendre_p(n, x); }); }
double bs_assoc_legendre_p(int n, int m, double x)    { return bs_wrap<double>([&] { return boost::math::legendre_p(n, m, x); }); }
double bs_legendre_p_prime(int n, double x)           { return bs_wrap<double>([&] { return boost::math::legendre_p_prime(n, x); }); }
void bs_legendre_p_zeros(int l, double* out) {
    if (!out) return;
    try {
        auto v = boost::math::legendre_p_zeros<double>(l);
        const std::size_t n = std::min<std::size_t>(static_cast<std::size_t>(l), v.size());
        std::copy(v.begin(), v.begin() + n, out);
    } catch (...) {
        // Leave output as-is on error
    }
}

float bs_legendre_p_f(int n, float x)                 { return bs_wrap<float>([&] { return boost::math::legendre_p(n, x); }); }
float bs_assoc_legendre_p_f(int n, int m, float x)    { return bs_wrap<float>([&] { return boost::math::legendre_p(n, m, x); }); }
float bs_legendre_p_prime_f(int n, float x)           { return bs_wrap<float>([&] { return boost::math::legendre_p_prime(n, x); }); }
void bs_legendre_p_zeros_f(int l, float* out) {
    if (!out) return;
    try {
        auto v = boost::math::legendre_p_zeros<float>(l);
        const std::size_t n = std::min<std::size_t>(static_cast<std::size_t>(l), v.size());
        std::copy(v.begin(), v.begin() + n, out);
    } catch (...) {
        // Leave output as-is on error
    }
}

long double bs_legendre_p_l(int n, long double x)              { return bs_wrap<long double>([&] { return boost::math::legendre_p(n, x); }); }
long double bs_assoc_legendre_p_l(int n, int m, long double x) { return bs_wrap<long double>([&] { return boost::math::legendre_p(n, m, x); }); }
long double bs_legendre_p_prime_l(int n, long double x)        { return bs_wrap<long double>([&] { return boost::math::legendre_p_prime(n, x); }); }
void bs_legendre_p_zeros_l(int l, long double* out) {
    if (!out) return;
    try {
        auto v = boost::math::legendre_p_zeros<long double>(l);
        const std::size_t n = std::min<std::size_t>(static_cast<std::size_t>(l), v.size());
        std::copy(v.begin(), v.begin() + n, out);
    } catch (...) {
        // Leave output as-is on error
    }
}

// Elliptic integrals (Legendre forms)
double bs_ellint_1_complete(double k)                 { return bs_wrap<double>([&] { return boost::math::ellint_1(k); }); }
double bs_ellint_1(double k, double phi)              { return bs_wrap<double>([&] { return boost::math::ellint_1(k, phi); }); }
double bs_ellint_2_complete(double k)                 { return bs_wrap<double>([&] { return boost::math::ellint_2(k); }); }
double bs_ellint_2(double k, double phi)              { return bs_wrap<double>([&] { return boost::math::ellint_2(k, phi); }); }
double bs_ellint_3(double k, double nu, double phi)   { return bs_wrap<double>([&] { return boost::math::ellint_3(k, nu, phi); }); } // incomplete elliptic integral of the third kind ∏(k, nu, phi)
double bs_ellint_3_complete(double k, double nu)      { return bs_wrap<double>([&] { return boost::math::ellint_3(k, nu); }); }


float bs_ellint_1_complete_f(float k)                 { return bs_wrap<float>([&] { return boost::math::ellint_1(k); }); }
float bs_ellint_1_f(float k, float phi)               { return bs_wrap<float>([&] { return boost::math::ellint_1(k, phi); }); }
float bs_ellint_2_complete_f(float k)                 { return bs_wrap<float>([&] { return boost::math::ellint_2(k); }); }
float bs_ellint_2_f(float k, float phi)               { return bs_wrap<float>([&] { return boost::math::ellint_2(k, phi); }); }
float bs_ellint_3_f(float k, float nu, float phi)     { return bs_wrap<float>([&] { return boost::math::ellint_3(k, nu, phi); }); } // incomplete elliptic integral of the third kind ∏(k, nu, phi)
float bs_ellint_3_complete_f(float k, float nu)      { return bs_wrap<float>([&] { return boost::math::ellint_3(k, nu); }); }

long double bs_ellint_1_complete_l(long double k)                 { return bs_wrap<long double>([&] { return boost::math::ellint_1(k); }); }
long double bs_ellint_1_l(long double k, long double phi)         { return bs_wrap<long double>([&] { return boost::math::ellint_1(k, phi); }); }
long double bs_ellint_2_complete_l(long double k)                 { return bs_wrap<long double>([&] { return boost::math::ellint_2(k); }); }
long double bs_ellint_2_l(long double k, long double phi)         { return bs_wrap<long double>([&] { return boost::math::ellint_2(k, phi); }); }
long double bs_ellint_3_l(long double k, long double nu,long double phi) { return bs_wrap<long double>([&] { return boost::math::ellint_3(k, nu, phi); }); } // incomplete elliptic integral of the third kind ∏(k, nu, phi)
long double bs_ellint_3_complete_l(long double k, long double nu)      { return bs_wrap<long double>([&] { return boost::math::ellint_3(k, nu); }); }

// Elliptic integrals (Carlson symmetric forms)
double bs_ellint_rc(double x, double y)               { return bs_wrap<double>([&] { return boost::math::ellint_rc(x, y); }); }
double bs_ellint_rf(double x, double y, double z)     { return bs_wrap<double>([&] { return boost::math::ellint_rf(x, y, z); }); }
double bs_ellint_rd(double x, double y, double z)     { return bs_wrap<double>([&] { return boost::math::ellint_rd(x, y, z); }); }
double bs_ellint_rj(double x, double y, double z, double p) { return bs_wrap<double>([&] { return boost::math::ellint_rj(x, y, z, p); }); }
double bs_ellint_rg(double x, double y, double z)     { return bs_wrap<double>([&] { return boost::math::ellint_rg(x, y, z); }); }

float bs_ellint_rc_f(float x, float y)                { return bs_wrap<float>([&] { return boost::math::ellint_rc(x, y); }); }
float bs_ellint_rf_f(float x, float y, float z)       { return bs_wrap<float>([&] { return boost::math::ellint_rf(x, y, z); }); }
float bs_ellint_rd_f(float x, float y, float z)       { return bs_wrap<float>([&] { return boost::math::ellint_rd(x, y, z); }); }
float bs_ellint_rj_f(float x, float y, float z, float p) { return bs_wrap<float>([&] { return boost::math::ellint_rj(x, y, z, p); }); }
float bs_ellint_rg_f(float x, float y, float z)       { return bs_wrap<float>([&] { return boost::math::ellint_rg(x, y, z); }); }

long double bs_ellint_rc_l(long double x, long double y)                { return bs_wrap<long double>([&] { return boost::math::ellint_rc(x, y); }); }
long double bs_ellint_rf_l(long double x, long double y, long double z) { return bs_wrap<long double>([&] { return boost::math::ellint_rf(x, y, z); }); }
long double bs_ellint_rd_l(long double x, long double y, long double z) { return bs_wrap<long double>([&] { return boost::math::ellint_rd(x, y, z); }); }
long double bs_ellint_rj_l(long double x, long double y, long double z, long double p) { return bs_wrap<long double>([&] { return boost::math::ellint_rj(x, y, z, p); }); }
long double bs_ellint_rg_l(long double x, long double y, long double z) { return bs_wrap<long double>([&] { return boost::math::ellint_rg(x, y, z); }); }

// Lambert W (real branches)
double bs_lambert_w0(double x)        { return bs_wrap<double>([&] { return boost::math::lambert_w0(x); }); }
double bs_lambert_wm1(double x)       { return bs_wrap<double>([&] { return boost::math::lambert_wm1(x); }); }

float bs_lambert_w0_f(float x)        { return bs_wrap<float>([&] { return boost::math::lambert_w0(x); }); }
float bs_lambert_wm1_f(float x)       { return bs_wrap<float>([&] { return boost::math::lambert_wm1(x); }); }

long double bs_lambert_w0_l(long double x)  { return bs_wrap<long double>([&] { return boost::math::lambert_w0(x); }); }
long double bs_lambert_wm1_l(long double x) { return bs_wrap<long double>([&] { return boost::math::lambert_wm1(x); }); }

// Γ(a) / Γ(b)
double bs_tgamma_ratio(double a, double b) { return bs_wrap<double>([&] { return boost::math::tgamma_ratio(a, b); }); }
float  bs_tgamma_ratio_f(float a, float b) { return bs_wrap<float>([&] { return boost::math::tgamma_ratio(a, b); }); }
long double bs_tgamma_ratio_l(long double a, long double b) { return bs_wrap<long double>([&] { return boost::math::tgamma_ratio(a, b); }); }

// Γ(a) / Γ(a + delta)
double bs_tgamma_delta_ratio(double a, double delta)  { return bs_wrap<double>([&] { return boost::math::tgamma_delta_ratio(a, delta); }); }
float  bs_tgamma_delta_ratio_f(float a, float delta)  { return bs_wrap<float>([&] { return boost::math::tgamma_delta_ratio(a, delta); }); }
long double bs_tgamma_delta_ratio_l(long double a, long double delta)  { return bs_wrap<long double>([&] { return boost::math::tgamma_delta_ratio(a, delta); }); }

// Hypergeometric functions
double bs_hypergeometric_1F0(double a, double z) { return bs_wrap<double>([&] { return boost::math::hypergeometric_1F0(a, z); }); }
float  bs_hypergeometric_1F0_f(float a, float z) { return bs_wrap<float>([&] { return boost::math::hypergeometric_1F0(a, z); }); }
long double bs_hypergeometric_1F0_l(long double a, long double z) { return bs_wrap<long double>([&] { return boost::math::hypergeometric_1F0(a, z); }); }

double bs_hypergeometric_0F1(double b, double z) { return bs_wrap<double>([&] { return boost::math::hypergeometric_0F1(b, z); }); }
float  bs_hypergeometric_0F1_f(float b, float z) { return bs_wrap<float>([&] { return boost::math::hypergeometric_0F1(b, z); }); }
long double bs_hypergeometric_0F1_l(long double b, long double z) { return bs_wrap<long double>([&] { return boost::math::hypergeometric_0F1(b, z); }); }

double bs_hypergeometric_2F0(double a, double b, double z) { return bs_wrap<double>([&] { return boost::math::hypergeometric_2F0(a, b, z); }); }
float  bs_hypergeometric_2F0_f(float a, float b, float z)  { return bs_wrap<float>([&] { return boost::math::hypergeometric_2F0(a, b, z); }); }
long double bs_hypergeometric_2F0_l(long double a, long double b, long double z) { return bs_wrap<long double>([&] { return boost::math::hypergeometric_2F0(a, b, z); }); }

double bs_hypergeometric_1F1(double a, double b, double z) { return bs_wrap<double>([&] { return boost::math::hypergeometric_1F1(a, b, z); }); }
float  bs_hypergeometric_1F1_f(float a, float b, float z)  { return bs_wrap<float>([&] { return boost::math::hypergeometric_1F1(a, b, z); }); }
long double bs_hypergeometric_1F1_l(long double a, long double b, long double z) { return bs_wrap<long double>([&] { return boost::math::hypergeometric_1F1(a, b, z); }); }

double bs_hypergeometric_pFq(const double* a, size_t p, const double* b, size_t q, double z) {
    return bs_wrap<double>([&] {
        std::vector<double> va(a, a + p);
        std::vector<double> vb(b, b + q);
        return boost::math::hypergeometric_pFq(va, vb, z);
    });
}
float bs_hypergeometric_pFq_f(const float* a, size_t p, const float* b, size_t q, float z) {
    return bs_wrap<float>([&] {
        std::vector<float> va(a, a + p);
        std::vector<float> vb(b, b + q);
        return boost::math::hypergeometric_pFq(va, vb, z);
    });
}
long double bs_hypergeometric_pFq_l(const long double* a, size_t p, const long double* b, size_t q, long double z) {
    return bs_wrap<long double>([&] {
        std::vector<long double> va(a, a + p);
        std::vector<long double> vb(b, b + q);
        return boost::math::hypergeometric_pFq(va, vb, z);
    });
}

// Bernoulli Numbers
double bs_bernoulli_b2n(const int n) {
    return bs_wrap<double>([&] { return boost::math::bernoulli_b2n<double>(n); });
}
float bs_bernoulli_b2n_f(const int n) {
    return bs_wrap<float>([&] { return boost::math::bernoulli_b2n<float>(n); });
}
long double bs_bernoulli_b2n_l(const int n) {
    return bs_wrap<long double>([&] { return boost::math::bernoulli_b2n<long double>(n); });
}

// Tangent Numbers (scalar)
double bs_tangent_t2n(const int n) {
    return bs_wrap<double>([&] { return boost::math::tangent_t2n<double>(n); });
}
float bs_tangent_t2n_f(const int n) {
    return bs_wrap<float>([&] { return boost::math::tangent_t2n<float>(n); });
}
long double bs_tangent_t2n_l(const int n) {
    return bs_wrap<long double>([&] { return boost::math::tangent_t2n<long double>(n); });
}

// Tangent Numbers (bulk sequence)
void bs_tangent_t2n_seq(int start_index, unsigned int count, double* out) {
    if (!out || count == 0) return;
    try {
        boost::math::tangent_t2n<double>(start_index, count, out);
    } catch (...) {
        // Leave output as-is on error
    }
}
void bs_tangent_t2n_seq_f(int start_index, unsigned int count, float* out) {
    if (!out || count == 0) return;
    try {
        boost::math::tangent_t2n<float>(start_index, count, out);
    } catch (...) {
        // Leave output as-is on error
    }
}
void bs_tangent_t2n_seq_l(int start_index, unsigned int count, long double* out) {
    if (!out || count == 0) return;
    try {
        boost::math::tangent_t2n<long double>(start_index, count, out);
    } catch (...) {
        // Leave output as-is on error
    }
}

// Fibonacci numbers
unsigned long long bs_fibonacci_ull(unsigned long long n) {
    return bs_wrap<unsigned long long>([&] { return boost::math::fibonacci<unsigned long long>(n); });
}

unsigned int bs_prime(unsigned int n) {
    return bs_wrap<unsigned int>([&] { return boost::math::prime(n); });
}

// Factorial
double bs_factorial(unsigned int i) {
    return bs_wrap<double>([&] { return boost::math::factorial<double>(i); });
}
float bs_factorial_f(unsigned int i) {
    return bs_wrap<float>([&] { return boost::math::factorial<float>(i); });
}
long double bs_factorial_l(unsigned int i) {
    return bs_wrap<long double>([&] { return boost::math::factorial<long double>(i); });
}

double bs_rising_factorial(double x, unsigned int i) {
    return bs_wrap<double>([&] { return boost::math::rising_factorial<double>(x, i); });
}
float bs_rising_factorial_f(float x, unsigned int i) {
    return bs_wrap<float>([&] { return boost::math::rising_factorial<float>(x, i); });
}
long double bs_rising_factorial_l(long double x, unsigned int i) {
    return bs_wrap<long double>([&] { return boost::math::rising_factorial<long double>(x, i); });
}

double bs_binomial_coefficient(unsigned int n, unsigned int k) {
    return bs_wrap<double>([&] { return boost::math::binomial_coefficient<double>(n, k); });
}
float bs_binomial_coefficient_f(unsigned int n, unsigned int k) {
    return bs_wrap<float>([&] { return boost::math::binomial_coefficient<float>(n, k); });
}
long double bs_binomial_coefficient_l(unsigned int n, unsigned int k) {
    return bs_wrap<long double>([&] { return boost::math::binomial_coefficient<long double>(n, k); });
}

double bs_double_factorial(unsigned int i) {
    return bs_wrap<double>([&] { return boost::math::double_factorial<double>(i); });
}
float bs_double_factorial_f(unsigned int i) {
    return bs_wrap<float>([&] { return boost::math::double_factorial<float>(i); });
}
long double bs_double_factorial_l(unsigned int i) {
    return bs_wrap<long double>([&] { return boost::math::double_factorial<long double>(i); });
}

// Laguerre
double bs_laguerre(unsigned int n, double x)    { return bs_wrap<double>([&] { return boost::math::laguerre(n, x); }); }
float bs_laguerre_f(unsigned int n, float x)    { return bs_wrap<float>([&] { return boost::math::laguerre(n, x); }); }
long double bs_laguerre_l(unsigned int n, long double x) { return bs_wrap<long double>([&] { return boost::math::laguerre(n, x); }); }

double bs_assoc_laguerre(unsigned int n, unsigned int m, double x) { return bs_wrap<double>([&] { return boost::math::laguerre(n, m, x); }); }
float bs_assoc_laguerre_f(unsigned int n, unsigned int m, float x) { return bs_wrap<float>([&] { return boost::math::laguerre(n, m, x); }); }
long double bs_assoc_laguerre_l(unsigned int n, unsigned int m, long double x) { return bs_wrap<long double>([&] { return boost::math::laguerre(n, m, x); }); }

// Chebyshev
double bs_chebyshev_T(unsigned int n, double x) { return bs_wrap<double>([&] { return boost::math::chebyshev_t(n, x); }); }
double bs_chebyshev_U(unsigned int n, double x) { return bs_wrap<double>([&] { return boost::math::chebyshev_u(n, x); }); }
float bs_chebyshev_T_f(unsigned int n, float x) { return bs_wrap<float>([&] { return boost::math::chebyshev_t(n, x); }); }
float bs_chebyshev_U_f(unsigned int n, float x) { return bs_wrap<float>([&] { return boost::math::chebyshev_u(n, x); }); }
long double bs_chebyshev_T_l(unsigned int n, long double x) { return bs_wrap<long double>([&] { return boost::math::chebyshev_t(n, x); }); }
long double bs_chebyshev_U_l(unsigned int n, long double x) { return bs_wrap<long double>([&] { return boost::math::chebyshev_u(n, x); }); }

// Chebyshev series evaluation via Clenshaw (first kind):
// Evaluates S(x) = c0/2 + sum_{k=1}^{count-1} c[k] * T_k(x)
// where c points to c_0..c_{count-1}.

double bs_chebyshev_clenshaw(const double* c, size_t count, double x) {
    return bs_wrap<double>([&] { return boost::math::chebyshev_clenshaw_recurrence(c, count, x);} );
}

float bs_chebyshev_clenshaw_f(const float* c, size_t count, float x) {
    return bs_wrap<float>([&] { return boost::math::chebyshev_clenshaw_recurrence(c, count, x);} );
}

long double bs_chebyshev_clenshaw_l(const long double* c, size_t count, long double x) {
    return bs_wrap<long double>([&] { return boost::math::chebyshev_clenshaw_recurrence(c, count, x);} );
}

bs_complex_d bs_spherical_harmonic(unsigned int n, int m, double theta, double phi) {
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


} //extern C

// MARK: Cardinal B-spline runtime dispatch (n in [0, 12])

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
/*        case 0:  return boost::math::cardinal_b_spline_prime<0,  Real>(x);
        case 1:  return boost::math::cardinal_b_spline_prime<1,  Real>(x);
        case 2:  return boost::math::cardinal_b_spline_prime<2,  Real>(x); */
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
/*        case 0:  return boost::math::cardinal_b_spline_double_prime<0,  Real>(x);
        case 1:  return boost::math::cardinal_b_spline_double_prime<1,  Real>(x);
        case 2:  return boost::math::cardinal_b_spline_double_prime<2,  Real>(x); */
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

// double
double bs_cardinal_b_spline(unsigned int n, double x) {
    return bs_wrap<double>([&]{ return dispatch_cardinal_b_spline<double>(n, x); });
}
double bs_cardinal_b_spline_prime(unsigned int n, double x) {
    return bs_wrap<double>([&]{ return dispatch_cardinal_b_spline_prime<double>(n, x); });
}
double bs_cardinal_b_spline_double_prime(unsigned int n, double x) {
    return bs_wrap<double>([&]{ return dispatch_cardinal_b_spline_double_prime<double>(n, x); });
}
double bs_forward_cardinal_b_spline(unsigned int n, double x) {
    return bs_wrap<double>([&]{ return dispatch_forward_cardinal_b_spline<double>(n, x); });
}

// float
float bs_cardinal_b_spline_f(unsigned int n, float x) {
    return bs_wrap<float>([&]{ return dispatch_cardinal_b_spline<float>(n, x); });
}
float bs_cardinal_b_spline_prime_f(unsigned int n, float x) {
    return bs_wrap<float>([&]{ return dispatch_cardinal_b_spline_prime<float>(n, x); });
}
float bs_cardinal_b_spline_double_prime_f(unsigned int n, float x) {
    return bs_wrap<float>([&]{ return dispatch_cardinal_b_spline_double_prime<float>(n, x); });
}
float bs_forward_cardinal_b_spline_f(unsigned int n, float x) {
    return bs_wrap<float>([&]{ return dispatch_forward_cardinal_b_spline<float>(n, x); });
}

// long double
long double bs_cardinal_b_spline_l(unsigned int n, long double x) {
    return bs_wrap<long double>([&]{ return dispatch_cardinal_b_spline<long double>(n, x); });
}
long double bs_cardinal_b_spline_prime_l(unsigned int n, long double x) {
    return bs_wrap<long double>([&]{ return dispatch_cardinal_b_spline_prime<long double>(n, x); });
}
long double bs_cardinal_b_spline_double_prime_l(unsigned int n, long double x) {
    return bs_wrap<long double>([&]{ return dispatch_cardinal_b_spline_double_prime<long double>(n, x); });
}
long double bs_forward_cardinal_b_spline_l(unsigned int n, long double x) {
    return bs_wrap<long double>([&]{ return dispatch_forward_cardinal_b_spline<long double>(n, x); });
}

double      bs_gegenbauer           (unsigned int n, double lambda, double x) {
    return bs_wrap<double>([&]{ return boost::math::gegenbauer(n, lambda, x) ; });
}
double      bs_gegenbauer_prime     (unsigned int n, double lambda, double x) {
    return bs_wrap<double>([&]{ return boost::math::gegenbauer_prime(n, lambda, x) ; });
}
double      bs_gegenbauer_derivative(unsigned int n, double lambda, double x, unsigned int k) {
    return bs_wrap<double>([&]{ return boost::math::gegenbauer_derivative(n, lambda, x, k) ; });
}


float      bs_gegenbauer_f           (unsigned int n, float lambda, float x) {
    return bs_wrap<float>([&]{ return boost::math::gegenbauer(n, lambda, x) ; });
}
float      bs_gegenbauer_prime_f     (unsigned int n, float lambda, float x) {
    return bs_wrap<float>([&]{ return boost::math::gegenbauer_prime(n, lambda, x) ; });
}
float      bs_gegenbauer_derivative_f(unsigned int n, float lambda, float x, unsigned int k) {
    return bs_wrap<float>([&]{ return boost::math::gegenbauer_derivative(n, lambda, x, k) ; });
}

long double      bs_gegenbauer_l           (unsigned int n, long double lambda, long double x) {
    return bs_wrap<long double>([&]{ return boost::math::gegenbauer(n, lambda, x) ; });
}
long double      bs_gegenbauer_prime_l     (unsigned int n, long double lambda, long double x) {
    return bs_wrap<long double>([&]{ return boost::math::gegenbauer_prime(n, lambda, x) ; });
}
long double      bs_gegenbauer_derivative_l(unsigned int n, long double lambda, long double x, unsigned int k) {
    return bs_wrap<long double>([&]{ return boost::math::gegenbauer_derivative(n, lambda, x, k) ; });
}


} // extern "C"
