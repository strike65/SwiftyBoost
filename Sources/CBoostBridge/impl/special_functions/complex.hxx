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
// Complex number helpers and elementary functions
#include <complex>
#include "../../internal/bs_internal.hpp"
#include "../include/special/complex.h"

// Converters for C POD <-> std::complex
static inline std::complex<double> to_std(complex_d z) noexcept { return { z.re, z.im }; }
static inline complex_d from_std(std::complex<double> z) noexcept { return { static_cast<double>(z.real()), static_cast<double>(z.imag()) }; }

static inline std::complex<float> to_std(complex_f z) noexcept { return { z.re, z.im }; }
static inline complex_f from_std(std::complex<float> z) noexcept { return { static_cast<float>(z.real()), static_cast<float>(z.imag()) }; }

static inline std::complex<long double> to_std(complex_l z) noexcept { return { z.re, z.im }; }
static inline complex_l from_std(std::complex<long double> z) noexcept { return { static_cast<long double>(z.real()), static_cast<long double>(z.imag()) }; }

#ifdef __cplusplus
extern "C" {
#endif

// Elementary arithmetic
complex_d bs_cadd_d(complex_d a, complex_d b) {
    auto r = bs_wrap_complex<double>([&]{ return to_std(a) + to_std(b); });
    return from_std(r);
}
complex_d bs_csub_d(complex_d a, complex_d b) {
    auto r = bs_wrap_complex<double>([&]{ return to_std(a) - to_std(b); });
    return from_std(r);
}
complex_d bs_cmul_d(complex_d a, complex_d b) {
    auto r = bs_wrap_complex<double>([&]{ return to_std(a) * to_std(b); });
    return from_std(r);
}
complex_d bs_cdiv_d(complex_d a, complex_d b) {
    auto r = bs_wrap_complex<double>([&]{ return to_std(a) / to_std(b); });
    return from_std(r);
}

complex_f bs_cadd_f(complex_f a, complex_f b) {
    auto r = bs_wrap_complex<float>([&]{ return to_std(a) + to_std(b); });
    return from_std(r);
}
complex_f bs_csub_f(complex_f a, complex_f b) {
    auto r = bs_wrap_complex<float>([&]{ return to_std(a) - to_std(b); });
    return from_std(r);
}
complex_f bs_cmul_f(complex_f a, complex_f b) {
    auto r = bs_wrap_complex<float>([&]{ return to_std(a) * to_std(b); });
    return from_std(r);
}
complex_f bs_cdiv_f(complex_f a, complex_f b) {
    auto r = bs_wrap_complex<float>([&]{ return to_std(a) / to_std(b); });
    return from_std(r);
}

complex_l bs_cadd_l(complex_l a, complex_l b) {
    auto r = bs_wrap_complex<long double>([&]{ return to_std(a) + to_std(b); });
    return from_std(r);
}
complex_l bs_csub_l(complex_l a, complex_l b) {
    auto r = bs_wrap_complex<long double>([&]{ return to_std(a) - to_std(b); });
    return from_std(r);
}
complex_l bs_cmul_l(complex_l a, complex_l b) {
    auto r = bs_wrap_complex<long double>([&]{ return to_std(a) * to_std(b); });
    return from_std(r);
}
complex_l bs_cdiv_l(complex_l a, complex_l b) {
    auto r = bs_wrap_complex<long double>([&]{ return to_std(a) / to_std(b); });
    return from_std(r);
}

// Elementary functions
complex_d bs_cexp_d(complex_d z) { auto r = bs_wrap_complex<double>([&]{ return std::exp(to_std(z)); }); return from_std(r); }
complex_d bs_clog_d(complex_d z) { auto r = bs_wrap_complex<double>([&]{ return std::log(to_std(z)); }); return from_std(r); }
complex_d bs_csqrt_d(complex_d z){ auto r = bs_wrap_complex<double>([&]{ return std::sqrt(to_std(z)); }); return from_std(r); }
complex_d bs_csin_d(complex_d z) { auto r = bs_wrap_complex<double>([&]{ return std::sin(to_std(z)); }); return from_std(r); }
complex_d bs_ccos_d(complex_d z) { auto r = bs_wrap_complex<double>([&]{ return std::cos(to_std(z)); }); return from_std(r); }
complex_d bs_ctan_d(complex_d z) { auto r = bs_wrap_complex<double>([&]{ return std::tan(to_std(z)); }); return from_std(r); }
complex_d bs_csinh_d(complex_d z){ auto r = bs_wrap_complex<double>([&]{ return std::sinh(to_std(z)); }); return from_std(r); }
complex_d bs_ccosh_d(complex_d z){ auto r = bs_wrap_complex<double>([&]{ return std::cosh(to_std(z)); }); return from_std(r); }
complex_d bs_ctanh_d(complex_d z){ auto r = bs_wrap_complex<double>([&]{ return std::tanh(to_std(z)); }); return from_std(r); }
complex_d bs_catan_d(complex_d z){ auto r = bs_wrap_complex<double>([&]{ return std::atan(to_std(z)); }); return from_std(r); }

complex_f bs_cexp_f(complex_f z) { auto r = bs_wrap_complex<float>([&]{ return std::exp(to_std(z)); }); return from_std(r); }
complex_f bs_clog_f(complex_f z) { auto r = bs_wrap_complex<float>([&]{ return std::log(to_std(z)); }); return from_std(r); }
complex_f bs_csqrt_f(complex_f z){ auto r = bs_wrap_complex<float>([&]{ return std::sqrt(to_std(z)); }); return from_std(r); }
complex_f bs_csin_f(complex_f z) { auto r = bs_wrap_complex<float>([&]{ return std::sin(to_std(z)); }); return from_std(r); }
complex_f bs_ccos_f(complex_f z) { auto r = bs_wrap_complex<float>([&]{ return std::cos(to_std(z)); }); return from_std(r); }
complex_f bs_ctan_f(complex_f z) { auto r = bs_wrap_complex<float>([&]{ return std::tan(to_std(z)); }); return from_std(r); }
complex_f bs_csinh_f(complex_f z){ auto r = bs_wrap_complex<float>([&]{ return std::sinh(to_std(z)); }); return from_std(r); }
complex_f bs_ccosh_f(complex_f z){ auto r = bs_wrap_complex<float>([&]{ return std::cosh(to_std(z)); }); return from_std(r); }
complex_f bs_ctanh_f(complex_f z){ auto r = bs_wrap_complex<float>([&]{ return std::tanh(to_std(z)); }); return from_std(r); }
complex_f bs_catan_f(complex_f z){ auto r = bs_wrap_complex<float>([&]{ return std::atan(to_std(z)); }); return from_std(r); }

complex_l bs_cexp_l(complex_l z) { auto r = bs_wrap_complex<long double>([&]{ return std::exp(to_std(z)); }); return from_std(r); }
complex_l bs_clog_l(complex_l z) { auto r = bs_wrap_complex<long double>([&]{ return std::log(to_std(z)); }); return from_std(r); }
complex_l bs_csqrt_l(complex_l z){ auto r = bs_wrap_complex<long double>([&]{ return std::sqrt(to_std(z)); }); return from_std(r); }
complex_l bs_csin_l(complex_l z) { auto r = bs_wrap_complex<long double>([&]{ return std::sin(to_std(z)); }); return from_std(r); }
complex_l bs_ccos_l(complex_l z) { auto r = bs_wrap_complex<long double>([&]{ return std::cos(to_std(z)); }); return from_std(r); }
complex_l bs_ctan_l(complex_l z) { auto r = bs_wrap_complex<long double>([&]{ return std::tan(to_std(z)); }); return from_std(r); }
complex_l bs_csinh_l(complex_l z){ auto r = bs_wrap_complex<long double>([&]{ return std::sinh(to_std(z)); }); return from_std(r); }
complex_l bs_ccosh_l(complex_l z){ auto r = bs_wrap_complex<long double>([&]{ return std::cosh(to_std(z)); }); return from_std(r); }
complex_l bs_ctanh_l(complex_l z){ auto r = bs_wrap_complex<long double>([&]{ return std::tanh(to_std(z)); }); return from_std(r); }
complex_l bs_catan_l(complex_l z){ auto r = bs_wrap_complex<long double>([&]{ return std::atan(to_std(z)); }); return from_std(r); }
#ifdef __cplusplus
}
#endif

