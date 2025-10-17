// Complex number helpers and elementary functions
#include <complex>
#include "../internal/bs_internal.hpp"
#include "../include/bs_complex.h"

// Converters for C POD <-> std::complex
static inline std::complex<double> to_std(bs_complex_d z) noexcept { return { z.re, z.im }; }
static inline bs_complex_d from_std(std::complex<double> z) noexcept { return { static_cast<double>(z.real()), static_cast<double>(z.imag()) }; }

static inline std::complex<float> to_std(bs_complex_f z) noexcept { return { z.re, z.im }; }
static inline bs_complex_f from_std(std::complex<float> z) noexcept { return { static_cast<float>(z.real()), static_cast<float>(z.imag()) }; }

static inline std::complex<long double> to_std(bs_complex_l z) noexcept { return { z.re, z.im }; }
static inline bs_complex_l from_std(std::complex<long double> z) noexcept { return { static_cast<long double>(z.real()), static_cast<long double>(z.imag()) }; }

extern "C" {

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

// Elementary functions
bs_complex_d bs_cexp(bs_complex_d z) { auto r = bs_wrap_complex<double>([&]{ return std::exp(to_std(z)); }); return from_std(r); }
bs_complex_d bs_clog(bs_complex_d z) { auto r = bs_wrap_complex<double>([&]{ return std::log(to_std(z)); }); return from_std(r); }
bs_complex_d bs_csqrt(bs_complex_d z){ auto r = bs_wrap_complex<double>([&]{ return std::sqrt(to_std(z)); }); return from_std(r); }
bs_complex_d bs_csin(bs_complex_d z) { auto r = bs_wrap_complex<double>([&]{ return std::sin(to_std(z)); }); return from_std(r); }
bs_complex_d bs_ccos(bs_complex_d z) { auto r = bs_wrap_complex<double>([&]{ return std::cos(to_std(z)); }); return from_std(r); }
bs_complex_d bs_ctan(bs_complex_d z) { auto r = bs_wrap_complex<double>([&]{ return std::tan(to_std(z)); }); return from_std(r); }
bs_complex_d bs_csinh(bs_complex_d z){ auto r = bs_wrap_complex<double>([&]{ return std::sinh(to_std(z)); }); return from_std(r); }
bs_complex_d bs_ccosh(bs_complex_d z){ auto r = bs_wrap_complex<double>([&]{ return std::cosh(to_std(z)); }); return from_std(r); }
bs_complex_d bs_ctanh(bs_complex_d z){ auto r = bs_wrap_complex<double>([&]{ return std::tanh(to_std(z)); }); return from_std(r); }
bs_complex_d bs_catan(bs_complex_d z){ auto r = bs_wrap_complex<double>([&]{ return std::atan(to_std(z)); }); return from_std(r); }

bs_complex_f bs_cexp_f(bs_complex_f z) { auto r = bs_wrap_complex<float>([&]{ return std::exp(to_std(z)); }); return from_std(r); }
bs_complex_f bs_clog_f(bs_complex_f z) { auto r = bs_wrap_complex<float>([&]{ return std::log(to_std(z)); }); return from_std(r); }
bs_complex_f bs_csqrt_f(bs_complex_f z){ auto r = bs_wrap_complex<float>([&]{ return std::sqrt(to_std(z)); }); return from_std(r); }
bs_complex_f bs_csin_f(bs_complex_f z) { auto r = bs_wrap_complex<float>([&]{ return std::sin(to_std(z)); }); return from_std(r); }
bs_complex_f bs_ccos_f(bs_complex_f z) { auto r = bs_wrap_complex<float>([&]{ return std::cos(to_std(z)); }); return from_std(r); }
bs_complex_f bs_ctan_f(bs_complex_f z) { auto r = bs_wrap_complex<float>([&]{ return std::tan(to_std(z)); }); return from_std(r); }
bs_complex_f bs_csinh_f(bs_complex_f z){ auto r = bs_wrap_complex<float>([&]{ return std::sinh(to_std(z)); }); return from_std(r); }
bs_complex_f bs_ccosh_f(bs_complex_f z){ auto r = bs_wrap_complex<float>([&]{ return std::cosh(to_std(z)); }); return from_std(r); }
bs_complex_f bs_ctanh_f(bs_complex_f z){ auto r = bs_wrap_complex<float>([&]{ return std::tanh(to_std(z)); }); return from_std(r); }
bs_complex_f bs_catan_f(bs_complex_f z){ auto r = bs_wrap_complex<float>([&]{ return std::atan(to_std(z)); }); return from_std(r); }

bs_complex_l bs_cexp_l(bs_complex_l z) { auto r = bs_wrap_complex<long double>([&]{ return std::exp(to_std(z)); }); return from_std(r); }
bs_complex_l bs_clog_l(bs_complex_l z) { auto r = bs_wrap_complex<long double>([&]{ return std::log(to_std(z)); }); return from_std(r); }
bs_complex_l bs_csqrt_l(bs_complex_l z){ auto r = bs_wrap_complex<long double>([&]{ return std::sqrt(to_std(z)); }); return from_std(r); }
bs_complex_l bs_csin_l(bs_complex_l z) { auto r = bs_wrap_complex<long double>([&]{ return std::sin(to_std(z)); }); return from_std(r); }
bs_complex_l bs_ccos_l(bs_complex_l z) { auto r = bs_wrap_complex<long double>([&]{ return std::cos(to_std(z)); }); return from_std(r); }
bs_complex_l bs_ctan_l(bs_complex_l z) { auto r = bs_wrap_complex<long double>([&]{ return std::tan(to_std(z)); }); return from_std(r); }
bs_complex_l bs_csinh_l(bs_complex_l z){ auto r = bs_wrap_complex<long double>([&]{ return std::sinh(to_std(z)); }); return from_std(r); }
bs_complex_l bs_ccosh_l(bs_complex_l z){ auto r = bs_wrap_complex<long double>([&]{ return std::cosh(to_std(z)); }); return from_std(r); }
bs_complex_l bs_ctanh_l(bs_complex_l z){ auto r = bs_wrap_complex<long double>([&]{ return std::tanh(to_std(z)); }); return from_std(r); }
bs_complex_l bs_catan_l(bs_complex_l z){ auto r = bs_wrap_complex<long double>([&]{ return std::atan(to_std(z)); }); return from_std(r); }

} // extern "C"

