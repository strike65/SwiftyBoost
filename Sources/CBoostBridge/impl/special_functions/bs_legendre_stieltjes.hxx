// Legendre-Stieltjes polynomial wrappers (Boost.Math special_functions/legendre_stieltjes.hpp)
#include <boost/math/special_functions/legendre_stieltjes.hpp>
#include <algorithm>
#include <cstddef>
#include <vector>
#include "../internal/bs_internal.hpp"

template <typename Real>
static inline Real legendre_stieltjes_eval(unsigned int m, Real x) noexcept {
    return bs_wrap<Real>([&] {
        boost::math::legendre_stieltjes<Real> poly(static_cast<std::size_t>(m));
        return poly(x);
    });
}

template <typename Real>
static inline Real legendre_stieltjes_prime(unsigned int m, Real x) noexcept {
    return bs_wrap<Real>([&] {
        boost::math::legendre_stieltjes<Real> poly(static_cast<std::size_t>(m));
        return poly.prime(x);
    });
}

template <typename Real>
static inline Real legendre_stieltjes_norm_sq(unsigned int m) noexcept {
    return bs_wrap<Real>([&] {
        boost::math::legendre_stieltjes<Real> poly(static_cast<std::size_t>(m));
        return poly.norm_sq();
    });
}

template <typename Real>
static inline std::size_t legendre_stieltjes_zeros(unsigned int m, Real *out, std::size_t count) noexcept {
    try {
        boost::math::legendre_stieltjes<Real> poly(static_cast<std::size_t>(m));
        auto zeros = poly.zeros();
        if (out && count > 0) {
            const std::size_t n = std::min<std::size_t>(count, zeros.size());
            std::copy(zeros.begin(), zeros.begin() + n, out);
        }
        return zeros.size();
    } catch (...) {
        return 0;
    }
}

extern "C" {

double bs_legendre_stieltjes_d(unsigned int m, double x) { return legendre_stieltjes_eval<double>(m, x); }
double bs_legendre_stieltjes_prime_d(unsigned int m, double x) { return legendre_stieltjes_prime<double>(m, x); }
double bs_legendre_stieltjes_norm_sq_d(unsigned int m) { return legendre_stieltjes_norm_sq<double>(m); }
std::size_t bs_legendre_stieltjes_zeros_d(unsigned int m, double *out, std::size_t count) {
    return legendre_stieltjes_zeros<double>(m, out, count);
}

float bs_legendre_stieltjes_f(unsigned int m, float x) { return legendre_stieltjes_eval<float>(m, x); }
float bs_legendre_stieltjes_prime_f(unsigned int m, float x) { return legendre_stieltjes_prime<float>(m, x); }
float bs_legendre_stieltjes_norm_sq_f(unsigned int m) { return legendre_stieltjes_norm_sq<float>(m); }
std::size_t bs_legendre_stieltjes_zeros_f(unsigned int m, float *out, std::size_t count) {
    return legendre_stieltjes_zeros<float>(m, out, count);
}

long double bs_legendre_stieltjes_l(unsigned int m, long double x) { return legendre_stieltjes_eval<long double>(m, x); }
long double bs_legendre_stieltjes_prime_l(unsigned int m, long double x) { return legendre_stieltjes_prime<long double>(m, x); }
long double bs_legendre_stieltjes_norm_sq_l(unsigned int m) { return legendre_stieltjes_norm_sq<long double>(m); }
std::size_t bs_legendre_stieltjes_zeros_l(unsigned int m, long double *out, std::size_t count) {
    return legendre_stieltjes_zeros<long double>(m, out, count);
}

} // extern "C"
