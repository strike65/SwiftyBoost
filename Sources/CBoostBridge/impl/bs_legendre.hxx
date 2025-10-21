// Legendre polynomials and helpers
#include <boost/math/special_functions/legendre.hpp>
#include <algorithm>
#include <vector>
#include "../internal/bs_internal.hpp"

#ifdef __cplusplus
extern "C" {
#endif

// Legendre
double bs_legendre_p_d(int n, double x)            { return bs_wrap<double>([&] { return boost::math::legendre_p(n, x); }); }
double bs_assoc_legendre_p_d(int n, int m, double x) { return bs_wrap<double>([&] { return boost::math::legendre_p(n, m, x); }); }
double bs_legendre_p_prime_d(int n, double x)      { return bs_wrap<double>([&] { return boost::math::legendre_p_prime(n, x); }); }
void bs_legendre_p_zeros_d(int l, double* out) {
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

#ifdef __cplusplus
}
#endif

