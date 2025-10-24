// Hermite polynomial wrappers (Boost.Math special_functions/hermite.hpp)
#include <boost/math/special_functions/hermite.hpp>
#include "../internal/bs_internal.hpp"

template <typename Real>
static inline Real hermite_eval(unsigned int n, Real x) noexcept {
    return bs_wrap<Real>([&] { return boost::math::hermite(n, x); });
}

template <typename Real, typename Value>
static inline Real hermite_next(unsigned int n, Real x, Value hn, Value hnm1) noexcept {
    return bs_wrap<Real>([&] { return boost::math::hermite_next(n, x, hn, hnm1); });
}

extern "C" {

double bs_hermite_d(unsigned int n, double x) { return hermite_eval<double>(n, x); }
double bs_hermite_next_d(unsigned int n, double x, double hn, double hnm1) {
    return hermite_next<double>(n, x, hn, hnm1);
}

float bs_hermite_f(unsigned int n, float x) { return hermite_eval<float>(n, x); }
float bs_hermite_next_f(unsigned int n, float x, float hn, float hnm1) {
    return hermite_next<float>(n, x, hn, hnm1);
}

long double bs_hermite_l(unsigned int n, long double x) { return hermite_eval<long double>(n, x); }
long double bs_hermite_next_l(unsigned int n, long double x, long double hn, long double hnm1) {
    return hermite_next<long double>(n, x, hn, hnm1);
}

} // extern "C"
