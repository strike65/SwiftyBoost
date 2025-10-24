// Cardinal B-spline wrappers
#include <boost/math/special_functions/cardinal_b_spline.hpp>
#include <limits>
#include "../internal/bs_internal.hpp"

// Internal dispatchers for cardinal B-spline (templated)
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

} // extern "C"

