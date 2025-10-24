// Jacobi polynomials (Boost.Math special_functions/jacobi.hpp)
#include <boost/math/special_functions/jacobi.hpp>
#include "../internal/bs_internal.hpp"

template <typename Real>
static inline Real jacobi_eval(unsigned int n, Real alpha, Real beta, Real x) noexcept {
    return bs_wrap<Real>([&] { return boost::math::jacobi<Real>(n, alpha, beta, x); });
}

template <typename Real>
static inline Real jacobi_prime(unsigned int n, Real alpha, Real beta, Real x) noexcept {
    return bs_wrap<Real>([&] { return boost::math::jacobi_prime<Real>(n, alpha, beta, x); });
}

template <typename Real>
static inline Real jacobi_double_prime(unsigned int n, Real alpha, Real beta, Real x) noexcept {
    return bs_wrap<Real>([&] { return boost::math::jacobi_double_prime<Real>(n, alpha, beta, x); });
}

template <typename Real>
static inline Real jacobi_derivative(unsigned int n, Real alpha, Real beta, Real x, unsigned int k) noexcept {
    return bs_wrap<Real>([&] { return boost::math::jacobi_derivative<Real>(n, alpha, beta, x, k); });
}

extern "C" {

double bs_jacobi_d(unsigned int n, double alpha, double beta, double x) {
    return jacobi_eval<double>(n, alpha, beta, x);
}
double bs_jacobi_prime_d(unsigned int n, double alpha, double beta, double x) {
    return jacobi_prime<double>(n, alpha, beta, x);
}
double bs_jacobi_double_prime_d(unsigned int n, double alpha, double beta, double x) {
    return jacobi_double_prime<double>(n, alpha, beta, x);
}
double bs_jacobi_derivative_d(unsigned int n, double alpha, double beta, double x, unsigned int k) {
    return jacobi_derivative<double>(n, alpha, beta, x, k);
}

float bs_jacobi_f(unsigned int n, float alpha, float beta, float x) {
    return jacobi_eval<float>(n, alpha, beta, x);
}
float bs_jacobi_prime_f(unsigned int n, float alpha, float beta, float x) {
    return jacobi_prime<float>(n, alpha, beta, x);
}
float bs_jacobi_double_prime_f(unsigned int n, float alpha, float beta, float x) {
    return jacobi_double_prime<float>(n, alpha, beta, x);
}
float bs_jacobi_derivative_f(unsigned int n, float alpha, float beta, float x, unsigned int k) {
    return jacobi_derivative<float>(n, alpha, beta, x, k);
}

long double bs_jacobi_l(unsigned int n, long double alpha, long double beta, long double x) {
    return jacobi_eval<long double>(n, alpha, beta, x);
}
long double bs_jacobi_prime_l(unsigned int n, long double alpha, long double beta, long double x) {
    return jacobi_prime<long double>(n, alpha, beta, x);
}
long double bs_jacobi_double_prime_l(unsigned int n, long double alpha, long double beta, long double x) {
    return jacobi_double_prime<long double>(n, alpha, beta, x);
}
long double bs_jacobi_derivative_l(unsigned int n, long double alpha, long double beta, long double x, unsigned int k) {
    return jacobi_derivative<long double>(n, alpha, beta, x, k);
}

} // extern "C"
