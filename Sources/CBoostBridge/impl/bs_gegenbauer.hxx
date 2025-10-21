// Gegenbauer polynomials
#include <boost/math/special_functions/gegenbauer.hpp>
#include "../internal/bs_internal.hpp"

extern "C" {

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

