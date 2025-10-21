// Factorials and combinatorics
#include <boost/math/special_functions/factorials.hpp>
#include "../internal/bs_internal.hpp"

extern "C" {

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

// Double factorial
double bs_double_factorial_d(unsigned int i) { return bs_wrap<double>([&] { return boost::math::double_factorial<double>(i); }); }
float bs_double_factorial_f(unsigned int i) { return bs_wrap<float>([&] { return boost::math::double_factorial<float>(i); }); }
long double bs_double_factorial_l(unsigned int i) { return bs_wrap<long double>([&] { return boost::math::double_factorial<long double>(i); }); }

} // extern "C"

