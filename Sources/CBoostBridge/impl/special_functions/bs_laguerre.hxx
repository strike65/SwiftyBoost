// Laguerre polynomials
#include <boost/math/special_functions/laguerre.hpp>
#include "../internal/bs_internal.hpp"

extern "C" {

double bs_laguerre_d(unsigned int n, double x) { return bs_wrap<double>([&] { return boost::math::laguerre(n, x); }); }
float bs_laguerre_f(unsigned int n, float x) { return bs_wrap<float>([&] { return boost::math::laguerre(n, x); }); }
long double bs_laguerre_l(unsigned int n, long double x) { return bs_wrap<long double>([&] { return boost::math::laguerre(n, x); }); }
double bs_assoc_laguerre_d(unsigned int n, unsigned int m, double x) { return bs_wrap<double>([&] { return boost::math::laguerre(n, m, x); }); }
float bs_assoc_laguerre_f(unsigned int n, unsigned int m, float x) { return bs_wrap<float>([&] { return boost::math::laguerre(n, m, x); }); }
long double bs_assoc_laguerre_l(unsigned int n, unsigned int m, long double x) { return bs_wrap<long double>([&] { return boost::math::laguerre(n, m, x); }); }

} // extern "C"

