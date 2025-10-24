// Digamma / Polygamma / Zeta wrappers
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/trigamma.hpp>
#include <boost/math/special_functions/polygamma.hpp>
#include <boost/math/special_functions/zeta.hpp>
#include "../internal/bs_internal.hpp"

#ifdef __cplusplus
extern "C" {
#endif

double bs_digamma_d(double x)            { return bs_wrap<double>([&] { return boost::math::digamma(x); }); }
double bs_trigamma_d(double x)           { return bs_wrap<double>([&] { return boost::math::trigamma(x); }); }
double bs_polygamma_d(int n, double x)   { return bs_wrap<double>([&] { return boost::math::polygamma(n, x); }); }
double bs_riemann_zeta_d(double x)       { return bs_wrap<double>([&] { return boost::math::zeta(x); }); }

float bs_digamma_f(float x)            { return bs_wrap<float>([&] { return boost::math::digamma(x); }); }
float bs_trigamma_f(float x)           { return bs_wrap<float>([&] { return boost::math::trigamma(x); }); }
float bs_polygamma_f(int n, float x)   { return bs_wrap<float>([&] { return boost::math::polygamma(n, x); }); }
float bs_riemann_zeta_f(float x)       { return bs_wrap<float>([&] { return boost::math::zeta(x); }); }

long double bs_digamma_l(long double x)          { return bs_wrap<long double>([&] { return boost::math::digamma(x); }); }
long double bs_trigamma_l(long double x)         { return bs_wrap<long double>([&] { return boost::math::trigamma(x); }); }
long double bs_polygamma_l(int n, long double x) { return bs_wrap<long double>([&] { return boost::math::polygamma(n, x); }); }
long double bs_riemann_zeta_l(long double x)     { return bs_wrap<long double>([&] { return boost::math::zeta(x); }); }

#ifdef __cplusplus
}
#endif

