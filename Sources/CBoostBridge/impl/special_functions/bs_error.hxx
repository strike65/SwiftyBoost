// Error function wrappers
#include <boost/math/special_functions/erf.hpp>
#include "../internal/bs_internal.hpp"

extern "C" {

double bs_erf_d(double x)                 { return bs_wrap<double>([&] { return boost::math::erf(x); }); }
double bs_erfc_d(double x)                { return bs_wrap<double>([&] { return boost::math::erfc(x); }); }
double bs_erf_inv_d(double p)             { return bs_wrap<double>([&] { return boost::math::erf_inv(p); }); }
double bs_erfc_inv_d(double p)            { return bs_wrap<double>([&] { return boost::math::erfc_inv(p); }); }

float bs_erf_f(float x)                 { return bs_wrap<float>([&] { return boost::math::erf(x); }); }
float bs_erfc_f(float x)                { return bs_wrap<float>([&] { return boost::math::erfc(x); }); }
float bs_erf_inv_f(float p)             { return bs_wrap<float>([&] { return boost::math::erf_inv(p); }); }
float bs_erfc_inv_f(float p)            { return bs_wrap<float>([&] { return boost::math::erfc_inv(p); }); }

long double bs_erf_l(long double x)     { return bs_wrap<long double>([&] { return boost::math::erf(x); }); }
long double bs_erfc_l(long double x)    { return bs_wrap<long double>([&] { return boost::math::erfc(x); }); }
long double bs_erf_inv_l(long double p) { return bs_wrap<long double>([&] { return boost::math::erf_inv(p); }); }
long double bs_erfc_inv_l(long double p){ return bs_wrap<long double>([&] { return boost::math::erfc_inv(p); }); }

} // extern "C"

