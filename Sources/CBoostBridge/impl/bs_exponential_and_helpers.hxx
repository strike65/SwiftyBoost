// Exponential integrals and related helpers
#include <boost/math/special_functions/expint.hpp>
#include <boost/math/special_functions/expm1.hpp>
#include <boost/math/special_functions/log1p.hpp>
#include <boost/math/special_functions/powm1.hpp>
#include <boost/math/special_functions/cbrt.hpp>
#include <boost/math/special_functions/sqrt1pm1.hpp>
#include <boost/math/special_functions/hypot.hpp>
#include <boost/math/special_functions/rsqrt.hpp>
#include "../internal/bs_internal.hpp"


extern "C" {

// Exponential integrals and related
double bs_expint_Ei_d(double x)        { return bs_wrap<double>([&] { return boost::math::expint(x); }); }
double bs_expint_En_d(int n, double x) { return bs_wrap<double>([&] { return boost::math::expint(n, x); }); }
double bs_expm1_d(double x)            { return bs_wrap<double>([&] { return boost::math::expm1(x); }); }
double bs_log1p_d(double x)            { return bs_wrap<double>([&] { return boost::math::log1p(x); }); }
double bs_log1pmx_d(double x)          { return bs_wrap<double>([&] { return boost::math::log1pmx(x); }); }
double bs_powm1_d(double x, double y)  { return bs_wrap<double>([&] { return boost::math::powm1(x, y); }); }
double bs_cbrt_d(double x)           { return bs_wrap<double>([&] { return boost::math::cbrt(x); }); }
double bs_sqrt1pm1_d(double x)         { return bs_wrap<double>([&] { return boost::math::sqrt1pm1(x); }); }
double bs_hypot_d(double x, double y)  { return bs_wrap<double>([&] { return boost::math::hypot(x, y); }); }
double bs_rsqrt_d(const double x)      { return bs_wrap<double>([&] { return boost::math::rsqrt(x); }); }

float bs_expint_Ei_f(float x)        { return bs_wrap<float>([&] { return boost::math::expint(x); }); }
float bs_expint_En_f(int n, float x) { return bs_wrap<float>([&] { return boost::math::expint(n, x); }); }
float bs_expm1_f(float x)            { return bs_wrap<float>([&] { return boost::math::expm1(x); }); }
float bs_log1p_f(float x)            { return bs_wrap<float>([&] { return boost::math::log1p(x); }); }
float bs_log1pmx_f(float x)          { return bs_wrap<float>([&] { return boost::math::log1pmx(x); }); }
float bs_powm1_f(float x, float y)   { return bs_wrap<float>([&] { return boost::math::powm1(x, y); }); }
float bs_cbrt_f(float x)             { return bs_wrap<float>([&] { return boost::math::cbrt(x); }); }
float bs_sqrt1pm1_f(float x)         { return bs_wrap<float>([&] { return boost::math::sqrt1pm1(x); }); }
float bs_hypot_f(float x, float y)   { return bs_wrap<float>([&] { return boost::math::hypot(x, y); }); }
float bs_rsqrt_f(const float x)      { return bs_wrap<float>([&] { return boost::math::rsqrt(x); }); }

long double bs_expint_Ei_l(long double x)          { return bs_wrap<long double>([&] { return boost::math::expint(x); }); }
long double bs_expint_En_l(int n, long double x)   { return bs_wrap<long double>([&] { return boost::math::expint(n, x); }); }
long double bs_expm1_l(long double x)              { return bs_wrap<long double>([&] { return boost::math::expm1(x); }); }
long double bs_log1p_l(long double x)              { return bs_wrap<long double>([&] { return boost::math::log1p(x); }); }
long double bs_log1pmx_l(long double x)            { return bs_wrap<long double>([&] { return boost::math::log1pmx(x); }); }
long double bs_powm1_l(long double x, long double y) { return bs_wrap<long double>([&] { return boost::math::powm1(x, y); }); }
long double bs_cbrt_l(long double x)               { return bs_wrap<long double>([&] { return boost::math::cbrt(x); }); }
long double bs_sqrt1pm1_l(long double x)           { return bs_wrap<long double>([&] { return boost::math::sqrt1pm1(x); }); }
long double bs_hypot_l(long double x, long double y)  { return bs_wrap<long double>([&] { return boost::math::hypot(x, y); }); }
long double bs_rsqrt_l(const long double x)        { return bs_wrap<long double>([&] { return boost::math::rsqrt(x); }); }


} // extern "C"

