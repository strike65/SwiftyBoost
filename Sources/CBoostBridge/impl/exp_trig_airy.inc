// Exponential integrals, special helpers, trig helpers, and Airy
#include <boost/math/special_functions/expint.hpp>
#include <boost/math/special_functions/expm1.hpp>
#include <boost/math/special_functions/log1p.hpp>
#include <boost/math/special_functions/log1p.hpp>
#include <boost/math/special_functions/log1p.hpp>
#include <boost/math/special_functions/powm1.hpp>
#include <boost/math/special_functions/cbrt.hpp>
#include <boost/math/special_functions/sqrt1pm1.hpp>
#include <boost/math/special_functions/hypot.hpp>
#include <boost/math/special_functions/sin_pi.hpp>
#include <boost/math/special_functions/cos_pi.hpp>
#include <boost/math/special_functions/airy.hpp>
#include "../internal/bs_internal.hpp"

extern "C" {

// Exponential integrals and related
double bs_expint_Ei(double x)        { return bs_wrap<double>([&] { return boost::math::expint(x); }); }
double bs_expint_En(int n, double x) { return bs_wrap<double>([&] { return boost::math::expint(n, x); }); }
double bs_expm1(double x)            { return bs_wrap<double>([&] { return boost::math::expm1(x); }); }
double bs_log1p(double x)            { return bs_wrap<double>([&] { return boost::math::log1p(x); }); }
double bs_log1pmx(double x)          { return bs_wrap<double>([&] { return boost::math::log1pmx(x); }); }
double bs_powm1(double x, double y)  { return bs_wrap<double>([&] { return boost::math::powm1(x, y); }); }
double bs_cbrt(double x)             { return bs_wrap<double>([&] { return boost::math::cbrt(x); }); }
double bs_sqrt1pm1(double x)         { return bs_wrap<double>([&] { return boost::math::sqrt1pm1(x); }); }
double bs_hypot(double x, double y)  { return bs_wrap<double>([&] { return boost::math::hypot(x, y); }); }

float bs_expint_Ei_f(float x)        { return bs_wrap<float>([&] { return boost::math::expint(x); }); }
float bs_expint_En_f(int n, float x) { return bs_wrap<float>([&] { return boost::math::expint(n, x); }); }
float bs_expm1_f(float x)            { return bs_wrap<float>([&] { return boost::math::expm1(x); }); }
float bs_log1p_f(float x)            { return bs_wrap<float>([&] { return boost::math::log1p(x); }); }
float bs_log1pmx_f(float x)          { return bs_wrap<float>([&] { return boost::math::log1pmx(x); }); }
float bs_powm1_f(float x, float y)   { return bs_wrap<float>([&] { return boost::math::powm1(x, y); }); }
float bs_cbrt_f(float x)             { return bs_wrap<float>([&] { return boost::math::cbrt(x); }); }
float bs_sqrt1pm1_f(float x)         { return bs_wrap<float>([&] { return boost::math::sqrt1pm1(x); }); }
float bs_hypot_f(float x, float y)   { return bs_wrap<float>([&] { return boost::math::hypot(x, y); }); }

long double bs_expint_Ei_l(long double x)          { return bs_wrap<long double>([&] { return boost::math::expint(x); }); }
long double bs_expint_En_l(int n, long double x)   { return bs_wrap<long double>([&] { return boost::math::expint(n, x); }); }
long double bs_expm1_l(long double x)              { return bs_wrap<long double>([&] { return boost::math::expm1(x); }); }
long double bs_log1p_l(long double x)              { return bs_wrap<long double>([&] { return boost::math::log1p(x); }); }
long double bs_log1pmx_l(long double x)            { return bs_wrap<long double>([&] { return boost::math::log1pmx(x); }); }
long double bs_powm1_l(long double x, long double y) { return bs_wrap<long double>([&] { return boost::math::powm1(x, y); }); }
long double bs_cbrt_l(long double x)               { return bs_wrap<long double>([&] { return boost::math::cbrt(x); }); }
long double bs_sqrt1pm1_l(long double x)           { return bs_wrap<long double>([&] { return boost::math::sqrt1pm1(x); }); }
long double bs_hypot_l(long double x, long double y)  { return bs_wrap<long double>([&] { return boost::math::hypot(x, y); }); }

// Trig helpers
double bs_sin_pi(double x)             { return bs_wrap<double>([&] { return boost::math::sin_pi(x); }); }
double bs_cos_pi(double x)             { return bs_wrap<double>([&] { return boost::math::cos_pi(x); }); }

float bs_sin_pi_f(float x)             { return bs_wrap<float>([&] { return boost::math::sin_pi(x); }); }
float bs_cos_pi_f(float x)             { return bs_wrap<float>([&] { return boost::math::cos_pi(x); }); }

long double bs_sin_pi_l(long double x) { return bs_wrap<long double>([&] { return boost::math::sin_pi(x); }); }
long double bs_cos_pi_l(long double x) { return bs_wrap<long double>([&] { return boost::math::cos_pi(x); }); }

// Airy
double bs_airy_ai(double x)            { return bs_wrap<double>([&] { return boost::math::airy_ai(x); }); }
double bs_airy_bi(double x)            { return bs_wrap<double>([&] { return boost::math::airy_bi(x); }); }
double bs_airy_ai_prime(double x)      { return bs_wrap<double>([&] { return boost::math::airy_ai_prime(x); }); }
double bs_airy_bi_prime(double x)      { return bs_wrap<double>([&] { return boost::math::airy_bi_prime(x); }); }

float bs_airy_ai_f(float x)            { return bs_wrap<float>([&] { return boost::math::airy_ai(x); }); }
float bs_airy_bi_f(float x)            { return bs_wrap<float>([&] { return boost::math::airy_bi(x); }); }
float bs_airy_ai_prime_f(float x)      { return bs_wrap<float>([&] { return boost::math::airy_ai_prime(x); }); }
float bs_airy_bi_prime_f(float x)      { return bs_wrap<float>([&] { return boost::math::airy_bi_prime(x); }); }

long double bs_airy_ai_l(long double x)       { return bs_wrap<long double>([&] { return boost::math::airy_ai(x); }); }
long double bs_airy_bi_l(long double x)       { return bs_wrap<long double>([&] { return boost::math::airy_bi(x); }); }
long double bs_airy_ai_prime_l(long double x) { return bs_wrap<long double>([&] { return boost::math::airy_ai_prime(x); }); }
long double bs_airy_bi_prime_l(long double x) { return bs_wrap<long double>([&] { return boost::math::airy_bi_prime(x); }); }

} // extern "C"

