// Airy functions
#include <boost/math/special_functions/airy.hpp>
#include "../internal/bs_internal.hpp"

#ifdef __cplusplus
extern "C" {
#endif

double bs_airy_ai_d(double x)            { return bs_wrap<double>([&] { return boost::math::airy_ai(x); }); }
double bs_airy_bi_d(double x)            { return bs_wrap<double>([&] { return boost::math::airy_bi(x); }); }
double bs_airy_ai_prime_d(double x)      { return bs_wrap<double>([&] { return boost::math::airy_ai_prime(x); }); }
double bs_airy_bi_prime_d(double x)      { return bs_wrap<double>([&] { return boost::math::airy_bi_prime(x); }); }

float bs_airy_ai_f(float x)            { return bs_wrap<float>([&] { return boost::math::airy_ai(x); }); }
float bs_airy_bi_f(float x)            { return bs_wrap<float>([&] { return boost::math::airy_bi(x); }); }
float bs_airy_ai_prime_f(float x)      { return bs_wrap<float>([&] { return boost::math::airy_ai_prime(x); }); }
float bs_airy_bi_prime_f(float x)      { return bs_wrap<float>([&] { return boost::math::airy_bi_prime(x); }); }

long double bs_airy_ai_l(long double x)       { return bs_wrap<long double>([&] { return boost::math::airy_ai(x); }); }
long double bs_airy_bi_l(long double x)       { return bs_wrap<long double>([&] { return boost::math::airy_bi(x); }); }
long double bs_airy_ai_prime_l(long double x) { return bs_wrap<long double>([&] { return boost::math::airy_ai_prime(x); }); }
long double bs_airy_bi_prime_l(long double x) { return bs_wrap<long double>([&] { return boost::math::airy_bi_prime(x); }); }

#ifdef __cplusplus
}
#endif

