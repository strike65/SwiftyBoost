// Trigonometric helpers (sin_pi, cos_pi)
#include <boost/math/special_functions/sin_pi.hpp>
#include <boost/math/special_functions/cos_pi.hpp>
#include "../internal/bs_internal.hpp"

#ifdef __cplusplus
extern "C" {
#endif

double bs_sin_pi_d(double x)             { return bs_wrap<double>([&] { return boost::math::sin_pi(x); }); }
double bs_cos_pi_d(double x)             { return bs_wrap<double>([&] { return boost::math::cos_pi(x); }); }

float bs_sin_pi_f(float x)             { return bs_wrap<float>([&] { return boost::math::sin_pi(x); }); }
float bs_cos_pi_f(float x)             { return bs_wrap<float>([&] { return boost::math::cos_pi(x); }); }

long double bs_sin_pi_l(long double x) { return bs_wrap<long double>([&] { return boost::math::sin_pi(x); }); }
long double bs_cos_pi_l(long double x) { return bs_wrap<long double>([&] { return boost::math::cos_pi(x); }); }

#ifdef __cplusplus
}
#endif

