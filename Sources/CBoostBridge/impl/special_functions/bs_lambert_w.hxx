// Lambert W wrappers
#include <boost/math/special_functions/lambert_w.hpp>
#include "../internal/bs_internal.hpp"

extern "C" {

double bs_lambert_w0(double x)        { return bs_wrap<double>([&] { return boost::math::lambert_w0(x); }); }
double bs_lambert_wm1(double x)       { return bs_wrap<double>([&] { return boost::math::lambert_wm1(x); }); }

float bs_lambert_w0_f(float x)        { return bs_wrap<float>([&] { return boost::math::lambert_w0(x); }); }
float bs_lambert_wm1_f(float x)       { return bs_wrap<float>([&] { return boost::math::lambert_wm1(x); }); }

long double bs_lambert_w0_l(long double x)  { return bs_wrap<long double>([&] { return boost::math::lambert_w0(x); }); }
long double bs_lambert_wm1_l(long double x) { return bs_wrap<long double>([&] { return boost::math::lambert_wm1(x); }); }

} // extern "C"

