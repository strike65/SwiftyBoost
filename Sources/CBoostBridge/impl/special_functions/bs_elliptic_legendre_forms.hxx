// Elliptic integrals (Legendre forms)
#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>
#include <boost/math/special_functions/ellint_3.hpp>
#include "../internal/bs_internal.hpp"

#ifdef __cplusplus
extern "C" {
#endif

double bs_ellint_1_complete(double k)                 { return bs_wrap<double>([&] { return boost::math::ellint_1(k); }); }
double bs_ellint_1(double k, double phi)              { return bs_wrap<double>([&] { return boost::math::ellint_1(k, phi); }); }
double bs_ellint_2_complete(double k)                 { return bs_wrap<double>([&] { return boost::math::ellint_2(k); }); }
double bs_ellint_2(double k, double phi)              { return bs_wrap<double>([&] { return boost::math::ellint_2(k, phi); }); }
double bs_ellint_3(double k, double nu, double phi)   { return bs_wrap<double>([&] { return boost::math::ellint_3(k, nu, phi); }); }
double bs_ellint_3_complete(double k, double nu)      { return bs_wrap<double>([&] { return boost::math::ellint_3(k, nu); }); }

float bs_ellint_1_complete_f(float k)                 { return bs_wrap<float>([&] { return boost::math::ellint_1(k); }); }
float bs_ellint_1_f(float k, float phi)               { return bs_wrap<float>([&] { return boost::math::ellint_1(k, phi); }); }
float bs_ellint_2_complete_f(float k)                 { return bs_wrap<float>([&] { return boost::math::ellint_2(k); }); }
float bs_ellint_2_f(float k, float phi)               { return bs_wrap<float>([&] { return boost::math::ellint_2(k, phi); }); }
float bs_ellint_3_f(float k, float nu, float phi)     { return bs_wrap<float>([&] { return boost::math::ellint_3(k, nu, phi); }); }
float bs_ellint_3_complete_f(float k, float nu)       { return bs_wrap<float>([&] { return boost::math::ellint_3(k, nu); }); }

long double bs_ellint_1_complete_l(long double k)                 { return bs_wrap<long double>([&] { return boost::math::ellint_1(k); }); }
long double bs_ellint_1_l(long double k, long double phi)         { return bs_wrap<long double>([&] { return boost::math::ellint_1(k, phi); }); }
long double bs_ellint_2_complete_l(long double k)                 { return bs_wrap<long double>([&] { return boost::math::ellint_2(k); }); }
long double bs_ellint_2_l(long double k, long double phi)         { return bs_wrap<long double>([&] { return boost::math::ellint_2(k, phi); }); }
long double bs_ellint_3_l(long double k, long double nu,long double phi) { return bs_wrap<long double>([&] { return boost::math::ellint_3(k, nu, phi); }); }
long double bs_ellint_3_complete_l(long double k, long double nu)      { return bs_wrap<long double>([&] { return boost::math::ellint_3(k, nu); }); }

#ifdef __cplusplus
}
#endif

