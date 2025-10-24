// Beta and incomplete beta wrappers
#include <boost/math/special_functions/beta.hpp>
#include "../internal/bs_internal.hpp"

#ifdef __cplusplus
extern "C" {
#endif

// Beta family
double bs_beta_d(double a, double b)                     { return bs_wrap<double>([&] { return boost::math::beta(a, b); }); }
double bs_fullBeta_d(double a, double b, double x)       { return bs_wrap<double>([&] { return boost::math::beta(a, b, x); }); }
double bs_ibeta_d(double a, double b, double x)          { return bs_wrap<double>([&] { return boost::math::ibeta(a, b, x); }); }
double bs_ibetac_d(double a, double b, double x)         { return bs_wrap<double>([&] { return boost::math::ibetac(a, b, x); }); }
double bs_ibeta_inv_d(double a, double b, double p)      { return bs_wrap<double>([&] { return boost::math::ibeta_inv(a, b, p); }); }
double bs_ibetac_inv_d(double a, double b, double p)     { return bs_wrap<double>([&] { return boost::math::ibetac_inv(a, b, p); }); }
double bs_ibeta_inva_d(double b, double x, double p)     { return bs_wrap<double>([&] { return boost::math::ibeta_inva(b, x, p); }); }
double bs_ibeta_invb_d(double a, double x, double p)     { return bs_wrap<double>([&] { return boost::math::ibeta_invb(a, x, p); }); }
double bs_ibeta_derivative_d(double a, double b, double x) { return bs_wrap<double>([&] { return boost::math::ibeta_derivative(a, b, x); }); }

float bs_beta_f(float a, float b)                      { return bs_wrap<float>([&] { return boost::math::beta(a, b); }); }
float bs_fullBeta_f(float a, float b, float x)         { return bs_wrap<float>([&] { return boost::math::beta(a, b, x); }); }
float bs_ibeta_f(float a, float b, float x)            { return bs_wrap<float>([&] { return boost::math::ibeta(a, b, x); }); }
float bs_ibetac_f(float a, float b, float x)           { return bs_wrap<float>([&] { return boost::math::ibetac(a, b, x); }); }
float bs_ibeta_inv_f(float a, float b, float p)        { return bs_wrap<float>([&] { return boost::math::ibeta_inv(a, b, p); }); }
float bs_ibetac_inv_f(float a, float b, float p)       { return bs_wrap<float>([&] { return boost::math::ibetac_inv(a, b, p); }); }
float bs_ibeta_inva_f(float b, float x, float p)       { return bs_wrap<float>([&] { return boost::math::ibeta_inva(b, x, p); }); }
float bs_ibeta_invb_f(float a, float x, float p)       { return bs_wrap<float>([&] { return boost::math::ibeta_invb(a, x, p); }); }
float bs_ibeta_derivative_f(float a, float b, float x) { return bs_wrap<float>([&] { return boost::math::ibeta_derivative(a, b, x); }); }

long double bs_beta_l(long double a, long double b)                      { return bs_wrap<long double>([&] { return boost::math::beta(a, b); }); }
long double bs_fullBeta_l(long double a, long double b, long double x)   { return bs_wrap<long double>([&] { return boost::math::beta(a, b, x); }); }
long double bs_ibeta_l(long double a, long double b, long double x)      { return bs_wrap<long double>([&] { return boost::math::ibeta(a, b, x); }); }
long double bs_ibetac_l(long double a, long double b, long double x)     { return bs_wrap<long double>([&] { return boost::math::ibetac(a, b, x); }); }
long double bs_ibeta_inv_l(long double a, long double b, long double p)  { return bs_wrap<long double>([&] { return boost::math::ibeta_inv(a, b, p); }); }
long double bs_ibetac_inv_l(long double a, long double b, long double p) { return bs_wrap<long double>([&] { return boost::math::ibetac_inv(a, b, p); }); }
long double bs_ibeta_inva_l(long double b, long double x, long double p) { return bs_wrap<long double>([&] { return boost::math::ibeta_inva(b, x, p); }); }
long double bs_ibeta_invb_l(long double a, long double x, long double p) { return bs_wrap<long double>([&] { return boost::math::ibeta_invb(a, x, p); }); }
long double bs_ibeta_derivative_l(long double a, long double b, long double x) { return bs_wrap<long double>([&] { return boost::math::ibeta_derivative(a, b, x); }); }

#ifdef __cplusplus
}
#endif

