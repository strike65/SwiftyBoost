// Gamma and incomplete gamma wrappers
#include <boost/math/special_functions/gamma.hpp>
#include "../internal/bs_internal.hpp"

extern "C" {

// Gamma
double bs_tgamma_d(double x)              { return bs_wrap<double>([&] { return boost::math::tgamma(x); }); }
double bs_lgamma_d(double x)              { return bs_wrap<double>([&] { return boost::math::lgamma(x); }); }
float  bs_tgamma_f(float x)               { return bs_wrap<float>([&] { return boost::math::tgamma(x); }); }
float  bs_lgamma_f(float x)               { return bs_wrap<float>([&] { return boost::math::lgamma(x); }); }
long double bs_tgamma_l(long double x)    { return bs_wrap<long double>([&] { return boost::math::tgamma(x); }); }
long double bs_lgamma_l(long double x)    { return bs_wrap<long double>([&] { return boost::math::lgamma(x); }); }

// Incomplete gamma (lower/upper, regularized, and inverses)
double bs_tgamma_lower_d(double a, double x)      { return bs_wrap<double>([&] { return boost::math::tgamma_lower(a, x); }); }
double bs_tgamma_upper_d(double a, double x)      { return bs_wrap<double>([&] { return boost::math::tgamma(a, x); }); }
double bs_gamma_p_d(double a, double x)           { return bs_wrap<double>([&] { return boost::math::gamma_p(a, x); }); }
double bs_gamma_q_d(double a, double x)           { return bs_wrap<double>([&] { return boost::math::gamma_q(a, x); }); }
double bs_gamma_p_inv_d(double a, double p)       { return bs_wrap<double>([&] { return boost::math::gamma_p_inv(a, p); }); }
double bs_gamma_q_inv_d(double a, double q)       { return bs_wrap<double>([&] { return boost::math::gamma_q_inv(a, q); }); }
double bs_gamma_p_derivative_d(double a, double x) { return bs_wrap<double>([&] { return boost::math::gamma_p_derivative(a, x); }); }

float bs_tgamma_lower_f(float a, float x)       { return bs_wrap<float>([&] { return boost::math::tgamma_lower(a, x); }); }
float bs_tgamma_upper_f(float a, float x)       { return bs_wrap<float>([&] { return boost::math::tgamma(a, x); }); }
float bs_gamma_p_f(float a, float x)            { return bs_wrap<float>([&] { return boost::math::gamma_p(a, x); }); }
float bs_gamma_q_f(float a, float x)            { return bs_wrap<float>([&] { return boost::math::gamma_q(a, x); }); }
float bs_gamma_p_inv_f(float a, float p)        { return bs_wrap<float>([&] { return boost::math::gamma_p_inv(a, p); }); }
float bs_gamma_q_inv_f(float a, float q)        { return bs_wrap<float>([&] { return boost::math::gamma_q_inv(a, q); }); }
float bs_gamma_p_derivative_f(float a, float x) { return bs_wrap<float>([&] { return boost::math::gamma_p_derivative(a, x); }); }

long double bs_tgamma_lower_l(long double a, long double x) { return bs_wrap<long double>([&] { return boost::math::tgamma_lower(a, x); }); }
long double bs_tgamma_upper_l(long double a, long double x) { return bs_wrap<long double>([&] { return boost::math::tgamma(a, x); }); }
long double bs_gamma_p_l(long double a, long double x)      { return bs_wrap<long double>([&] { return boost::math::gamma_p(a, x); }); }
long double bs_gamma_q_l(long double a, long double x)      { return bs_wrap<long double>([&] { return boost::math::gamma_q(a, x); }); }
long double bs_gamma_p_inv_l(long double a, long double p)  { return bs_wrap<long double>([&] { return boost::math::gamma_p_inv(a, p); }); }
long double bs_gamma_q_inv_l(long double a, long double q)  { return bs_wrap<long double>([&] { return boost::math::gamma_q_inv(a, q); }); }
long double bs_gamma_p_derivative_l(long double a, long double x) { return bs_wrap<long double>([&] { return boost::math::gamma_p_derivative(a, x); }); }

// Γ(a) / Γ(b) and Γ(a) / Γ(a+delta)
double bs_tgamma_ratio_d(double a, double b) { return bs_wrap<double>([&] { return boost::math::tgamma_ratio(a, b); }); }
float  bs_tgamma_ratio_f(float a, float b) { return bs_wrap<float>([&] { return boost::math::tgamma_ratio(a, b); }); }
long double bs_tgamma_ratio_l(long double a, long double b) { return bs_wrap<long double>([&] { return boost::math::tgamma_ratio(a, b); }); }
double bs_tgamma_delta_ratio_d(double a, double delta)  { return bs_wrap<double>([&] { return boost::math::tgamma_delta_ratio(a, delta); }); }
float  bs_tgamma_delta_ratio_f(float a, float delta)  { return bs_wrap<float>([&] { return boost::math::tgamma_delta_ratio(a, delta); }); }
long double bs_tgamma_delta_ratio_l(long double a, long double delta)  { return bs_wrap<long double>([&] { return boost::math::tgamma_delta_ratio(a, delta); }); }

} // extern "C"

