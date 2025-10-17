// Gamma, Error functions, and incomplete gamma helpers
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/erf.hpp>
#include "../internal/bs_internal.hpp"

extern "C" {

// Gamma / Error
double bs_tgamma(double x)              { return bs_wrap<double>([&] { return boost::math::tgamma(x); }); }
double bs_lgamma(double x)              { return bs_wrap<double>([&] { return boost::math::lgamma(x); }); }
double bs_erf(double x)                 { return bs_wrap<double>([&] { return boost::math::erf(x); }); }
double bs_erfc(double x)                { return bs_wrap<double>([&] { return boost::math::erfc(x); }); }
double bs_erf_inv(double p)             { return bs_wrap<double>([&] { return boost::math::erf_inv(p); }); }
double bs_erfc_inv(double p)            { return bs_wrap<double>([&] { return boost::math::erfc_inv(p); }); }

float bs_tgamma_f(float x)              { return bs_wrap<float>([&] { return boost::math::tgamma(x); }); }
float bs_lgamma_f(float x)              { return bs_wrap<float>([&] { return boost::math::lgamma(x); }); }
float bs_erf_f(float x)                 { return bs_wrap<float>([&] { return boost::math::erf(x); }); }
float bs_erfc_f(float x)                { return bs_wrap<float>([&] { return boost::math::erfc(x); }); }
float bs_erf_inv_f(float p)             { return bs_wrap<float>([&] { return boost::math::erf_inv(p); }); }
float bs_erfc_inv_f(float p)            { return bs_wrap<float>([&] { return boost::math::erfc_inv(p); }); }

long double bs_tgamma_l(long double x)  { return bs_wrap<long double>([&] { return boost::math::tgamma(x); }); }
long double bs_lgamma_l(long double x)  { return bs_wrap<long double>([&] { return boost::math::lgamma(x); }); }
long double bs_erf_l(long double x)     { return bs_wrap<long double>([&] { return boost::math::erf(x); }); }
long double bs_erfc_l(long double x)    { return bs_wrap<long double>([&] { return boost::math::erfc(x); }); }
long double bs_erf_inv_l(long double p) { return bs_wrap<long double>([&] { return boost::math::erf_inv(p); }); }
long double bs_erfc_inv_l(long double p){ return bs_wrap<long double>([&] { return boost::math::erfc_inv(p); }); }

// Incomplete gamma (lower/upper, regularized, and inverses)
double bs_tgamma_lower(double a, double x)      { return bs_wrap<double>([&] { return boost::math::tgamma_lower(a, x); }); }
double bs_tgamma_upper(double a, double x)      { return bs_wrap<double>([&] { return boost::math::tgamma(a, x); }); }
double bs_gamma_p(double a, double x)           { return bs_wrap<double>([&] { return boost::math::gamma_p(a, x); }); }
double bs_gamma_q(double a, double x)           { return bs_wrap<double>([&] { return boost::math::gamma_q(a, x); }); }
double bs_gamma_p_inv(double a, double p)       { return bs_wrap<double>([&] { return boost::math::gamma_p_inv(a, p); }); }
double bs_gamma_q_inv(double a, double q)       { return bs_wrap<double>([&] { return boost::math::gamma_q_inv(a, q); }); }

// Derivatives of regularized incomplete gamma functions
double bs_gamma_p_derivative(double a, double x) { return bs_wrap<double>([&] { return boost::math::gamma_p_derivative(a, x); }); }

float bs_tgamma_lower_f(float a, float x)       { return bs_wrap<float>([&] { return boost::math::tgamma_lower(a, x); }); }
float bs_tgamma_upper_f(float a, float x)       { return bs_wrap<float>([&] { return boost::math::tgamma(a, x); }); }
float bs_gamma_p_f(float a, float x)            { return bs_wrap<float>([&] { return boost::math::gamma_p(a, x); }); }
float bs_gamma_q_f(float a, float x)            { return bs_wrap<float>([&] { return boost::math::gamma_q(a, x); }); }
float bs_gamma_p_inv_f(float a, float p)        { return bs_wrap<float>([&] { return boost::math::gamma_p_inv(a, p); }); }
float bs_gamma_q_inv_f(float a, float q)        { return bs_wrap<float>([&] { return boost::math::gamma_q_inv(a, q); }); }

// Derivatives (float)
float bs_gamma_p_derivative_f(float a, float x) { return bs_wrap<float>([&] { return boost::math::gamma_p_derivative(a, x); }); }

long double bs_tgamma_lower_l(long double a, long double x) { return bs_wrap<long double>([&] { return boost::math::tgamma_lower(a, x); }); }
long double bs_tgamma_upper_l(long double a, long double x) { return bs_wrap<long double>([&] { return boost::math::tgamma(a, x); }); }
long double bs_gamma_p_l(long double a, long double x)      { return bs_wrap<long double>([&] { return boost::math::gamma_p(a, x); }); }
long double bs_gamma_q_l(long double a, long double x)      { return bs_wrap<long double>([&] { return boost::math::gamma_q(a, x); }); }
long double bs_gamma_p_inv_l(long double a, long double p)  { return bs_wrap<long double>([&] { return boost::math::gamma_p_inv(a, p); }); }
long double bs_gamma_q_inv_l(long double a, long double q)  { return bs_wrap<long double>([&] { return boost::math::gamma_q_inv(a, q); }); }

// Derivatives (long double)
long double bs_gamma_p_derivative_l(long double a, long double x) { return bs_wrap<long double>([&] { return boost::math::gamma_p_derivative(a, x); }); }

// Γ(a) / Γ(b) and Γ(a) / Γ(a+delta)
double bs_tgamma_ratio(double a, double b) { return bs_wrap<double>([&] { return boost::math::tgamma_ratio(a, b); }); }
float  bs_tgamma_ratio_f(float a, float b) { return bs_wrap<float>([&] { return boost::math::tgamma_ratio(a, b); }); }
long double bs_tgamma_ratio_l(long double a, long double b) { return bs_wrap<long double>([&] { return boost::math::tgamma_ratio(a, b); }); }

double bs_tgamma_delta_ratio(double a, double delta)  { return bs_wrap<double>([&] { return boost::math::tgamma_delta_ratio(a, delta); }); }
float  bs_tgamma_delta_ratio_f(float a, float delta)  { return bs_wrap<float>([&] { return boost::math::tgamma_delta_ratio(a, delta); }); }
long double bs_tgamma_delta_ratio_l(long double a, long double delta)  { return bs_wrap<long double>([&] { return boost::math::tgamma_delta_ratio(a, delta); }); }

} // extern "C"

