#pragma once

#ifdef __cplusplus
extern "C" {
#endif

// Gamma and log-gamma
double bs_tgamma(double x);
double bs_lgamma(double x);
double bs_tgamma_ratio(double a, double b);
double bs_tgamma_delta_ratio(double a, double delta);

float bs_tgamma_f(float x);
float bs_lgamma_f(float x);
float bs_tgamma_ratio_f(float a, float b);
float bs_tgamma_delta_ratio_f(float a, float delta);

long double bs_tgamma_l(long double x);
long double bs_lgamma_l(long double x);
long double bs_tgamma_ratio_l(long double a, long double b);
long double bs_tgamma_delta_ratio_l(long double a, long double delta);

// Incomplete gamma (lower/upper, regularized, and inverses)
double bs_tgamma_lower(double a, double x);
double bs_tgamma_upper(double a, double x);
double bs_gamma_p(double a, double x);
double bs_gamma_q(double a, double x);
double bs_gamma_p_inv(double a, double p);
double bs_gamma_q_inv(double a, double q);

// Derivatives of regularized incomplete gamma (w.r.t. x)
double bs_gamma_p_derivative(double a, double x);

float bs_tgamma_lower_f(float a, float x);
float bs_tgamma_upper_f(float a, float x);
float bs_gamma_p_f(float a, float x);
float bs_gamma_q_f(float a, float x);
float bs_gamma_p_inv_f(float a, float p);
float bs_gamma_q_inv_f(float a, float q);

// Derivatives (float)
float bs_gamma_p_derivative_f(float a, float x);

long double bs_tgamma_lower_l(long double a, long double x);
long double bs_tgamma_upper_l(long double a, long double x);
long double bs_gamma_p_l(long double a, long double x);
long double bs_gamma_q_l(long double a, long double x);
long double bs_gamma_p_inv_l(long double a, long double p);
long double bs_gamma_q_inv_l(long double a, long double q);

// Derivatives (long double)
long double bs_gamma_p_derivative_l(long double a, long double x);

#ifdef __cplusplus
}
#endif

