#pragma once

#ifdef __cplusplus
extern "C" {
#endif

// Cardinal B-spline (Boost.Math special_functions/cardinal_b_spline.hpp)
double bs_cardinal_b_spline(unsigned int n, double x);
double bs_cardinal_b_spline_prime(unsigned int n, double x);
double bs_cardinal_b_spline_double_prime(unsigned int n, double x);
double bs_forward_cardinal_b_spline(unsigned int n, double x);

float bs_cardinal_b_spline_f(unsigned int n, float x);
float bs_cardinal_b_spline_prime_f(unsigned int n, float x);
float bs_cardinal_b_spline_double_prime_f(unsigned int n, float x);
float bs_forward_cardinal_b_spline_f(unsigned int n, float x);

long double bs_cardinal_b_spline_l(unsigned int n, long double x);
long double bs_cardinal_b_spline_prime_l(unsigned int n, long double x);
long double bs_cardinal_b_spline_double_prime_l(unsigned int n, long double x);
long double bs_forward_cardinal_b_spline_l(unsigned int n, long double x);

#ifdef __cplusplus
}
#endif

