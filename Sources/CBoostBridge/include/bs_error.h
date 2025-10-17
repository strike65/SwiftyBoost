#pragma once

#ifdef __cplusplus
extern "C" {
#endif

// Error functions and inverses
double bs_erf(double x);
double bs_erfc(double x);
double bs_erf_inv(double p);
double bs_erfc_inv(double p);

float bs_erf_f(float x);
float bs_erfc_f(float x);
float bs_erf_inv_f(float p);
float bs_erfc_inv_f(float p);

long double bs_erf_l(long double x);
long double bs_erfc_l(long double x);
long double bs_erf_inv_l(long double p);
long double bs_erfc_inv_l(long double p);

#ifdef __cplusplus
}
#endif

