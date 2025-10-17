#pragma once

#ifdef __cplusplus
extern "C" {
#endif

// Laguerre
double bs_laguerre(unsigned int n, double x);
float bs_laguerre_f(unsigned int n, float x);
long double bs_laguerre_l(unsigned int n, long double x);

double bs_assoc_laguerre(unsigned int n, unsigned int m, double x);
float bs_assoc_laguerre_f(unsigned int n, unsigned int m, float x);
long double bs_assoc_laguerre_l(unsigned int n, unsigned int m, long double x);

#ifdef __cplusplus
}
#endif

