#pragma once

#ifdef __cplusplus
extern "C" {
#endif

// Digamma / Polygamma / Zeta
double bs_digamma(double x);
double bs_trigamma(double x);
double bs_polygamma(int n, double x);
double bs_riemann_zeta(double x);

float bs_digamma_f(float x);
float bs_trigamma_f(float x);
float bs_polygamma_f(int n, float x);
float bs_riemann_zeta_f(float x);

long double bs_digamma_l(long double x);
long double bs_trigamma_l(long double x);
long double bs_polygamma_l(int n, long double x);
long double bs_riemann_zeta_l(long double x);

#ifdef __cplusplus
}
#endif

