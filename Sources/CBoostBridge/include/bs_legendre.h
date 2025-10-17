#pragma once

#ifdef __cplusplus
extern "C" {
#endif

// Legendre
double bs_legendre_p(int n, double x);
double bs_assoc_legendre_p(int n, int m, double x);
double bs_legendre_p_prime(int n, double x);
void bs_legendre_p_zeros(int l, double *out);

float bs_legendre_p_f(int n, float x);
float bs_assoc_legendre_p_f(int n, int m, float x);
float bs_legendre_p_prime_f(int n, float x);
void bs_legendre_p_zeros_f(int l, float *out);

long double bs_legendre_p_l(int n, long double x);
long double bs_assoc_legendre_p_l(int n, int m, long double x);
long double bs_legendre_p_prime_l(int n, long double x);
void bs_legendre_p_zeros_l(int l, long double *out);

#ifdef __cplusplus
}
#endif

