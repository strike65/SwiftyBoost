#pragma once

#ifdef __cplusplus
extern "C" {
#endif

// Gegenbauer polynomials
double bs_gegenbauer(unsigned int n, double lambda, double x);
double bs_gegenbauer_prime(unsigned int n, double lambda, double x);
double bs_gegenbauer_derivative(unsigned int n, double lambda, double x, unsigned int k);
float bs_gegenbauer_f(unsigned int n, float lambda, float x);
float bs_gegenbauer_prime_f(unsigned int n, float lambda, float x);
float bs_gegenbauer_derivative_f(unsigned int n, float lambda, float x, unsigned int k);
long double bs_gegenbauer_l(unsigned int n, long double lambda, long double x);
long double bs_gegenbauer_prime_l(unsigned int n, long double lambda, long double x);
long double bs_gegenbauer_derivative_l(unsigned int n, long double lambda, long double x, unsigned int k);

#ifdef __cplusplus
}
#endif

