#pragma once

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

// Bernoulli Numbers
double bs_bernoulli_b2n(const int n);
float bs_bernoulli_b2n_f(const int n);
long double bs_bernoulli_b2n_l(const int n);

// Tangent Numbers (scalar)
double bs_tangent_t2n(const int n);
float bs_tangent_t2n_f(const int n);
long double bs_tangent_t2n_l(const int n);

// Tangent Numbers (bulk sequence)
void bs_tangent_t2n_seq(int start_index, unsigned int count, double *out);
void bs_tangent_t2n_seq_f(int start_index, unsigned int count, float *out);
void bs_tangent_t2n_seq_l(int start_index, unsigned int count, long double *out);

// prime numbers up to 10000th
unsigned int bs_prime(unsigned int n);

// Fibonacci numbers
unsigned long long bs_fibonacci_ull(unsigned long long n);

#ifdef __cplusplus
}
#endif

