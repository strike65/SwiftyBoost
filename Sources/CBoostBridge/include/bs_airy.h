#pragma once

#ifdef __cplusplus
extern "C" {
#endif

// Airy functions
double bs_airy_ai(double x);
double bs_airy_bi(double x);
double bs_airy_ai_prime(double x);
double bs_airy_bi_prime(double x);

float bs_airy_ai_f(float x);
float bs_airy_bi_f(float x);
float bs_airy_ai_prime_f(float x);
float bs_airy_bi_prime_f(float x);

long double bs_airy_ai_l(long double x);
long double bs_airy_bi_l(long double x);
long double bs_airy_ai_prime_l(long double x);
long double bs_airy_bi_prime_l(long double x);

#ifdef __cplusplus
}
#endif

