#pragma once

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

// Hypergeometric functions (Gauss/confluent/general)
double bs_hypergeometric_1F0(double a, double z);
float bs_hypergeometric_1F0_f(float a, float z);
long double bs_hypergeometric_1F0_l(long double a, long double z);

double bs_hypergeometric_0F1(double b, double z);
float bs_hypergeometric_0F1_f(float b, float z);
long double bs_hypergeometric_0F1_l(long double b, long double z);

double bs_hypergeometric_2F0(double a, double b, double z);
float bs_hypergeometric_2F0_f(float a, float b, float z);
long double bs_hypergeometric_2F0_l(long double a, long double b, long double z);

double bs_hypergeometric_1F1(double a, double b, double z);
float bs_hypergeometric_1F1_f(float a, float b, float z);
long double bs_hypergeometric_1F1_l(long double a, long double b, long double z);

// General pFq with arrays of parameters
double bs_hypergeometric_pFq(const double *a, size_t p, const double *b, size_t q, double z);
float bs_hypergeometric_pFq_f(const float *a, size_t p, const float *b, size_t q, float z);
long double bs_hypergeometric_pFq_l(const long double *a, size_t p, const long double *b, size_t q, long double z);

#ifdef __cplusplus
}
#endif

