#pragma once

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

// Exponential integrals and related helpers
double bs_expint_Ei(double x);
double bs_expint_En(int n, double x);
double bs_expm1(double x);
double bs_log1p(double x);
double bs_log1pmx(double x);
double bs_powm1(double x, double y);
double bs_cbrt(double x);
double bs_sqrt1pm1(double x);
double bs_hypot(double x, double y);

float bs_expint_Ei_f(float x);
float bs_expint_En_f(int n, float x);
float bs_expm1_f(float x);
float bs_log1p_f(float x);
float bs_log1pmx_f(float x);
float bs_powm1_f(float x, float y);
float bs_cbrt_f(float x);
float bs_sqrt1pm1_f(float x);
float bs_hypot_f(float x, float y);

long double bs_expint_Ei_l(long double x);
long double bs_expint_En_l(int n, long double x);
long double bs_expm1_l(long double x);
long double bs_log1p_l(long double x);
long double bs_log1pmx_l(long double x);
long double bs_powm1_l(long double x, long double y);
long double bs_cbrt_l(long double x);
long double bs_sqrt1pm1_l(long double x);
long double bs_hypot_l(long double x, long double y);

#ifdef __cplusplus
}
#endif

