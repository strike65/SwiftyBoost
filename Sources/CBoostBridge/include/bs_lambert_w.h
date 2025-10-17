#pragma once

#ifdef __cplusplus
extern "C" {
#endif

// Lambert W
double bs_lambert_w0(double x);
double bs_lambert_wm1(double x);

float bs_lambert_w0_f(float x);
float bs_lambert_wm1_f(float x);

long double bs_lambert_w0_l(long double x);
long double bs_lambert_wm1_l(long double x);

#ifdef __cplusplus
}
#endif

