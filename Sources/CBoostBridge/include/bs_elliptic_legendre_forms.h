#pragma once

#ifdef __cplusplus
extern "C" {
#endif

// Elliptic integrals (Legendre forms)
double bs_ellint_1_complete(double k);
double bs_ellint_1(double k, double phi);
double bs_ellint_2_complete(double k);
double bs_ellint_2(double k, double phi);
double bs_ellint_3(double k, double nu, double phi);
double bs_ellint_3_complete(double k, double nu);

float bs_ellint_1_complete_f(float k);
float bs_ellint_1_f(float k, float phi);
float bs_ellint_2_complete_f(float k);
float bs_ellint_2_f(float k, float phi);
float bs_ellint_3_f(float k, float nu, float phi);
float bs_ellint_3_complete_f(float k, float nu);

long double bs_ellint_1_complete_l(long double k);
long double bs_ellint_1_l(long double k, long double phi);
long double bs_ellint_2_complete_l(long double k);
long double bs_ellint_2_l(long double k, long double phi);
long double bs_ellint_3_l(long double k, long double nu, long double phi);
long double bs_ellint_3_complete_l(long double k, long double nu);

#ifdef __cplusplus
}
#endif

