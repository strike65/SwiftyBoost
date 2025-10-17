#pragma once

#ifdef __cplusplus
extern "C" {
#endif

// Elliptic integrals (Carlson symmetric forms)
double bs_ellint_rc(double x, double y);
double bs_ellint_rf(double x, double y, double z);
double bs_ellint_rd(double x, double y, double z);
double bs_ellint_rj(double x, double y, double z, double p);
double bs_ellint_rg(double x, double y, double z);

float bs_ellint_rc_f(float x, float y);
float bs_ellint_rf_f(float x, float y, float z);
float bs_ellint_rd_f(float x, float y, float z);
float bs_ellint_rj_f(float x, float y, float z, float p);
float bs_ellint_rg_f(float x, float y, float z);

long double bs_ellint_rc_l(long double x, long double y);
long double bs_ellint_rf_l(long double x, long double y, long double z);
long double bs_ellint_rd_l(long double x, long double y, long double z);
long double bs_ellint_rj_l(long double x, long double y, long double z, long double p);
long double bs_ellint_rg_l(long double x, long double y, long double z);

#ifdef __cplusplus
}
#endif

