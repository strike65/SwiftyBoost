//
//  Created by Volker Thieme 2025.
//  Gamma distribution bridge (Boost.Math) — C ABI
//
#pragma once

#ifdef __cplusplus
extern "C" {
#endif

// Parameterisierung: shape k (alpha) > 0, scale theta > 0
// Unterstützte Operationen für float/double/long double.

// PDF
float       bs_gamma_pdf_f(float k, float theta, float x);
double      bs_gamma_pdf(double k, double theta, double x);
long double bs_gamma_pdf_l(long double k, long double theta, long double x);

// CDF (lower tail)
float       bs_gamma_cdf_f(float k, float theta, float x);
double      bs_gamma_cdf(double k, double theta, double x);
long double bs_gamma_cdf_l(long double k, long double theta, long double x);

// CCDF (upper tail, complement)
float       bs_gamma_ccdf_f(float k, float theta, float x);
double      bs_gamma_ccdf(double k, double theta, double x);
long double bs_gamma_ccdf_l(long double k, long double theta, long double x);

// Quantile (lower tail)
float       bs_gamma_quantile_f(float k, float theta, float p);
double      bs_gamma_quantile(double k, double theta, double p);
long double bs_gamma_quantile_l(long double k, long double theta, long double p);

// Quantile complement (upper tail)
float       bs_gamma_quantile_complement_f(float k, float theta, float q);
double      bs_gamma_quantile_complement(double k, double theta, double q);
long double bs_gamma_quantile_complement_l(long double k, long double theta, long double q);

// Momente und Kenngrößen
float       bs_gamma_mean_f(float k, float theta);
double      bs_gamma_mean(double k, double theta);
long double bs_gamma_mean_l(long double k, long double theta);

float       bs_gamma_variance_f(float k, float theta);
double      bs_gamma_variance(double k, double theta);
long double bs_gamma_variance_l(long double k, long double theta);

float       bs_gamma_mode_f(float k, float theta);
double      bs_gamma_mode(double k, double theta);
long double bs_gamma_mode_l(long double k, long double theta);

float       bs_gamma_skewness_f(float k, float theta);
double      bs_gamma_skewness(double k, double theta);
long double bs_gamma_skewness_l(long double k, long double theta);

float       bs_gamma_kurtosis_excess_f(float k, float theta);
double      bs_gamma_kurtosis_excess(double k, double theta);
long double bs_gamma_kurtosis_excess_l(long double k, long double theta);

#ifdef __cplusplus
}
#endif

