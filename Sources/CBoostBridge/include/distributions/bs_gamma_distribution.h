//
//  Created by Volker Thieme 2025.
//  Gamma distribution bridge (Boost.Math) — C ABI
//
#pragma once
#ifdef __cplusplus
extern "C" {
#endif

// C-friendly range PODs for [lower, upper]
typedef struct { float lower, upper; } bs_range_f;
typedef struct { double lower, upper; } bs_range_d;
typedef struct { long double lower, upper; } bs_range_l;

// Opaque-handle API (construct once, reuse)
// Lifetime
void* bs_gamma_make_f(float k, float theta);
void* bs_gamma_make(double k, double theta);
void* bs_gamma_make_l(long double k, long double theta);
void  bs_gamma_free_f(void* handle);
void  bs_gamma_free(void* handle);
void  bs_gamma_free_l(void* handle);

// Operations via handle
float       bs_gamma_cdf_f(const void* handle, float x);
double      bs_gamma_cdf(const void* handle, double x);
long double bs_gamma_cdf_l(const void* handle, long double x);

float       bs_gamma_ccdf_f(const void* handle, float x);
double      bs_gamma_ccdf(const void* handle, double x);
long double bs_gamma_ccdf_l(const void* handle, long double x);

/// MARK: hazard function

float       bs_gamma_hazard_f(const void* handle, float x);
double      bs_gamma_hazard(const void* handle, double x);
long double bs_gamma_hazard_l(const void* handle, long double x);

/// MARK: cumulative hazard function (CHF)
float       bs_gamma_chf_f(const void* handle, float x);
double      bs_gamma_chf(const void* handle, double x);
long double bs_gamma_chf_l(const void* handle, long double x);

float       bs_gamma_kurtosis_f(const void* handle);
double      bs_gamma_kurtosis(const void* handle);
long double bs_gamma_kurtosis_l(const void* handle);

float       bs_gamma_kurtosis_excess_f(const void* handle);
double      bs_gamma_kurtosis_excess(const void* handle);
long double bs_gamma_kurtosis_excess_l(const void* handle);

float       bs_gamma_mean_f(const void* handle);
double      bs_gamma_mean (const void* handle);
long double bs_gamma_mean_l(const void* handle);

float       bs_gamma_median_f(const void* handle);
double      bs_gamma_median (const void* handle);
long double bs_gamma_median_l(const void* handle);

float       bs_gamma_mode_f(const void* handle);
double      bs_gamma_mode(const void* handle);
long double bs_gamma_mode_l(const void* handle);

float       bs_gamma_pdf_f(const void* handle, float x);
double      bs_gamma_pdf(const void* handle, double x);
long double bs_gamma_pdf_l(const void* handle, long double x);

bs_range_f  bs_gamma_range_f(const void* handle);
bs_range_d  bs_gamma_range (const void* handle);
bs_range_l  bs_gamma_range_l(const void* handle);

float       bs_gamma_quantile_f(const void* handle, float p);
double      bs_gamma_quantile(const void* handle, double p);
long double bs_gamma_quantile_l(const void* handle, long double p);

float       bs_gamma_quantile_complement_f(const void* handle, float q);
double      bs_gamma_quantile_complement(const void* handle, double q);
long double bs_gamma_quantile_complement_l(const void* handle, long double q);

float       bs_gamma_skewness_f(const void* handle);
double      bs_gamma_skewness (const void* handle);
long double bs_gamma_skewness_l(const void* handle);

float       bs_gamma_standard_deviation_f(const void* handle);
double      bs_gamma_standard_deviation (const void* handle);
long double bs_gamma_standard_deviation_l(const void* handle);

bs_range_f  bs_gamma_support_f(const void* handle);
bs_range_d  bs_gamma_support (const void* handle);
bs_range_l  bs_gamma_support_l(const void* handle);

float       bs_gamma_variance_f(const void* handle);
double      bs_gamma_variance(const void* handle);
long double bs_gamma_variance_l(const void* handle);

float       bs_gamma_entropy_f(const void* handle);
double      bs_gamma_entropy(const void* handle);
long double bs_gamma_entropy_l(const void* handle);

#ifdef __cplusplus
}
#endif
