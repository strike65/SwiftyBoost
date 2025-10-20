//
//  Created by Volker Thieme 2025.
//  Students t distribution bridge (Boost.Math) — C ABI
//
#pragma once

#ifdef __cplusplus
extern "C" {
#endif

// Opaque-handle API (construct once, reuse)
// Lifetime
void* bs_arcsine_make_f(float x_min, float x_max);
void* bs_arcsine_make_d(double x_min, double x_max);
void* bs_arcsine_make_l(long double x_min, long double x_max);
void  bs_arcsine_free_f(void* handle);
void  bs_arcsine_free_d(void* handle);
void  bs_arcsine_free_l(void* handle);

// Operations via handle
float       bs_arcsine_cdf_f(const void* handle, float x);
double      bs_arcsine_cdf_d(const void* handle, double x);
long double bs_arcsine_cdf_l(const void* handle, long double x);

float       bs_arcsine_ccdf_f(const void* handle, float x);
double      bs_arcsine_ccdf_d(const void* handle, double x);
long double bs_arcsine_ccdf_l(const void* handle, long double x);

/// MARK: hazard function

float       bs_arcsine_hazard_f(const void* handle, float x);
double      bs_arcsine_hazard_d(const void* handle, double x);
long double bs_arcsine_hazard_l(const void* handle, long double x);

/// MARK: cumulative hazard function (CHF)
float       bs_arcsine_chf_f(const void* handle, float x);
double      bs_arcsine_chf_d(const void* handle, double x);
long double bs_arcsine_chf_l(const void* handle, long double x);

float       bs_arcsine_kurtosis_f(const void* handle);
double      bs_arcsine_kurtosis_d(const void* handle);
long double bs_arcsine_kurtosis_l(const void* handle);

float       bs_arcsine_kurtosis_excess_f(const void* handle);
double      bs_arcsine_kurtosis_excess_d(const void* handle);
long double bs_arcsine_kurtosis_excess_l(const void* handle);

float       bs_arcsine_mean_f(const void* handle);
double      bs_arcsine_mean_d(const void* handle);
long double bs_arcsine_mean_l(const void* handle);

float       bs_arcsine_median_f(const void* handle);
double      bs_arcsine_median_d(const void* handle);
long double bs_arcsine_median_l(const void* handle);

float       bs_arcsine_mode_f(const void* handle);
double      bs_arcsine_mode_d(const void* handle);
long double bs_arcsine_mode_l(const void* handle);

float       bs_arcsine_pdf_f(const void* handle, float x);
double      bs_arcsine_pdf_d(const void* handle, double x);
long double bs_arcsine_pdf_l(const void* handle, long double x);

bs_range_f  bs_arcsine_range_f(const void* handle);
bs_range_d  bs_arcsine_range_d(const void* handle);
bs_range_l  bs_arcsine_range_l(const void* handle);

float       bs_arcsine_quantile_f(const void* handle, float p);
double      bs_arcsine_quantile_d(const void* handle, double p);
long double bs_arcsine_quantile_l(const void* handle, long double p);

float       bs_arcsine_quantile_complement_f(const void* handle, float q);
double      bs_arcsine_quantile_complement_d(const void* handle, double q);
long double bs_arcsine_quantile_complement_l(const void* handle, long double q);

float       bs_arcsine_skewness_f(const void* handle);
double      bs_arcsine_skewness_d(const void* handle);
long double bs_arcsine_skewness_l(const void* handle);

float       bs_arcsine_standard_deviation_f(const void* handle);
double      bs_arcsine_standard_deviation_d(const void* handle);
long double bs_arcsine_standard_deviation_l(const void* handle);

bs_range_f  bs_arcsine_support_f(const void* handle);
bs_range_d  bs_arcsine_support_d(const void* handle);
bs_range_l  bs_arcsine_support_l(const void* handle);

float       bs_arcsine_variance_f(const void* handle);
double      bs_arcsine_variance_d(const void* handle);
long double bs_arcsine_variance_l(const void* handle);

float       bs_arcsine_entropy_f(const void* handle);
double      bs_arcsine_entropy_d(const void* handle);
long double bs_arcsine_entropy_l(const void* handle);

#ifdef __cplusplus
}
#endif
