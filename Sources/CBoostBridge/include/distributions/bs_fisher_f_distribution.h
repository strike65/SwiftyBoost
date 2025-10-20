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
void* bs_fisher_f_make_f(float df1, float df2);
void* bs_fisher_f_make_d(double df1, double df2);
void* bs_fisher_f_make_l(long double df1, long double df2);
void  bs_fisher_f_free_f(void* handle);
void  bs_fisher_f_free_d(void* handle);
void  bs_fisher_f_free_l(void* handle);

// Operations via handle
float       bs_fisher_f_cdf_f(const void* handle, float x);
double      bs_fisher_f_cdf_d(const void* handle, double x);
long double bs_fisher_f_cdf_l(const void* handle, long double x);

float       bs_fisher_f_ccdf_f(const void* handle, float x);
double      bs_fisher_f_ccdf_d(const void* handle, double x);
long double bs_fisher_f_ccdf_l(const void* handle, long double x);

/// MARK: hazard function

float       bs_fisher_f_hazard_f(const void* handle, float x);
double      bs_fisher_f_hazard_d(const void* handle, double x);
long double bs_fisher_f_hazard_l(const void* handle, long double x);

/// MARK: cumulative hazard function (CHF)
float       bs_fisher_f_chf_f(const void* handle, float x);
double      bs_fisher_f_chf_d(const void* handle, double x);
long double bs_fisher_f_chf_l(const void* handle, long double x);

float       bs_fisher_f_kurtosis_f(const void* handle);
double      bs_fisher_f_kurtosis_d(const void* handle);
long double bs_fisher_f_kurtosis_l(const void* handle);

float       bs_fisher_f_kurtosis_excess_f(const void* handle);
double      bs_fisher_f_kurtosis_excess_d(const void* handle);
long double bs_fisher_f_kurtosis_excess_l(const void* handle);

float       bs_fisher_f_mean_f(const void* handle);
double      bs_fisher_f_mean_d(const void* handle);
long double bs_fisher_f_mean_l(const void* handle);

float       bs_fisher_f_median_f(const void* handle);
double      bs_fisher_f_median_d(const void* handle);
long double bs_fisher_f_median_l(const void* handle);

float       bs_fisher_f_mode_f(const void* handle);
double      bs_fisher_f_mode_d(const void* handle);
long double bs_fisher_f_mode_l(const void* handle);

float       bs_fisher_f_pdf_f(const void* handle, float x);
double      bs_fisher_f_pdf_d(const void* handle, double x);
long double bs_fisher_f_pdf_l(const void* handle, long double x);

bs_range_f  bs_fisher_f_range_f(const void* handle);
bs_range_d  bs_fisher_f_range_d(const void* handle);
bs_range_l  bs_fisher_f_range_l(const void* handle);

float       bs_fisher_f_quantile_f(const void* handle, float p);
double      bs_fisher_f_quantile_d(const void* handle, double p);
long double bs_fisher_f_quantile_l(const void* handle, long double p);

float       bs_fisher_f_quantile_complement_f(const void* handle, float q);
double      bs_fisher_f_quantile_complement_d(const void* handle, double q);
long double bs_fisher_f_quantile_complement_l(const void* handle, long double q);

float       bs_fisher_f_skewness_f(const void* handle);
double      bs_fisher_f_skewness_d(const void* handle);
long double bs_fisher_f_skewness_l(const void* handle);

float       bs_fisher_f_standard_deviation_f(const void* handle);
double      bs_fisher_f_standard_deviation_d(const void* handle);
long double bs_fisher_f_standard_deviation_l(const void* handle);

bs_range_f  bs_fisher_f_support_f(const void* handle);
bs_range_d  bs_fisher_f_support_d(const void* handle);
bs_range_l  bs_fisher_f_support_l(const void* handle);

float       bs_fisher_f_variance_f(const void* handle);
double      bs_fisher_f_variance_d(const void* handle);
long double bs_fisher_f_variance_l(const void* handle);

#ifdef __cplusplus
}
#endif
