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
void* bs_student_t_make_f(float v);
void* bs_student_t_make(double v);
void* bs_student_t_make_l(long double v);
void  bs_student_t_free_f(void* handle);
void  bs_student_t_free(void* handle);
void  bs_student_t_free_l(void* handle);

// Operations via handle
float       bs_student_t_cdf_f(const void* handle, float x);
double      bs_student_t_cdf(const void* handle, double x);
long double bs_student_t_cdf_l(const void* handle, long double x);

float       bs_student_t_ccdf_f(const void* handle, float x);
double      bs_student_t_ccdf(const void* handle, double x);
long double bs_student_t_ccdf_l(const void* handle, long double x);

/// MARK: hazard function

float       bs_student_t_hazard_f(const void* handle, float x);
double      bs_student_t_hazard(const void* handle, double x);
long double bs_student_t_hazard_l(const void* handle, long double x);

/// MARK: cumulative hazard function (CHF)
float       bs_student_t_chf_f(const void* handle, float x);
double      bs_student_t_chf(const void* handle, double x);
long double bs_student_t_chf_l(const void* handle, long double x);

float       bs_student_t_kurtosis_f(const void* handle);
double      bs_student_t_kurtosis(const void* handle);
long double bs_student_t_kurtosis_l(const void* handle);

float       bs_student_t_kurtosis_excess_f(const void* handle);
double      bs_student_t_kurtosis_excess(const void* handle);
long double bs_student_t_kurtosis_excess_l(const void* handle);

float       bs_student_t_mean_f(const void* handle);
double      bs_student_t_mean (const void* handle);
long double bs_student_t_mean_l(const void* handle);

float       bs_student_t_median_f(const void* handle);
double      bs_student_t_median (const void* handle);
long double bs_student_t_median_l(const void* handle);

float       bs_student_t_mode_f(const void* handle);
double      bs_student_t_mode(const void* handle);
long double bs_student_t_mode_l(const void* handle);

float       bs_student_t_pdf_f(const void* handle, float x);
double      bs_student_t_pdf(const void* handle, double x);
long double bs_student_t_pdf_l(const void* handle, long double x);

bs_range_f  bs_student_t_range_f(const void* handle);
bs_range_d  bs_student_t_range (const void* handle);
bs_range_l  bs_student_t_range_l(const void* handle);

float       bs_student_t_quantile_f(const void* handle, float p);
double      bs_student_t_quantile(const void* handle, double p);
long double bs_student_t_quantile_l(const void* handle, long double p);

float       bs_student_t_quantile_complement_f(const void* handle, float q);
double      bs_student_t_quantile_complement(const void* handle, double q);
long double bs_student_t_quantile_complement_l(const void* handle, long double q);

float       bs_student_t_skewness_f(const void* handle);
double      bs_student_t_skewness (const void* handle);
long double bs_student_t_skewness_l(const void* handle);

float       bs_student_t_standard_deviation_f(const void* handle);
double      bs_student_t_standard_deviation (const void* handle);
long double bs_student_t_standard_deviation_l(const void* handle);

bs_range_f  bs_student_t_support_f(const void* handle);
bs_range_d  bs_student_t_support (const void* handle);
bs_range_l  bs_student_t_support_l(const void* handle);

float       bs_student_t_variance_f(const void* handle);
double      bs_student_t_variance(const void* handle);
long double bs_student_t_variance_l(const void* handle);

float       bs_student_t_entropy_f(const void* handle);
double      bs_student_t_entropy(const void* handle);
long double bs_student_t_entropy_l(const void* handle);


// Planning utility: find required degrees of freedom.
// Mirrors boost::math::students_t_distribution<>::find_degrees_of_freedom.
// Parameters:
//  - difference_from_mean: effect size (difference in means)
//  - alpha: significance level (Type I error)
//  - beta: Type II error (1 - power)
//  - sd: standard deviation
//  - hint: initial guess for df (<= 0 -> default 1)
float       bs_student_t_find_degrees_of_freedom_f(float difference_from_mean, float alpha, float beta, float sd, float hint);
double      bs_student_t_find_degrees_of_freedom(double difference_from_mean, double alpha, double beta, double sd, double hint);
long double bs_student_t_find_degrees_of_freedom_l(long double difference_from_mean, long double alpha, long double beta, long double sd, long double hint);

#ifdef __cplusplus
}
#endif
