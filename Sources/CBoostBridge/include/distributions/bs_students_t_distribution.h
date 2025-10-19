//
//  Created by Volker Thieme 2025.
//  Students t distribution bridge (Boost.Math) — C ABI
//
#pragma once

#ifdef __cplusplus
extern "C" {
#endif

// Students t distribution with degrees of freedom v > 0.
// Supported operations for float/double/long double.

// PDF
float       bs_students_t_pdf_f(float v, float x);
double      bs_students_t_pdf(double v, double x);
long double bs_students_t_pdf_l(long double v, long double x);

// CDF (lower tail)
float       bs_students_t_cdf_f(float v, float x);
double      bs_students_t_cdf(double v, double x);
long double bs_students_t_cdf_l(long double v, long double x);

// CCDF (upper tail, complement)
float       bs_students_t_ccdf_f(float v, float x);
double      bs_students_t_ccdf(double v, double x);
long double bs_students_t_ccdf_l(long double v, long double x);

// Quantile (lower tail)
float       bs_students_t_quantile_f(float v, float p);
double      bs_students_t_quantile(double v, double p);
long double bs_students_t_quantile_l(long double v, long double p);

// Quantile complement (upper tail)
float       bs_students_t_quantile_complement_f(float v, float q);
double      bs_students_t_quantile_complement(double v, double q);
long double bs_students_t_quantile_complement_l(long double v, long double q);

// Summary stats
float       bs_students_t_mean_f(float v);
double      bs_students_t_mean(double v);
long double bs_students_t_mean_l(long double v);

float       bs_students_t_variance_f(float v);
double      bs_students_t_variance(double v);
long double bs_students_t_variance_l(long double v);

float       bs_students_t_mode_f(float v);
double      bs_students_t_mode(double v);
long double bs_students_t_mode_l(long double v);

float       bs_students_t_kurtosis_excess_f(float v);
double      bs_students_t_kurtosis_excess(double v);
long double bs_students_t_kurtosis_excess_l(long double v);

// Planning utility: find required degrees of freedom.
// Mirrors boost::math::students_t_distribution<>::find_degrees_of_freedom.
// Parameters:
//  - difference_from_mean: effect size (difference in means)
//  - alpha: significance level (Type I error)
//  - beta: Type II error (1 - power)
//  - sd: standard deviation
//  - hint: initial guess for df (<= 0 -> default 1)
float       bs_students_t_find_degrees_of_freedom_f(float difference_from_mean, float alpha, float beta, float sd, float hint);
double      bs_students_t_find_degrees_of_freedom(double difference_from_mean, double alpha, double beta, double sd, double hint);
long double bs_students_t_find_degrees_of_freedom_l(long double difference_from_mean, long double alpha, long double beta, long double sd, long double hint);

#ifdef __cplusplus
}
#endif

