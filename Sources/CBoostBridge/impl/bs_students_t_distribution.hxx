//
//  Created by Volker Thieme 2025.
//  Students t distribution bridge (Boost.Math) — implementation TU include
//
#include "../internal/bs_internal.hpp"
#include "../include/distributions/bs_students_t_distribution.h"

#include <boost/math/distributions/students_t.hpp>
#include <boost/math/distributions/complement.hpp>

using boost::math::students_t_distribution;
using boost::math::complement;

template <class Real>
static inline students_t_distribution<Real> make_students_t(Real v)
{
    return students_t_distribution<Real>(v);
}

// PDF
float bs_students_t_pdf_f(float v, float x) {
    return bs_wrap<float>([&]{ return pdf(make_students_t<float>(v), x); });
}
double bs_students_t_pdf(double v, double x) {
    return bs_wrap<double>([&]{ return pdf(make_students_t<double>(v), x); });
}
long double bs_students_t_pdf_l(long double v, long double x) {
    return bs_wrap<long double>([&]{ return pdf(make_students_t<long double>(v), x); });
}

// CDF (lower)
float bs_students_t_cdf_f(float v, float x) {
    return bs_wrap<float>([&]{ return cdf(make_students_t<float>(v), x); });
}
double bs_students_t_cdf(double v, double x) {
    return bs_wrap<double>([&]{ return cdf(make_students_t<double>(v), x); });
}
long double bs_students_t_cdf_l(long double v, long double x) {
    return bs_wrap<long double>([&]{ return cdf(make_students_t<long double>(v), x); });
}

// CCDF (upper)
float bs_students_t_ccdf_f(float v, float x) {
    return bs_wrap<float>([&]{ return cdf(complement(make_students_t<float>(v), x)); });
}
double bs_students_t_ccdf(double v, double x) {
    return bs_wrap<double>([&]{ return cdf(complement(make_students_t<double>(v), x)); });
}
long double bs_students_t_ccdf_l(long double v, long double x) {
    return bs_wrap<long double>([&]{ return cdf(complement(make_students_t<long double>(v), x)); });
}

// Quantiles (lower)
float bs_students_t_quantile_f(float v, float p) {
    return bs_wrap<float>([&]{ return quantile(make_students_t<float>(v), p); });
}
double bs_students_t_quantile(double v, double p) {
    return bs_wrap<double>([&]{ return quantile(make_students_t<double>(v), p); });
}
long double bs_students_t_quantile_l(long double v, long double p) {
    return bs_wrap<long double>([&]{ return quantile(make_students_t<long double>(v), p); });
}

// Quantiles (upper, complement)
float bs_students_t_quantile_complement_f(float v, float q) {
    return bs_wrap<float>([&]{ return quantile(complement(make_students_t<float>(v), q)); });
}
double bs_students_t_quantile_complement(double v, double q) {
    return bs_wrap<double>([&]{ return quantile(complement(make_students_t<double>(v), q)); });
}
long double bs_students_t_quantile_complement_l(long double v, long double q) {
    return bs_wrap<long double>([&]{ return quantile(complement(make_students_t<long double>(v), q)); });
}

// Moments / characteristics
float bs_students_t_mean_f(float v) {
    return bs_wrap<float>([&]{ return mean(make_students_t<float>(v)); });
}
double bs_students_t_mean(double v) {
    return bs_wrap<double>([&]{ return mean(make_students_t<double>(v)); });
}
long double bs_students_t_mean_l(long double v) {
    return bs_wrap<long double>([&]{ return mean(make_students_t<long double>(v)); });
}

float bs_students_t_variance_f(float v) {
    return bs_wrap<float>([&]{ return variance(make_students_t<float>(v)); });
}
double bs_students_t_variance(double v) {
    return bs_wrap<double>([&]{ return variance(make_students_t<double>(v)); });
}
long double bs_students_t_variance_l(long double v) {
    return bs_wrap<long double>([&]{ return variance(make_students_t<long double>(v)); });
}

float bs_students_t_mode_f(float v) {
    (void)v;
    return 0.0f;
}
double bs_students_t_mode(double v) {
    (void)v;
    return 0.0;
}
long double bs_students_t_mode_l(long double v) {
    (void)v;
    return static_cast<long double>(0.0L);
}

float bs_students_t_kurtosis_excess_f(float v) {
    return bs_wrap<float>([&]{ return kurtosis_excess(make_students_t<float>(v)); });
}
double bs_students_t_kurtosis_excess(double v) {
    return bs_wrap<double>([&]{ return kurtosis_excess(make_students_t<double>(v)); });
}
long double bs_students_t_kurtosis_excess_l(long double v) {
    return bs_wrap<long double>([&]{ return kurtosis_excess(make_students_t<long double>(v)); });
}

// Static utility: required degrees of freedom
float bs_students_t_find_degrees_of_freedom_f(float difference_from_mean, float alpha, float beta, float sd, float hint) {
    return bs_wrap<float>([&]{ return students_t_distribution<float>::find_degrees_of_freedom(difference_from_mean, alpha, beta, sd, hint); });
}
double bs_students_t_find_degrees_of_freedom(double difference_from_mean, double alpha, double beta, double sd, double hint) {
    return bs_wrap<double>([&]{ return students_t_distribution<double>::find_degrees_of_freedom(difference_from_mean, alpha, beta, sd, hint); });
}
long double bs_students_t_find_degrees_of_freedom_l(long double difference_from_mean, long double alpha, long double beta, long double sd, long double hint) {
    return bs_wrap<long double>([&]{ return students_t_distribution<long double>::find_degrees_of_freedom(difference_from_mean, alpha, beta, sd, hint); });
}

