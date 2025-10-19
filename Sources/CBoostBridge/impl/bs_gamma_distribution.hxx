//
//  Created by Volker Thieme 2025.
//  Gamma distribution bridge (Boost.Math) — implementation TU include
//
#include "../internal/bs_internal.hpp"
#include "../include/distributions/bs_gamma_distribution.h"

#include <boost/math/distributions/gamma.hpp>
#include <boost/math/distributions/complement.hpp>

using boost::math::gamma_distribution;
using boost::math::complement;

template <class Real>
static inline gamma_distribution<Real> make_gamma(Real k, Real theta)
{
    return gamma_distribution<Real>(k, theta);
}

// PDF
float bs_gamma_pdf_f(float k, float theta, float x) {
    return bs_wrap<float>([&]{ return pdf(make_gamma<float>(k, theta), x); });
}
double bs_gamma_pdf(double k, double theta, double x) {
    return bs_wrap<double>([&]{ return pdf(make_gamma<double>(k, theta), x); });
}
long double bs_gamma_pdf_l(long double k, long double theta, long double x) {
    return bs_wrap<long double>([&]{ return pdf(make_gamma<long double>(k, theta), x); });
}

// CDF (lower)
float bs_gamma_cdf_f(float k, float theta, float x) {
    return bs_wrap<float>([&]{ return cdf(make_gamma<float>(k, theta), x); });
}
double bs_gamma_cdf(double k, double theta, double x) {
    return bs_wrap<double>([&]{ return cdf(make_gamma<double>(k, theta), x); });
}
long double bs_gamma_cdf_l(long double k, long double theta, long double x) {
    return bs_wrap<long double>([&]{ return cdf(make_gamma<long double>(k, theta), x); });
}

// CCDF (upper)
float bs_gamma_ccdf_f(float k, float theta, float x) {
    return bs_wrap<float>([&]{ return cdf(complement(make_gamma<float>(k, theta), x)); });
}
double bs_gamma_ccdf(double k, double theta, double x) {
    return bs_wrap<double>([&]{ return cdf(complement(make_gamma<double>(k, theta), x)); });
}
long double bs_gamma_ccdf_l(long double k, long double theta, long double x) {
    return bs_wrap<long double>([&]{ return cdf(complement(make_gamma<long double>(k, theta), x)); });
}

// Quantiles (lower)
float bs_gamma_quantile_f(float k, float theta, float p) {
    return bs_wrap<float>([&]{ return quantile(make_gamma<float>(k, theta), p); });
}
double bs_gamma_quantile(double k, double theta, double p) {
    return bs_wrap<double>([&]{ return quantile(make_gamma<double>(k, theta), p); });
}
long double bs_gamma_quantile_l(long double k, long double theta, long double p) {
    return bs_wrap<long double>([&]{ return quantile(make_gamma<long double>(k, theta), p); });
}

// Quantiles (upper, complement)
float bs_gamma_quantile_complement_f(float k, float theta, float q) {
    return bs_wrap<float>([&]{ return quantile(complement(make_gamma<float>(k, theta), q)); });
}
double bs_gamma_quantile_complement(double k, double theta, double q) {
    return bs_wrap<double>([&]{ return quantile(complement(make_gamma<double>(k, theta), q)); });
}
long double bs_gamma_quantile_complement_l(long double k, long double theta, long double q) {
    return bs_wrap<long double>([&]{ return quantile(complement(make_gamma<long double>(k, theta), q)); });
}

// Moments/characteristics
float bs_gamma_mean_f(float k, float theta) {
    return bs_wrap<float>([&]{ return mean(make_gamma<float>(k, theta)); });
}
double bs_gamma_mean(double k, double theta) {
    return bs_wrap<double>([&]{ return mean(make_gamma<double>(k, theta)); });
}
long double bs_gamma_mean_l(long double k, long double theta) {
    return bs_wrap<long double>([&]{ return mean(make_gamma<long double>(k, theta)); });
}

float bs_gamma_variance_f(float k, float theta) {
    return bs_wrap<float>([&]{ return variance(make_gamma<float>(k, theta)); });
}
double bs_gamma_variance(double k, double theta) {
    return bs_wrap<double>([&]{ return variance(make_gamma<double>(k, theta)); });
}
long double bs_gamma_variance_l(long double k, long double theta) {
    return bs_wrap<long double>([&]{ return variance(make_gamma<long double>(k, theta)); });
}

float bs_gamma_mode_f(float k, float theta) {
    return bs_wrap<float>([&]{ return mode(make_gamma<float>(k, theta)); });
}
double bs_gamma_mode(double k, double theta) {
    return bs_wrap<double>([&]{ return mode(make_gamma<double>(k, theta)); });
}
long double bs_gamma_mode_l(long double k, long double theta) {
    return bs_wrap<long double>([&]{ return mode(make_gamma<long double>(k, theta)); });
}

float bs_gamma_skewness_f(float k, float theta) {
    (void)theta;
    return bs_wrap<float>([&]{ return skewness(make_gamma<float>(k, theta)); });
}
double bs_gamma_skewness(double k, double theta) {
    (void)theta;
    return bs_wrap<double>([&]{ return skewness(make_gamma<double>(k, theta)); });
}
long double bs_gamma_skewness_l(long double k, long double theta) {
    (void)theta;
    return bs_wrap<long double>([&]{ return skewness(make_gamma<long double>(k, theta)); });
}

float bs_gamma_kurtosis_excess_f(float k, float theta) {
    (void)theta;
    return bs_wrap<float>([&]{ return kurtosis_excess(make_gamma<float>(k, theta)); });
}
double bs_gamma_kurtosis_excess(double k, double theta) {
    (void)theta;
    return bs_wrap<double>([&]{ return kurtosis_excess(make_gamma<double>(k, theta)); });
}
long double bs_gamma_kurtosis_excess_l(long double k, long double theta) {
    (void)theta;
    return bs_wrap<long double>([&]{ return kurtosis_excess(make_gamma<long double>(k, theta)); });
}

