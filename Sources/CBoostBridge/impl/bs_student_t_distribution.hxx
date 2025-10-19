//
//  Created by Volker Thieme 2025.
//  Students t distribution bridge (Boost.Math) — implementation TU include
//
#include "../internal/bs_internal.hpp"
#include "../include/distributions/bs_student_t_distribution.h"

#include <boost/math/distributions/students_t.hpp>
#include <boost/math/distributions/complement.hpp>

using boost::math::students_t_distribution;
using boost::math::complement;

// -----------------------
// Handle-based API
// -----------------------

struct bs_student_t_f_handle { students_t_distribution<float> dist; explicit bs_student_t_f_handle(float v) : dist(v) {} };
struct bs_student_t_d_handle { students_t_distribution<double> dist; explicit bs_student_t_d_handle(double v) : dist(v) {} };
struct bs_student_t_l_handle { students_t_distribution<long double> dist; explicit bs_student_t_l_handle(long double v) : dist(v) {} };

extern "C" {
void* bs_student_t_make_f(float v) {
    try { return new bs_student_t_f_handle(v); } catch (...) { return nullptr; }
}
void* bs_student_t_make(double v) {
    try { return new bs_student_t_d_handle(v); } catch (...) { return nullptr; }
}
void* bs_student_t_make_l(long double v) {
    try { return new bs_student_t_l_handle(v); } catch (...) { return nullptr; }
}
void bs_student_t_free_f(void* handle) { delete static_cast<bs_student_t_f_handle*>(handle); }
void bs_student_t_free(void* handle) { delete static_cast<bs_student_t_d_handle*>(handle); }
void bs_student_t_free_l(void* handle) { delete static_cast<bs_student_t_l_handle*>(handle); }

float bs_student_t_cdf_f(const void* handle, float x) {
    auto h = static_cast<const bs_student_t_f_handle*>(handle);
    if (!h) return std::numeric_limits<float>::quiet_NaN();
    return bs_wrap<float>([&]{ return cdf(h->dist, x); });
}
double bs_student_t_cdf(const void* handle, double x) {
    auto h = static_cast<const bs_student_t_d_handle*>(handle);
    if (!h) return std::numeric_limits<double>::quiet_NaN();
    return bs_wrap<double>([&]{ return cdf(h->dist, x); });
}
long double bs_student_t_cdf_l(const void* handle, long double x) {
    auto h = static_cast<const bs_student_t_l_handle*>(handle);
    if (!h) return std::numeric_limits<long double>::quiet_NaN();
    return bs_wrap<long double>([&]{ return cdf(h->dist, x); });
}

float bs_student_t_ccdf_f(const void* handle, float x) {
    auto h = static_cast<const bs_student_t_f_handle*>(handle);
    if (!h) return std::numeric_limits<float>::quiet_NaN();
    return bs_wrap<float>([&]{ return cdf(complement(h->dist, x)); });
}
double bs_student_t_ccdf(const void* handle, double x) {
    auto h = static_cast<const bs_student_t_d_handle*>(handle);
    if (!h) return std::numeric_limits<double>::quiet_NaN();
    return bs_wrap<double>([&]{ return cdf(complement(h->dist, x)); });
}
long double bs_student_t_ccdf_l(const void* handle, long double x) {
    auto h = static_cast<const bs_student_t_l_handle*>(handle);
    if (!h) return std::numeric_limits<long double>::quiet_NaN();
    return bs_wrap<long double>([&]{ return cdf(complement(h->dist, x)); });
}
/// MARK: hazard function

float bs_student_t_hazard_f(const void* handle, float x) {
    auto h = static_cast<const bs_student_t_f_handle*>(handle);
    if (!h) return std::numeric_limits<float>::quiet_NaN();
    float p = bs_wrap<float>([&]{ return cdf(complement(h->dist, x)); });
    float d = bs_wrap<float>([&]{ return pdf(h -> dist, x); });
    if (d > p * boost::math::tools::max_value<float>()) {
        return std::numeric_limits<float>::quiet_NaN();
    }
    if (d == 0) {
        return 0;
    }
    return d / p;
}
double bs_student_t_hazard(const void* handle, double x) {
    auto h = static_cast<const bs_student_t_d_handle*>(handle);
    if (!h) return std::numeric_limits<double>::quiet_NaN();
    double p = bs_wrap<double>([&]{ return cdf(complement(h->dist, x)); });
    double d = bs_wrap<double>([&]{ return pdf(h -> dist, x); });
    if (d > p * boost::math::tools::max_value<double>()) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (d == 0) {
        return 0;
    }
    return d / p;
}
long double bs_student_t_hazard_l(const void* handle, long double x)  {
    auto h = static_cast<const bs_student_t_l_handle*>(handle);
    if (!h) return std::numeric_limits<long double>::quiet_NaN();
    long double p = bs_wrap<long double>([&]{ return cdf(complement(h->dist, x)); });
    long double d = bs_wrap<long double>([&]{ return pdf(h -> dist, x); });
    if (d > p * boost::math::tools::max_value<long double>()) {
        return std::numeric_limits<long double>::quiet_NaN();
    }
    if (d == 0) {
        return 0;
    }
    return d / p;
}

/// MARK: cumulative hazard function
float bs_student_t_chf_f(const void* handle, float x) {
    auto h = static_cast<const bs_student_t_f_handle*>(handle);
    if (!h) return std::numeric_limits<float>::quiet_NaN();
    float p = bs_wrap<float>([&]{ return cdf(complement(h->dist, x)); });
    return -log(p);
}
double     bs_student_t_chf(const void* handle, double x) {
    auto h = static_cast<const bs_student_t_d_handle*>(handle);
    if (!h) return std::numeric_limits<double>::quiet_NaN();
    double p = bs_wrap<double>([&]{ return cdf(complement(h->dist, x)); });
    return -log(p);
}

long double bs_student_t_chf_l(const void* handle, long double x)  {
    auto h = static_cast<const bs_student_t_l_handle*>(handle);
    if (!h) return std::numeric_limits<long double>::quiet_NaN();
    long double p = bs_wrap<long double>([&]{ return cdf(complement(h->dist, x)); });
    return -log(p);
}

float bs_student_t_kurtosis_f(const void* handle) {
    auto h = static_cast<const bs_student_t_f_handle*>(handle);
    if (!h) return std::numeric_limits<float>::quiet_NaN();
    return bs_wrap<float>([&]{ return kurtosis(h->dist); });
}
double bs_student_t_kurtosis(const void* handle) {
    auto h = static_cast<const bs_student_t_d_handle*>(handle);
    if (!h) return std::numeric_limits<double>::quiet_NaN();
    return bs_wrap<float>([&]{ return kurtosis(h->dist); });
}
long double bs_student_t_kurtosis_l(const void* handle) {
    auto h = static_cast<const bs_student_t_l_handle*>(handle);
    if (!h) return std::numeric_limits<long double>::quiet_NaN();
    return bs_wrap<float>([&]{ return kurtosis(h->dist); });
}

float bs_student_t_kurtosis_excess_f(const void* handle) {
    auto h = static_cast<const bs_student_t_f_handle*>(handle);
    if (!h) return std::numeric_limits<float>::quiet_NaN();
    return bs_wrap<float>([&]{ return kurtosis_excess(h->dist); });
}
double bs_student_t_kurtosis_excess(const void* handle) {
    auto h = static_cast<const bs_student_t_d_handle*>(handle);
    if (!h) return std::numeric_limits<double>::quiet_NaN();
    return bs_wrap<float>([&]{ return kurtosis_excess(h->dist); });
}
long double bs_student_t_kurtosis_excess_l(const void* handle) {
    auto h = static_cast<const bs_student_t_l_handle*>(handle);
    if (!h) return std::numeric_limits<long double>::quiet_NaN();
    return bs_wrap<float>([&]{ return kurtosis_excess(h->dist); });
}

float bs_student_t_mean_f(const void* handle) {
    auto h = static_cast<const bs_student_t_f_handle*>(handle);
    if (!h) return std::numeric_limits<float>::quiet_NaN();
    return bs_wrap<float>([&]{ return mean(h->dist); });
}
double bs_student_t_mean(const void* handle) {
    auto h = static_cast<const bs_student_t_d_handle*>(handle);
    if (!h) return std::numeric_limits<double>::quiet_NaN();
    return bs_wrap<double>([&]{ return mean(h->dist); });
}
long double bs_student_t_mean_l(const void* handle) {
    auto h = static_cast<const bs_student_t_l_handle*>(handle);
    if (!h) return std::numeric_limits<long double>::quiet_NaN();
    return bs_wrap<long double>([&]{ return mean(h->dist); });
}

float bs_student_t_median_f(const void* handle) {
    auto h = static_cast<const bs_student_t_f_handle*>(handle);
    if (!h) return std::numeric_limits<float>::quiet_NaN();
    return bs_wrap<float>([&]{ return median(h->dist); });
}
double bs_student_t_median (const void* handle) {
    auto h = static_cast<const bs_student_t_d_handle*>(handle);
    if (!h) return std::numeric_limits<double>::quiet_NaN();
    return bs_wrap<long double>([&]{ return median(h->dist); });
}
long double bs_student_t_median_l(const void* handle) {
    auto h = static_cast<const bs_student_t_l_handle*>(handle);
    if (!h) return std::numeric_limits<long double>::quiet_NaN();
    return bs_wrap<long double>([&]{ return median(h->dist); });
}

float bs_student_t_mode_f(const void* handle) {
    auto h = static_cast<const bs_student_t_f_handle*>(handle);
    if (!h) return std::numeric_limits<float>::quiet_NaN();
    return bs_wrap<float>([&]{ return mode(h->dist); });
}
double bs_student_t_mode(const void* handle) {
    auto h = static_cast<const bs_student_t_d_handle*>(handle);
    if (!h) return std::numeric_limits<double>::quiet_NaN();
    return bs_wrap<double>([&]{ return mode(h->dist); });
}
long double bs_student_t_mode_l(const void* handle) {
    auto h = static_cast<const bs_student_t_l_handle*>(handle);
    if (!h) return std::numeric_limits<long double>::quiet_NaN();
    return bs_wrap<long double>([&]{ return mode(h->dist); });
}

float bs_student_t_pdf_f(const void* handle, float x) {
    auto h = static_cast<const bs_student_t_f_handle*>(handle);
    if (!h) return std::numeric_limits<float>::quiet_NaN();
    return bs_wrap<float>([&]{ return pdf(h->dist, x); });
}
double bs_student_t_pdf(const void* handle, double x) {
    auto h = static_cast<const bs_student_t_d_handle*>(handle);
    if (!h) return std::numeric_limits<double>::quiet_NaN();
    return bs_wrap<double>([&]{ return pdf(h->dist, x); });
}
long double bs_student_t_pdf_l(const void* handle, long double x) {
    auto h = static_cast<const bs_student_t_l_handle*>(handle);
    if (!h) return std::numeric_limits<long double>::quiet_NaN();
    return bs_wrap<long double>([&]{ return pdf(h->dist, x); });
}

bs_range_f bs_student_t_range_f(const void* handle) {
    auto h = static_cast<const bs_student_t_f_handle*>(handle);
    if (!h) {
        float nan = std::numeric_limits<float>::quiet_NaN();
        return bs_range_f{ nan, nan };
    }
    auto pr = range(h->dist);
    return bs_range_f{ static_cast<float>(pr.first), static_cast<float>(pr.second) };
}

bs_range_d bs_student_t_range(const void* handle) {
    auto h = static_cast<const bs_student_t_d_handle*>(handle);
    if (!h) {
        double nan = std::numeric_limits<double>::quiet_NaN();
        return bs_range_d{ nan, nan };
    }
    auto pr = range(h->dist);
    return bs_range_d{ static_cast<double>(pr.first), static_cast<double>(pr.second) };
}

bs_range_l bs_student_t_range_l(const void* handle) {
    auto h = static_cast<const bs_student_t_l_handle*>(handle);
    if (!h) {
        long double nan = std::numeric_limits<long double>::quiet_NaN();
        return bs_range_l{ nan, nan };
    }
    auto pr = range(h->dist);
    return bs_range_l{ static_cast<long double>(pr.first), static_cast<long double>(pr.second) };
}

bs_range_f bs_student_t_support_f(const void* handle) {
    auto h = static_cast<const bs_student_t_f_handle*>(handle);
    if (!h) {
        float nan = std::numeric_limits<float>::quiet_NaN();
        return bs_range_f{ nan, nan };
    }
    auto pr = support(h->dist);
    return bs_range_f{ static_cast<float>(pr.first), static_cast<float>(pr.second) };
}

bs_range_d bs_student_t_support(const void* handle) {
    auto h = static_cast<const bs_student_t_d_handle*>(handle);
    if (!h) {
        double nan = std::numeric_limits<double>::quiet_NaN();
        return bs_range_d{ nan, nan };
    }
    auto pr = support(h->dist);
    return bs_range_d{ static_cast<double>(pr.first), static_cast<double>(pr.second) };
}

bs_range_l bs_student_t_support_l(const void* handle) {
    auto h = static_cast<const bs_student_t_l_handle*>(handle);
    if (!h) {
        long double nan = std::numeric_limits<long double>::quiet_NaN();
        return bs_range_l{ nan, nan };
    }
    auto pr = support(h->dist);
    return bs_range_l{ static_cast<long double>(pr.first), static_cast<long double>(pr.second) };
}

float bs_student_t_quantile_f(const void* handle, float p) {
    auto h = static_cast<const bs_student_t_f_handle*>(handle);
    if (!h) return std::numeric_limits<float>::quiet_NaN();
    return bs_wrap<float>([&]{ return quantile(h->dist, p); });
}
double bs_student_t_quantile(const void* handle, double p) {
    auto h = static_cast<const bs_student_t_d_handle*>(handle);
    if (!h) return std::numeric_limits<double>::quiet_NaN();
    return bs_wrap<double>([&]{ return quantile(h->dist, p); });
}
long double bs_student_t_quantile_l(const void* handle, long double p) {
    auto h = static_cast<const bs_student_t_l_handle*>(handle);
    if (!h) return std::numeric_limits<long double>::quiet_NaN();
    return bs_wrap<long double>([&]{ return quantile(h->dist, p); });
}

float bs_student_t_quantile_complement_f(const void* handle, float q) {
    auto h = static_cast<const bs_student_t_f_handle*>(handle);
    if (!h) return std::numeric_limits<float>::quiet_NaN();
    return bs_wrap<float>([&]{ return quantile(complement(h->dist, q)); });
}
double bs_student_t_quantile_complement(const void* handle, double q) {
    auto h = static_cast<const bs_student_t_d_handle*>(handle);
    if (!h) return std::numeric_limits<double>::quiet_NaN();
    return bs_wrap<double>([&]{ return quantile(complement(h->dist, q)); });
}
long double bs_student_t_quantile_complement_l(const void* handle, long double q) {
    auto h = static_cast<const bs_student_t_l_handle*>(handle);
    if (!h) return std::numeric_limits<long double>::quiet_NaN();
    return bs_wrap<long double>([&]{ return quantile(complement(h->dist, q)); });
}

float bs_student_t_skewness_f(const void* handle) {
    auto h = static_cast<const bs_student_t_f_handle*>(handle);
    if (!h) return std::numeric_limits<float>::quiet_NaN();
    return bs_wrap<float>([&]{ return skewness(h->dist); });
}
double      bs_student_t_skewness (const void* handle) {
    auto h = static_cast<const bs_student_t_d_handle*>(handle);
    if (!h) return std::numeric_limits<double>::quiet_NaN();
    return bs_wrap<double>([&]{ return skewness(h->dist); });
}
long double bs_student_t_skewness_l(const void* handle) {
    auto h = static_cast<const bs_student_t_l_handle*>(handle);
    if (!h) return std::numeric_limits<long double>::quiet_NaN();
    return bs_wrap<long double>([&]{ return skewness(h->dist); });
}

float       bs_student_t_standard_deviation_f(const void* handle) {
    auto h = static_cast<const bs_student_t_f_handle*>(handle);
    if (!h) return std::numeric_limits<float>::quiet_NaN();
    return bs_wrap<float>([&]{ return standard_deviation(h->dist); });
}
double      bs_student_t_standard_deviation (const void* handle) {
    auto h = static_cast<const bs_student_t_l_handle*>(handle);
    if (!h) return std::numeric_limits<double>::quiet_NaN();
    return bs_wrap<double>([&]{ return standard_deviation(h->dist); });
}
long double bs_student_t_standard_deviation_l(const void* handle) {
    auto h = static_cast<const bs_student_t_l_handle*>(handle);
    if (!h) return std::numeric_limits<long double>::quiet_NaN();
    return bs_wrap<long double>([&]{ return standard_deviation(h->dist); });
}

float bs_student_t_variance_f(const void* handle) {
    auto h = static_cast<const bs_student_t_f_handle*>(handle);
    if (!h) return std::numeric_limits<float>::quiet_NaN();
    return bs_wrap<float>([&]{ return variance(h->dist); });
}
double bs_student_t_variance(const void* handle) {
    auto h = static_cast<const bs_student_t_d_handle*>(handle);
    if (!h) return std::numeric_limits<double>::quiet_NaN();
    return bs_wrap<double>([&]{ return variance(h->dist); });
}
long double bs_student_t_variance_l(const void* handle) {
    auto h = static_cast<const bs_student_t_l_handle*>(handle);
    if (!h) return std::numeric_limits<long double>::quiet_NaN();
    return bs_wrap<long double>([&]{ return variance(h->dist); });
}

float bs_student_t_entropy_f(const void* handle) {
    auto h = static_cast<const bs_student_t_f_handle*>(handle);
    if (!h) return std::numeric_limits<float>::quiet_NaN();
    return bs_wrap<float>([&]{ return entropy(h->dist); });
}
double bs_student_t_entropy(const void* handle) {
    auto h = static_cast<const bs_student_t_d_handle*>(handle);
    if (!h) return std::numeric_limits<double>::quiet_NaN();
    return bs_wrap<double>([&]{ return entropy(h->dist); });
}
long double bs_student_t_entropy_l(const void* handle) {
    auto h = static_cast<const bs_student_t_l_handle*>(handle);
    if (!h) return std::numeric_limits<long double>::quiet_NaN();
    return bs_wrap<long double>([&]{ return entropy(h->dist); });
}

// Static utility: required degrees of freedom
float bs_student_t_find_degrees_of_freedom_f(float difference_from_mean, float alpha, float beta, float sd, float hint) {
    return bs_wrap<float>([&]{ return students_t_distribution<float>::find_degrees_of_freedom(difference_from_mean, alpha, beta, sd, hint); });
}
double bs_student_t_find_degrees_of_freedom(double difference_from_mean, double alpha, double beta, double sd, double hint) {
    return bs_wrap<double>([&]{ return students_t_distribution<double>::find_degrees_of_freedom(difference_from_mean, alpha, beta, sd, hint); });
}
long double bs_student_t_find_degrees_of_freedom_l(long double difference_from_mean, long double alpha, long double beta, long double sd, long double hint) {
    return bs_wrap<long double>([&]{ return students_t_distribution<long double>::find_degrees_of_freedom(difference_from_mean, alpha, beta, sd, hint); });
}


} // extern "C"
