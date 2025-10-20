//
//  Created by Volker Thieme 2025.
//  Students t distribution bridge (Boost.Math) — implementation TU include
//
#include "../internal/bs_internal.hpp"
#include "../include/distributions/bs_fisher_f_distribution.h"

#include <boost/math/distributions/fisher_f.hpp>
#include <boost/math/distributions/complement.hpp>
#include <boost/math/distributions/detail/derived_accessors.hpp>

using boost::math::fisher_f_distribution;
using boost::math::complement;

// -----------------------
// Handle-based API
// -----------------------

struct bs_fisher_f_f_handle { fisher_f_distribution<float> dist; explicit bs_fisher_f_f_handle(float df1, float df2) : dist(df1, df2) {} };
struct bs_fisher_f_d_handle { fisher_f_distribution<double> dist; explicit bs_fisher_f_d_handle(double df1, double df2) : dist(df1, df2) {} };
struct bs_fisher_f_l_handle { fisher_f_distribution<long double> dist; explicit bs_fisher_f_l_handle(long double df1, long double df2) : dist(df1, df2) {} };

extern "C" {
void* bs_fisher_f_make_f(float df1, float df2) {
    try { return new bs_fisher_f_f_handle(df1, df2); } catch (...) { return nullptr; }
}
void* bs_fisher_f_make_d(double df1, double df2) {
    try { return new bs_fisher_f_d_handle(df1, df2); } catch (...) { return nullptr; }
}
void* bs_fisher_f_make_l(long double df1, long double df2) {
    try { return new bs_fisher_f_l_handle(df1, df2); } catch (...) { return nullptr; }
}
void bs_fisher_f_free_f(void* handle) { delete static_cast<bs_fisher_f_f_handle*>(handle); }
void bs_fisher_f_free_d(void* handle) { delete static_cast<bs_fisher_f_d_handle*>(handle); }
void bs_fisher_f_free_l(void* handle) { delete static_cast<bs_fisher_f_l_handle*>(handle); }

float bs_fisher_f_cdf_f(const void* handle, float x) {
    auto h = static_cast<const bs_fisher_f_f_handle*>(handle);
    if (!h) return std::numeric_limits<float>::quiet_NaN();
    return bs_wrap<float>([&]{ return cdf(h->dist, x); });
}
double bs_fisher_f_cdf_d(const void* handle, double x) {
    auto h = static_cast<const bs_fisher_f_d_handle*>(handle);
    if (!h) return std::numeric_limits<double>::quiet_NaN();
    return bs_wrap<double>([&]{ return cdf(h->dist, x); });
}
long double bs_fisher_f_cdf_l(const void* handle, long double x) {
    auto h = static_cast<const bs_fisher_f_l_handle*>(handle);
    if (!h) return std::numeric_limits<long double>::quiet_NaN();
    return bs_wrap<long double>([&]{ return cdf(h->dist, x); });
}

float bs_fisher_f_ccdf_f(const void* handle, float x) {
    auto h = static_cast<const bs_fisher_f_f_handle*>(handle);
    if (!h) return std::numeric_limits<float>::quiet_NaN();
    return bs_wrap<float>([&]{ return cdf(complement(h->dist, x)); });
}
double bs_fisher_f_ccdf_d(const void* handle, double x) {
    auto h = static_cast<const bs_fisher_f_d_handle*>(handle);
    if (!h) return std::numeric_limits<double>::quiet_NaN();
    return bs_wrap<double>([&]{ return cdf(complement(h->dist, x)); });
}
long double bs_fisher_f_ccdf_l(const void* handle, long double x) {
    auto h = static_cast<const bs_fisher_f_l_handle*>(handle);
    if (!h) return std::numeric_limits<long double>::quiet_NaN();
    return bs_wrap<long double>([&]{ return cdf(complement(h->dist, x)); });
}
/// MARK: hazard function

float bs_fisher_f_hazard_f(const void* handle, float x) {
    auto h = static_cast<const bs_fisher_f_f_handle*>(handle);
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
double bs_fisher_f_hazard_d(const void* handle, double x) {
    auto h = static_cast<const bs_fisher_f_d_handle*>(handle);
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
long double bs_fisher_f_hazard_l(const void* handle, long double x)  {
    auto h = static_cast<const bs_fisher_f_l_handle*>(handle);
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
float bs_fisher_f_chf_f(const void* handle, float x) {
    auto h = static_cast<const bs_fisher_f_f_handle*>(handle);
    if (!h) return std::numeric_limits<float>::quiet_NaN();
    float p = bs_wrap<float>([&]{ return cdf(complement(h->dist, x)); });
    return -log(p);
}
double     bs_fisher_f_chf_d(const void* handle, double x) {
    auto h = static_cast<const bs_fisher_f_d_handle*>(handle);
    if (!h) return std::numeric_limits<double>::quiet_NaN();
    double p = bs_wrap<double>([&]{ return cdf(complement(h->dist, x)); });
    return -log(p);
}

long double bs_fisher_f_chf_l(const void* handle, long double x)  {
    auto h = static_cast<const bs_fisher_f_d_handle*>(handle);
    if (!h) return std::numeric_limits<long double>::quiet_NaN();
    long double p = bs_wrap<long double>([&]{ return cdf(complement(h->dist, x)); });
    return -log(p);
}

float bs_fisher_f_kurtosis_f(const void* handle) {
    auto h = static_cast<const bs_fisher_f_f_handle*>(handle);
    if (!h) return std::numeric_limits<float>::quiet_NaN();
    return bs_wrap<float>([&]{ return kurtosis(h->dist); });
}
double bs_fisher_f_kurtosis_d(const void* handle) {
    auto h = static_cast<const bs_fisher_f_d_handle*>(handle);
    if (!h) return std::numeric_limits<double>::quiet_NaN();
    return bs_wrap<float>([&]{ return kurtosis(h->dist); });
}
long double bs_fisher_f_kurtosis_l(const void* handle) {
    auto h = static_cast<const bs_fisher_f_l_handle*>(handle);
    if (!h) return std::numeric_limits<long double>::quiet_NaN();
    return bs_wrap<float>([&]{ return kurtosis(h->dist); });
}

float bs_fisher_f_kurtosis_excess_f(const void* handle) {
    auto h = static_cast<const bs_fisher_f_f_handle*>(handle);
    if (!h) return std::numeric_limits<float>::quiet_NaN();
    return bs_wrap<float>([&]{ return kurtosis_excess(h->dist); });
}
double bs_fisher_f_kurtosis_excess_d(const void* handle) {
    auto h = static_cast<const bs_fisher_f_d_handle*>(handle);
    if (!h) return std::numeric_limits<double>::quiet_NaN();
    return bs_wrap<float>([&]{ return kurtosis_excess(h->dist); });
}
long double bs_fisher_f_kurtosis_excess_l(const void* handle) {
    auto h = static_cast<const bs_fisher_f_l_handle*>(handle);
    if (!h) return std::numeric_limits<long double>::quiet_NaN();
    return bs_wrap<float>([&]{ return kurtosis_excess(h->dist); });
}

float bs_fisher_f_mean_f(const void* handle) {
    auto h = static_cast<const bs_fisher_f_f_handle*>(handle);
    if (!h) return std::numeric_limits<float>::quiet_NaN();
    return bs_wrap<float>([&]{ return mean(h->dist); });
}
double bs_fisher_f_mean_d(const void* handle) {
    auto h = static_cast<const bs_fisher_f_d_handle*>(handle);
    if (!h) return std::numeric_limits<double>::quiet_NaN();
    return bs_wrap<double>([&]{ return mean(h->dist); });
}
long double bs_fisher_f_mean_l(const void* handle) {
    auto h = static_cast<const bs_fisher_f_l_handle*>(handle);
    if (!h) return std::numeric_limits<long double>::quiet_NaN();
    return bs_wrap<long double>([&]{ return mean(h->dist); });
}

float bs_fisher_f_median_f(const void* handle) {
    auto h = static_cast<const bs_fisher_f_f_handle*>(handle);
    if (!h) return std::numeric_limits<float>::quiet_NaN();
    return bs_wrap<float>([&]{ return median(h->dist); });
}
double bs_fisher_f_median_d(const void* handle) {
    auto h = static_cast<const bs_fisher_f_d_handle*>(handle);
    if (!h) return std::numeric_limits<double>::quiet_NaN();
    return bs_wrap<long double>([&]{ return median(h->dist); });
}
long double bs_fisher_f_median_l(const void* handle) {
    auto h = static_cast<const bs_fisher_f_l_handle*>(handle);
    if (!h) return std::numeric_limits<long double>::quiet_NaN();
    return bs_wrap<long double>([&]{ return median(h->dist); });
}

float bs_fisher_f_mode_f(const void* handle) {
    auto h = static_cast<const bs_fisher_f_f_handle*>(handle);
    if (!h) return std::numeric_limits<float>::quiet_NaN();
    return bs_wrap<float>([&]{ return mode(h->dist); });
}
double bs_fisher_f_mode_d(const void* handle) {
    auto h = static_cast<const bs_fisher_f_d_handle*>(handle);
    if (!h) return std::numeric_limits<double>::quiet_NaN();
    return bs_wrap<double>([&]{ return mode(h->dist); });
}
long double bs_fisher_f_mode_l(const void* handle) {
    auto h = static_cast<const bs_fisher_f_l_handle*>(handle);
    if (!h) return std::numeric_limits<long double>::quiet_NaN();
    return bs_wrap<long double>([&]{ return mode(h->dist); });
}

float bs_fisher_f_pdf_f(const void* handle, float x) {
    auto h = static_cast<const bs_fisher_f_f_handle*>(handle);
    if (!h) return std::numeric_limits<float>::quiet_NaN();
    return bs_wrap<float>([&]{ return pdf(h->dist, x); });
}
double bs_fisher_f_pdf_d(const void* handle, double x) {
    auto h = static_cast<const bs_fisher_f_d_handle*>(handle);
    if (!h) return std::numeric_limits<double>::quiet_NaN();
    return bs_wrap<double>([&]{ return pdf(h->dist, x); });
}
long double bs_fisher_f_pdf_l(const void* handle, long double x) {
    auto h = static_cast<const bs_fisher_f_l_handle*>(handle);
    if (!h) return std::numeric_limits<long double>::quiet_NaN();
    return bs_wrap<long double>([&]{ return pdf(h->dist, x); });
}

bs_range_f bs_fisher_f_range_f(const void* handle) {
    auto h = static_cast<const bs_fisher_f_f_handle*>(handle);
    if (!h) {
        float nan = std::numeric_limits<float>::quiet_NaN();
        return bs_range_f{ nan, nan };
    }
    auto pr = range(h->dist);
    return bs_range_f{ static_cast<float>(pr.first), static_cast<float>(pr.second) };
}

bs_range_d bs_fisher_f_range_d(const void* handle) {
    auto h = static_cast<const bs_fisher_f_d_handle*>(handle);
    if (!h) {
        double nan = std::numeric_limits<double>::quiet_NaN();
        return bs_range_d{ nan, nan };
    }
    auto pr = range(h->dist);
    return bs_range_d{ static_cast<double>(pr.first), static_cast<double>(pr.second) };
}

bs_range_l bs_fisher_f_range_l(const void* handle) {
    auto h = static_cast<const bs_fisher_f_l_handle*>(handle);
    if (!h) {
        long double nan = std::numeric_limits<long double>::quiet_NaN();
        return bs_range_l{ nan, nan };
    }
    auto pr = range(h->dist);
    return bs_range_l{ static_cast<long double>(pr.first), static_cast<long double>(pr.second) };
}

bs_range_f bs_fisher_f_support_f(const void* handle) {
    auto h = static_cast<const bs_fisher_f_f_handle*>(handle);
    if (!h) {
        float nan = std::numeric_limits<float>::quiet_NaN();
        return bs_range_f{ nan, nan };
    }
    auto pr = support(h->dist);
    return bs_range_f{ static_cast<float>(pr.first), static_cast<float>(pr.second) };
}

bs_range_d bs_fisher_f_support_d(const void* handle) {
    auto h = static_cast<const bs_fisher_f_d_handle*>(handle);
    if (!h) {
        double nan = std::numeric_limits<double>::quiet_NaN();
        return bs_range_d{ nan, nan };
    }
    auto pr = support(h->dist);
    return bs_range_d{ static_cast<double>(pr.first), static_cast<double>(pr.second) };
}

bs_range_l bs_fisher_f_support_l(const void* handle) {
    auto h = static_cast<const bs_fisher_f_l_handle*>(handle);
    if (!h) {
        long double nan = std::numeric_limits<long double>::quiet_NaN();
        return bs_range_l{ nan, nan };
    }
    auto pr = support(h->dist);
    return bs_range_l{ static_cast<long double>(pr.first), static_cast<long double>(pr.second) };
}

float bs_fisher_f_quantile_f(const void* handle, float p) {
    auto h = static_cast<const bs_fisher_f_f_handle*>(handle);
    if (!h) return std::numeric_limits<float>::quiet_NaN();
    return bs_wrap<float>([&]{ return quantile(h->dist, p); });
}
double bs_fisher_f_quantile_d(const void* handle, double p) {
    auto h = static_cast<const bs_fisher_f_d_handle*>(handle);
    if (!h) return std::numeric_limits<double>::quiet_NaN();
    return bs_wrap<double>([&]{ return quantile(h->dist, p); });
}
long double bs_fisher_f_quantile_l(const void* handle, long double p) {
    auto h = static_cast<const bs_fisher_f_l_handle*>(handle);
    if (!h) return std::numeric_limits<long double>::quiet_NaN();
    return bs_wrap<long double>([&]{ return quantile(h->dist, p); });
}

float bs_fisher_f_quantile_complement_f(const void* handle, float q) {
    auto h = static_cast<const bs_fisher_f_f_handle*>(handle);
    if (!h) return std::numeric_limits<float>::quiet_NaN();
    return bs_wrap<float>([&]{ return quantile(complement(h->dist, q)); });
}
double bs_fisher_f_quantile_complement_d(const void* handle, double q) {
    auto h = static_cast<const bs_fisher_f_d_handle*>(handle);
    if (!h) return std::numeric_limits<double>::quiet_NaN();
    return bs_wrap<double>([&]{ return quantile(complement(h->dist, q)); });
}
long double bs_fisher_f_quantile_complement_l(const void* handle, long double q) {
    auto h = static_cast<const bs_fisher_f_l_handle*>(handle);
    if (!h) return std::numeric_limits<long double>::quiet_NaN();
    return bs_wrap<long double>([&]{ return quantile(complement(h->dist, q)); });
}

float bs_fisher_f_skewness_f(const void* handle) {
    auto h = static_cast<const bs_fisher_f_f_handle*>(handle);
    if (!h) return std::numeric_limits<float>::quiet_NaN();
    return bs_wrap<float>([&]{ return skewness(h->dist); });
}
double      bs_fisher_f_skewness_d(const void* handle) {
    auto h = static_cast<const bs_fisher_f_d_handle*>(handle);
    if (!h) return std::numeric_limits<double>::quiet_NaN();
    return bs_wrap<double>([&]{ return skewness(h->dist); });
}
long double bs_fisher_f_skewness_l(const void* handle) {
    auto h = static_cast<const bs_fisher_f_l_handle*>(handle);
    if (!h) return std::numeric_limits<long double>::quiet_NaN();
    return bs_wrap<long double>([&]{ return skewness(h->dist); });
}

float       bs_fisher_f_standard_deviation_f(const void* handle) {
    auto h = static_cast<const bs_fisher_f_f_handle*>(handle);
    if (!h) return std::numeric_limits<float>::quiet_NaN();
    return bs_wrap<float>([&]{ return standard_deviation(h->dist); });
}
double      bs_fisher_f_standard_deviation_d(const void* handle) {
    auto h = static_cast<const bs_fisher_f_l_handle*>(handle);
    if (!h) return std::numeric_limits<double>::quiet_NaN();
    return bs_wrap<double>([&]{ return standard_deviation(h->dist); });
}
long double bs_fisher_f_standard_deviation_l(const void* handle) {
    auto h = static_cast<const bs_fisher_f_l_handle*>(handle);
    if (!h) return std::numeric_limits<long double>::quiet_NaN();
    return bs_wrap<long double>([&]{ return standard_deviation(h->dist); });
}

float bs_fisher_f_variance_f(const void* handle) {
    auto h = static_cast<const bs_fisher_f_f_handle*>(handle);
    if (!h) return std::numeric_limits<float>::quiet_NaN();
    return bs_wrap<float>([&]{ return variance(h->dist); });
}
double bs_fisher_f_variance_d(const void* handle) {
    auto h = static_cast<const bs_fisher_f_d_handle*>(handle);
    if (!h) return std::numeric_limits<double>::quiet_NaN();
    return bs_wrap<double>([&]{ return variance(h->dist); });
}
long double bs_fisher_f_variance_l(const void* handle) {
    auto h = static_cast<const bs_fisher_f_l_handle*>(handle);
    if (!h) return std::numeric_limits<long double>::quiet_NaN();
    return bs_wrap<long double>([&]{ return variance(h->dist); });
}

} // extern "C"
