//
//  Generic distribution vtable (Boost.Math) — implementation TU include
//
#include "../internal/bs_internal.hpp"
#include "../include/distributions/bs_generic_distribution.h"

#include <string>
#include <algorithm>
#include <cctype>

#include <boost/math/distributions/gamma.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/distributions/fisher_f.hpp>
#include <boost/math/distributions/arcsine.hpp>
#include <boost/math/distributions/beta.hpp>
#include <boost/math/distributions/complement.hpp>
#include <boost/math/tools/numeric_limits.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/bernoulli.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions/negative_binomial.hpp>
#include <boost/math/distributions/cauchy.hpp>
#include <boost/math/distributions/exponential.hpp>
#include <boost/math/distributions/extreme_value.hpp>
#include <boost/math/distributions/geometric.hpp>
#include <boost/math/distributions/holtsmark.hpp>
#include <boost/math/distributions/hyperexponential.hpp>
#include <boost/math/distributions/hypergeometric.hpp>
#include <boost/math/distributions/inverse_chi_squared.hpp>

using boost::math::gamma_distribution;
using boost::math::students_t_distribution;
using boost::math::fisher_f_distribution;
using boost::math::arcsine_distribution;
using boost::math::beta_distribution;
using boost::math::chi_squared_distribution;
using boost::math::complement;
using boost::math::bernoulli_distribution;
using boost::math::binomial_distribution;
using boost::math::negative_binomial_distribution;
using boost::math::cauchy_distribution;
using boost::math::exponential_distribution;
using boost::math::extreme_value_distribution;
using boost::math::geometric_distribution;
using boost::math::holtsmark_distribution;
using boost::math::hyperexponential_distribution;
using boost::math::hypergeometric_distribution;
using boost::math::inverse_chi_squared_distribution;

namespace {
// Lowercase ASCII helper
inline std::string lower_ascii(const char* s) {
    std::string r = s ? std::string(s) : std::string();
    std::transform(r.begin(), r.end(), r.begin(), [](unsigned char c){ return static_cast<char>(std::tolower(c)); });
    return r;
}

template <class Param, class ValueT>
inline bool find_param(const Param* params, size_t n, const char* const* keys, size_t keyCount, ValueT* out) {
    for (size_t i = 0; i < n; ++i) {
        std::string k = lower_ascii(params[i].key);
        for (size_t j = 0; j < keyCount; ++j) {
            if (k == lower_ascii(keys[j])) { *out = params[i].value; return true; }
        }
    }
    return false;
}

// Templated delete function
template <class Handle>
static void free_fn(void* p) { delete static_cast<Handle*>(p); }

// Thunks for real-valued functions
template <class Handle, class T>
static T pdf_fn(const void* p, T x) {
    auto h = static_cast<const Handle*>(p);
    if (!h) return std::numeric_limits<T>::quiet_NaN();
    return bs_wrap<T>([&]{ return pdf(h->dist, x); });
}
template <class Handle, class T>
static T logpdf_fn(const void* p, T x) {
    auto h = static_cast<const Handle*>(p);
    if (!h) return std::numeric_limits<T>::quiet_NaN();
    return bs_wrap<T>([&]{ return std::log(pdf(h->dist, x)); });
}
template <class Handle, class T>
static T cdf_fn(const void* p, T x) {
    auto h = static_cast<const Handle*>(p);
    if (!h) return std::numeric_limits<T>::quiet_NaN();
    return bs_wrap<T>([&]{ return cdf(h->dist, x); });
}
template <class Handle, class T>
static T sf_fn(const void* p, T x) {
    auto h = static_cast<const Handle*>(p);
    if (!h) return std::numeric_limits<T>::quiet_NaN();
    return bs_wrap<T>([&]{ return cdf(complement(h->dist, x)); });
}
template <class Handle, class T>
static T hazard_fn(const void* p, T x) {
    auto h = static_cast<const Handle*>(p);
    if (!h) return std::numeric_limits<T>::quiet_NaN();
    T s = bs_wrap<T>([&]{ return cdf(complement(h->dist, x)); });
    T d = bs_wrap<T>([&]{ return pdf(h->dist, x); });
    if (d == 0) return static_cast<T>(0);
    if (d > s * boost::math::tools::max_value<T>()) return std::numeric_limits<T>::quiet_NaN();
    return d / s;
}
template <class Handle, class T>
static T chf_fn(const void* p, T x) {
    auto h = static_cast<const Handle*>(p);
    if (!h) return std::numeric_limits<T>::quiet_NaN();
    T s = bs_wrap<T>([&]{ return cdf(complement(h->dist, x)); });
    return static_cast<T>(-std::log(s));
}
template <class Handle, class T>
static T quantile_fn(const void* p, T q) {
    auto h = static_cast<const Handle*>(p);
    if (!h) return std::numeric_limits<T>::quiet_NaN();
    return bs_wrap<T>([&]{ return quantile(h->dist, q); });
}
template <class Handle, class T>
static T quantile_complement_fn(const void* p, T q) {
    auto h = static_cast<const Handle*>(p);
    if (!h) return std::numeric_limits<T>::quiet_NaN();
    return bs_wrap<T>([&]{ return quantile(complement(h->dist, q)); });
}

// Range and stats
template <class Handle, class Range, class T>
static Range range_fn(const void* p) {
    auto h = static_cast<const Handle*>(p);
    if (!h) {
        T nan = std::numeric_limits<T>::quiet_NaN();
        return Range{ nan, nan };
    }
    auto pr = range(h->dist);
    return Range{ static_cast<T>(pr.first), static_cast<T>(pr.second) };
}
template <class Handle, class T>
static T mean_fn(const void* p) {
    auto h = static_cast<const Handle*>(p);
    if (!h) return std::numeric_limits<T>::quiet_NaN();
    return bs_wrap<T>([&]{ return mean(h->dist); });
}
template <class Handle, class T>
static T variance_fn(const void* p) {
    auto h = static_cast<const Handle*>(p);
    if (!h) return std::numeric_limits<T>::quiet_NaN();
    return bs_wrap<T>([&]{ return variance(h->dist); });
}
template <class Handle, class T>
static T skewness_fn(const void* p) {
    auto h = static_cast<const Handle*>(p);
    if (!h) return std::numeric_limits<T>::quiet_NaN();
    return bs_wrap<T>([&]{ return skewness(h->dist); });
}
template <class Handle, class T>
static T kurtosis_fn(const void* p) {
    auto h = static_cast<const Handle*>(p);
    if (!h) return std::numeric_limits<T>::quiet_NaN();
    return bs_wrap<T>([&]{ return kurtosis(h->dist); });
}
template <class Handle, class T>
static T kurtosis_excess_fn(const void* p) {
    auto h = static_cast<const Handle*>(p);
    if (!h) return std::numeric_limits<T>::quiet_NaN();
    return bs_wrap<T>([&]{ return kurtosis_excess(h->dist); });
}
template <class Handle, class T>
static T mode_fn(const void* p) {
    auto h = static_cast<const Handle*>(p);
    if (!h) return std::numeric_limits<T>::quiet_NaN();
    return bs_wrap<T>([&]{ return mode(h->dist); });
}
template <class Handle, class T>
static T median_fn(const void* p) {
    auto h = static_cast<const Handle*>(p);
    if (!h) return std::numeric_limits<T>::quiet_NaN();
    return bs_wrap<T>([&]{ return median(h->dist); });
}
template <class Handle, class T>
static T entropy_fn(const void* p) {
    auto h = static_cast<const Handle*>(p);
    if (!h) return std::numeric_limits<T>::quiet_NaN();
    return bs_wrap<T>([&]{ return entropy(h->dist); });
}

// Handle types for current distributions
struct gamma_f_handle { gamma_distribution<float> dist; gamma_f_handle(float k, float th): dist(k, th) {} };
struct gamma_d_handle { gamma_distribution<double> dist; gamma_d_handle(double k, double th): dist(k, th) {} };
struct gamma_l_handle { gamma_distribution<long double> dist; gamma_l_handle(long double k, long double th): dist(k, th) {} };

struct student_t_f_handle { students_t_distribution<float> dist; explicit student_t_f_handle(float v): dist(v) {} };
struct student_t_d_handle { students_t_distribution<double> dist; explicit student_t_d_handle(double v): dist(v) {} };
struct student_t_l_handle { students_t_distribution<long double> dist; explicit student_t_l_handle(long double v): dist(v) {} };

// Fisher–Snedecor F
struct fisher_f_f_handle { fisher_f_distribution<float> dist; fisher_f_f_handle(float df1, float df2): dist(df1, df2) {} };
struct fisher_f_d_handle { fisher_f_distribution<double> dist; fisher_f_d_handle(double df1, double df2): dist(df1, df2) {} };
struct fisher_f_l_handle { fisher_f_distribution<long double> dist; fisher_f_l_handle(long double df1, long double df2): dist(df1, df2) {} };

// Arcsine
struct arcsine_f_handle { arcsine_distribution<float> dist; arcsine_f_handle(float a, float b): dist(a, b) {} };
struct arcsine_d_handle { arcsine_distribution<double> dist; arcsine_d_handle(double a, double b): dist(a, b) {} };
struct arcsine_l_handle { arcsine_distribution<long double> dist; arcsine_l_handle(long double a, long double b): dist(a, b) {} };

// Beta
struct beta_f_handle { beta_distribution<float> dist; beta_f_handle(float a, float b): dist(a, b) {} };
struct beta_d_handle { beta_distribution<double> dist; beta_d_handle(double a, double b): dist(a, b) {} };
struct beta_l_handle { beta_distribution<long double> dist; beta_l_handle(long double a, long double b): dist(a, b) {} };

// Chi-Squared
struct chi_squared_f_handle { chi_squared_distribution<float> dist; explicit chi_squared_f_handle(float v): dist(v) {} };
struct chi_squared_d_handle { chi_squared_distribution<double> dist; explicit chi_squared_d_handle(double v): dist(v) {} };
struct chi_squared_l_handle { chi_squared_distribution<long double> dist; explicit chi_squared_l_handle(long double v): dist(v) {} };

// Bernoulli
struct bernoulli_f_handle { bernoulli_distribution<float> dist; explicit bernoulli_f_handle(float p): dist(p) {} };
struct bernoulli_d_handle { bernoulli_distribution<double> dist; explicit bernoulli_d_handle(double p): dist(p) {} };
struct bernoulli_l_handle { bernoulli_distribution<long double> dist; explicit bernoulli_l_handle(long double p): dist(p) {} };

// Binomial
struct binomial_f_handle { binomial_distribution<float> dist; binomial_f_handle(unsigned n, float p): dist(n, p) {} };
struct binomial_d_handle { binomial_distribution<double> dist; binomial_d_handle(unsigned n, double p): dist(n, p) {} };
struct binomial_l_handle { binomial_distribution<long double> dist; binomial_l_handle(unsigned n, long double p): dist(n, p) {} };

// Cauchy
struct cauchy_f_handle { cauchy_distribution<float> dist; cauchy_f_handle(float loc, float scale): dist(loc, scale) {} };
struct cauchy_d_handle { cauchy_distribution<double> dist; cauchy_d_handle(double loc, double scale): dist(loc, scale) {} };
struct cauchy_l_handle { cauchy_distribution<long double> dist; cauchy_l_handle(long double loc, long double scale): dist(loc, scale) {} };

// Exponential
struct exponential_f_handle { exponential_distribution<float> dist; explicit exponential_f_handle(float lambda): dist(lambda) {} };
struct exponential_d_handle { exponential_distribution<double> dist; explicit exponential_d_handle(double lambda): dist(lambda) {} };
struct exponential_l_handle { exponential_distribution<long double> dist; explicit exponential_l_handle(long double lambda): dist(lambda) {} };

// Extreme Value (Gumbel)
struct extreme_value_f_handle { extreme_value_distribution<float> dist; extreme_value_f_handle(float loc, float scale): dist(loc, scale) {} };
struct extreme_value_d_handle { extreme_value_distribution<double> dist; extreme_value_d_handle(double loc, double scale): dist(loc, scale) {} };
struct extreme_value_l_handle { extreme_value_distribution<long double> dist; extreme_value_l_handle(long double loc, long double scale): dist(loc, scale) {} };

// Geometric
struct geometric_f_handle { geometric_distribution<float> dist; explicit geometric_f_handle(float p): dist(p) {} };
struct geometric_d_handle { geometric_distribution<double> dist; explicit geometric_d_handle(double p): dist(p) {} };
struct geometric_l_handle { geometric_distribution<long double> dist; explicit geometric_l_handle(long double p): dist(p) {} };

// Holtsmark
struct holtsmark_f_handle { holtsmark_distribution<float> dist; holtsmark_f_handle(float loc, float scale): dist(loc, scale) {} };
struct holtsmark_d_handle { holtsmark_distribution<double> dist; holtsmark_d_handle(double loc, double scale): dist(loc, scale) {} };
struct holtsmark_l_handle { holtsmark_distribution<long double> dist; holtsmark_l_handle(long double loc, long double scale): dist(loc, scale) {} };

} // namespace

extern "C" {

bool bs_dist_make_d(const char* name, const bs_param_d* params, size_t count, bs_dist_d* out) {
    if (!name || !out) return false;
    try {
        std::string n = lower_ascii(name);
    if (n == "gamma" || n == "gamma_distribution") {
        // Gamma(k, theta)
        double k = 0, th = 1;
        const char* kKeys[] = { "shape", "k" };
        const char* thKeys[] = { "scale", "theta" };
        if (!find_param(params, count, kKeys, 2, &k)) return false;
        find_param(params, count, thKeys, 2, &th);
        auto* h = new (std::nothrow) gamma_d_handle(k, th);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<gamma_d_handle, double>;
        out->logpdf = &logpdf_fn<gamma_d_handle, double>;
        out->cdf = &cdf_fn<gamma_d_handle, double>;
        out->sf = &sf_fn<gamma_d_handle, double>;
        out->hazard = &hazard_fn<gamma_d_handle, double>;
        out->chf = &chf_fn<gamma_d_handle, double>;
        out->quantile = &quantile_fn<gamma_d_handle, double>;
        out->quantile_complement = &quantile_complement_fn<gamma_d_handle, double>;
        out->range = &range_fn<gamma_d_handle, bs_range_d, double>;
        out->mean = &mean_fn<gamma_d_handle, double>;
        out->variance = &variance_fn<gamma_d_handle, double>;
        out->skewness = &skewness_fn<gamma_d_handle, double>;
        out->kurtosis = &kurtosis_fn<gamma_d_handle, double>;
        out->kurtosis_excess = &kurtosis_excess_fn<gamma_d_handle, double>;
        out->mode = &mode_fn<gamma_d_handle, double>;
        out->median = &median_fn<gamma_d_handle, double>;
        out->entropy = &entropy_fn<gamma_d_handle, double>;
        out->free = &free_fn<gamma_d_handle>;
        return true;
    } else if (n == "studentt" || n == "student_t" || n == "students_t" || n == "t" || n == "t_distribution") {
        double v = 0;
        const char* dfKeys[] = { "df", "nu", "degreesoffreedom" };
        if (!find_param(params, count, dfKeys, 3, &v)) return false;
        auto* h = new (std::nothrow) student_t_d_handle(v);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<student_t_d_handle, double>;
        out->logpdf = &logpdf_fn<student_t_d_handle, double>;
        out->cdf = &cdf_fn<student_t_d_handle, double>;
        out->sf = &sf_fn<student_t_d_handle, double>;
        out->hazard = &hazard_fn<student_t_d_handle, double>;
        out->chf = &chf_fn<student_t_d_handle, double>;
        out->quantile = &quantile_fn<student_t_d_handle, double>;
        out->quantile_complement = &quantile_complement_fn<student_t_d_handle, double>;
        out->range = &range_fn<student_t_d_handle, bs_range_d, double>;
        out->mean = &mean_fn<student_t_d_handle, double>;
        out->variance = &variance_fn<student_t_d_handle, double>;
        out->skewness = &skewness_fn<student_t_d_handle, double>;
        out->kurtosis = &kurtosis_fn<student_t_d_handle, double>;
        out->kurtosis_excess = &kurtosis_excess_fn<student_t_d_handle, double>;
        out->mode = &mode_fn<student_t_d_handle, double>;
        out->median = &median_fn<student_t_d_handle, double>;
        out->entropy = &entropy_fn<student_t_d_handle, double>;
        out->free = &free_fn<student_t_d_handle>;
        return true;
    } else if (n == "fisherf" || n == "f" || n == "f_distribution") {
        double df1 = 0, df2 = 0;
        const char* df1Keys[] = { "df1", "d1", "m", "degreesoffreedom1" };
        const char* df2Keys[] = { "df2", "d2", "n", "degreesoffreedom2" };
        if (!find_param(params, count, df1Keys, 4, &df1)) return false;
        if (!find_param(params, count, df2Keys, 4, &df2)) return false;
        auto* h = new (std::nothrow) fisher_f_d_handle(df1, df2);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<fisher_f_d_handle, double>;
        out->logpdf = &logpdf_fn<fisher_f_d_handle, double>;
        out->cdf = &cdf_fn<fisher_f_d_handle, double>;
        out->sf = &sf_fn<fisher_f_d_handle, double>;
        out->hazard = &hazard_fn<fisher_f_d_handle, double>;
        out->chf = &chf_fn<fisher_f_d_handle, double>;
        out->quantile = &quantile_fn<fisher_f_d_handle, double>;
        out->quantile_complement = &quantile_complement_fn<fisher_f_d_handle, double>;
        out->range = &range_fn<fisher_f_d_handle, bs_range_d, double>;
        out->mean = &mean_fn<fisher_f_d_handle, double>;
        out->variance = &variance_fn<fisher_f_d_handle, double>;
        out->skewness = &skewness_fn<fisher_f_d_handle, double>;
        out->kurtosis = &kurtosis_fn<fisher_f_d_handle, double>;
        out->kurtosis_excess = &kurtosis_excess_fn<fisher_f_d_handle, double>;
        out->mode = &mode_fn<fisher_f_d_handle, double>;
        out->median = &median_fn<fisher_f_d_handle, double>;
        out->entropy = nullptr; // Boost fisher_f doesn’t expose entropy accessor; leave null
        out->free = &free_fn<fisher_f_d_handle>;
        return true;
    } else if (n == "arcsine" || n == "arcsine_distribution") {
        double a = 0, b = 0;
        const char* aKeys[] = { "minx", "min", "a", "lower" };
        const char* bKeys[] = { "maxx", "max", "b", "upper" };
        if (!find_param(params, count, aKeys, 4, &a)) return false;
        if (!find_param(params, count, bKeys, 4, &b)) return false;
        auto* h = new (std::nothrow) arcsine_d_handle(a, b);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<arcsine_d_handle, double>;
        out->logpdf = &logpdf_fn<arcsine_d_handle, double>;
        out->cdf = &cdf_fn<arcsine_d_handle, double>;
        out->sf = &sf_fn<arcsine_d_handle, double>;
        out->hazard = &hazard_fn<arcsine_d_handle, double>;
        out->chf = &chf_fn<arcsine_d_handle, double>;
        out->quantile = &quantile_fn<arcsine_d_handle, double>;
        out->quantile_complement = &quantile_complement_fn<arcsine_d_handle, double>;
        out->range = &range_fn<arcsine_d_handle, bs_range_d, double>;
        out->mean = &mean_fn<arcsine_d_handle, double>;
        out->variance = &variance_fn<arcsine_d_handle, double>;
        out->skewness = &skewness_fn<arcsine_d_handle, double>;
        out->kurtosis = &kurtosis_fn<arcsine_d_handle, double>;
        out->kurtosis_excess = &kurtosis_excess_fn<arcsine_d_handle, double>;
        out->mode = &mode_fn<arcsine_d_handle, double>;
        out->median = &median_fn<arcsine_d_handle, double>;
        out->entropy = nullptr; // not provided in Boost
        out->free = &free_fn<arcsine_d_handle>;
        return true;
    } else if (n == "beta" || n == "beta_distribution") {
        double a = 0, b = 0;
        const char* aKeys[] = { "alpha", "a", "p", "shape1" };
        const char* bKeys[] = { "beta", "b", "q", "shape2" };
        if (!find_param(params, count, aKeys, 4, &a)) return false;
        if (!find_param(params, count, bKeys, 4, &b)) return false;
        auto* h = new (std::nothrow) beta_d_handle(a, b);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<beta_d_handle, double>;
        out->logpdf = &logpdf_fn<beta_d_handle, double>;
        out->cdf = &cdf_fn<beta_d_handle, double>;
        out->sf = &sf_fn<beta_d_handle, double>;
        out->hazard = &hazard_fn<beta_d_handle, double>;
        out->chf = &chf_fn<beta_d_handle, double>;
        out->quantile = &quantile_fn<beta_d_handle, double>;
        out->quantile_complement = &quantile_complement_fn<beta_d_handle, double>;
        out->range = &range_fn<beta_d_handle, bs_range_d, double>;
        out->mean = &mean_fn<beta_d_handle, double>;
        out->variance = &variance_fn<beta_d_handle, double>;
        out->skewness = &skewness_fn<beta_d_handle, double>;
        out->kurtosis = &kurtosis_fn<beta_d_handle, double>;
        out->kurtosis_excess = &kurtosis_excess_fn<beta_d_handle, double>;
        out->mode = &mode_fn<beta_d_handle, double>;
        out->median = &median_fn<beta_d_handle, double>;
        out->entropy = nullptr; // not provided in Boost
        out->free = &free_fn<beta_d_handle>;
        return true;
    } else if (n == "chisquared" || n == "chi_squared" || n == "chi2" || n == "chi-squared" || n == "chisquare") {
        double v = 0;
        const char* dfKeys[] = { "df", "nu", "degreesoffreedom" };
        if (!find_param(params, count, dfKeys, 3, &v)) return false;
        auto* h = new (std::nothrow) chi_squared_d_handle(v);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<chi_squared_d_handle, double>;
        out->logpdf = &logpdf_fn<chi_squared_d_handle, double>;
        out->cdf = &cdf_fn<chi_squared_d_handle, double>;
        out->sf = &sf_fn<chi_squared_d_handle, double>;
        out->hazard = &hazard_fn<chi_squared_d_handle, double>;
        out->chf = &chf_fn<chi_squared_d_handle, double>;
        out->quantile = &quantile_fn<chi_squared_d_handle, double>;
        out->quantile_complement = &quantile_complement_fn<chi_squared_d_handle, double>;
        out->range = &range_fn<chi_squared_d_handle, bs_range_d, double>;
        out->mean = &mean_fn<chi_squared_d_handle, double>;
        out->variance = &variance_fn<chi_squared_d_handle, double>;
        out->skewness = &skewness_fn<chi_squared_d_handle, double>;
        out->kurtosis = &kurtosis_fn<chi_squared_d_handle, double>;
        out->kurtosis_excess = &kurtosis_excess_fn<chi_squared_d_handle, double>;
        out->mode = &mode_fn<chi_squared_d_handle, double>;
        out->median = &median_fn<chi_squared_d_handle, double>;
        out->entropy = nullptr; // not provided in Boost for chi-squared
        out->free = &free_fn<chi_squared_d_handle>;
        return true;
        
    } else if (n == "bernoulli" || n == "bernoulli_distribution") {
        double p = 0;
        const char* pKeys[] = { "p", "prob", "probability", "success", "theta" };
        if (!find_param(params, count, pKeys, 5, &p)) return false;
        auto* h = new (std::nothrow) bernoulli_d_handle(p);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<bernoulli_d_handle, double>;
        out->logpdf = &logpdf_fn<bernoulli_d_handle, double>;
        out->cdf = &cdf_fn<bernoulli_d_handle, double>;
        out->sf = &sf_fn<bernoulli_d_handle, double>;
        out->hazard = &hazard_fn<bernoulli_d_handle, double>;
        out->chf = &chf_fn<bernoulli_d_handle, double>;
        out->quantile = &quantile_fn<bernoulli_d_handle, double>;
        out->quantile_complement = &quantile_complement_fn<bernoulli_d_handle, double>;
        out->range = &range_fn<bernoulli_d_handle, bs_range_d, double>;
        out->mean = &mean_fn<bernoulli_d_handle, double>;
        out->variance = &variance_fn<bernoulli_d_handle, double>;
        out->skewness = &skewness_fn<bernoulli_d_handle, double>;
        out->kurtosis = &kurtosis_fn<bernoulli_d_handle, double>;
        out->kurtosis_excess = &kurtosis_excess_fn<bernoulli_d_handle, double>;
        out->mode = &mode_fn<bernoulli_d_handle, double>;
        out->median = &median_fn<bernoulli_d_handle, double>;
        out->entropy = nullptr;
        out->free = &free_fn<bernoulli_d_handle>;
        return true;
    } else if (n == "binomial" || n == "binomial_distribution") {
        double nTrials = 0, p = 0;
        const char* nKeys[] = { "n", "trials" };
        const char* pKeys[] = { "p", "prob", "probability", "success" };
        if (!find_param(params, count, nKeys, 2, &nTrials)) return false;
        if (!find_param(params, count, pKeys, 4, &p)) return false;
        if (nTrials < 0) return false;
        auto* h = new (std::nothrow) binomial_d_handle(static_cast<unsigned>(nTrials), p);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<binomial_d_handle, double>;
        out->logpdf = &logpdf_fn<binomial_d_handle, double>;
        out->cdf = &cdf_fn<binomial_d_handle, double>;
        out->sf = &sf_fn<binomial_d_handle, double>;
        out->hazard = &hazard_fn<binomial_d_handle, double>;
        out->chf = &chf_fn<binomial_d_handle, double>;
        out->quantile = &quantile_fn<binomial_d_handle, double>;
        out->quantile_complement = &quantile_complement_fn<binomial_d_handle, double>;
        out->range = &range_fn<binomial_d_handle, bs_range_d, double>;
        out->mean = &mean_fn<binomial_d_handle, double>;
        out->variance = &variance_fn<binomial_d_handle, double>;
        out->skewness = &skewness_fn<binomial_d_handle, double>;
        out->kurtosis = &kurtosis_fn<binomial_d_handle, double>;
        out->kurtosis_excess = &kurtosis_excess_fn<binomial_d_handle, double>;
        out->mode = &mode_fn<binomial_d_handle, double>;
        out->median = &median_fn<binomial_d_handle, double>;
        out->entropy = nullptr;
        out->free = &free_fn<binomial_d_handle>;
        return true;
    } else if (n == "cauchy" || n == "cauchy_distribution") {
        double loc = 0, scale = 0;
        const char* locKeys[] = { "location", "loc", "mu", "median", "x0" };
        const char* scaleKeys[] = { "scale", "gamma", "sigma", "b" };
        find_param(params, count, locKeys, 5, &loc);
        if (!find_param(params, count, scaleKeys, 4, &scale)) return false;
        auto* h = new (std::nothrow) cauchy_d_handle(loc, scale);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<cauchy_d_handle, double>;
        out->logpdf = &logpdf_fn<cauchy_d_handle, double>;
        out->cdf = &cdf_fn<cauchy_d_handle, double>;
        out->sf = &sf_fn<cauchy_d_handle, double>;
        out->hazard = &hazard_fn<cauchy_d_handle, double>;
        out->chf = &chf_fn<cauchy_d_handle, double>;
        out->quantile = &quantile_fn<cauchy_d_handle, double>;
        out->quantile_complement = &quantile_complement_fn<cauchy_d_handle, double>;
        out->range = &range_fn<cauchy_d_handle, bs_range_d, double>;
        out->mean = nullptr;
        out->variance = nullptr;
        out->skewness = nullptr;
        out->kurtosis = nullptr;
        out->kurtosis_excess = nullptr;
        out->mode = &mode_fn<cauchy_d_handle, double>;
        out->median = &median_fn<cauchy_d_handle, double>;
        out->entropy = &entropy_fn<cauchy_d_handle, double>;
        out->free = &free_fn<cauchy_d_handle>;
        return true;
    } else if (n == "exponential" || n == "exponential_distribution" || n == "exp") {
        double lambda = 0, scale = 0;
        const char* lamKeys[] = { "lambda", "rate" };
        const char* scKeys[] = { "scale", "theta" };
        bool haveLambda = find_param(params, count, lamKeys, 2, &lambda);
        bool haveScale = find_param(params, count, scKeys, 2, &scale);
        if (!haveLambda) {
            if (!haveScale) return false;
            if (scale == 0) return false;
            lambda = 1.0 / scale;
        }
        auto* h = new (std::nothrow) exponential_d_handle(lambda);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<exponential_d_handle, double>;
        out->logpdf = &logpdf_fn<exponential_d_handle, double>;
        out->cdf = &cdf_fn<exponential_d_handle, double>;
        out->sf = &sf_fn<exponential_d_handle, double>;
        out->hazard = &hazard_fn<exponential_d_handle, double>;
        out->chf = &chf_fn<exponential_d_handle, double>;
        out->quantile = &quantile_fn<exponential_d_handle, double>;
        out->quantile_complement = &quantile_complement_fn<exponential_d_handle, double>;
        out->range = &range_fn<exponential_d_handle, bs_range_d, double>;
        out->mean = &mean_fn<exponential_d_handle, double>;
        out->variance = &variance_fn<exponential_d_handle, double>;
        out->skewness = &skewness_fn<exponential_d_handle, double>;
        out->kurtosis = &kurtosis_fn<exponential_d_handle, double>;
        out->kurtosis_excess = &kurtosis_excess_fn<exponential_d_handle, double>;
        out->mode = &mode_fn<exponential_d_handle, double>;
        out->median = &median_fn<exponential_d_handle, double>;
        out->entropy = &entropy_fn<exponential_d_handle, double>;
        out->free = &free_fn<exponential_d_handle>;
        return true;
    } else if (n == "extremevalue" || n == "extreme_value" || n == "gumbel" || n == "extreme_value_distribution") {
        double loc = 0, scale = 0;
        const char* locKeys[] = { "location", "loc", "mu" };
        const char* scaleKeys[] = { "scale", "gamma", "sigma", "b" };
        find_param(params, count, locKeys, 3, &loc);
        if (!find_param(params, count, scaleKeys, 4, &scale)) return false;
        auto* h = new (std::nothrow) extreme_value_d_handle(loc, scale);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<extreme_value_d_handle, double>;
        out->logpdf = &logpdf_fn<extreme_value_d_handle, double>;
        out->cdf = &cdf_fn<extreme_value_d_handle, double>;
        out->sf = &sf_fn<extreme_value_d_handle, double>;
        out->hazard = &hazard_fn<extreme_value_d_handle, double>;
        out->chf = &chf_fn<extreme_value_d_handle, double>;
        out->quantile = &quantile_fn<extreme_value_d_handle, double>;
        out->quantile_complement = &quantile_complement_fn<extreme_value_d_handle, double>;
        out->range = &range_fn<extreme_value_d_handle, bs_range_d, double>;
        out->mean = &mean_fn<extreme_value_d_handle, double>;
        out->variance = &variance_fn<extreme_value_d_handle, double>;
        out->skewness = &skewness_fn<extreme_value_d_handle, double>;
        out->kurtosis = &kurtosis_fn<extreme_value_d_handle, double>;
        out->kurtosis_excess = &kurtosis_excess_fn<extreme_value_d_handle, double>;
        out->mode = &mode_fn<extreme_value_d_handle, double>;
        out->median = &median_fn<extreme_value_d_handle, double>;
        out->entropy = nullptr;
        out->free = &free_fn<extreme_value_d_handle>;
        return true;
    } else if (n == "geometric" || n == "geometric_distribution") {
        double p = 0;
        const char* pKeys[] = { "p", "prob", "probability", "success", "theta" };
        if (!find_param(params, count, pKeys, 5, &p)) return false;
        auto* h = new (std::nothrow) geometric_d_handle(p);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<geometric_d_handle, double>;
        out->logpdf = &logpdf_fn<geometric_d_handle, double>;
        out->cdf = &cdf_fn<geometric_d_handle, double>;
        out->sf = &sf_fn<geometric_d_handle, double>;
        out->hazard = &hazard_fn<geometric_d_handle, double>;
        out->chf = &chf_fn<geometric_d_handle, double>;
        out->quantile = &quantile_fn<geometric_d_handle, double>;
        out->quantile_complement = &quantile_complement_fn<geometric_d_handle, double>;
        out->range = &range_fn<geometric_d_handle, bs_range_d, double>;
        out->mean = &mean_fn<geometric_d_handle, double>;
        out->variance = &variance_fn<geometric_d_handle, double>;
        out->skewness = &skewness_fn<geometric_d_handle, double>;
        out->kurtosis = &kurtosis_fn<geometric_d_handle, double>;
        out->kurtosis_excess = &kurtosis_excess_fn<geometric_d_handle, double>;
        out->mode = &mode_fn<geometric_d_handle, double>;
        out->median = &median_fn<geometric_d_handle, double>;
        out->entropy = nullptr;
        out->free = &free_fn<geometric_d_handle>;
        return true;
    } else if (n == "holtsmark" || n == "holtsmark_distribution") {
        double loc = 0, scale = 0;
        const char* locKeys[] = { "location", "loc", "mu", "median", "x0" };
        const char* scaleKeys[] = { "scale", "gamma", "sigma", "b" };
        find_param(params, count, locKeys, 5, &loc);
        if (!find_param(params, count, scaleKeys, 4, &scale)) return false;
        auto* h = new (std::nothrow) holtsmark_d_handle(loc, scale);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<holtsmark_d_handle, double>;
        out->logpdf = &logpdf_fn<holtsmark_d_handle, double>;
        out->cdf = &cdf_fn<holtsmark_d_handle, double>;
        out->sf = &sf_fn<holtsmark_d_handle, double>;
        out->hazard = &hazard_fn<holtsmark_d_handle, double>;
        out->chf = &chf_fn<holtsmark_d_handle, double>;
        out->quantile = &quantile_fn<holtsmark_d_handle, double>;
        out->quantile_complement = &quantile_complement_fn<holtsmark_d_handle, double>;
        out->range = &range_fn<holtsmark_d_handle, bs_range_d, double>;
        out->mean = &mean_fn<holtsmark_d_handle, double>;
        out->variance = &variance_fn<holtsmark_d_handle, double>;
        out->skewness = nullptr;
        out->kurtosis = nullptr;
        out->kurtosis_excess = nullptr;
        out->mode = &mode_fn<holtsmark_d_handle, double>;
        out->median = &median_fn<holtsmark_d_handle, double>;
        out->entropy = nullptr;
        out->free = &free_fn<holtsmark_d_handle>;
        return true;
    }
        return false;
    } catch (...) {
        return false;
    }
}

bool bs_dist_make_f(const char* name, const bs_param_f* params, size_t count, bs_dist_f* out) {
    if (!name || !out) return false;
    try {
        std::string n = lower_ascii(name);
    if (n == "gamma" || n == "gamma_distribution") {
        float k = 0, th = 1;
        const char* kKeys[] = { "shape", "k" };
        const char* thKeys[] = { "scale", "theta" };
        if (!find_param(params, count, kKeys, 2, &k)) return false;
        find_param(params, count, thKeys, 2, &th);
        auto* h = new (std::nothrow) gamma_f_handle(k, th);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<gamma_f_handle, float>;
        out->logpdf = &logpdf_fn<gamma_f_handle, float>;
        out->cdf = &cdf_fn<gamma_f_handle, float>;
        out->sf = &sf_fn<gamma_f_handle, float>;
        out->hazard = &hazard_fn<gamma_f_handle, float>;
        out->chf = &chf_fn<gamma_f_handle, float>;
        out->quantile = &quantile_fn<gamma_f_handle, float>;
        out->quantile_complement = &quantile_complement_fn<gamma_f_handle, float>;
        out->range = &range_fn<gamma_f_handle, bs_range_f, float>;
        out->mean = &mean_fn<gamma_f_handle, float>;
        out->variance = &variance_fn<gamma_f_handle, float>;
        out->skewness = &skewness_fn<gamma_f_handle, float>;
        out->kurtosis = &kurtosis_fn<gamma_f_handle, float>;
        out->kurtosis_excess = &kurtosis_excess_fn<gamma_f_handle, float>;
        out->mode = &mode_fn<gamma_f_handle, float>;
        out->median = &median_fn<gamma_f_handle, float>;
        out->entropy = &entropy_fn<gamma_f_handle, float>;
        out->free = &free_fn<gamma_f_handle>;
        return true;
    } else if (n == "studentt" || n == "student_t" || n == "students_t" || n == "t" || n == "t_distribution") {
        float v = 0;
        const char* dfKeys[] = { "df", "nu", "degreesoffreedom" };
        if (!find_param(params, count, dfKeys, 3, &v)) return false;
        auto* h = new (std::nothrow) student_t_f_handle(v);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<student_t_f_handle, float>;
        out->logpdf = &logpdf_fn<student_t_f_handle, float>;
        out->cdf = &cdf_fn<student_t_f_handle, float>;
        out->sf = &sf_fn<student_t_f_handle, float>;
        out->hazard = &hazard_fn<student_t_f_handle, float>;
        out->chf = &chf_fn<student_t_f_handle, float>;
        out->quantile = &quantile_fn<student_t_f_handle, float>;
        out->quantile_complement = &quantile_complement_fn<student_t_f_handle, float>;
        out->range = &range_fn<student_t_f_handle, bs_range_f, float>;
        out->mean = &mean_fn<student_t_f_handle, float>;
        out->variance = &variance_fn<student_t_f_handle, float>;
        out->skewness = &skewness_fn<student_t_f_handle, float>;
        out->kurtosis = &kurtosis_fn<student_t_f_handle, float>;
        out->kurtosis_excess = &kurtosis_excess_fn<student_t_f_handle, float>;
        out->mode = &mode_fn<student_t_f_handle, float>;
        out->median = &median_fn<student_t_f_handle, float>;
        out->entropy = &entropy_fn<student_t_f_handle, float>;
        out->free = &free_fn<student_t_f_handle>;
        return true;
    } else if (n == "fisherf" || n == "f" || n == "f_distribution") {
        float df1 = 0, df2 = 0;
        const char* df1Keys[] = { "df1", "d1", "m", "degreesoffreedom1" };
        const char* df2Keys[] = { "df2", "d2", "n", "degreesoffreedom2" };
        if (!find_param(params, count, df1Keys, 4, &df1)) return false;
        if (!find_param(params, count, df2Keys, 4, &df2)) return false;
        auto* h = new (std::nothrow) fisher_f_f_handle(df1, df2);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<fisher_f_f_handle, float>;
        out->logpdf = &logpdf_fn<fisher_f_f_handle, float>;
        out->cdf = &cdf_fn<fisher_f_f_handle, float>;
        out->sf = &sf_fn<fisher_f_f_handle, float>;
        out->hazard = &hazard_fn<fisher_f_f_handle, float>;
        out->chf = &chf_fn<fisher_f_f_handle, float>;
        out->quantile = &quantile_fn<fisher_f_f_handle, float>;
        out->quantile_complement = &quantile_complement_fn<fisher_f_f_handle, float>;
        out->range = &range_fn<fisher_f_f_handle, bs_range_f, float>;
        out->mean = &mean_fn<fisher_f_f_handle, float>;
        out->variance = &variance_fn<fisher_f_f_handle, float>;
        out->skewness = &skewness_fn<fisher_f_f_handle, float>;
        out->kurtosis = &kurtosis_fn<fisher_f_f_handle, float>;
        out->kurtosis_excess = &kurtosis_excess_fn<fisher_f_f_handle, float>;
        out->mode = &mode_fn<fisher_f_f_handle, float>;
        out->median = &median_fn<fisher_f_f_handle, float>;
        out->entropy = nullptr;
        out->free = &free_fn<fisher_f_f_handle>;
        return true;
    } else if (n == "arcsine" || n == "arcsine_distribution") {
        float a = 0, b = 0;
        const char* aKeys[] = { "minx", "min", "a", "lower" };
        const char* bKeys[] = { "maxx", "max", "b", "upper" };
        if (!find_param(params, count, aKeys, 4, &a)) return false;
        if (!find_param(params, count, bKeys, 4, &b)) return false;
        auto* h = new (std::nothrow) arcsine_f_handle(a, b);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<arcsine_f_handle, float>;
        out->logpdf = &logpdf_fn<arcsine_f_handle, float>;
        out->cdf = &cdf_fn<arcsine_f_handle, float>;
        out->sf = &sf_fn<arcsine_f_handle, float>;
        out->hazard = &hazard_fn<arcsine_f_handle, float>;
        out->chf = &chf_fn<arcsine_f_handle, float>;
        out->quantile = &quantile_fn<arcsine_f_handle, float>;
        out->quantile_complement = &quantile_complement_fn<arcsine_f_handle, float>;
        out->range = &range_fn<arcsine_f_handle, bs_range_f, float>;
        out->mean = &mean_fn<arcsine_f_handle, float>;
        out->variance = &variance_fn<arcsine_f_handle, float>;
        out->skewness = &skewness_fn<arcsine_f_handle, float>;
        out->kurtosis = &kurtosis_fn<arcsine_f_handle, float>;
        out->kurtosis_excess = &kurtosis_excess_fn<arcsine_f_handle, float>;
        out->mode = &mode_fn<arcsine_f_handle, float>;
        out->median = &median_fn<arcsine_f_handle, float>;
        out->entropy = nullptr;
        out->free = &free_fn<arcsine_f_handle>;
        return true;
    } else if (n == "beta" || n == "beta_distribution") {
        double a = 0, b = 0;
        const char* aKeys[] = { "alpha", "a", "p", "shape1" };
        const char* bKeys[] = { "beta", "b", "q", "shape2" };
        if (!find_param(params, count, aKeys, 4, &a)) return false;
        if (!find_param(params, count, bKeys, 4, &b)) return false;
        auto* h = new (std::nothrow) beta_d_handle(a, b);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<beta_f_handle, float>;
        out->logpdf = &logpdf_fn<beta_f_handle, float>;
        out->cdf = &cdf_fn<beta_f_handle, float>;
        out->sf = &sf_fn<beta_f_handle, float>;
        out->hazard = &hazard_fn<beta_f_handle, float>;
        out->chf = &chf_fn<beta_f_handle, float>;
        out->quantile = &quantile_fn<beta_f_handle, float>;
        out->quantile_complement = &quantile_complement_fn<beta_f_handle, float>;
        out->range = &range_fn<beta_f_handle, bs_range_f, float>;
        out->mean = &mean_fn<beta_f_handle, float>;
        out->variance = &variance_fn<beta_f_handle, float>;
        out->skewness = &skewness_fn<beta_f_handle, float>;
        out->kurtosis = &kurtosis_fn<beta_f_handle, float>;
        out->kurtosis_excess = &kurtosis_excess_fn<beta_f_handle, float>;
        out->mode = &mode_fn<beta_f_handle, float>;
        out->median = &median_fn<beta_f_handle, float>;
        out->entropy = nullptr; // not provided in Boost
        out->free = &free_fn<beta_f_handle>;
        return true;
    } else if (n == "chisquared" || n == "chi_squared" || n == "chi2" || n == "chi-squared" || n == "chisquare") {
        float v = 0;
        const char* dfKeys[] = { "df", "nu", "degreesoffreedom" };
        if (!find_param(params, count, dfKeys, 3, &v)) return false;
        auto* h = new (std::nothrow) chi_squared_f_handle(v);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<chi_squared_f_handle, float>;
        out->logpdf = &logpdf_fn<chi_squared_f_handle, float>;
        out->cdf = &cdf_fn<chi_squared_f_handle, float>;
        out->sf = &sf_fn<chi_squared_f_handle, float>;
        out->hazard = &hazard_fn<chi_squared_f_handle, float>;
        out->chf = &chf_fn<chi_squared_f_handle, float>;
        out->quantile = &quantile_fn<chi_squared_f_handle, float>;
        out->quantile_complement = &quantile_complement_fn<chi_squared_f_handle, float>;
        out->range = &range_fn<chi_squared_f_handle, bs_range_f, float>;
        out->mean = &mean_fn<chi_squared_f_handle, float>;
        out->variance = &variance_fn<chi_squared_f_handle, float>;
        out->skewness = &skewness_fn<chi_squared_f_handle, float>;
        out->kurtosis = &kurtosis_fn<chi_squared_f_handle, float>;
        out->kurtosis_excess = &kurtosis_excess_fn<chi_squared_f_handle, float>;
        out->mode = &mode_fn<chi_squared_f_handle, float>;
        out->median = &median_fn<chi_squared_f_handle, float>;
        out->entropy = nullptr; // not provided in Boost for chi-squared
        out->free = &free_fn<chi_squared_f_handle>;
        return true;
    } else if (n == "bernoulli" || n == "bernoulli_distribution") {
        float p = 0;
        const char* pKeys[] = { "p", "prob", "probability", "success", "theta" };
        if (!find_param(params, count, pKeys, 5, &p)) return false;
        auto* h = new (std::nothrow) bernoulli_f_handle(p);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<bernoulli_f_handle, float>;
        out->logpdf = &logpdf_fn<bernoulli_f_handle, float>;
        out->cdf = &cdf_fn<bernoulli_f_handle, float>;
        out->sf = &sf_fn<bernoulli_f_handle, float>;
        out->hazard = &hazard_fn<bernoulli_f_handle, float>;
        out->chf = &chf_fn<bernoulli_f_handle, float>;
        out->quantile = &quantile_fn<bernoulli_f_handle, float>;
        out->quantile_complement = &quantile_complement_fn<bernoulli_f_handle, float>;
        out->range = &range_fn<bernoulli_f_handle, bs_range_f, float>;
        out->mean = &mean_fn<bernoulli_f_handle, float>;
        out->variance = &variance_fn<bernoulli_f_handle, float>;
        out->skewness = &skewness_fn<bernoulli_f_handle, float>;
        out->kurtosis = &kurtosis_fn<bernoulli_f_handle, float>;
        out->kurtosis_excess = &kurtosis_excess_fn<bernoulli_f_handle, float>;
        out->mode = &mode_fn<bernoulli_f_handle, float>;
        out->median = &median_fn<bernoulli_f_handle, float>;
        out->entropy = nullptr;
        out->free = &free_fn<bernoulli_f_handle>;
        return true;
    } else if (n == "binomial" || n == "binomial_distribution") {
        float nTrials = 0, p = 0;
        const char* nKeys[] = { "n", "trials" };
        const char* pKeys[] = { "p", "prob", "probability", "success" };
        if (!find_param(params, count, nKeys, 2, &nTrials)) return false;
        if (!find_param(params, count, pKeys, 4, &p)) return false;
        if (nTrials < 0) return false;
        auto* h = new (std::nothrow) binomial_f_handle(static_cast<unsigned>(nTrials), p);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<binomial_f_handle, float>;
        out->logpdf = &logpdf_fn<binomial_f_handle, float>;
        out->cdf = &cdf_fn<binomial_f_handle, float>;
        out->sf = &sf_fn<binomial_f_handle, float>;
        out->hazard = &hazard_fn<binomial_f_handle, float>;
        out->chf = &chf_fn<binomial_f_handle, float>;
        out->quantile = &quantile_fn<binomial_f_handle, float>;
        out->quantile_complement = &quantile_complement_fn<binomial_f_handle, float>;
        out->range = &range_fn<binomial_f_handle, bs_range_f, float>;
        out->mean = &mean_fn<binomial_f_handle, float>;
        out->variance = &variance_fn<binomial_f_handle, float>;
        out->skewness = &skewness_fn<binomial_f_handle, float>;
        out->kurtosis = &kurtosis_fn<binomial_f_handle, float>;
        out->kurtosis_excess = &kurtosis_excess_fn<binomial_f_handle, float>;
        out->mode = &mode_fn<binomial_f_handle, float>;
        out->median = &median_fn<binomial_f_handle, float>;
        out->entropy = nullptr;
        out->free = &free_fn<binomial_f_handle>;
        return true;
    } else if (n == "cauchy" || n == "cauchy_distribution") {
        float loc = 0, scale = 0;
        const char* locKeys[] = { "location", "loc", "mu", "median", "x0" };
        const char* scaleKeys[] = { "scale", "gamma", "sigma", "b" };
        find_param(params, count, locKeys, 5, &loc);
        if (!find_param(params, count, scaleKeys, 4, &scale)) return false;
        auto* h = new (std::nothrow) cauchy_f_handle(loc, scale);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<cauchy_f_handle, float>;
        out->logpdf = &logpdf_fn<cauchy_f_handle, float>;
        out->cdf = &cdf_fn<cauchy_f_handle, float>;
        out->sf = &sf_fn<cauchy_f_handle, float>;
        out->hazard = &hazard_fn<cauchy_f_handle, float>;
        out->chf = &chf_fn<cauchy_f_handle, float>;
        out->quantile = &quantile_fn<cauchy_f_handle, float>;
        out->quantile_complement = &quantile_complement_fn<cauchy_f_handle, float>;
        out->range = &range_fn<cauchy_f_handle, bs_range_f, float>;
        out->mean = nullptr;
        out->variance = nullptr;
        out->skewness = nullptr;
        out->kurtosis = nullptr;
        out->kurtosis_excess = nullptr;
        out->mode = &mode_fn<cauchy_f_handle, float>;
        out->median = &median_fn<cauchy_f_handle, float>;
        out->entropy = &entropy_fn<cauchy_f_handle, float>;
        out->free = &free_fn<cauchy_f_handle>;
        return true;
    } else if (n == "exponential" || n == "exponential_distribution" || n == "exp") {
        float lambda = 0, scale = 0;
        const char* lamKeys[] = { "lambda", "rate" };
        const char* scKeys[] = { "scale", "theta" };
        bool haveLambda = find_param(params, count, lamKeys, 2, &lambda);
        bool haveScale = find_param(params, count, scKeys, 2, &scale);
        if (!haveLambda) {
            if (!haveScale) return false;
            if (scale == 0) return false;
            lambda = 1.0f / scale;
        }
        auto* h = new (std::nothrow) exponential_f_handle(lambda);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<exponential_f_handle, float>;
        out->logpdf = &logpdf_fn<exponential_f_handle, float>;
        out->cdf = &cdf_fn<exponential_f_handle, float>;
        out->sf = &sf_fn<exponential_f_handle, float>;
        out->hazard = &hazard_fn<exponential_f_handle, float>;
        out->chf = &chf_fn<exponential_f_handle, float>;
        out->quantile = &quantile_fn<exponential_f_handle, float>;
        out->quantile_complement = &quantile_complement_fn<exponential_f_handle, float>;
        out->range = &range_fn<exponential_f_handle, bs_range_f, float>;
        out->mean = &mean_fn<exponential_f_handle, float>;
        out->variance = &variance_fn<exponential_f_handle, float>;
        out->skewness = &skewness_fn<exponential_f_handle, float>;
        out->kurtosis = &kurtosis_fn<exponential_f_handle, float>;
        out->kurtosis_excess = &kurtosis_excess_fn<exponential_f_handle, float>;
        out->mode = &mode_fn<exponential_f_handle, float>;
        out->median = &median_fn<exponential_f_handle, float>;
        out->entropy = &entropy_fn<exponential_f_handle, float>;
        out->free = &free_fn<exponential_f_handle>;
        return true;
    } else if (n == "extremevalue" || n == "extreme_value" || n == "gumbel" || n == "extreme_value_distribution") {
        float loc = 0, scale = 0;
        const char* locKeys[] = { "location", "loc", "mu" };
        const char* scaleKeys[] = { "scale", "gamma", "sigma", "b" };
        find_param(params, count, locKeys, 3, &loc);
        if (!find_param(params, count, scaleKeys, 4, &scale)) return false;
        auto* h = new (std::nothrow) extreme_value_f_handle(loc, scale);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<extreme_value_f_handle, float>;
        out->logpdf = &logpdf_fn<extreme_value_f_handle, float>;
        out->cdf = &cdf_fn<extreme_value_f_handle, float>;
        out->sf = &sf_fn<extreme_value_f_handle, float>;
        out->hazard = &hazard_fn<extreme_value_f_handle, float>;
        out->chf = &chf_fn<extreme_value_f_handle, float>;
        out->quantile = &quantile_fn<extreme_value_f_handle, float>;
        out->quantile_complement = &quantile_complement_fn<extreme_value_f_handle, float>;
        out->range = &range_fn<extreme_value_f_handle, bs_range_f, float>;
        out->mean = &mean_fn<extreme_value_f_handle, float>;
        out->variance = &variance_fn<extreme_value_f_handle, float>;
        out->skewness = &skewness_fn<extreme_value_f_handle, float>;
        out->kurtosis = &kurtosis_fn<extreme_value_f_handle, float>;
        out->kurtosis_excess = &kurtosis_excess_fn<extreme_value_f_handle, float>;
        out->mode = &mode_fn<extreme_value_f_handle, float>;
        out->median = &median_fn<extreme_value_f_handle, float>;
        out->entropy = nullptr;
        out->free = &free_fn<extreme_value_f_handle>;
        return true;
    } else if (n == "geometric" || n == "geometric_distribution") {
        float p = 0;
        const char* pKeys[] = { "p", "prob", "probability", "success", "theta" };
        if (!find_param(params, count, pKeys, 5, &p)) return false;
        auto* h = new (std::nothrow) geometric_f_handle(p);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<geometric_f_handle, float>;
        out->logpdf = &logpdf_fn<geometric_f_handle, float>;
        out->cdf = &cdf_fn<geometric_f_handle, float>;
        out->sf = &sf_fn<geometric_f_handle, float>;
        out->hazard = &hazard_fn<geometric_f_handle, float>;
        out->chf = &chf_fn<geometric_f_handle, float>;
        out->quantile = &quantile_fn<geometric_f_handle, float>;
        out->quantile_complement = &quantile_complement_fn<geometric_f_handle, float>;
        out->range = &range_fn<geometric_f_handle, bs_range_f, float>;
        out->mean = &mean_fn<geometric_f_handle, float>;
        out->variance = &variance_fn<geometric_f_handle, float>;
        out->skewness = &skewness_fn<geometric_f_handle, float>;
        out->kurtosis = &kurtosis_fn<geometric_f_handle, float>;
        out->kurtosis_excess = &kurtosis_excess_fn<geometric_f_handle, float>;
        out->mode = &mode_fn<geometric_f_handle, float>;
        out->median = &median_fn<geometric_f_handle, float>;
        out->entropy = nullptr;
        out->free = &free_fn<geometric_f_handle>;
        return true;
    } else if (n == "holtsmark" || n == "holtsmark_distribution") {
        float loc = 0, scale = 0;
        const char* locKeys[] = { "location", "loc", "mu", "median", "x0" };
        const char* scaleKeys[] = { "scale", "gamma", "sigma", "b" };
        find_param(params, count, locKeys, 5, &loc);
        if (!find_param(params, count, scaleKeys, 4, &scale)) return false;
        auto* h = new (std::nothrow) holtsmark_f_handle(loc, scale);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<holtsmark_f_handle, float>;
        out->logpdf = &logpdf_fn<holtsmark_f_handle, float>;
        out->cdf = &cdf_fn<holtsmark_f_handle, float>;
        out->sf = &sf_fn<holtsmark_f_handle, float>;
        out->hazard = &hazard_fn<holtsmark_f_handle, float>;
        out->chf = &chf_fn<holtsmark_f_handle, float>;
        out->quantile = &quantile_fn<holtsmark_f_handle, float>;
        out->quantile_complement = &quantile_complement_fn<holtsmark_f_handle, float>;
        out->range = &range_fn<holtsmark_f_handle, bs_range_f, float>;
        out->mean = &mean_fn<holtsmark_f_handle, float>;
        out->variance = &variance_fn<holtsmark_f_handle, float>;
        out->skewness = nullptr;
        out->kurtosis = nullptr;
        out->kurtosis_excess = nullptr;
        out->mode = &mode_fn<holtsmark_f_handle, float>;
        out->median = &median_fn<holtsmark_f_handle, float>;
        out->entropy = nullptr;
        out->free = &free_fn<holtsmark_f_handle>;
        return true;
    }

        return false;
    } catch (...) {
        return false;
    }
}

bool bs_dist_make_l(const char* name, const bs_param_l* params, size_t count, bs_dist_l* out) {
    if (!name || !out) return false;
    try {
        std::string n = lower_ascii(name);
    if (n == "gamma" || n == "gamma_distribution") {
        long double k = 0, th = 1;
        const char* kKeys[] = { "shape", "k" };
        const char* thKeys[] = { "scale", "theta" };
        if (!find_param(params, count, kKeys, 2, &k)) return false;
        find_param(params, count, thKeys, 2, &th);
        auto* h = new (std::nothrow) gamma_l_handle(k, th);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<gamma_l_handle, long double>;
        out->logpdf = &logpdf_fn<gamma_l_handle, long double>;
        out->cdf = &cdf_fn<gamma_l_handle, long double>;
        out->sf = &sf_fn<gamma_l_handle, long double>;
        out->hazard = &hazard_fn<gamma_l_handle, long double>;
        out->chf = &chf_fn<gamma_l_handle, long double>;
        out->quantile = &quantile_fn<gamma_l_handle, long double>;
        out->quantile_complement = &quantile_complement_fn<gamma_l_handle, long double>;
        out->range = &range_fn<gamma_l_handle, bs_range_l, long double>;
        out->mean = &mean_fn<gamma_l_handle, long double>;
        out->variance = &variance_fn<gamma_l_handle, long double>;
        out->skewness = &skewness_fn<gamma_l_handle, long double>;
        out->kurtosis = &kurtosis_fn<gamma_l_handle, long double>;
        out->kurtosis_excess = &kurtosis_excess_fn<gamma_l_handle, long double>;
        out->mode = &mode_fn<gamma_l_handle, long double>;
        out->median = &median_fn<gamma_l_handle, long double>;
        out->entropy = &entropy_fn<gamma_l_handle, long double>;
        out->free = &free_fn<gamma_l_handle>;
        return true;
    } else if (n == "studentt" || n == "student_t" || n == "students_t" || n == "t" || n == "t_distribution") {
        long double v = 0;
        const char* dfKeys[] = { "df", "nu", "degreesoffreedom" };
        if (!find_param(params, count, dfKeys, 3, &v)) return false;
        auto* h = new (std::nothrow) student_t_l_handle(v);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<student_t_l_handle, long double>;
        out->logpdf = &logpdf_fn<student_t_l_handle, long double>;
        out->cdf = &cdf_fn<student_t_l_handle, long double>;
        out->sf = &sf_fn<student_t_l_handle, long double>;
        out->hazard = &hazard_fn<student_t_l_handle, long double>;
        out->chf = &chf_fn<student_t_l_handle, long double>;
        out->quantile = &quantile_fn<student_t_l_handle, long double>;
        out->quantile_complement = &quantile_complement_fn<student_t_l_handle, long double>;
        out->range = &range_fn<student_t_l_handle, bs_range_l, long double>;
        out->mean = &mean_fn<student_t_l_handle, long double>;
        out->variance = &variance_fn<student_t_l_handle, long double>;
        out->skewness = &skewness_fn<student_t_l_handle, long double>;
        out->kurtosis = &kurtosis_fn<student_t_l_handle, long double>;
        out->kurtosis_excess = &kurtosis_excess_fn<student_t_l_handle, long double>;
        out->mode = &mode_fn<student_t_l_handle, long double>;
        out->median = &median_fn<student_t_l_handle, long double>;
        out->entropy = &entropy_fn<student_t_l_handle, long double>;
        out->free = &free_fn<student_t_l_handle>;
        return true;
    } else if (n == "fisherf" || n == "f" || n == "f_distribution") {
        long double df1 = 0, df2 = 0;
        const char* df1Keys[] = { "df1", "d1", "m", "degreesoffreedom1" };
        const char* df2Keys[] = { "df2", "d2", "n", "degreesoffreedom2" };
        if (!find_param(params, count, df1Keys, 4, &df1)) return false;
        if (!find_param(params, count, df2Keys, 4, &df2)) return false;
        auto* h = new (std::nothrow) fisher_f_l_handle(df1, df2);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<fisher_f_l_handle, long double>;
        out->logpdf = &logpdf_fn<fisher_f_l_handle, long double>;
        out->cdf = &cdf_fn<fisher_f_l_handle, long double>;
        out->sf = &sf_fn<fisher_f_l_handle, long double>;
        out->hazard = &hazard_fn<fisher_f_l_handle, long double>;
        out->chf = &chf_fn<fisher_f_l_handle, long double>;
        out->quantile = &quantile_fn<fisher_f_l_handle, long double>;
        out->quantile_complement = &quantile_complement_fn<fisher_f_l_handle, long double>;
        out->range = &range_fn<fisher_f_l_handle, bs_range_l, long double>;
        out->mean = &mean_fn<fisher_f_l_handle, long double>;
        out->variance = &variance_fn<fisher_f_l_handle, long double>;
        out->skewness = &skewness_fn<fisher_f_l_handle, long double>;
        out->kurtosis = &kurtosis_fn<fisher_f_l_handle, long double>;
        out->kurtosis_excess = &kurtosis_excess_fn<fisher_f_l_handle, long double>;
        out->mode = &mode_fn<fisher_f_l_handle, long double>;
        out->median = &median_fn<fisher_f_l_handle, long double>;
        out->entropy = nullptr;
        out->free = &free_fn<fisher_f_l_handle>;
        return true;
    } else if (n == "arcsine" || n == "arcsine_distribution") {
        long double a = 0, b = 0;
        const char* aKeys[] = { "minx", "min", "a", "lower" };
        const char* bKeys[] = { "maxx", "max", "b", "upper" };
        if (!find_param(params, count, aKeys, 4, &a)) return false;
        if (!find_param(params, count, bKeys, 4, &b)) return false;
        auto* h = new (std::nothrow) arcsine_l_handle(a, b);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<arcsine_l_handle, long double>;
        out->logpdf = &logpdf_fn<arcsine_l_handle, long double>;
        out->cdf = &cdf_fn<arcsine_l_handle, long double>;
        out->sf = &sf_fn<arcsine_l_handle, long double>;
        out->hazard = &hazard_fn<arcsine_l_handle, long double>;
        out->chf = &chf_fn<arcsine_l_handle, long double>;
        out->quantile = &quantile_fn<arcsine_l_handle, long double>;
        out->quantile_complement = &quantile_complement_fn<arcsine_l_handle, long double>;
        out->range = &range_fn<arcsine_l_handle, bs_range_l, long double>;
        out->mean = &mean_fn<arcsine_l_handle, long double>;
        out->variance = &variance_fn<arcsine_l_handle, long double>;
        out->skewness = &skewness_fn<arcsine_l_handle, long double>;
        out->kurtosis = &kurtosis_fn<arcsine_l_handle, long double>;
        out->kurtosis_excess = &kurtosis_excess_fn<arcsine_l_handle, long double>;
        out->mode = &mode_fn<arcsine_l_handle, long double>;
        out->median = &median_fn<arcsine_l_handle, long double>;
        out->entropy = nullptr;
        out->free = &free_fn<arcsine_l_handle>;
        return true;
    } else if (n == "beta" || n == "beta_distribution") {
        double a = 0, b = 0;
        const char* aKeys[] = { "alpha", "a", "p", "shape1" };
        const char* bKeys[] = { "beta", "b", "q", "shape2" };
        if (!find_param(params, count, aKeys, 4, &a)) return false;
        if (!find_param(params, count, bKeys, 4, &b)) return false;
        auto* h = new (std::nothrow) beta_d_handle(a, b);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<beta_l_handle, long double>;
        out->logpdf = &logpdf_fn<beta_l_handle, long double>;
        out->cdf = &cdf_fn<beta_l_handle, long double>;
        out->sf = &sf_fn<beta_l_handle, long double>;
        out->hazard = &hazard_fn<beta_l_handle, long double>;
        out->chf = &chf_fn<beta_l_handle, long double>;
        out->quantile = &quantile_fn<beta_l_handle, long double>;
        out->quantile_complement = &quantile_complement_fn<beta_l_handle, long double>;
        out->range = &range_fn<beta_d_handle, bs_range_l, long double>;
        out->mean = &mean_fn<beta_l_handle, long double>;
        out->variance = &variance_fn<beta_l_handle, long double>;
        out->skewness = &skewness_fn<beta_l_handle, long double>;
        out->kurtosis = &kurtosis_fn<beta_l_handle, long double>;
        out->kurtosis_excess = &kurtosis_excess_fn<beta_l_handle, long double>;
        out->mode = &mode_fn<beta_l_handle, long double>;
        out->median = &median_fn<beta_l_handle, long double>;
        out->entropy = nullptr; // not provided in Boost
        out->free = &free_fn<beta_d_handle>;
        return true;
    } else if (n == "chisquared" || n == "chi_squared" || n == "chi2" || n == "chi-squared" || n == "chisquare") {
        long double v = 0;
        const char* dfKeys[] = { "df", "nu", "degreesoffreedom" };
        if (!find_param(params, count, dfKeys, 3, &v)) return false;
        auto* h = new (std::nothrow) chi_squared_l_handle(v);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<chi_squared_l_handle, long double>;
        out->logpdf = &logpdf_fn<chi_squared_l_handle, long double>;
        out->cdf = &cdf_fn<chi_squared_l_handle, long double>;
        out->sf = &sf_fn<chi_squared_l_handle, long double>;
        out->hazard = &hazard_fn<chi_squared_l_handle, long double>;
        out->chf = &chf_fn<chi_squared_l_handle, long double>;
        out->quantile = &quantile_fn<chi_squared_l_handle, long double>;
        out->quantile_complement = &quantile_complement_fn<chi_squared_l_handle, long double>;
        out->range = &range_fn<chi_squared_l_handle, bs_range_l, long double>;
        out->mean = &mean_fn<chi_squared_l_handle, long double>;
        out->variance = &variance_fn<chi_squared_l_handle, long double>;
        out->skewness = &skewness_fn<chi_squared_l_handle, long double>;
        out->kurtosis = &kurtosis_fn<chi_squared_l_handle, long double>;
        out->kurtosis_excess = &kurtosis_excess_fn<chi_squared_l_handle, long double>;
        out->mode = &mode_fn<chi_squared_l_handle, long double>;
        out->median = &median_fn<chi_squared_l_handle, long double>;
        out->entropy = nullptr; // not provided in Boost for chi-squared
        out->free = &free_fn<chi_squared_l_handle>;
        return true;
    } else if (n == "bernoulli" || n == "bernoulli_distribution") {
        long double p = 0;
        const char* pKeys[] = { "p", "prob", "probability", "success", "theta" };
        if (!find_param(params, count, pKeys, 5, &p)) return false;
        auto* h = new (std::nothrow) bernoulli_l_handle(p);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<bernoulli_l_handle, long double>;
        out->logpdf = &logpdf_fn<bernoulli_l_handle, long double>;
        out->cdf = &cdf_fn<bernoulli_l_handle, long double>;
        out->sf = &sf_fn<bernoulli_l_handle, long double>;
        out->hazard = &hazard_fn<bernoulli_l_handle, long double>;
        out->chf = &chf_fn<bernoulli_l_handle, long double>;
        out->quantile = &quantile_fn<bernoulli_l_handle, long double>;
        out->quantile_complement = &quantile_complement_fn<bernoulli_l_handle, long double>;
        out->range = &range_fn<bernoulli_l_handle, bs_range_l, long double>;
        out->mean = &mean_fn<bernoulli_l_handle, long double>;
        out->variance = &variance_fn<bernoulli_l_handle, long double>;
        out->skewness = &skewness_fn<bernoulli_l_handle, long double>;
        out->kurtosis = &kurtosis_fn<bernoulli_l_handle, long double>;
        out->kurtosis_excess = &kurtosis_excess_fn<bernoulli_l_handle, long double>;
        out->mode = &mode_fn<bernoulli_l_handle, long double>;
        out->median = &median_fn<bernoulli_l_handle, long double>;
        out->entropy = nullptr;
        out->free = &free_fn<bernoulli_l_handle>;
        return true;
    } else if (n == "binomial" || n == "binomial_distribution") {
        long double nTrials = 0, p = 0;
        const char* nKeys[] = { "n", "trials" };
        const char* pKeys[] = { "p", "prob", "probability", "success" };
        if (!find_param(params, count, nKeys, 2, &nTrials)) return false;
        if (!find_param(params, count, pKeys, 4, &p)) return false;
        if (nTrials < 0) return false;
        auto* h = new (std::nothrow) binomial_l_handle(static_cast<unsigned>(nTrials), p);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<binomial_l_handle, long double>;
        out->logpdf = &logpdf_fn<binomial_l_handle, long double>;
        out->cdf = &cdf_fn<binomial_l_handle, long double>;
        out->sf = &sf_fn<binomial_l_handle, long double>;
        out->hazard = &hazard_fn<binomial_l_handle, long double>;
        out->chf = &chf_fn<binomial_l_handle, long double>;
        out->quantile = &quantile_fn<binomial_l_handle, long double>;
        out->quantile_complement = &quantile_complement_fn<binomial_l_handle, long double>;
        out->range = &range_fn<binomial_l_handle, bs_range_l, long double>;
        out->mean = &mean_fn<binomial_l_handle, long double>;
        out->variance = &variance_fn<binomial_l_handle, long double>;
        out->skewness = &skewness_fn<binomial_l_handle, long double>;
        out->kurtosis = &kurtosis_fn<binomial_l_handle, long double>;
        out->kurtosis_excess = &kurtosis_excess_fn<binomial_l_handle, long double>;
        out->mode = &mode_fn<binomial_l_handle, long double>;
        out->median = &median_fn<binomial_l_handle, long double>;
        out->entropy = nullptr;
        out->free = &free_fn<binomial_l_handle>;
        return true;
    } else if (n == "cauchy" || n == "cauchy_distribution") {
        long double loc = 0, scale = 0;
        const char* locKeys[] = { "location", "loc", "mu", "median", "x0" };
        const char* scaleKeys[] = { "scale", "gamma", "sigma", "b" };
        find_param(params, count, locKeys, 5, &loc);
        if (!find_param(params, count, scaleKeys, 4, &scale)) return false;
        auto* h = new (std::nothrow) cauchy_l_handle(loc, scale);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<cauchy_l_handle, long double>;
        out->logpdf = &logpdf_fn<cauchy_l_handle, long double>;
        out->cdf = &cdf_fn<cauchy_l_handle, long double>;
        out->sf = &sf_fn<cauchy_l_handle, long double>;
        out->hazard = &hazard_fn<cauchy_l_handle, long double>;
        out->chf = &chf_fn<cauchy_l_handle, long double>;
        out->quantile = &quantile_fn<cauchy_l_handle, long double>;
        out->quantile_complement = &quantile_complement_fn<cauchy_l_handle, long double>;
        out->range = &range_fn<cauchy_l_handle, bs_range_l, long double>;
        out->mean = nullptr;
        out->variance = nullptr;
        out->skewness = nullptr;
        out->kurtosis = nullptr;
        out->kurtosis_excess = nullptr;
        out->mode = &mode_fn<cauchy_l_handle, long double>;
        out->median = &median_fn<cauchy_l_handle, long double>;
        out->entropy = &entropy_fn<cauchy_l_handle, long double>;
        out->free = &free_fn<cauchy_l_handle>;
        return true;
    } else if (n == "exponential" || n == "exponential_distribution" || n == "exp") {
        long double lambda = 0, scale = 0;
        const char* lamKeys[] = { "lambda", "rate" };
        const char* scKeys[] = { "scale", "theta" };
        bool haveLambda = find_param(params, count, lamKeys, 2, &lambda);
        bool haveScale = find_param(params, count, scKeys, 2, &scale);
        if (!haveLambda) {
            if (!haveScale) return false;
            if (scale == 0) return false;
            lambda = 1.0L / scale;
        }
        auto* h = new (std::nothrow) exponential_l_handle(lambda);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<exponential_l_handle, long double>;
        out->logpdf = &logpdf_fn<exponential_l_handle, long double>;
        out->cdf = &cdf_fn<exponential_l_handle, long double>;
        out->sf = &sf_fn<exponential_l_handle, long double>;
        out->hazard = &hazard_fn<exponential_l_handle, long double>;
        out->chf = &chf_fn<exponential_l_handle, long double>;
        out->quantile = &quantile_fn<exponential_l_handle, long double>;
        out->quantile_complement = &quantile_complement_fn<exponential_l_handle, long double>;
        out->range = &range_fn<exponential_l_handle, bs_range_l, long double>;
        out->mean = &mean_fn<exponential_l_handle, long double>;
        out->variance = &variance_fn<exponential_l_handle, long double>;
        out->skewness = &skewness_fn<exponential_l_handle, long double>;
        out->kurtosis = &kurtosis_fn<exponential_l_handle, long double>;
        out->kurtosis_excess = &kurtosis_excess_fn<exponential_l_handle, long double>;
        out->mode = &mode_fn<exponential_l_handle, long double>;
        out->median = &median_fn<exponential_l_handle, long double>;
        out->entropy = &entropy_fn<exponential_l_handle, long double>;
        out->free = &free_fn<exponential_l_handle>;
        return true;
    } else if (n == "extremevalue" || n == "extreme_value" || n == "gumbel" || n == "extreme_value_distribution") {
        long double loc = 0, scale = 0;
        const char* locKeys[] = { "location", "loc", "mu" };
        const char* scaleKeys[] = { "scale", "gamma", "sigma", "b" };
        find_param(params, count, locKeys, 3, &loc);
        if (!find_param(params, count, scaleKeys, 4, &scale)) return false;
        auto* h = new (std::nothrow) extreme_value_l_handle(loc, scale);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<extreme_value_l_handle, long double>;
        out->logpdf = &logpdf_fn<extreme_value_l_handle, long double>;
        out->cdf = &cdf_fn<extreme_value_l_handle, long double>;
        out->sf = &sf_fn<extreme_value_l_handle, long double>;
        out->hazard = &hazard_fn<extreme_value_l_handle, long double>;
        out->chf = &chf_fn<extreme_value_l_handle, long double>;
        out->quantile = &quantile_fn<extreme_value_l_handle, long double>;
        out->quantile_complement = &quantile_complement_fn<extreme_value_l_handle, long double>;
        out->range = &range_fn<extreme_value_l_handle, bs_range_l, long double>;
        out->mean = &mean_fn<extreme_value_l_handle, long double>;
        out->variance = &variance_fn<extreme_value_l_handle, long double>;
        out->skewness = &skewness_fn<extreme_value_l_handle, long double>;
        out->kurtosis = &kurtosis_fn<extreme_value_l_handle, long double>;
        out->kurtosis_excess = &kurtosis_excess_fn<extreme_value_l_handle, long double>;
        out->mode = &mode_fn<extreme_value_l_handle, long double>;
        out->median = &median_fn<extreme_value_l_handle, long double>;
        out->entropy = nullptr;
        out->free = &free_fn<extreme_value_l_handle>;
        return true;
    } else if (n == "geometric" || n == "geometric_distribution") {
        long double p = 0;
        const char* pKeys[] = { "p", "prob", "probability", "success", "theta" };
        if (!find_param(params, count, pKeys, 5, &p)) return false;
        auto* h = new (std::nothrow) geometric_l_handle(p);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<geometric_l_handle, long double>;
        out->logpdf = &logpdf_fn<geometric_l_handle, long double>;
        out->cdf = &cdf_fn<geometric_l_handle, long double>;
        out->sf = &sf_fn<geometric_l_handle, long double>;
        out->hazard = &hazard_fn<geometric_l_handle, long double>;
        out->chf = &chf_fn<geometric_l_handle, long double>;
        out->quantile = &quantile_fn<geometric_l_handle, long double>;
        out->quantile_complement = &quantile_complement_fn<geometric_l_handle, long double>;
        out->range = &range_fn<geometric_l_handle, bs_range_l, long double>;
        out->mean = &mean_fn<geometric_l_handle, long double>;
        out->variance = &variance_fn<geometric_l_handle, long double>;
        out->skewness = &skewness_fn<geometric_l_handle, long double>;
        out->kurtosis = &kurtosis_fn<geometric_l_handle, long double>;
        out->kurtosis_excess = &kurtosis_excess_fn<geometric_l_handle, long double>;
        out->mode = &mode_fn<geometric_l_handle, long double>;
        out->median = &median_fn<geometric_l_handle, long double>;
        out->entropy = nullptr;
        out->free = &free_fn<geometric_l_handle>;
        return true;
    } else if (n == "holtsmark" || n == "holtsmark_distribution") {
        long double loc = 0, scale = 0;
        const char* locKeys[] = { "location", "loc", "mu", "median", "x0" };
        const char* scaleKeys[] = { "scale", "gamma", "sigma", "b" };
        find_param(params, count, locKeys, 5, &loc);
        if (!find_param(params, count, scaleKeys, 4, &scale)) return false;
        auto* h = new (std::nothrow) holtsmark_l_handle(loc, scale);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<holtsmark_l_handle, long double>;
        out->logpdf = &logpdf_fn<holtsmark_l_handle, long double>;
        out->cdf = &cdf_fn<holtsmark_l_handle, long double>;
        out->sf = &sf_fn<holtsmark_l_handle, long double>;
        out->hazard = &hazard_fn<holtsmark_l_handle, long double>;
        out->chf = &chf_fn<holtsmark_l_handle, long double>;
        out->quantile = &quantile_fn<holtsmark_l_handle, long double>;
        out->quantile_complement = &quantile_complement_fn<holtsmark_l_handle, long double>;
        out->range = &range_fn<holtsmark_l_handle, bs_range_l, long double>;
        out->mean = &mean_fn<holtsmark_l_handle, long double>;
        out->variance = &variance_fn<holtsmark_l_handle, long double>;
        out->skewness = nullptr;
        out->kurtosis = nullptr;
        out->kurtosis_excess = nullptr;
        out->mode = &mode_fn<holtsmark_l_handle, long double>;
        out->median = &median_fn<holtsmark_l_handle, long double>;
        out->entropy = nullptr;
        out->free = &free_fn<holtsmark_l_handle>;
        return true;
    }
        return false;
    } catch (...) {
        return false;
    }
}

} // extern "C"

