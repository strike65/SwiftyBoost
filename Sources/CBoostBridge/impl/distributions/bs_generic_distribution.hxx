//
//  Generic distribution vtable (Boost.Math) — implementation TU include
//
#include "../internal/bs_internal.hpp"
#include "../include/distributions/bs_generic_distribution.h"

#include <string>
#include <algorithm>
#include <cctype>
#include <vector>
#include <initializer_list>
#include <limits>

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
#include <boost/math/distributions/non_central_chi_squared.hpp>
#include <boost/math/distributions/non_central_f.hpp>
#include <boost/math/distributions/non_central_t.hpp>
#include <boost/math/distributions/pareto.hpp>
#include <boost/math/distributions/poisson.hpp>
#include <boost/math/distributions/non_central_f.hpp>
#include <boost/math/distributions/non_central_t.hpp>
#include <boost/math/distributions/inverse_gamma.hpp>
#include <boost/math/distributions/inverse_gaussian.hpp>
#include <boost/math/distributions/kolmogorov_smirnov.hpp>
#include <boost/math/distributions/landau.hpp>
#include <boost/math/distributions/laplace.hpp>
#include <boost/math/distributions/logistic.hpp>
#include <boost/math/distributions/lognormal.hpp>
#include <boost/math/distributions/mapairy.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/negative_binomial.hpp>
#include <boost/math/distributions/non_central_beta.hpp>
#include <boost/math/distributions/rayleigh.hpp>
#include <boost/math/distributions/saspoint5.hpp>
#include <boost/math/distributions/skew_normal.hpp>
#include <boost/math/distributions/triangular.hpp>
#include <boost/math/distributions/uniform.hpp>
#include <boost/math/distributions/weibull.hpp>

using boost::math::non_central_beta_distribution;
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
using boost::math::inverse_gamma_distribution;
using boost::math::inverse_gaussian_distribution;
using boost::math::kolmogorov_smirnov_distribution;
using boost::math::landau_distribution;
using boost::math::laplace_distribution;
using boost::math::logistic_distribution;
using boost::math::lognormal_distribution;
using boost::math::mapairy_distribution;
using boost::math::normal_distribution;
using boost::math::negative_binomial_distribution;
using boost::math::non_central_chi_squared_distribution;
using boost::math::non_central_f_distribution;
using boost::math::non_central_t_distribution;
using boost::math::pareto_distribution;
using boost::math::poisson_distribution;
using boost::math::non_central_f_distribution;
using boost::math::non_central_t_distribution;
using boost::math::rayleigh_distribution;
using boost::math::saspoint5_distribution;
using boost::math::skew_normal_distribution;
using boost::math::triangular_distribution;
using boost::math::uniform_distribution;
using boost::math::weibull_distribution;

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

template <class Value>
static bool build_indexed_vector(const std::vector<std::pair<size_t, Value>>& entries, std::vector<Value>& out) {
    if (entries.empty()) {
        out.clear();
        return true;
    }
    size_t maxIndex = 0;
    for (const auto& item : entries) {
        if (item.first > maxIndex) {
            maxIndex = item.first;
        }
    }
    const size_t size = maxIndex + 1;
    out.assign(size, Value(0));
    std::vector<bool> seen(size, false);
    for (const auto& item : entries) {
        const size_t idx = item.first;
        if (idx >= size) {
            return false;
        }
        if (seen[idx]) {
            return false;
        }
        out[idx] = item.second;
        seen[idx] = true;
    }
    for (bool flag : seen) {
        if (!flag) {
            return false;
        }
    }
    return true;
}

inline std::string trim_trailing_delimiters(std::string prefix) {
    while (!prefix.empty() && (prefix.back() == '_' || prefix.back() == '-' || prefix.back() == ' ')) {
        prefix.pop_back();
    }
    return prefix;
}

template <class Param, class Value>
static bool collect_hyperex_parameters(const Param* params, size_t count, std::vector<Value>& probs, std::vector<Value>& rates) {
    std::vector<std::pair<size_t, Value>> probEntries;
    std::vector<std::pair<size_t, Value>> rateEntries;
    size_t nextProbIndex = 0;
    size_t nextRateIndex = 0;

    for (size_t i = 0; i < count; ++i) {
        std::string key = lower_ascii(params[i].key);
        Value value = static_cast<Value>(params[i].value);

        bool hadDigits = false;
        size_t pos = key.size();
        while (pos > 0 && std::isdigit(static_cast<unsigned char>(key[pos - 1]))) {
            hadDigits = true;
            --pos;
        }
        std::string prefix = key.substr(0, pos);
        std::string digits = key.substr(pos);
        prefix = trim_trailing_delimiters(prefix);

        auto parseIndex = [&](const std::string& s) -> size_t {
            size_t result = 0;
            for (char c : s) {
                if (!std::isdigit(static_cast<unsigned char>(c))) {
                    return std::numeric_limits<size_t>::max();
                }
                result = result * 10 + static_cast<size_t>(c - '0');
            }
            return result;
        };

        auto matchAny = [&](const std::string& candidate, std::initializer_list<const char*> variants) -> bool {
            for (const char* option : variants) {
                if (candidate == lower_ascii(option)) {
                    return true;
                }
            }
            return false;
        };

        if (matchAny(prefix, { "rate", "rates", "lambda", "lam", "lambdaphase" })) {
            size_t index = hadDigits ? parseIndex(digits) : nextRateIndex++;
            if (index == std::numeric_limits<size_t>::max()) {
                return false;
            }
            rateEntries.emplace_back(index, value);
        } else if (matchAny(prefix, { "prob", "probability", "p", "weight", "w" })) {
            size_t index = hadDigits ? parseIndex(digits) : nextProbIndex++;
            if (index == std::numeric_limits<size_t>::max()) {
                return false;
            }
            probEntries.emplace_back(index, value);
        }
    }

    if (rateEntries.empty()) {
        return false;
    }
    if (!build_indexed_vector(rateEntries, rates)) {
        return false;
    }

    if (!probEntries.empty()) {
        if (!build_indexed_vector(probEntries, probs)) {
            return false;
        }
        if (probs.size() != rates.size()) {
            return false;
        }
    } else {
        probs.clear();
    }
    return true;
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
static T median_via_quantile_fn(const void* p) {
    auto h = static_cast<const Handle*>(p);
    if (!h) return std::numeric_limits<T>::quiet_NaN();
    return bs_wrap<T>([&]{ return quantile(h->dist, static_cast<T>(0.5)); });
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

// noncentral Beta
struct non_central_beta_f_handle { non_central_beta_distribution<float> dist; non_central_beta_f_handle(float a, float b, float lambda): dist(a, b, lambda) {} };
struct non_central_beta_d_handle { non_central_beta_distribution<double> dist; non_central_beta_d_handle(double a, double b, double lambda): dist(a, b, lambda) {} };
struct non_central_beta_l_handle { non_central_beta_distribution<long double> dist; non_central_beta_l_handle(long double a, long double b, long double lambda): dist(a, b,lambda) {} };

// Chi-Squared
struct chi_squared_f_handle { chi_squared_distribution<float> dist; explicit chi_squared_f_handle(float v): dist(v) {} };
struct chi_squared_d_handle { chi_squared_distribution<double> dist; explicit chi_squared_d_handle(double v): dist(v) {} };
struct chi_squared_l_handle { chi_squared_distribution<long double> dist; explicit chi_squared_l_handle(long double v): dist(v) {} };

// noncentral Chi-Squared
struct non_central_chi_squared_f_handle { non_central_chi_squared_distribution<float> dist; explicit non_central_chi_squared_f_handle(float v, float l): dist(v,l) {} };
struct non_central_chi_squared_d_handle { non_central_chi_squared_distribution<double> dist; explicit non_central_chi_squared_d_handle(double v,double l): dist(v,l) {} };
struct non_central_chi_squared_l_handle { non_central_chi_squared_distribution<long double> dist; explicit non_central_chi_squared_l_handle(long double v, long double l): dist(v, l) {} };
// noncentral F
struct non_central_f_f_handle { non_central_f_distribution<float> dist; non_central_f_f_handle(float df1, float df2, float lambda): dist(df1, df2, lambda) {} };
struct non_central_f_d_handle { non_central_f_distribution<double> dist; non_central_f_d_handle(double df1, double df2, double lambda): dist(df1, df2, lambda) {} };
struct non_central_f_l_handle { non_central_f_distribution<long double> dist; non_central_f_l_handle(long double df1, long double df2, long double lambda): dist(df1, df2, lambda) {} };
// noncentral Student T
struct non_central_t_f_handle { non_central_t_distribution<float> dist; non_central_t_f_handle(float v, float lambda): dist(v, lambda) {} };
struct non_central_t_d_handle { non_central_t_distribution<double> dist; non_central_t_d_handle(double v, double lambda): dist(v, lambda) {} };
struct non_central_t_l_handle { non_central_t_distribution<long double> dist; non_central_t_l_handle(long double v, long double lambda): dist(v, lambda) {} };
// Pareto
struct pareto_f_handle { pareto_distribution<float> dist; pareto_f_handle(float xm, float alpha): dist(xm, alpha) {} };
struct pareto_d_handle { pareto_distribution<double> dist; pareto_d_handle(double xm, double alpha): dist(xm, alpha) {} };
struct pareto_l_handle { pareto_distribution<long double> dist; pareto_l_handle(long double xm, long double alpha): dist(xm, alpha) {} };
// Poisson
struct poisson_f_handle { poisson_distribution<float> dist; explicit poisson_f_handle(float mean): dist(mean) {} };
struct poisson_d_handle { poisson_distribution<double> dist; explicit poisson_d_handle(double mean): dist(mean) {} };
struct poisson_l_handle { poisson_distribution<long double> dist; explicit poisson_l_handle(long double mean): dist(mean) {} };
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

// Gaussian
struct normal_f_handle { normal_distribution<float> dist; normal_f_handle(float loc, float scale): dist(loc, scale) {} };
struct normal_d_handle { normal_distribution<double> dist; normal_d_handle(double loc, double scale): dist(loc, scale) {} };
struct normal_l_handle { normal_distribution<long double> dist; normal_l_handle(long double loc, long double scale): dist(loc, scale) {} };

// Kolmogorov–Smirnov
struct kolmogorov_smirnov_f_handle { kolmogorov_smirnov_distribution<float> dist; explicit kolmogorov_smirnov_f_handle(float n): dist(n) {} };
struct kolmogorov_smirnov_d_handle { kolmogorov_smirnov_distribution<double> dist; explicit kolmogorov_smirnov_d_handle(double n): dist(n) {} };
struct kolmogorov_smirnov_l_handle { kolmogorov_smirnov_distribution<long double> dist; explicit kolmogorov_smirnov_l_handle(long double n): dist(n) {} };

// Landau
struct landau_f_handle { landau_distribution<float> dist; landau_f_handle(float loc, float scale): dist(loc, scale) {} };
struct landau_d_handle { landau_distribution<double> dist; landau_d_handle(double loc, double scale): dist(loc, scale) {} };
struct landau_l_handle { landau_distribution<long double> dist; landau_l_handle(long double loc, long double scale): dist(loc, scale) {} };

// Laplace (double exponential)
struct laplace_f_handle { laplace_distribution<float> dist; laplace_f_handle(float loc, float scale): dist(loc, scale) {} };
struct laplace_d_handle { laplace_distribution<double> dist; laplace_d_handle(double loc, double scale): dist(loc, scale) {} };
struct laplace_l_handle { laplace_distribution<long double> dist; laplace_l_handle(long double loc, long double scale): dist(loc, scale) {} };

// Hyperexponential (mixture of exponentials)
struct hyperexponential_f_handle {
    hyperexponential_distribution<float> dist;
    hyperexponential_f_handle(const std::vector<float>& probs, const std::vector<float>& rates)
        : dist(probs.begin(), probs.end(), rates.begin(), rates.end()) {}
    explicit hyperexponential_f_handle(const std::vector<float>& rates)
        : dist(rates.begin(), rates.end()) {}
};
struct hyperexponential_d_handle {
    hyperexponential_distribution<double> dist;
    hyperexponential_d_handle(const std::vector<double>& probs, const std::vector<double>& rates)
        : dist(probs.begin(), probs.end(), rates.begin(), rates.end()) {}
    explicit hyperexponential_d_handle(const std::vector<double>& rates)
        : dist(rates.begin(), rates.end()) {}
};
struct hyperexponential_l_handle {
    hyperexponential_distribution<long double> dist;
    hyperexponential_l_handle(const std::vector<long double>& probs, const std::vector<long double>& rates)
        : dist(probs.begin(), probs.end(), rates.begin(), rates.end()) {}
    explicit hyperexponential_l_handle(const std::vector<long double>& rates)
        : dist(rates.begin(), rates.end()) {}
};

// Logistic
struct logistic_f_handle { logistic_distribution<float> dist; logistic_f_handle(float loc, float scale): dist(loc, scale) {} };
struct logistic_d_handle { logistic_distribution<double> dist; logistic_d_handle(double loc, double scale): dist(loc, scale) {} };
struct logistic_l_handle { logistic_distribution<long double> dist; logistic_l_handle(long double loc, long double scale): dist(loc, scale) {} };

// Lognormal
struct lognormal_f_handle { lognormal_distribution<float> dist; lognormal_f_handle(float mu, float sigma): dist(mu, sigma) {} };
struct lognormal_d_handle { lognormal_distribution<double> dist; lognormal_d_handle(double mu, double sigma): dist(mu, sigma) {} };
struct lognormal_l_handle { lognormal_distribution<long double> dist; lognormal_l_handle(long double mu, long double sigma): dist(mu, sigma) {} };

// Map-Airy
struct mapairy_f_handle { mapairy_distribution<float> dist; mapairy_f_handle(float loc, float scale): dist(loc, scale) {} };
struct mapairy_d_handle { mapairy_distribution<double> dist; mapairy_d_handle(double loc, double scale): dist(loc, scale) {} };
struct mapairy_l_handle { mapairy_distribution<long double> dist; mapairy_l_handle(long double loc, long double scale): dist(loc, scale) {} };

// Rayleigh
struct rayleigh_f_handle { rayleigh_distribution<float> dist; explicit rayleigh_f_handle(float scale): dist(scale) {} };
struct rayleigh_d_handle { rayleigh_distribution<double> dist; explicit rayleigh_d_handle(double scale): dist(scale) {} };
struct rayleigh_l_handle { rayleigh_distribution<long double> dist; explicit rayleigh_l_handle(long double scale): dist(scale) {} };

// SaS point 5 (symmetric alpha-stable α = 0.5)
struct saspoint5_f_handle { saspoint5_distribution<float> dist; saspoint5_f_handle(float location, float scale): dist(location, scale) {} };
struct saspoint5_d_handle { saspoint5_distribution<double> dist; saspoint5_d_handle(double location, double scale): dist(location, scale) {} };
struct saspoint5_l_handle { saspoint5_distribution<long double> dist; saspoint5_l_handle(long double location, long double scale): dist(location, scale) {} };

// Skew normal
struct skew_normal_f_handle { skew_normal_distribution<float> dist; skew_normal_f_handle(float location, float scale, float shape): dist(location, scale, shape) {} };
struct skew_normal_d_handle { skew_normal_distribution<double> dist; skew_normal_d_handle(double location, double scale, double shape): dist(location, scale, shape) {} };
struct skew_normal_l_handle { skew_normal_distribution<long double> dist; skew_normal_l_handle(long double location, long double scale, long double shape): dist(location, scale, shape) {} };

// Triangular
struct triangular_f_handle { triangular_distribution<float> dist; triangular_f_handle(float lower, float mode, float upper): dist(lower, mode, upper) {} };
struct triangular_d_handle { triangular_distribution<double> dist; triangular_d_handle(double lower, double mode, double upper): dist(lower, mode, upper) {} };
struct triangular_l_handle { triangular_distribution<long double> dist; triangular_l_handle(long double lower, long double mode, long double upper): dist(lower, mode, upper) {} };

// Uniform
struct uniform_f_handle { uniform_distribution<float> dist; uniform_f_handle(float lower, float upper): dist(lower, upper) {} };
struct uniform_d_handle { uniform_distribution<double> dist; uniform_d_handle(double lower, double upper): dist(lower, upper) {} };
struct uniform_l_handle { uniform_distribution<long double> dist; uniform_l_handle(long double lower, long double upper): dist(lower, upper) {} };

// Weibull
struct weibull_f_handle { weibull_distribution<float> dist; weibull_f_handle(float shape, float scale): dist(shape, scale) {} };
struct weibull_d_handle { weibull_distribution<double> dist; weibull_d_handle(double shape, double scale): dist(shape, scale) {} };
struct weibull_l_handle { weibull_distribution<long double> dist; weibull_l_handle(long double shape, long double scale): dist(shape, scale) {} };


// inverse Chi-Squared
struct inverse_chi_squared_f_handle { inverse_chi_squared_distribution<float> dist; explicit inverse_chi_squared_f_handle(float v, float s): dist(v,s) {} };
struct inverse_chi_squared_d_handle { inverse_chi_squared_distribution<double> dist; explicit inverse_chi_squared_d_handle(double v, double s): dist(v,s) {} };
struct inverse_chi_squared_l_handle { inverse_chi_squared_distribution<long double> dist; explicit inverse_chi_squared_l_handle(long double v, long double s): dist(v,s) {} };

// inverse Gamma
struct inverse_gamma_f_handle { inverse_gamma_distribution<float> dist; inverse_gamma_f_handle(float shape, float scale): dist(shape, scale) {} };
struct inverse_gamma_d_handle { inverse_gamma_distribution<double> dist; inverse_gamma_d_handle(double shape, double scale): dist(shape, scale) {} };
struct inverse_gamma_l_handle { inverse_gamma_distribution<long double> dist; inverse_gamma_l_handle(long double shape, long double scale): dist(shape, scale) {} };

// inverse Gaussian (a.k.a. inverse Normal, Wald)
struct inverse_gaussian_f_handle { inverse_gaussian_distribution<float> dist; inverse_gaussian_f_handle(float mean, float scale): dist(mean, scale) {} };
struct inverse_gaussian_d_handle { inverse_gaussian_distribution<double> dist; inverse_gaussian_d_handle(double mean, double scale): dist(mean, scale) {} };
struct inverse_gaussian_l_handle { inverse_gaussian_distribution<long double> dist; inverse_gaussian_l_handle(long double mean, long double scale): dist(mean, scale) {} };

// negative Binomial
struct negative_binomial_f_handle { negative_binomial_distribution<float> dist; negative_binomial_f_handle(float r, float p): dist(r, p) {} };
struct negative_binomial_d_handle { negative_binomial_distribution<double> dist; negative_binomial_d_handle(double r, double p): dist(r, p) {} };
struct negative_binomial_l_handle { negative_binomial_distribution<long double> dist; negative_binomial_l_handle(long double r, long double p): dist(r, p) {} };


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
    } else if (n == "non_central_beta" || n == "non_central_beta_distribution") {
        double a = 1, b = 1, l = 0;
        const char* aKeys[] = { "alpha", "a", "p", "shape1" };
        const char* bKeys[] = { "beta", "b", "q", "shape2" };
        const char* lKeys[] = {"lambda", "l"};
        if (!find_param(params, count, aKeys, 4, &a)) return false;
        if (!find_param(params, count, bKeys, 4, &b)) return false;
        if (!find_param(params, count, lKeys, 2, &l)) return false;
        auto* h = new (std::nothrow) non_central_beta_d_handle(a, b, l);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<non_central_beta_d_handle, double>;
        out->logpdf = &logpdf_fn<non_central_beta_d_handle, double>;
        out->cdf = &cdf_fn<non_central_beta_d_handle, double>;
        out->sf = &sf_fn<non_central_beta_d_handle, double>;
        out->hazard = &hazard_fn<non_central_beta_d_handle, double>;
        out->chf = &chf_fn<non_central_beta_d_handle, double>;
        out->quantile = &quantile_fn<non_central_beta_d_handle, double>;
        out->quantile_complement = &quantile_complement_fn<non_central_beta_d_handle, double>;
        out->range = &range_fn<non_central_beta_d_handle, bs_range_d, double>;
        out->mean = &mean_fn<non_central_beta_d_handle, double>;
        out->variance = &variance_fn<non_central_beta_d_handle, double>;
        out->skewness = nullptr;
        out->kurtosis = nullptr; //&kurtosis_fn<non_central_beta_d_handle, double>;
        out->kurtosis_excess = nullptr;
        out->mode = &mode_fn<non_central_beta_d_handle, double>;
        out->median = &median_fn<non_central_beta_d_handle, double>;
        out->entropy = nullptr; // not provided in Boost
        out->free = &free_fn<non_central_beta_d_handle>;
        return true;
    } else if (n == "non_central_chi_squared" || n == "noncentral_chi_squared" || n == "noncentralchisquared" || n == "non_central_chi2" || n == "noncentral_chi2" || n == "nc_chi_squared" || n == "ncchisquared") {
        double v = 0;
        double lambda = 0;
        const char* dfKeys[] = { "df", "nu", "degreesoffreedom" };
        const char* lambdaKeys[] = { "lambda", "noncentrality", "delta", "nc" };
        if (!find_param(params, count, dfKeys, 3, &v)) return false;
        if (!find_param(params, count, lambdaKeys, 4, &lambda)) return false;
        if (!(v > 0)) return false;
        if (!(lambda >= 0)) return false;
        auto* h = new (std::nothrow) non_central_chi_squared_d_handle(v, lambda);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<non_central_chi_squared_d_handle, double>;
        out->logpdf = &logpdf_fn<non_central_chi_squared_d_handle, double>;
        out->cdf = &cdf_fn<non_central_chi_squared_d_handle, double>;
        out->sf = &sf_fn<non_central_chi_squared_d_handle, double>;
        out->hazard = &hazard_fn<non_central_chi_squared_d_handle, double>;
        out->chf = &chf_fn<non_central_chi_squared_d_handle, double>;
        out->quantile = &quantile_fn<non_central_chi_squared_d_handle, double>;
        out->quantile_complement = &quantile_complement_fn<non_central_chi_squared_d_handle, double>;
        out->range = &range_fn<non_central_chi_squared_d_handle, bs_range_d, double>;
        out->mean = &mean_fn<non_central_chi_squared_d_handle, double>;
        out->variance = &variance_fn<non_central_chi_squared_d_handle, double>;
        out->skewness = &skewness_fn<non_central_chi_squared_d_handle, double>;
        out->kurtosis = &kurtosis_fn<non_central_chi_squared_d_handle, double>;
        out->kurtosis_excess = &kurtosis_excess_fn<non_central_chi_squared_d_handle, double>;
        out->mode = &mode_fn<non_central_chi_squared_d_handle, double>;
        out->median = &median_fn<non_central_chi_squared_d_handle, double>;
        out->entropy = nullptr;
        out->free = &free_fn<non_central_chi_squared_d_handle>;
        return true;
    } else if (n == "non_central_f" || n == "noncentral_f" || n == "noncentralf" || n == "non_central_f_ratio" || n == "noncentral_f_ratio" || n == "nc_f" || n == "ncf") {
        double df1 = 0;
        double df2 = 0;
        double lambda = 0;
        const char* df1Keys[] = { "df1", "d1", "m", "degreesoffreedom1" };
        const char* df2Keys[] = { "df2", "d2", "n", "degreesoffreedom2" };
        const char* lambdaKeys[] = { "lambda", "noncentrality", "delta", "nc" };
        if (!find_param(params, count, df1Keys, 4, &df1)) return false;
        if (!find_param(params, count, df2Keys, 4, &df2)) return false;
        if (!find_param(params, count, lambdaKeys, 4, &lambda)) return false;
        if (!(df1 > 0) || !(df2 > 0) || !(lambda >= 0)) return false;
        auto* h = new (std::nothrow) non_central_f_d_handle(df1, df2, lambda);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<non_central_f_d_handle, double>;
        out->logpdf = &logpdf_fn<non_central_f_d_handle, double>;
        out->cdf = &cdf_fn<non_central_f_d_handle, double>;
        out->sf = &sf_fn<non_central_f_d_handle, double>;
        out->hazard = &hazard_fn<non_central_f_d_handle, double>;
        out->chf = &chf_fn<non_central_f_d_handle, double>;
        out->quantile = &quantile_fn<non_central_f_d_handle, double>;
        out->quantile_complement = &quantile_complement_fn<non_central_f_d_handle, double>;
        out->range = &range_fn<non_central_f_d_handle, bs_range_d, double>;
        out->mean = &mean_fn<non_central_f_d_handle, double>;
        out->variance = &variance_fn<non_central_f_d_handle, double>;
        out->skewness = &skewness_fn<non_central_f_d_handle, double>;
        out->kurtosis = &kurtosis_fn<non_central_f_d_handle, double>;
        out->kurtosis_excess = &kurtosis_excess_fn<non_central_f_d_handle, double>;
        out->mode = &mode_fn<non_central_f_d_handle, double>;
        out->median = &median_fn<non_central_f_d_handle, double>;
        out->entropy = nullptr;
        out->free = &free_fn<non_central_f_d_handle>;
        return true;
    } else if (n == "non_central_t" || n == "noncentral_t" || n == "noncentralt" || n == "non_central_student_t" || n == "noncentral_student_t" || n == "noncentralstudentt" || n == "nc_t") {
        double v = 0;
        double lambda = 0;
        const char* dfKeys[] = { "df", "nu", "degreesoffreedom" };
        const char* lambdaKeys[] = { "lambda", "noncentrality", "delta", "nc" };
        if (!find_param(params, count, dfKeys, 3, &v)) return false;
        if (!find_param(params, count, lambdaKeys, 4, &lambda)) return false;
        if (!(v > 0)) return false;
        if (!std::isfinite(lambda)) return false;
        auto* h = new (std::nothrow) non_central_t_d_handle(v, lambda);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<non_central_t_d_handle, double>;
        out->logpdf = &logpdf_fn<non_central_t_d_handle, double>;
        out->cdf = &cdf_fn<non_central_t_d_handle, double>;
        out->sf = &sf_fn<non_central_t_d_handle, double>;
        out->hazard = &hazard_fn<non_central_t_d_handle, double>;
        out->chf = &chf_fn<non_central_t_d_handle, double>;
        out->quantile = &quantile_fn<non_central_t_d_handle, double>;
        out->quantile_complement = &quantile_complement_fn<non_central_t_d_handle, double>;
        out->range = &range_fn<non_central_t_d_handle, bs_range_d, double>;
        out->mean = &mean_fn<non_central_t_d_handle, double>;
        out->variance = &variance_fn<non_central_t_d_handle, double>;
        out->skewness = &skewness_fn<non_central_t_d_handle, double>;
        out->kurtosis = &kurtosis_fn<non_central_t_d_handle, double>;
        out->kurtosis_excess = &kurtosis_excess_fn<non_central_t_d_handle, double>;
        out->mode = &mode_fn<non_central_t_d_handle, double>;
        out->median = &median_fn<non_central_t_d_handle, double>;
        out->entropy = nullptr;
        out->free = &free_fn<non_central_t_d_handle>;
        return true;
    } else if (n == "pareto" || n == "pareto_distribution") {
        double scale = 0;
        double shape = 0;
        const char* scaleKeys[] = { "scale", "xm", "minimum", "lower", "x0" };
        const char* shapeKeys[] = { "shape", "alpha" };
        if (!find_param(params, count, scaleKeys, 5, &scale)) return false;
        if (!find_param(params, count, shapeKeys, 2, &shape)) return false;
        if (!(scale > 0) || !(shape > 0)) return false;
        auto* h = new (std::nothrow) pareto_d_handle(scale, shape);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<pareto_d_handle, double>;
        out->logpdf = &logpdf_fn<pareto_d_handle, double>;
        out->cdf = &cdf_fn<pareto_d_handle, double>;
        out->sf = &sf_fn<pareto_d_handle, double>;
        out->hazard = &hazard_fn<pareto_d_handle, double>;
        out->chf = &chf_fn<pareto_d_handle, double>;
        out->quantile = &quantile_fn<pareto_d_handle, double>;
        out->quantile_complement = &quantile_complement_fn<pareto_d_handle, double>;
        out->range = &range_fn<pareto_d_handle, bs_range_d, double>;
        out->mean = &mean_fn<pareto_d_handle, double>;
        out->variance = &variance_fn<pareto_d_handle, double>;
        out->skewness = &skewness_fn<pareto_d_handle, double>;
        out->kurtosis = &kurtosis_fn<pareto_d_handle, double>;
        out->kurtosis_excess = &kurtosis_excess_fn<pareto_d_handle, double>;
        out->mode = &mode_fn<pareto_d_handle, double>;
        out->median = &median_fn<pareto_d_handle, double>;
        out->entropy = nullptr;
        out->free = &free_fn<pareto_d_handle>;
        return true;
    } else if (n == "poisson" || n == "poisson_distribution") {
        double mean = 0;
        const char* meanKeys[] = { "mean", "lambda", "mu" };
        if (!find_param(params, count, meanKeys, 3, &mean)) return false;
        if (!(mean >= 0)) return false;
        auto* h = new (std::nothrow) poisson_d_handle(mean);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<poisson_d_handle, double>;
        out->logpdf = &logpdf_fn<poisson_d_handle, double>;
        out->cdf = &cdf_fn<poisson_d_handle, double>;
        out->sf = &sf_fn<poisson_d_handle, double>;
        out->hazard = &hazard_fn<poisson_d_handle, double>;
        out->chf = &chf_fn<poisson_d_handle, double>;
        out->quantile = &quantile_fn<poisson_d_handle, double>;
        out->quantile_complement = &quantile_complement_fn<poisson_d_handle, double>;
        out->range = &range_fn<poisson_d_handle, bs_range_d, double>;
        out->mean = &mean_fn<poisson_d_handle, double>;
        out->variance = &variance_fn<poisson_d_handle, double>;
        out->skewness = &skewness_fn<poisson_d_handle, double>;
        out->kurtosis = &kurtosis_fn<poisson_d_handle, double>;
        out->kurtosis_excess = &kurtosis_excess_fn<poisson_d_handle, double>;
        out->mode = &mode_fn<poisson_d_handle, double>;
        out->median = &median_fn<poisson_d_handle, double>;
        out->entropy = nullptr;
        out->free = &free_fn<poisson_d_handle>;
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
    } else if (n == "inverse_chi_squared" || n == "inversechisquared" || n == "inv_chi_squared" || n == "invchisquared" || n == "inverse_chi2" || n == "inv_chi2") {
        double df = 0;
        const char* dfKeys[] = { "df", "nu", "degreesoffreedom", "v" };
        if (!find_param(params, count, dfKeys, 4, &df)) return false;
        double scale = 0;
        const char* scaleKeys[] = { "scale", "sigma2", "xi" };
        bool hasScale = find_param(params, count, scaleKeys, 3, &scale);
        if (!hasScale) {
            if (!(df > 0)) return false;
            scale = 1.0 / df;
        }
        auto* h = new (std::nothrow) inverse_chi_squared_d_handle(df, scale);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<inverse_chi_squared_d_handle, double>;
        out->logpdf = &logpdf_fn<inverse_chi_squared_d_handle, double>;
        out->cdf = &cdf_fn<inverse_chi_squared_d_handle, double>;
        out->sf = &sf_fn<inverse_chi_squared_d_handle, double>;
        out->hazard = &hazard_fn<inverse_chi_squared_d_handle, double>;
        out->chf = &chf_fn<inverse_chi_squared_d_handle, double>;
        out->quantile = &quantile_fn<inverse_chi_squared_d_handle, double>;
        out->quantile_complement = &quantile_complement_fn<inverse_chi_squared_d_handle, double>;
        out->range = &range_fn<inverse_chi_squared_d_handle, bs_range_d, double>;
        out->mean = &mean_fn<inverse_chi_squared_d_handle, double>;
        out->variance = &variance_fn<inverse_chi_squared_d_handle, double>;
        out->skewness = &skewness_fn<inverse_chi_squared_d_handle, double>;
        out->kurtosis = &kurtosis_fn<inverse_chi_squared_d_handle, double>;
        out->kurtosis_excess = &kurtosis_excess_fn<inverse_chi_squared_d_handle, double>;
        out->mode = &mode_fn<inverse_chi_squared_d_handle, double>;
        out->median = &median_fn<inverse_chi_squared_d_handle, double>;
        out->entropy = nullptr;
        out->free = &free_fn<inverse_chi_squared_d_handle>;
        return true;
    } else if (n == "inverse_gamma" || n == "inversegamma" || n == "inv_gamma" || n == "invgamma") {
        double shape = 0;
        const char* shapeKeys[] = { "shape", "alpha", "k" };
        if (!find_param(params, count, shapeKeys, 3, &shape)) return false;
        double scale = 1;
        const char* scaleKeys[] = { "scale", "theta", "beta" };
        find_param(params, count, scaleKeys, 3, &scale);
        auto* h = new (std::nothrow) inverse_gamma_d_handle(shape, scale);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<inverse_gamma_d_handle, double>;
        out->logpdf = &logpdf_fn<inverse_gamma_d_handle, double>;
        out->cdf = &cdf_fn<inverse_gamma_d_handle, double>;
        out->sf = &sf_fn<inverse_gamma_d_handle, double>;
        out->hazard = &hazard_fn<inverse_gamma_d_handle, double>;
        out->chf = &chf_fn<inverse_gamma_d_handle, double>;
        out->quantile = &quantile_fn<inverse_gamma_d_handle, double>;
        out->quantile_complement = &quantile_complement_fn<inverse_gamma_d_handle, double>;
        out->range = &range_fn<inverse_gamma_d_handle, bs_range_d, double>;
        out->mean = &mean_fn<inverse_gamma_d_handle, double>;
        out->variance = &variance_fn<inverse_gamma_d_handle, double>;
        out->skewness = &skewness_fn<inverse_gamma_d_handle, double>;
        out->kurtosis = &kurtosis_fn<inverse_gamma_d_handle, double>;
        out->kurtosis_excess = &kurtosis_excess_fn<inverse_gamma_d_handle, double>;
        out->mode = &mode_fn<inverse_gamma_d_handle, double>;
        out->median = &median_fn<inverse_gamma_d_handle, double>;
        out->entropy = nullptr;
        out->free = &free_fn<inverse_gamma_d_handle>;
        return true;
    } else if (n == "inverse_gaussian" || n == "inversegaussian" || n == "inverse_normal" || n == "inversenormal" || n == "wald") {
        double mean = 0;
        const char* meanKeys[] = { "mean", "mu", "location" };
        if (!find_param(params, count, meanKeys, 3, &mean)) return false;
        double scale = 1;
        const char* scaleKeys[] = { "scale", "lambda", "shape" };
        find_param(params, count, scaleKeys, 3, &scale);
        auto* h = new (std::nothrow) inverse_gaussian_d_handle(mean, scale);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<inverse_gaussian_d_handle, double>;
        out->logpdf = &logpdf_fn<inverse_gaussian_d_handle, double>;
        out->cdf = &cdf_fn<inverse_gaussian_d_handle, double>;
        out->sf = &sf_fn<inverse_gaussian_d_handle, double>;
        out->hazard = &hazard_fn<inverse_gaussian_d_handle, double>;
        out->chf = &chf_fn<inverse_gaussian_d_handle, double>;
        out->quantile = &quantile_fn<inverse_gaussian_d_handle, double>;
        out->quantile_complement = &quantile_complement_fn<inverse_gaussian_d_handle, double>;
        out->range = &range_fn<inverse_gaussian_d_handle, bs_range_d, double>;
        out->mean = &mean_fn<inverse_gaussian_d_handle, double>;
        out->variance = &variance_fn<inverse_gaussian_d_handle, double>;
        out->skewness = &skewness_fn<inverse_gaussian_d_handle, double>;
        out->kurtosis = &kurtosis_fn<inverse_gaussian_d_handle, double>;
        out->kurtosis_excess = &kurtosis_excess_fn<inverse_gaussian_d_handle, double>;
        out->mode = &mode_fn<inverse_gaussian_d_handle, double>;
        out->median = &median_fn<inverse_gaussian_d_handle, double>;
        out->entropy = nullptr;
        out->free = &free_fn<inverse_gaussian_d_handle>;
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
    } else if (n == "negative_binomial" || n == "negativebinomial" || n == "neg_binomial" || n == "negative_binomial_distribution" || n == "nbinom") {
        double successes = 0, p = 0;
        const char* rKeys[] = { "r", "successes", "target", "count" };
        const char* pKeys[] = { "p", "prob", "probability", "success" };
        if (!find_param(params, count, rKeys, 4, &successes)) return false;
        if (!find_param(params, count, pKeys, 4, &p)) return false;
        if (!(successes > 0)) return false;
        if (!(p > 0 && p <= 1)) return false;
        auto* h = new (std::nothrow) negative_binomial_d_handle(successes, p);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<negative_binomial_d_handle, double>;
        out->logpdf = &logpdf_fn<negative_binomial_d_handle, double>;
        out->cdf = &cdf_fn<negative_binomial_d_handle, double>;
        out->sf = &sf_fn<negative_binomial_d_handle, double>;
        out->hazard = &hazard_fn<negative_binomial_d_handle, double>;
        out->chf = &chf_fn<negative_binomial_d_handle, double>;
        out->quantile = &quantile_fn<negative_binomial_d_handle, double>;
        out->quantile_complement = &quantile_complement_fn<negative_binomial_d_handle, double>;
        out->range = &range_fn<negative_binomial_d_handle, bs_range_d, double>;
        out->mean = &mean_fn<negative_binomial_d_handle, double>;
        out->variance = &variance_fn<negative_binomial_d_handle, double>;
        out->skewness = &skewness_fn<negative_binomial_d_handle, double>;
        out->kurtosis = &kurtosis_fn<negative_binomial_d_handle, double>;
        out->kurtosis_excess = &kurtosis_excess_fn<negative_binomial_d_handle, double>;
        out->mode = &mode_fn<negative_binomial_d_handle, double>;
        out->median = &median_via_quantile_fn<negative_binomial_d_handle, double>;
        out->entropy = nullptr;
        out->free = &free_fn<negative_binomial_d_handle>;
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
        find_param(params, count, lamKeys, 2, &lambda);
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
        double loc = 0, scale = 1;
        const char* locKeys[] = { "location", "loc", "mu", "median", "x0" };
        const char* scaleKeys[] = { "scale", "gamma", "sigma", "b" };
        find_param(params, count, locKeys, 5, &loc);
        find_param(params, count, scaleKeys, 4, &scale);
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
    } else if (n == "hyperexponential" || n == "hyper_exponential" || n == "hyperexp" || n == "hyperexponential_distribution") {
        std::vector<double> probs;
        std::vector<double> rates;
        if (!collect_hyperex_parameters(params, count, probs, rates)) return false;
        hyperexponential_d_handle* h = nullptr;
        if (probs.empty()) {
            h = new (std::nothrow) hyperexponential_d_handle(rates);
        } else {
            h = new (std::nothrow) hyperexponential_d_handle(probs, rates);
        }
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<hyperexponential_d_handle, double>;
        out->logpdf = &logpdf_fn<hyperexponential_d_handle, double>;
        out->cdf = &cdf_fn<hyperexponential_d_handle, double>;
        out->sf = &sf_fn<hyperexponential_d_handle, double>;
        out->hazard = &hazard_fn<hyperexponential_d_handle, double>;
        out->chf = &chf_fn<hyperexponential_d_handle, double>;
        out->quantile = &quantile_fn<hyperexponential_d_handle, double>;
        out->quantile_complement = &quantile_complement_fn<hyperexponential_d_handle, double>;
        out->range = &range_fn<hyperexponential_d_handle, bs_range_d, double>;
        out->mean = &mean_fn<hyperexponential_d_handle, double>;
        out->variance = &variance_fn<hyperexponential_d_handle, double>;
        out->skewness = &skewness_fn<hyperexponential_d_handle, double>;
        out->kurtosis = &kurtosis_fn<hyperexponential_d_handle, double>;
        out->kurtosis_excess = &kurtosis_excess_fn<hyperexponential_d_handle, double>;
        out->mode = &mode_fn<hyperexponential_d_handle, double>;
        out->median = &median_fn<hyperexponential_d_handle, double>;
        out->entropy = nullptr;
        out->free = &free_fn<hyperexponential_d_handle>;
        return true;
    } else if (n == "normal" || n == "normal_distribution" || n == "gauss" || n == "gaussian" || n == "gaussian_distribution" || n == "gauss_distribution") {
        double loc = 0, scale = 1;
        const char* locKeys[] = { "location", "loc", "mu", "mean" };
        const char* scaleKeys[] = { "sd", "standard_deviation", "sigma" };
        find_param(params, count, locKeys, 4, &loc);
        find_param(params, count, scaleKeys, 3, &scale);
        auto* h = new (std::nothrow) normal_d_handle(loc, scale);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<normal_d_handle, double>;
        out->logpdf = &logpdf_fn<normal_d_handle, double>;
        out->cdf = &cdf_fn<normal_d_handle, double>;
        out->sf = &sf_fn<normal_d_handle, double>;
        out->hazard = &hazard_fn<normal_d_handle, double>;
        out->chf = &chf_fn<normal_d_handle, double>;
        out->quantile = &quantile_fn<normal_d_handle, double>;
        out->quantile_complement = &quantile_complement_fn<normal_d_handle, double>;
        out->range = &range_fn<normal_d_handle, bs_range_d, double>;
        out->mean = &mean_fn<normal_d_handle, double>;
        out->variance = &variance_fn<normal_d_handle, double>;
        out->skewness = &skewness_fn<normal_d_handle, double>;
        out->kurtosis = &kurtosis_fn<normal_d_handle, double>;
        out->kurtosis_excess = &kurtosis_excess_fn<normal_d_handle, double>;;
        out->mode = &mode_fn<normal_d_handle, double>;
        out->median = &median_fn<normal_d_handle, double>;
        out->entropy = &entropy_fn<normal_d_handle, double>;;
        out->free = &free_fn<normal_d_handle>;
        return true;
    } else if (n == "logistic" || n == "logistic_distribution") {
        double loc = 0, scale = 1;
        const char* locKeys[] = { "location", "loc", "mu", "median" };
        const char* scaleKeys[] = { "scale", "s", "sigma", "diversity" };
        find_param(params, count, locKeys, 4, &loc);
        find_param(params, count, scaleKeys, 4, &scale);
        if (!(scale > 0)) return false;
        auto* h = new (std::nothrow) logistic_d_handle(loc, scale);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<logistic_d_handle, double>;
        out->logpdf = &logpdf_fn<logistic_d_handle, double>;
        out->cdf = &cdf_fn<logistic_d_handle, double>;
        out->sf = &sf_fn<logistic_d_handle, double>;
        out->hazard = &hazard_fn<logistic_d_handle, double>;
        out->chf = &chf_fn<logistic_d_handle, double>;
        out->quantile = &quantile_fn<logistic_d_handle, double>;
        out->quantile_complement = &quantile_complement_fn<logistic_d_handle, double>;
        out->range = &range_fn<logistic_d_handle, bs_range_d, double>;
        out->mean = &mean_fn<logistic_d_handle, double>;
        out->variance = &variance_fn<logistic_d_handle, double>;
        out->skewness = &skewness_fn<logistic_d_handle, double>;
        out->kurtosis = &kurtosis_fn<logistic_d_handle, double>;
        out->kurtosis_excess = &kurtosis_excess_fn<logistic_d_handle, double>;
        out->mode = &mode_fn<logistic_d_handle, double>;
        out->median = &median_fn<logistic_d_handle, double>;
        out->entropy = &entropy_fn<logistic_d_handle, double>;
        out->free = &free_fn<logistic_d_handle>;
        return true;
    } else if (n == "lognormal" || n == "log_normal" || n == "lognormal_distribution") {
        double mu = 0, sigma = 1;
        const char* muKeys[] = { "location", "loc", "mu", "meanlog" };
        const char* sigmaKeys[] = { "scale", "sigma", "sd", "standard_deviation" };
        find_param(params, count, muKeys, 4, &mu);
        find_param(params, count, sigmaKeys, 4, &sigma);
        if (!(sigma > 0)) return false;
        auto* h = new (std::nothrow) lognormal_d_handle(mu, sigma);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<lognormal_d_handle, double>;
        out->logpdf = &logpdf_fn<lognormal_d_handle, double>;
        out->cdf = &cdf_fn<lognormal_d_handle, double>;
        out->sf = &sf_fn<lognormal_d_handle, double>;
        out->hazard = &hazard_fn<lognormal_d_handle, double>;
        out->chf = &chf_fn<lognormal_d_handle, double>;
        out->quantile = &quantile_fn<lognormal_d_handle, double>;
        out->quantile_complement = &quantile_complement_fn<lognormal_d_handle, double>;
        out->range = &range_fn<lognormal_d_handle, bs_range_d, double>;
        out->mean = &mean_fn<lognormal_d_handle, double>;
        out->variance = &variance_fn<lognormal_d_handle, double>;
        out->skewness = &skewness_fn<lognormal_d_handle, double>;
        out->kurtosis = &kurtosis_fn<lognormal_d_handle, double>;
        out->kurtosis_excess = &kurtosis_excess_fn<lognormal_d_handle, double>;
        out->mode = &mode_fn<lognormal_d_handle, double>;
        out->median = &median_fn<lognormal_d_handle, double>;
        out->entropy = &entropy_fn<lognormal_d_handle, double>;
        out->free = &free_fn<lognormal_d_handle>;
        return true;
    } else if (n == "mapairy" || n == "map_airy" || n == "mapairy_distribution") {
        double loc = 0, scale = 1;
        const char* locKeys[] = { "location", "loc", "mu" };
        const char* scaleKeys[] = { "scale", "c", "sigma" };
        find_param(params, count, locKeys, 3, &loc);
        find_param(params, count, scaleKeys, 3, &scale);
        if (!(scale > 0)) return false;
        auto* h = new (std::nothrow) mapairy_d_handle(loc, scale);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<mapairy_d_handle, double>;
        out->logpdf = &logpdf_fn<mapairy_d_handle, double>;
        out->cdf = &cdf_fn<mapairy_d_handle, double>;
        out->sf = &sf_fn<mapairy_d_handle, double>;
        out->hazard = &hazard_fn<mapairy_d_handle, double>;
        out->chf = &chf_fn<mapairy_d_handle, double>;
        out->quantile = &quantile_fn<mapairy_d_handle, double>;
        out->quantile_complement = &quantile_complement_fn<mapairy_d_handle, double>;
        out->range = &range_fn<mapairy_d_handle, bs_range_d, double>;
        out->mean = &mean_fn<mapairy_d_handle, double>;
        out->variance = &variance_fn<mapairy_d_handle, double>;
        out->skewness = nullptr;
        out->kurtosis = nullptr;
        out->kurtosis_excess = nullptr;
        out->mode = &mode_fn<mapairy_d_handle, double>;
        out->median = &median_fn<mapairy_d_handle, double>;
        out->entropy = &entropy_fn<mapairy_d_handle, double>;
        out->free = &free_fn<mapairy_d_handle>;
        return true;
    } else if (n == "rayleigh" || n == "rayleigh_distribution") {
        double scale = 0;
        const char* scaleKeys[] = { "scale", "sigma", "beta" };
        if (!find_param(params, count, scaleKeys, 3, &scale)) return false;
        if (!(scale > 0)) return false;
        auto* h = new (std::nothrow) rayleigh_d_handle(scale);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<rayleigh_d_handle, double>;
        out->logpdf = &logpdf_fn<rayleigh_d_handle, double>;
        out->cdf = &cdf_fn<rayleigh_d_handle, double>;
        out->sf = &sf_fn<rayleigh_d_handle, double>;
        out->hazard = &hazard_fn<rayleigh_d_handle, double>;
        out->chf = &chf_fn<rayleigh_d_handle, double>;
        out->quantile = &quantile_fn<rayleigh_d_handle, double>;
        out->quantile_complement = &quantile_complement_fn<rayleigh_d_handle, double>;
        out->range = &range_fn<rayleigh_d_handle, bs_range_d, double>;
        out->mean = &mean_fn<rayleigh_d_handle, double>;
        out->variance = &variance_fn<rayleigh_d_handle, double>;
        out->skewness = &skewness_fn<rayleigh_d_handle, double>;
        out->kurtosis = &kurtosis_fn<rayleigh_d_handle, double>;
        out->kurtosis_excess = &kurtosis_excess_fn<rayleigh_d_handle, double>;
        out->mode = &mode_fn<rayleigh_d_handle, double>;
        out->median = &median_fn<rayleigh_d_handle, double>;
        out->entropy = &entropy_fn<rayleigh_d_handle, double>;
        out->free = &free_fn<rayleigh_d_handle>;
        return true;
    } else if (
        n == "saspoint5" || n == "sas_point5" || n == "sas_point_5" ||
        n == "saspointfive" || n == "sas0.5" || n == "sas_alpha_half" ||
        n == "stable_point5" || n == "stable_alpha_half"
    ) {
        double location = 0;
        double scale = 1;
        const char* locKeys[] = { "location", "loc", "mu" };
        const char* scaleKeys[] = { "scale", "sigma", "c" };
        find_param(params, count, locKeys, 3, &location);
        find_param(params, count, scaleKeys, 3, &scale);
        if (!(scale > 0)) return false;
        auto* h = new (std::nothrow) saspoint5_d_handle(location, scale);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<saspoint5_d_handle, double>;
        out->logpdf = &logpdf_fn<saspoint5_d_handle, double>;
        out->cdf = &cdf_fn<saspoint5_d_handle, double>;
        out->sf = &sf_fn<saspoint5_d_handle, double>;
        out->hazard = &hazard_fn<saspoint5_d_handle, double>;
        out->chf = &chf_fn<saspoint5_d_handle, double>;
        out->quantile = &quantile_fn<saspoint5_d_handle, double>;
        out->quantile_complement = &quantile_complement_fn<saspoint5_d_handle, double>;
        out->range = &range_fn<saspoint5_d_handle, bs_range_d, double>;
        out->mean = nullptr;
        out->variance = nullptr;
        out->skewness = nullptr;
        out->kurtosis = nullptr;
        out->kurtosis_excess = nullptr;
        out->mode = &mode_fn<saspoint5_d_handle, double>;
        out->median = &median_fn<saspoint5_d_handle, double>;
        out->entropy = &entropy_fn<saspoint5_d_handle, double>;
        out->free = &free_fn<saspoint5_d_handle>;
        return true;
    } else if (n == "skew_normal" || n == "skewnormal" || n == "skew_normal_distribution") {
        double location = 0;
        double scale = 1;
        double shape = 0;
        const char* locKeys[] = { "location", "loc", "mu" };
        const char* scaleKeys[] = { "scale", "sigma", "omega" };
        const char* shapeKeys[] = { "shape", "alpha", "skew" };
        find_param(params, count, locKeys, 3, &location);
        find_param(params, count, scaleKeys, 3, &scale);
        find_param(params, count, shapeKeys, 3, &shape);
        if (!(scale > 0)) return false;
        auto* h = new (std::nothrow) skew_normal_d_handle(location, scale, shape);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<skew_normal_d_handle, double>;
        out->logpdf = &logpdf_fn<skew_normal_d_handle, double>;
        out->cdf = &cdf_fn<skew_normal_d_handle, double>;
        out->sf = &sf_fn<skew_normal_d_handle, double>;
        out->hazard = &hazard_fn<skew_normal_d_handle, double>;
        out->chf = &chf_fn<skew_normal_d_handle, double>;
        out->quantile = &quantile_fn<skew_normal_d_handle, double>;
        out->quantile_complement = &quantile_complement_fn<skew_normal_d_handle, double>;
        out->range = &range_fn<skew_normal_d_handle, bs_range_d, double>;
        out->mean = &mean_fn<skew_normal_d_handle, double>;
        out->variance = &variance_fn<skew_normal_d_handle, double>;
        out->skewness = &skewness_fn<skew_normal_d_handle, double>;
        out->kurtosis = &kurtosis_fn<skew_normal_d_handle, double>;
        out->kurtosis_excess = &kurtosis_excess_fn<skew_normal_d_handle, double>;
        out->mode = &mode_fn<skew_normal_d_handle, double>;
        out->median = &median_via_quantile_fn<skew_normal_d_handle, double>;
        out->entropy = nullptr;
        out->free = &free_fn<skew_normal_d_handle>;
        return true;
    } else if (n == "triangular" || n == "triangular_distribution" || n == "triangle") {
        double lower = 0;
        double upper = 0;
        double mode = 0;
        const char* lowerKeys[] = { "lower", "min", "minimum", "a" };
        const char* upperKeys[] = { "upper", "max", "maximum", "b" };
        const char* modeKeys[] = { "mode", "peak", "c" };
        if (!find_param(params, count, lowerKeys, 4, &lower)) return false;
        if (!find_param(params, count, upperKeys, 4, &upper)) return false;
        if (!find_param(params, count, modeKeys, 3, &mode)) return false;
        if (!(upper > lower)) return false;
        if (!(mode >= lower && mode <= upper)) return false;
        auto* h = new (std::nothrow) triangular_d_handle(lower, mode, upper);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<triangular_d_handle, double>;
        out->logpdf = &logpdf_fn<triangular_d_handle, double>;
        out->cdf = &cdf_fn<triangular_d_handle, double>;
        out->sf = &sf_fn<triangular_d_handle, double>;
        out->hazard = &hazard_fn<triangular_d_handle, double>;
        out->chf = &chf_fn<triangular_d_handle, double>;
        out->quantile = &quantile_fn<triangular_d_handle, double>;
        out->quantile_complement = &quantile_complement_fn<triangular_d_handle, double>;
        out->range = &range_fn<triangular_d_handle, bs_range_d, double>;
        out->mean = &mean_fn<triangular_d_handle, double>;
        out->variance = &variance_fn<triangular_d_handle, double>;
        out->skewness = &skewness_fn<triangular_d_handle, double>;
        out->kurtosis = &kurtosis_fn<triangular_d_handle, double>;
        out->kurtosis_excess = &kurtosis_excess_fn<triangular_d_handle, double>;
        out->mode = &mode_fn<triangular_d_handle, double>;
        out->median = &median_fn<triangular_d_handle, double>;
        out->entropy = &entropy_fn<triangular_d_handle, double>;
        out->free = &free_fn<triangular_d_handle>;
        return true;
    } else if (
        n == "uniform" || n == "uniform_distribution" || n == "uniform_real" ||
        n == "rectangular" || n == "rectangular_distribution"
    ) {
        double lower = 0;
        double upper = 1;
        const char* lowerKeys[] = { "lower", "min", "minimum", "a" };
        const char* upperKeys[] = { "upper", "max", "maximum", "b" };
        find_param(params, count, lowerKeys, 4, &lower);
        find_param(params, count, upperKeys, 4, &upper);
        if (!(upper > lower)) return false;
        auto* h = new (std::nothrow) uniform_d_handle(lower, upper);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<uniform_d_handle, double>;
        out->logpdf = &logpdf_fn<uniform_d_handle, double>;
        out->cdf = &cdf_fn<uniform_d_handle, double>;
        out->sf = &sf_fn<uniform_d_handle, double>;
        out->hazard = &hazard_fn<uniform_d_handle, double>;
        out->chf = &chf_fn<uniform_d_handle, double>;
        out->quantile = &quantile_fn<uniform_d_handle, double>;
        out->quantile_complement = &quantile_complement_fn<uniform_d_handle, double>;
        out->range = &range_fn<uniform_d_handle, bs_range_d, double>;
        out->mean = &mean_fn<uniform_d_handle, double>;
        out->variance = &variance_fn<uniform_d_handle, double>;
        out->skewness = &skewness_fn<uniform_d_handle, double>;
        out->kurtosis = &kurtosis_fn<uniform_d_handle, double>;
        out->kurtosis_excess = &kurtosis_excess_fn<uniform_d_handle, double>;
        out->mode = &mode_fn<uniform_d_handle, double>;
        out->median = &median_fn<uniform_d_handle, double>;
        out->entropy = &entropy_fn<uniform_d_handle, double>;
        out->free = &free_fn<uniform_d_handle>;
        return true;
    } else if (n == "weibull" || n == "weibull_distribution") {
        double shape = 0;
        double scale = 1;
        const char* shapeKeys[] = { "shape", "k", "alpha" };
        const char* scaleKeys[] = { "scale", "lambda", "beta" };
        if (!find_param(params, count, shapeKeys, 3, &shape)) return false;
        find_param(params, count, scaleKeys, 3, &scale);
        if (!(shape > 0) || !(scale > 0)) return false;
        auto* h = new (std::nothrow) weibull_d_handle(shape, scale);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<weibull_d_handle, double>;
        out->logpdf = &logpdf_fn<weibull_d_handle, double>;
        out->cdf = &cdf_fn<weibull_d_handle, double>;
        out->sf = &sf_fn<weibull_d_handle, double>;
        out->hazard = &hazard_fn<weibull_d_handle, double>;
        out->chf = &chf_fn<weibull_d_handle, double>;
        out->quantile = &quantile_fn<weibull_d_handle, double>;
        out->quantile_complement = &quantile_complement_fn<weibull_d_handle, double>;
        out->range = &range_fn<weibull_d_handle, bs_range_d, double>;
        out->mean = &mean_fn<weibull_d_handle, double>;
        out->variance = &variance_fn<weibull_d_handle, double>;
        out->skewness = &skewness_fn<weibull_d_handle, double>;
        out->kurtosis = &kurtosis_fn<weibull_d_handle, double>;
        out->kurtosis_excess = &kurtosis_excess_fn<weibull_d_handle, double>;
        out->mode = &mode_fn<weibull_d_handle, double>;
        out->median = &median_fn<weibull_d_handle, double>;
        out->entropy = &entropy_fn<weibull_d_handle, double>;
        out->free = &free_fn<weibull_d_handle>;
        return true;
    } else if (n == "kolmogorov_smirnov" || n == "kolmogorov-smirnov" || n == "kolmogorovsmirnov" || n == "kolmogorov_smirnov_distribution" || n == "ks" || n == "ks_distribution") {
        double observations = 0;
        const char* nKeys[] = { "n", "sample_count", "samplecount", "sample_size", "samples", "observations" };
        if (!find_param(params, count, nKeys, 6, &observations)) return false;
        auto* h = new (std::nothrow) kolmogorov_smirnov_d_handle(observations);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<kolmogorov_smirnov_d_handle, double>;
        out->logpdf = &logpdf_fn<kolmogorov_smirnov_d_handle, double>;
        out->cdf = &cdf_fn<kolmogorov_smirnov_d_handle, double>;
        out->sf = &sf_fn<kolmogorov_smirnov_d_handle, double>;
        out->hazard = &hazard_fn<kolmogorov_smirnov_d_handle, double>;
        out->chf = &chf_fn<kolmogorov_smirnov_d_handle, double>;
        out->quantile = &quantile_fn<kolmogorov_smirnov_d_handle, double>;
        out->quantile_complement = &quantile_complement_fn<kolmogorov_smirnov_d_handle, double>;
        out->range = &range_fn<kolmogorov_smirnov_d_handle, bs_range_d, double>;
        out->mean = &mean_fn<kolmogorov_smirnov_d_handle, double>;
        out->variance = &variance_fn<kolmogorov_smirnov_d_handle, double>;
        out->skewness = &skewness_fn<kolmogorov_smirnov_d_handle, double>;
        out->kurtosis = &kurtosis_fn<kolmogorov_smirnov_d_handle, double>;
        out->kurtosis_excess = &kurtosis_excess_fn<kolmogorov_smirnov_d_handle, double>;
        out->mode = &mode_fn<kolmogorov_smirnov_d_handle, double>;
        out->median = &median_via_quantile_fn<kolmogorov_smirnov_d_handle, double>;
        out->entropy = nullptr;
        out->free = &free_fn<kolmogorov_smirnov_d_handle>;
        return true;
    } else if (n == "landau" || n == "landau_distribution") {
        double loc = 0, scale = 1;
        const char* locKeys[] = { "location", "loc", "mu" };
        const char* scaleKeys[] = { "scale", "c", "sigma" };
        find_param(params, count, locKeys, 3, &loc);
        find_param(params, count, scaleKeys, 3, &scale);
        auto* h = new (std::nothrow) landau_d_handle(loc, scale);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<landau_d_handle, double>;
        out->logpdf = &logpdf_fn<landau_d_handle, double>;
        out->cdf = &cdf_fn<landau_d_handle, double>;
        out->sf = &sf_fn<landau_d_handle, double>;
        out->hazard = &hazard_fn<landau_d_handle, double>;
        out->chf = &chf_fn<landau_d_handle, double>;
        out->quantile = &quantile_fn<landau_d_handle, double>;
        out->quantile_complement = &quantile_complement_fn<landau_d_handle, double>;
        out->range = &range_fn<landau_d_handle, bs_range_d, double>;
        out->mean = nullptr;
        out->variance = nullptr;
        out->skewness = nullptr;
        out->kurtosis = nullptr;
        out->kurtosis_excess = nullptr;
        out->mode = &mode_fn<landau_d_handle, double>;
        out->median = &median_fn<landau_d_handle, double>;
        out->entropy = &entropy_fn<landau_d_handle, double>;
        out->free = &free_fn<landau_d_handle>;
        return true;
    } else if (n == "laplace" || n == "laplace_distribution" || n == "double_exponential" || n == "doubleexponential") {
        double loc = 0, scale = 1;
        const char* locKeys[] = { "location", "loc", "mu", "mean" };
        const char* scaleKeys[] = { "scale", "diversity", "b" };
        find_param(params, count, locKeys, 4, &loc);
        find_param(params, count, scaleKeys, 3, &scale);
        auto* h = new (std::nothrow) laplace_d_handle(loc, scale);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<laplace_d_handle, double>;
        out->logpdf = &logpdf_fn<laplace_d_handle, double>;
        out->cdf = &cdf_fn<laplace_d_handle, double>;
        out->sf = &sf_fn<laplace_d_handle, double>;
        out->hazard = &hazard_fn<laplace_d_handle, double>;
        out->chf = &chf_fn<laplace_d_handle, double>;
        out->quantile = &quantile_fn<laplace_d_handle, double>;
        out->quantile_complement = &quantile_complement_fn<laplace_d_handle, double>;
        out->range = &range_fn<laplace_d_handle, bs_range_d, double>;
        out->mean = &mean_fn<laplace_d_handle, double>;
        out->variance = &variance_fn<laplace_d_handle, double>;
        out->skewness = &skewness_fn<laplace_d_handle, double>;
        out->kurtosis = &kurtosis_fn<laplace_d_handle, double>;
        out->kurtosis_excess = &kurtosis_excess_fn<laplace_d_handle, double>;
        out->mode = &mode_fn<laplace_d_handle, double>;
        out->median = &median_fn<laplace_d_handle, double>;
        out->entropy = &entropy_fn<laplace_d_handle, double>;
        out->free = &free_fn<laplace_d_handle>;
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
    } else if (n == "non_central_beta" || n == "non_central_beta_distribution") {
        float a = 0, b = 0, l = 0;
        const char* aKeys[] = { "alpha", "a", "p", "shape1" };
        const char* bKeys[] = { "beta", "b", "q", "shape2" };
        const char* lKeys[] = {"lambda", "l"};
        if (!find_param(params, count, aKeys, 4, &a)) return false;
        if (!find_param(params, count, bKeys, 4, &b)) return false;
        if (!find_param(params, count, lKeys, 2, &l)) return false;
        auto* h = new (std::nothrow) non_central_beta_f_handle(a, b, l);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<non_central_beta_f_handle, float>;
        out->logpdf = &logpdf_fn<non_central_beta_f_handle, float>;
        out->cdf = &cdf_fn<non_central_beta_f_handle, float>;
        out->sf = &sf_fn<non_central_beta_f_handle, float>;
        out->hazard = &hazard_fn<non_central_beta_f_handle, float>;
        out->chf = &chf_fn<non_central_beta_f_handle, float>;
        out->quantile = &quantile_fn<non_central_beta_f_handle, float>;
        out->quantile_complement = &quantile_complement_fn<non_central_beta_f_handle, float>;
        out->range = &range_fn<non_central_beta_f_handle, bs_range_f, float>;
        out->mean = &mean_fn<non_central_beta_f_handle, float>;
        out->variance = &variance_fn<non_central_beta_f_handle, float>;
        out->skewness = nullptr;;
        out->kurtosis = nullptr;;
        out->kurtosis_excess = nullptr;;
        out->mode = &mode_fn<non_central_beta_f_handle, float>;
        out->median = &median_fn<non_central_beta_f_handle, float>;
        out->entropy = nullptr; // not provided in Boost
        out->free = &free_fn<non_central_beta_f_handle>;
        return true;
    } else if (n == "non_central_chi_squared" || n == "noncentral_chi_squared" || n == "noncentralchisquared" || n == "non_central_chi2" || n == "noncentral_chi2" || n == "nc_chi_squared" || n == "ncchisquared") {
        float v = 0;
        float lambda = 0;
        const char* dfKeys[] = { "df", "nu", "degreesoffreedom" };
        const char* lambdaKeys[] = { "lambda", "noncentrality", "delta", "nc" };
        if (!find_param(params, count, dfKeys, 3, &v)) return false;
        if (!find_param(params, count, lambdaKeys, 4, &lambda)) return false;
        if (!(v > 0) || !(lambda >= 0)) return false;
        auto* h = new (std::nothrow) non_central_chi_squared_f_handle(v, lambda);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<non_central_chi_squared_f_handle, float>;
        out->logpdf = &logpdf_fn<non_central_chi_squared_f_handle, float>;
        out->cdf = &cdf_fn<non_central_chi_squared_f_handle, float>;
        out->sf = &sf_fn<non_central_chi_squared_f_handle, float>;
        out->hazard = &hazard_fn<non_central_chi_squared_f_handle, float>;
        out->chf = &chf_fn<non_central_chi_squared_f_handle, float>;
        out->quantile = &quantile_fn<non_central_chi_squared_f_handle, float>;
        out->quantile_complement = &quantile_complement_fn<non_central_chi_squared_f_handle, float>;
        out->range = &range_fn<non_central_chi_squared_f_handle, bs_range_f, float>;
        out->mean = &mean_fn<non_central_chi_squared_f_handle, float>;
        out->variance = &variance_fn<non_central_chi_squared_f_handle, float>;
        out->skewness = &skewness_fn<non_central_chi_squared_f_handle, float>;
        out->kurtosis = &kurtosis_fn<non_central_chi_squared_f_handle, float>;
        out->kurtosis_excess = &kurtosis_excess_fn<non_central_chi_squared_f_handle, float>;
        out->mode = &mode_fn<non_central_chi_squared_f_handle, float>;
        out->median = &median_fn<non_central_chi_squared_f_handle, float>;
        out->entropy = nullptr;
        out->free = &free_fn<non_central_chi_squared_f_handle>;
        return true;
    } else if (n == "non_central_f" || n == "noncentral_f" || n == "noncentralf" || n == "non_central_f_ratio" || n == "noncentral_f_ratio" || n == "nc_f" || n == "ncf") {
        float df1 = 0;
        float df2 = 0;
        float lambda = 0;
        const char* df1Keys[] = { "df1", "d1", "m", "degreesoffreedom1" };
        const char* df2Keys[] = { "df2", "d2", "n", "degreesoffreedom2" };
        const char* lambdaKeys[] = { "lambda", "noncentrality", "delta", "nc" };
        if (!find_param(params, count, df1Keys, 4, &df1)) return false;
        if (!find_param(params, count, df2Keys, 4, &df2)) return false;
        if (!find_param(params, count, lambdaKeys, 4, &lambda)) return false;
        if (!(df1 > 0) || !(df2 > 0) || !(lambda >= 0)) return false;
        auto* h = new (std::nothrow) non_central_f_f_handle(df1, df2, lambda);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<non_central_f_f_handle, float>;
        out->logpdf = &logpdf_fn<non_central_f_f_handle, float>;
        out->cdf = &cdf_fn<non_central_f_f_handle, float>;
        out->sf = &sf_fn<non_central_f_f_handle, float>;
        out->hazard = &hazard_fn<non_central_f_f_handle, float>;
        out->chf = &chf_fn<non_central_f_f_handle, float>;
        out->quantile = &quantile_fn<non_central_f_f_handle, float>;
        out->quantile_complement = &quantile_complement_fn<non_central_f_f_handle, float>;
        out->range = &range_fn<non_central_f_f_handle, bs_range_f, float>;
        out->mean = &mean_fn<non_central_f_f_handle, float>;
        out->variance = &variance_fn<non_central_f_f_handle, float>;
        out->skewness = &skewness_fn<non_central_f_f_handle, float>;
        out->kurtosis = &kurtosis_fn<non_central_f_f_handle, float>;
        out->kurtosis_excess = &kurtosis_excess_fn<non_central_f_f_handle, float>;
        out->mode = &mode_fn<non_central_f_f_handle, float>;
        out->median = &median_fn<non_central_f_f_handle, float>;
        out->entropy = nullptr;
        out->free = &free_fn<non_central_f_f_handle>;
        return true;
    } else if (n == "non_central_t" || n == "noncentral_t" || n == "noncentralt" || n == "non_central_student_t" || n == "noncentral_student_t" || n == "noncentralstudentt" || n == "nc_t") {
        float v = 0;
        float lambda = 0;
        const char* dfKeys[] = { "df", "nu", "degreesoffreedom" };
        const char* lambdaKeys[] = { "lambda", "noncentrality", "delta", "nc" };
        if (!find_param(params, count, dfKeys, 3, &v)) return false;
        if (!find_param(params, count, lambdaKeys, 4, &lambda)) return false;
        if (!(v > 0)) return false;
        if (!std::isfinite(lambda)) return false;
        auto* h = new (std::nothrow) non_central_t_f_handle(v, lambda);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<non_central_t_f_handle, float>;
        out->logpdf = &logpdf_fn<non_central_t_f_handle, float>;
        out->cdf = &cdf_fn<non_central_t_f_handle, float>;
        out->sf = &sf_fn<non_central_t_f_handle, float>;
        out->hazard = &hazard_fn<non_central_t_f_handle, float>;
        out->chf = &chf_fn<non_central_t_f_handle, float>;
        out->quantile = &quantile_fn<non_central_t_f_handle, float>;
        out->quantile_complement = &quantile_complement_fn<non_central_t_f_handle, float>;
        out->range = &range_fn<non_central_t_f_handle, bs_range_f, float>;
        out->mean = &mean_fn<non_central_t_f_handle, float>;
        out->variance = &variance_fn<non_central_t_f_handle, float>;
        out->skewness = &skewness_fn<non_central_t_f_handle, float>;
        out->kurtosis = &kurtosis_fn<non_central_t_f_handle, float>;
        out->kurtosis_excess = &kurtosis_excess_fn<non_central_t_f_handle, float>;
        out->mode = &mode_fn<non_central_t_f_handle, float>;
        out->median = &median_fn<non_central_t_f_handle, float>;
        out->entropy = nullptr;
        out->free = &free_fn<non_central_t_f_handle>;
        return true;
    } else if (n == "pareto" || n == "pareto_distribution") {
        float scale = 0;
        float shape = 0;
        const char* scaleKeys[] = { "scale", "xm", "minimum", "lower", "x0" };
        const char* shapeKeys[] = { "shape", "alpha" };
        if (!find_param(params, count, scaleKeys, 5, &scale)) return false;
        if (!find_param(params, count, shapeKeys, 2, &shape)) return false;
        if (!(scale > 0) || !(shape > 0)) return false;
        auto* h = new (std::nothrow) pareto_f_handle(scale, shape);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<pareto_f_handle, float>;
        out->logpdf = &logpdf_fn<pareto_f_handle, float>;
        out->cdf = &cdf_fn<pareto_f_handle, float>;
        out->sf = &sf_fn<pareto_f_handle, float>;
        out->hazard = &hazard_fn<pareto_f_handle, float>;
        out->chf = &chf_fn<pareto_f_handle, float>;
        out->quantile = &quantile_fn<pareto_f_handle, float>;
        out->quantile_complement = &quantile_complement_fn<pareto_f_handle, float>;
        out->range = &range_fn<pareto_f_handle, bs_range_f, float>;
        out->mean = &mean_fn<pareto_f_handle, float>;
        out->variance = &variance_fn<pareto_f_handle, float>;
        out->skewness = &skewness_fn<pareto_f_handle, float>;
        out->kurtosis = &kurtosis_fn<pareto_f_handle, float>;
        out->kurtosis_excess = &kurtosis_excess_fn<pareto_f_handle, float>;
        out->mode = &mode_fn<pareto_f_handle, float>;
        out->median = &median_fn<pareto_f_handle, float>;
        out->entropy = nullptr;
        out->free = &free_fn<pareto_f_handle>;
        return true;
    } else if (n == "poisson" || n == "poisson_distribution") {
        float mean = 0;
        const char* meanKeys[] = { "mean", "lambda", "mu" };
        if (!find_param(params, count, meanKeys, 3, &mean)) return false;
        if (!(mean >= 0)) return false;
        auto* h = new (std::nothrow) poisson_f_handle(mean);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<poisson_f_handle, float>;
        out->logpdf = &logpdf_fn<poisson_f_handle, float>;
        out->cdf = &cdf_fn<poisson_f_handle, float>;
        out->sf = &sf_fn<poisson_f_handle, float>;
        out->hazard = &hazard_fn<poisson_f_handle, float>;
        out->chf = &chf_fn<poisson_f_handle, float>;
        out->quantile = &quantile_fn<poisson_f_handle, float>;
        out->quantile_complement = &quantile_complement_fn<poisson_f_handle, float>;
        out->range = &range_fn<poisson_f_handle, bs_range_f, float>;
        out->mean = &mean_fn<poisson_f_handle, float>;
        out->variance = &variance_fn<poisson_f_handle, float>;
        out->skewness = &skewness_fn<poisson_f_handle, float>;
        out->kurtosis = &kurtosis_fn<poisson_f_handle, float>;
        out->kurtosis_excess = &kurtosis_excess_fn<poisson_f_handle, float>;
        out->mode = &mode_fn<poisson_f_handle, float>;
        out->median = &median_fn<poisson_f_handle, float>;
        out->entropy = nullptr;
        out->free = &free_fn<poisson_f_handle>;
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
    } else if (n == "non_central_chi_squared" || n == "noncentral_chi_squared" || n == "noncentralchisquared" || n == "non_central_chi2" || n == "noncentral_chi2" || n == "nc_chi_squared" || n == "ncchisquared") {
        float v = 0;
        float lambda = 0;
        const char* dfKeys[] = { "df", "nu", "degreesoffreedom" };
        const char* lambdaKeys[] = { "lambda", "noncentrality", "delta", "nc" };
        if (!find_param(params, count, dfKeys, 3, &v)) return false;
        if (!find_param(params, count, lambdaKeys, 4, &lambda)) return false;
        if (!(v > 0)) return false;
        if (!(lambda >= 0)) return false;
        auto* h = new (std::nothrow) non_central_chi_squared_f_handle(v, lambda);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<non_central_chi_squared_f_handle, float>;
        out->logpdf = &logpdf_fn<non_central_chi_squared_f_handle, float>;
        out->cdf = &cdf_fn<non_central_chi_squared_f_handle, float>;
        out->sf = &sf_fn<non_central_chi_squared_f_handle, float>;
        out->hazard = &hazard_fn<non_central_chi_squared_f_handle, float>;
        out->chf = &chf_fn<non_central_chi_squared_f_handle, float>;
        out->quantile = &quantile_fn<non_central_chi_squared_f_handle, float>;
        out->quantile_complement = &quantile_complement_fn<non_central_chi_squared_f_handle, float>;
        out->range = &range_fn<non_central_chi_squared_f_handle, bs_range_f, float>;
        out->mean = &mean_fn<non_central_chi_squared_f_handle, float>;
        out->variance = &variance_fn<non_central_chi_squared_f_handle, float>;
        out->skewness = &skewness_fn<non_central_chi_squared_f_handle, float>;
        out->kurtosis = &kurtosis_fn<non_central_chi_squared_f_handle, float>;
        out->kurtosis_excess = &kurtosis_excess_fn<non_central_chi_squared_f_handle, float>;
        out->mode = &mode_fn<non_central_chi_squared_f_handle, float>;
        out->median = &median_fn<non_central_chi_squared_f_handle, float>;
        out->entropy = nullptr;
        out->free = &free_fn<non_central_chi_squared_f_handle>;
        return true;
    } else if (n == "inverse_chi_squared" || n == "inversechisquared" || n == "inv_chi_squared" || n == "invchisquared" || n == "inverse_chi2" || n == "inv_chi2") {
        float df = 0;
        const char* dfKeys[] = { "df", "nu", "degreesoffreedom", "v" };
        if (!find_param(params, count, dfKeys, 4, &df)) return false;
        float scale = 0;
        const char* scaleKeys[] = { "scale", "sigma2", "xi" };
        bool hasScale = find_param(params, count, scaleKeys, 3, &scale);
        if (!hasScale) {
            if (!(df > 0)) return false;
            scale = 1.0f / df;
        }
        auto* h = new (std::nothrow) inverse_chi_squared_f_handle(df, scale);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<inverse_chi_squared_f_handle, float>;
        out->logpdf = &logpdf_fn<inverse_chi_squared_f_handle, float>;
        out->cdf = &cdf_fn<inverse_chi_squared_f_handle, float>;
        out->sf = &sf_fn<inverse_chi_squared_f_handle, float>;
        out->hazard = &hazard_fn<inverse_chi_squared_f_handle, float>;
        out->chf = &chf_fn<inverse_chi_squared_f_handle, float>;
        out->quantile = &quantile_fn<inverse_chi_squared_f_handle, float>;
        out->quantile_complement = &quantile_complement_fn<inverse_chi_squared_f_handle, float>;
        out->range = &range_fn<inverse_chi_squared_f_handle, bs_range_f, float>;
        out->mean = &mean_fn<inverse_chi_squared_f_handle, float>;
        out->variance = &variance_fn<inverse_chi_squared_f_handle, float>;
        out->skewness = &skewness_fn<inverse_chi_squared_f_handle, float>;
        out->kurtosis = &kurtosis_fn<inverse_chi_squared_f_handle, float>;
        out->kurtosis_excess = &kurtosis_excess_fn<inverse_chi_squared_f_handle, float>;
        out->mode = &mode_fn<inverse_chi_squared_f_handle, float>;
        out->median = &median_fn<inverse_chi_squared_f_handle, float>;
        out->entropy = nullptr;
        out->free = &free_fn<inverse_chi_squared_f_handle>;
        return true;
    } else if (n == "inverse_gamma" || n == "inversegamma" || n == "inv_gamma" || n == "invgamma") {
        float shape = 0;
        const char* shapeKeys[] = { "shape", "alpha", "k" };
        if (!find_param(params, count, shapeKeys, 3, &shape)) return false;
        float scale = 1;
        const char* scaleKeys[] = { "scale", "theta", "beta" };
        find_param(params, count, scaleKeys, 3, &scale);
        auto* h = new (std::nothrow) inverse_gamma_f_handle(shape, scale);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<inverse_gamma_f_handle, float>;
        out->logpdf = &logpdf_fn<inverse_gamma_f_handle, float>;
        out->cdf = &cdf_fn<inverse_gamma_f_handle, float>;
        out->sf = &sf_fn<inverse_gamma_f_handle, float>;
        out->hazard = &hazard_fn<inverse_gamma_f_handle, float>;
        out->chf = &chf_fn<inverse_gamma_f_handle, float>;
        out->quantile = &quantile_fn<inverse_gamma_f_handle, float>;
        out->quantile_complement = &quantile_complement_fn<inverse_gamma_f_handle, float>;
        out->range = &range_fn<inverse_gamma_f_handle, bs_range_f, float>;
        out->mean = &mean_fn<inverse_gamma_f_handle, float>;
        out->variance = &variance_fn<inverse_gamma_f_handle, float>;
        out->skewness = &skewness_fn<inverse_gamma_f_handle, float>;
        out->kurtosis = &kurtosis_fn<inverse_gamma_f_handle, float>;
        out->kurtosis_excess = &kurtosis_excess_fn<inverse_gamma_f_handle, float>;
        out->mode = &mode_fn<inverse_gamma_f_handle, float>;
        out->median = &median_fn<inverse_gamma_f_handle, float>;
        out->entropy = nullptr;
        out->free = &free_fn<inverse_gamma_f_handle>;
        return true;
    } else if (n == "inverse_gaussian" || n == "inversegaussian" || n == "inverse_normal" || n == "inversenormal" || n == "wald") {
        float mean = 0;
        const char* meanKeys[] = { "mean", "mu", "location" };
        if (!find_param(params, count, meanKeys, 3, &mean)) return false;
        float scale = 1;
        const char* scaleKeys[] = { "scale", "lambda", "shape" };
        find_param(params, count, scaleKeys, 3, &scale);
        auto* h = new (std::nothrow) inverse_gaussian_f_handle(mean, scale);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<inverse_gaussian_f_handle, float>;
        out->logpdf = &logpdf_fn<inverse_gaussian_f_handle, float>;
        out->cdf = &cdf_fn<inverse_gaussian_f_handle, float>;
        out->sf = &sf_fn<inverse_gaussian_f_handle, float>;
        out->hazard = &hazard_fn<inverse_gaussian_f_handle, float>;
        out->chf = &chf_fn<inverse_gaussian_f_handle, float>;
        out->quantile = &quantile_fn<inverse_gaussian_f_handle, float>;
        out->quantile_complement = &quantile_complement_fn<inverse_gaussian_f_handle, float>;
        out->range = &range_fn<inverse_gaussian_f_handle, bs_range_f, float>;
        out->mean = &mean_fn<inverse_gaussian_f_handle, float>;
        out->variance = &variance_fn<inverse_gaussian_f_handle, float>;
        out->skewness = &skewness_fn<inverse_gaussian_f_handle, float>;
        out->kurtosis = &kurtosis_fn<inverse_gaussian_f_handle, float>;
        out->kurtosis_excess = &kurtosis_excess_fn<inverse_gaussian_f_handle, float>;
        out->mode = &mode_fn<inverse_gaussian_f_handle, float>;
        out->median = &median_fn<inverse_gaussian_f_handle, float>;
        out->entropy = nullptr;
        out->free = &free_fn<inverse_gaussian_f_handle>;
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
    } else if (n == "negative_binomial" || n == "negativebinomial" || n == "neg_binomial" || n == "negative_binomial_distribution" || n == "nbinom") {
        double successes = 0, p = 0;
        const char* rKeys[] = { "r", "successes", "target", "count" };
        const char* pKeys[] = { "p", "prob", "probability", "success" };
        if (!find_param(params, count, rKeys, 4, &successes)) return false;
        if (!find_param(params, count, pKeys, 4, &p)) return false;
        if (!(successes > 0)) return false;
        if (!(p > 0 && p <= 1)) return false;
        auto* h = new (std::nothrow) negative_binomial_f_handle(static_cast<float>(successes), static_cast<float>(p));
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<negative_binomial_f_handle, float>;
        out->logpdf = &logpdf_fn<negative_binomial_f_handle, float>;
        out->cdf = &cdf_fn<negative_binomial_f_handle, float>;
        out->sf = &sf_fn<negative_binomial_f_handle, float>;
        out->hazard = &hazard_fn<negative_binomial_f_handle, float>;
        out->chf = &chf_fn<negative_binomial_f_handle, float>;
        out->quantile = &quantile_fn<negative_binomial_f_handle, float>;
        out->quantile_complement = &quantile_complement_fn<negative_binomial_f_handle, float>;
        out->range = &range_fn<negative_binomial_f_handle, bs_range_f, float>;
        out->mean = &mean_fn<negative_binomial_f_handle, float>;
        out->variance = &variance_fn<negative_binomial_f_handle, float>;
        out->skewness = &skewness_fn<negative_binomial_f_handle, float>;
        out->kurtosis = &kurtosis_fn<negative_binomial_f_handle, float>;
        out->kurtosis_excess = &kurtosis_excess_fn<negative_binomial_f_handle, float>;
        out->mode = &mode_fn<negative_binomial_f_handle, float>;
        out->median = &median_via_quantile_fn<negative_binomial_f_handle, float>;
        out->entropy = nullptr;
        out->free = &free_fn<negative_binomial_f_handle>;
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
    } else if (n == "hyperexponential" || n == "hyper_exponential" || n == "hyperexp" || n == "hyperexponential_distribution") {
        std::vector<float> probs;
        std::vector<float> rates;
        if (!collect_hyperex_parameters(params, count, probs, rates)) return false;
        hyperexponential_f_handle* h = nullptr;
        if (probs.empty()) {
            h = new (std::nothrow) hyperexponential_f_handle(rates);
        } else {
            h = new (std::nothrow) hyperexponential_f_handle(probs, rates);
        }
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<hyperexponential_f_handle, float>;
        out->logpdf = &logpdf_fn<hyperexponential_f_handle, float>;
        out->cdf = &cdf_fn<hyperexponential_f_handle, float>;
        out->sf = &sf_fn<hyperexponential_f_handle, float>;
        out->hazard = &hazard_fn<hyperexponential_f_handle, float>;
        out->chf = &chf_fn<hyperexponential_f_handle, float>;
        out->quantile = &quantile_fn<hyperexponential_f_handle, float>;
        out->quantile_complement = &quantile_complement_fn<hyperexponential_f_handle, float>;
        out->range = &range_fn<hyperexponential_f_handle, bs_range_f, float>;
        out->mean = &mean_fn<hyperexponential_f_handle, float>;
        out->variance = &variance_fn<hyperexponential_f_handle, float>;
        out->skewness = &skewness_fn<hyperexponential_f_handle, float>;
        out->kurtosis = &kurtosis_fn<hyperexponential_f_handle, float>;
        out->kurtosis_excess = &kurtosis_excess_fn<hyperexponential_f_handle, float>;
        out->mode = &mode_fn<hyperexponential_f_handle, float>;
        out->median = &median_fn<hyperexponential_f_handle, float>;
        out->entropy = nullptr;
        out->free = &free_fn<hyperexponential_f_handle>;
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
    } else if (n == "normal" || n == "normal_distribution" || n == "gauss" || n == "gaussian" || n == "gaussian_distribution" || n == "gauss_distribution") {
        double loc = 0, scale = 1;
        const char* locKeys[] = { "location", "loc", "mu", "mean" };
        const char* scaleKeys[] = { "sd", "standard_deviation", "sigma" };
        find_param(params, count, locKeys, 4, &loc);
        find_param(params, count, scaleKeys, 3, &scale);
        auto* h = new (std::nothrow) normal_f_handle(loc, scale);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<normal_f_handle, float>;
        out->logpdf = &logpdf_fn<normal_f_handle, float>;
        out->cdf = &cdf_fn<normal_f_handle, float>;
        out->sf = &sf_fn<normal_f_handle, float>;
        out->hazard = &hazard_fn<normal_f_handle, float>;
        out->chf = &chf_fn<normal_f_handle, float>;
        out->quantile = &quantile_fn<normal_f_handle, float>;
        out->quantile_complement = &quantile_complement_fn<normal_f_handle, float>;
        out->range = &range_fn<normal_f_handle, bs_range_f, float>;
        out->mean = &mean_fn<normal_f_handle, float>;
        out->variance = &variance_fn<normal_f_handle, float>;
        out->skewness = &skewness_fn<normal_f_handle, float>;
        out->kurtosis = &kurtosis_fn<normal_f_handle, float>;
        out->kurtosis_excess = &kurtosis_excess_fn<normal_f_handle, float>;;
        out->mode = &mode_fn<normal_f_handle, float>;
        out->median = &median_fn<normal_f_handle, float>;
        out->entropy = &entropy_fn<normal_f_handle, float>;;
        out->free = &free_fn<normal_f_handle>;
        return true;
    } else if (n == "logistic" || n == "logistic_distribution") {
        double loc = 0, scale = 1;
        const char* locKeys[] = { "location", "loc", "mu", "median" };
        const char* scaleKeys[] = { "scale", "s", "sigma", "diversity" };
        find_param(params, count, locKeys, 4, &loc);
        find_param(params, count, scaleKeys, 4, &scale);
        if (!(scale > 0)) return false;
        auto* h = new (std::nothrow) logistic_f_handle(static_cast<float>(loc), static_cast<float>(scale));
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<logistic_f_handle, float>;
        out->logpdf = &logpdf_fn<logistic_f_handle, float>;
        out->cdf = &cdf_fn<logistic_f_handle, float>;
        out->sf = &sf_fn<logistic_f_handle, float>;
        out->hazard = &hazard_fn<logistic_f_handle, float>;
        out->chf = &chf_fn<logistic_f_handle, float>;
        out->quantile = &quantile_fn<logistic_f_handle, float>;
        out->quantile_complement = &quantile_complement_fn<logistic_f_handle, float>;
        out->range = &range_fn<logistic_f_handle, bs_range_f, float>;
        out->mean = &mean_fn<logistic_f_handle, float>;
        out->variance = &variance_fn<logistic_f_handle, float>;
        out->skewness = &skewness_fn<logistic_f_handle, float>;
        out->kurtosis = &kurtosis_fn<logistic_f_handle, float>;
        out->kurtosis_excess = &kurtosis_excess_fn<logistic_f_handle, float>;
        out->mode = &mode_fn<logistic_f_handle, float>;
        out->median = &median_fn<logistic_f_handle, float>;
        out->entropy = &entropy_fn<logistic_f_handle, float>;
        out->free = &free_fn<logistic_f_handle>;
        return true;
    } else if (n == "lognormal" || n == "log_normal" || n == "lognormal_distribution") {
        double mu = 0, sigma = 1;
        const char* muKeys[] = { "location", "loc", "mu", "meanlog" };
        const char* sigmaKeys[] = { "scale", "sigma", "sd", "standard_deviation" };
        find_param(params, count, muKeys, 4, &mu);
        find_param(params, count, sigmaKeys, 4, &sigma);
        if (!(sigma > 0)) return false;
        auto* h = new (std::nothrow) lognormal_f_handle(static_cast<float>(mu), static_cast<float>(sigma));
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<lognormal_f_handle, float>;
        out->logpdf = &logpdf_fn<lognormal_f_handle, float>;
        out->cdf = &cdf_fn<lognormal_f_handle, float>;
        out->sf = &sf_fn<lognormal_f_handle, float>;
        out->hazard = &hazard_fn<lognormal_f_handle, float>;
        out->chf = &chf_fn<lognormal_f_handle, float>;
        out->quantile = &quantile_fn<lognormal_f_handle, float>;
        out->quantile_complement = &quantile_complement_fn<lognormal_f_handle, float>;
        out->range = &range_fn<lognormal_f_handle, bs_range_f, float>;
        out->mean = &mean_fn<lognormal_f_handle, float>;
        out->variance = &variance_fn<lognormal_f_handle, float>;
        out->skewness = &skewness_fn<lognormal_f_handle, float>;
        out->kurtosis = &kurtosis_fn<lognormal_f_handle, float>;
        out->kurtosis_excess = &kurtosis_excess_fn<lognormal_f_handle, float>;
        out->mode = &mode_fn<lognormal_f_handle, float>;
        out->median = &median_fn<lognormal_f_handle, float>;
        out->entropy = &entropy_fn<lognormal_f_handle, float>;
        out->free = &free_fn<lognormal_f_handle>;
        return true;
    } else if (n == "mapairy" || n == "map_airy" || n == "mapairy_distribution") {
        double loc = 0, scale = 1;
        const char* locKeys[] = { "location", "loc", "mu" };
        const char* scaleKeys[] = { "scale", "c", "sigma" };
        find_param(params, count, locKeys, 3, &loc);
        find_param(params, count, scaleKeys, 3, &scale);
        if (!(scale > 0)) return false;
        auto* h = new (std::nothrow) mapairy_f_handle(static_cast<float>(loc), static_cast<float>(scale));
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<mapairy_f_handle, float>;
        out->logpdf = &logpdf_fn<mapairy_f_handle, float>;
        out->cdf = &cdf_fn<mapairy_f_handle, float>;
        out->sf = &sf_fn<mapairy_f_handle, float>;
        out->hazard = &hazard_fn<mapairy_f_handle, float>;
        out->chf = &chf_fn<mapairy_f_handle, float>;
        out->quantile = &quantile_fn<mapairy_f_handle, float>;
        out->quantile_complement = &quantile_complement_fn<mapairy_f_handle, float>;
        out->range = &range_fn<mapairy_f_handle, bs_range_f, float>;
        out->mean = &mean_fn<mapairy_f_handle, float>;
        out->variance = &variance_fn<mapairy_f_handle, float>;
        out->skewness = nullptr;
        out->kurtosis = nullptr;
        out->kurtosis_excess = nullptr;
        out->mode = &mode_fn<mapairy_f_handle, float>;
        out->median = &median_fn<mapairy_f_handle, float>;
        out->entropy = &entropy_fn<mapairy_f_handle, float>;
        out->free = &free_fn<mapairy_f_handle>;
        return true;
    } else if (n == "rayleigh" || n == "rayleigh_distribution") {
        double scale = 0;
        const char* scaleKeys[] = { "scale", "sigma", "beta" };
        if (!find_param(params, count, scaleKeys, 3, &scale)) return false;
        if (!(scale > 0)) return false;
        auto* h = new (std::nothrow) rayleigh_f_handle(static_cast<float>(scale));
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<rayleigh_f_handle, float>;
        out->logpdf = &logpdf_fn<rayleigh_f_handle, float>;
        out->cdf = &cdf_fn<rayleigh_f_handle, float>;
        out->sf = &sf_fn<rayleigh_f_handle, float>;
        out->hazard = &hazard_fn<rayleigh_f_handle, float>;
        out->chf = &chf_fn<rayleigh_f_handle, float>;
        out->quantile = &quantile_fn<rayleigh_f_handle, float>;
        out->quantile_complement = &quantile_complement_fn<rayleigh_f_handle, float>;
        out->range = &range_fn<rayleigh_f_handle, bs_range_f, float>;
        out->mean = &mean_fn<rayleigh_f_handle, float>;
        out->variance = &variance_fn<rayleigh_f_handle, float>;
        out->skewness = &skewness_fn<rayleigh_f_handle, float>;
        out->kurtosis = &kurtosis_fn<rayleigh_f_handle, float>;
        out->kurtosis_excess = &kurtosis_excess_fn<rayleigh_f_handle, float>;
        out->mode = &mode_fn<rayleigh_f_handle, float>;
        out->median = &median_fn<rayleigh_f_handle, float>;
        out->entropy = &entropy_fn<rayleigh_f_handle, float>;
        out->free = &free_fn<rayleigh_f_handle>;
        return true;
    } else if (
        n == "saspoint5" || n == "sas_point5" || n == "sas_point_5" ||
        n == "saspointfive" || n == "sas0.5" || n == "sas_alpha_half" ||
        n == "stable_point5" || n == "stable_alpha_half"
    ) {
        double location = 0;
        double scale = 1;
        const char* locKeys[] = { "location", "loc", "mu" };
        const char* scaleKeys[] = { "scale", "sigma", "c" };
        find_param(params, count, locKeys, 3, &location);
        find_param(params, count, scaleKeys, 3, &scale);
        if (!(scale > 0)) return false;
        auto* h = new (std::nothrow) saspoint5_f_handle(static_cast<float>(location), static_cast<float>(scale));
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<saspoint5_f_handle, float>;
        out->logpdf = &logpdf_fn<saspoint5_f_handle, float>;
        out->cdf = &cdf_fn<saspoint5_f_handle, float>;
        out->sf = &sf_fn<saspoint5_f_handle, float>;
        out->hazard = &hazard_fn<saspoint5_f_handle, float>;
        out->chf = &chf_fn<saspoint5_f_handle, float>;
        out->quantile = &quantile_fn<saspoint5_f_handle, float>;
        out->quantile_complement = &quantile_complement_fn<saspoint5_f_handle, float>;
        out->range = &range_fn<saspoint5_f_handle, bs_range_f, float>;
        out->mean = nullptr;
        out->variance = nullptr;
        out->skewness = nullptr;
        out->kurtosis = nullptr;
        out->kurtosis_excess = nullptr;
        out->mode = &mode_fn<saspoint5_f_handle, float>;
        out->median = &median_fn<saspoint5_f_handle, float>;
        out->entropy = &entropy_fn<saspoint5_f_handle, float>;
        out->free = &free_fn<saspoint5_f_handle>;
        return true;
    } else if (n == "skew_normal" || n == "skewnormal" || n == "skew_normal_distribution") {
        double location = 0;
        double scale = 1;
        double shape = 0;
        const char* locKeys[] = { "location", "loc", "mu" };
        const char* scaleKeys[] = { "scale", "sigma", "omega" };
        const char* shapeKeys[] = { "shape", "alpha", "skew" };
        find_param(params, count, locKeys, 3, &location);
        find_param(params, count, scaleKeys, 3, &scale);
        find_param(params, count, shapeKeys, 3, &shape);
        if (!(scale > 0)) return false;
        auto* h = new (std::nothrow) skew_normal_f_handle(static_cast<float>(location), static_cast<float>(scale), static_cast<float>(shape));
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<skew_normal_f_handle, float>;
        out->logpdf = &logpdf_fn<skew_normal_f_handle, float>;
        out->cdf = &cdf_fn<skew_normal_f_handle, float>;
        out->sf = &sf_fn<skew_normal_f_handle, float>;
        out->hazard = &hazard_fn<skew_normal_f_handle, float>;
        out->chf = &chf_fn<skew_normal_f_handle, float>;
        out->quantile = &quantile_fn<skew_normal_f_handle, float>;
        out->quantile_complement = &quantile_complement_fn<skew_normal_f_handle, float>;
        out->range = &range_fn<skew_normal_f_handle, bs_range_f, float>;
        out->mean = &mean_fn<skew_normal_f_handle, float>;
        out->variance = &variance_fn<skew_normal_f_handle, float>;
        out->skewness = &skewness_fn<skew_normal_f_handle, float>;
        out->kurtosis = &kurtosis_fn<skew_normal_f_handle, float>;
        out->kurtosis_excess = &kurtosis_excess_fn<skew_normal_f_handle, float>;
        out->mode = &mode_fn<skew_normal_f_handle, float>;
        out->median = &median_via_quantile_fn<skew_normal_f_handle, float>;
        out->entropy = nullptr;
        out->free = &free_fn<skew_normal_f_handle>;
        return true;
    } else if (n == "triangular" || n == "triangular_distribution" || n == "triangle") {
        double lower = 0;
        double upper = 0;
        double mode = 0;
        const char* lowerKeys[] = { "lower", "min", "minimum", "a" };
        const char* upperKeys[] = { "upper", "max", "maximum", "b" };
        const char* modeKeys[] = { "mode", "peak", "c" };
        if (!find_param(params, count, lowerKeys, 4, &lower)) return false;
        if (!find_param(params, count, upperKeys, 4, &upper)) return false;
        if (!find_param(params, count, modeKeys, 3, &mode)) return false;
        if (!(upper > lower)) return false;
        if (!(mode >= lower && mode <= upper)) return false;
        auto* h = new (std::nothrow) triangular_f_handle(static_cast<float>(lower), static_cast<float>(mode), static_cast<float>(upper));
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<triangular_f_handle, float>;
        out->logpdf = &logpdf_fn<triangular_f_handle, float>;
        out->cdf = &cdf_fn<triangular_f_handle, float>;
        out->sf = &sf_fn<triangular_f_handle, float>;
        out->hazard = &hazard_fn<triangular_f_handle, float>;
        out->chf = &chf_fn<triangular_f_handle, float>;
        out->quantile = &quantile_fn<triangular_f_handle, float>;
        out->quantile_complement = &quantile_complement_fn<triangular_f_handle, float>;
        out->range = &range_fn<triangular_f_handle, bs_range_f, float>;
        out->mean = &mean_fn<triangular_f_handle, float>;
        out->variance = &variance_fn<triangular_f_handle, float>;
        out->skewness = &skewness_fn<triangular_f_handle, float>;
        out->kurtosis = &kurtosis_fn<triangular_f_handle, float>;
        out->kurtosis_excess = &kurtosis_excess_fn<triangular_f_handle, float>;
        out->mode = &mode_fn<triangular_f_handle, float>;
        out->median = &median_fn<triangular_f_handle, float>;
        out->entropy = &entropy_fn<triangular_f_handle, float>;
        out->free = &free_fn<triangular_f_handle>;
        return true;
    } else if (
        n == "uniform" || n == "uniform_distribution" || n == "uniform_real" ||
        n == "rectangular" || n == "rectangular_distribution"
    ) {
        double lower = 0;
        double upper = 1;
        const char* lowerKeys[] = { "lower", "min", "minimum", "a" };
        const char* upperKeys[] = { "upper", "max", "maximum", "b" };
        find_param(params, count, lowerKeys, 4, &lower);
        find_param(params, count, upperKeys, 4, &upper);
        if (!(upper > lower)) return false;
        auto* h = new (std::nothrow) uniform_f_handle(static_cast<float>(lower), static_cast<float>(upper));
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<uniform_f_handle, float>;
        out->logpdf = &logpdf_fn<uniform_f_handle, float>;
        out->cdf = &cdf_fn<uniform_f_handle, float>;
        out->sf = &sf_fn<uniform_f_handle, float>;
        out->hazard = &hazard_fn<uniform_f_handle, float>;
        out->chf = &chf_fn<uniform_f_handle, float>;
        out->quantile = &quantile_fn<uniform_f_handle, float>;
        out->quantile_complement = &quantile_complement_fn<uniform_f_handle, float>;
        out->range = &range_fn<uniform_f_handle, bs_range_f, float>;
        out->mean = &mean_fn<uniform_f_handle, float>;
        out->variance = &variance_fn<uniform_f_handle, float>;
        out->skewness = &skewness_fn<uniform_f_handle, float>;
        out->kurtosis = &kurtosis_fn<uniform_f_handle, float>;
        out->kurtosis_excess = &kurtosis_excess_fn<uniform_f_handle, float>;
        out->mode = &mode_fn<uniform_f_handle, float>;
        out->median = &median_fn<uniform_f_handle, float>;
        out->entropy = &entropy_fn<uniform_f_handle, float>;
        out->free = &free_fn<uniform_f_handle>;
        return true;
    } else if (n == "weibull" || n == "weibull_distribution") {
        double shape = 0;
        double scale = 1;
        const char* shapeKeys[] = { "shape", "k", "alpha" };
        const char* scaleKeys[] = { "scale", "lambda", "beta" };
        if (!find_param(params, count, shapeKeys, 3, &shape)) return false;
        find_param(params, count, scaleKeys, 3, &scale);
        if (!(shape > 0) || !(scale > 0)) return false;
        auto* h = new (std::nothrow) weibull_f_handle(static_cast<float>(shape), static_cast<float>(scale));
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<weibull_f_handle, float>;
        out->logpdf = &logpdf_fn<weibull_f_handle, float>;
        out->cdf = &cdf_fn<weibull_f_handle, float>;
        out->sf = &sf_fn<weibull_f_handle, float>;
        out->hazard = &hazard_fn<weibull_f_handle, float>;
        out->chf = &chf_fn<weibull_f_handle, float>;
        out->quantile = &quantile_fn<weibull_f_handle, float>;
        out->quantile_complement = &quantile_complement_fn<weibull_f_handle, float>;
        out->range = &range_fn<weibull_f_handle, bs_range_f, float>;
        out->mean = &mean_fn<weibull_f_handle, float>;
        out->variance = &variance_fn<weibull_f_handle, float>;
        out->skewness = &skewness_fn<weibull_f_handle, float>;
        out->kurtosis = &kurtosis_fn<weibull_f_handle, float>;
        out->kurtosis_excess = &kurtosis_excess_fn<weibull_f_handle, float>;
        out->mode = &mode_fn<weibull_f_handle, float>;
        out->median = &median_fn<weibull_f_handle, float>;
        out->entropy = &entropy_fn<weibull_f_handle, float>;
        out->free = &free_fn<weibull_f_handle>;
        return true;
    } else if (n == "kolmogorov_smirnov" || n == "kolmogorov-smirnov" || n == "kolmogorovsmirnov" || n == "kolmogorov_smirnov_distribution" || n == "ks" || n == "ks_distribution") {
        float observations = 0;
        const char* nKeys[] = { "n", "sample_count", "samplecount", "sample_size", "samples", "observations" };
        if (!find_param(params, count, nKeys, 6, &observations)) return false;
        auto* h = new (std::nothrow) kolmogorov_smirnov_f_handle(observations);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<kolmogorov_smirnov_f_handle, float>;
        out->logpdf = &logpdf_fn<kolmogorov_smirnov_f_handle, float>;
        out->cdf = &cdf_fn<kolmogorov_smirnov_f_handle, float>;
        out->sf = &sf_fn<kolmogorov_smirnov_f_handle, float>;
        out->hazard = &hazard_fn<kolmogorov_smirnov_f_handle, float>;
        out->chf = &chf_fn<kolmogorov_smirnov_f_handle, float>;
        out->quantile = &quantile_fn<kolmogorov_smirnov_f_handle, float>;
        out->quantile_complement = &quantile_complement_fn<kolmogorov_smirnov_f_handle, float>;
        out->range = &range_fn<kolmogorov_smirnov_f_handle, bs_range_f, float>;
        out->mean = &mean_fn<kolmogorov_smirnov_f_handle, float>;
        out->variance = &variance_fn<kolmogorov_smirnov_f_handle, float>;
        out->skewness = &skewness_fn<kolmogorov_smirnov_f_handle, float>;
        out->kurtosis = &kurtosis_fn<kolmogorov_smirnov_f_handle, float>;
        out->kurtosis_excess = &kurtosis_excess_fn<kolmogorov_smirnov_f_handle, float>;
        out->mode = &mode_fn<kolmogorov_smirnov_f_handle, float>;
        out->median = &median_via_quantile_fn<kolmogorov_smirnov_f_handle, float>;
        out->entropy = nullptr;
        out->free = &free_fn<kolmogorov_smirnov_f_handle>;
        return true;
    } else if (n == "landau" || n == "landau_distribution") {
        float loc = 0, scale = 1;
        const char* locKeys[] = { "location", "loc", "mu" };
        const char* scaleKeys[] = { "scale", "c", "sigma" };
        find_param(params, count, locKeys, 3, &loc);
        find_param(params, count, scaleKeys, 3, &scale);
        auto* h = new (std::nothrow) landau_f_handle(loc, scale);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<landau_f_handle, float>;
        out->logpdf = &logpdf_fn<landau_f_handle, float>;
        out->cdf = &cdf_fn<landau_f_handle, float>;
        out->sf = &sf_fn<landau_f_handle, float>;
        out->hazard = &hazard_fn<landau_f_handle, float>;
        out->chf = &chf_fn<landau_f_handle, float>;
        out->quantile = &quantile_fn<landau_f_handle, float>;
        out->quantile_complement = &quantile_complement_fn<landau_f_handle, float>;
        out->range = &range_fn<landau_f_handle, bs_range_f, float>;
        out->mean = nullptr;
        out->variance = nullptr;
        out->skewness = nullptr;
        out->kurtosis = nullptr;
        out->kurtosis_excess = nullptr;
        out->mode = &mode_fn<landau_f_handle, float>;
        out->median = &median_fn<landau_f_handle, float>;
        out->entropy = &entropy_fn<landau_f_handle, float>;
        out->free = &free_fn<landau_f_handle>;
        return true;
    } else if (n == "laplace" || n == "laplace_distribution" || n == "double_exponential" || n == "doubleexponential") {
        float loc = 0, scale = 1;
        const char* locKeys[] = { "location", "loc", "mu", "mean" };
        const char* scaleKeys[] = { "scale", "diversity", "b" };
        find_param(params, count, locKeys, 4, &loc);
        find_param(params, count, scaleKeys, 3, &scale);
        auto* h = new (std::nothrow) laplace_f_handle(loc, scale);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<laplace_f_handle, float>;
        out->logpdf = &logpdf_fn<laplace_f_handle, float>;
        out->cdf = &cdf_fn<laplace_f_handle, float>;
        out->sf = &sf_fn<laplace_f_handle, float>;
        out->hazard = &hazard_fn<laplace_f_handle, float>;
        out->chf = &chf_fn<laplace_f_handle, float>;
        out->quantile = &quantile_fn<laplace_f_handle, float>;
        out->quantile_complement = &quantile_complement_fn<laplace_f_handle, float>;
        out->range = &range_fn<laplace_f_handle, bs_range_f, float>;
        out->mean = &mean_fn<laplace_f_handle, float>;
        out->variance = &variance_fn<laplace_f_handle, float>;
        out->skewness = &skewness_fn<laplace_f_handle, float>;
        out->kurtosis = &kurtosis_fn<laplace_f_handle, float>;
        out->kurtosis_excess = &kurtosis_excess_fn<laplace_f_handle, float>;
        out->mode = &mode_fn<laplace_f_handle, float>;
        out->median = &median_fn<laplace_f_handle, float>;
        out->entropy = &entropy_fn<laplace_f_handle, float>;
        out->free = &free_fn<laplace_f_handle>;
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
    } else if (n == "non_central_beta" || n == "non_central_beta_distribution") {
        long double a = 0, b = 0, l = 0;
        const char* aKeys[] = { "alpha", "a", "p", "shape1" };
        const char* bKeys[] = { "beta", "b", "q", "shape2" };
        const char* lKeys[] = {"lambda", "l"};
        if (!find_param(params, count, aKeys, 4, &a)) return false;
        if (!find_param(params, count, bKeys, 4, &b)) return false;
        if (!find_param(params, count, lKeys, 2, &l)) return false;
        auto* h = new (std::nothrow) non_central_beta_l_handle(a, b, l);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<non_central_beta_l_handle, long double>;
        out->logpdf = &logpdf_fn<non_central_beta_l_handle, long double>;
        out->cdf = &cdf_fn<non_central_beta_l_handle, long double>;
        out->sf = &sf_fn<non_central_beta_l_handle, long double>;
        out->hazard = &hazard_fn<non_central_beta_l_handle, long double>;
        out->chf = &chf_fn<non_central_beta_l_handle, long double>;
        out->quantile = &quantile_fn<non_central_beta_l_handle, long double>;
        out->quantile_complement = &quantile_complement_fn<non_central_beta_l_handle, long double>;
        out->range = &range_fn<non_central_beta_l_handle, bs_range_l, long double>;
        out->mean = &mean_fn<non_central_beta_l_handle, long double>;
        out->variance = &variance_fn<non_central_beta_l_handle, long double>;
        out->skewness = nullptr;
        out->kurtosis = nullptr;
        out->kurtosis_excess =nullptr;
        out->mode = &mode_fn<non_central_beta_l_handle, long double>;
        out->median = &median_fn<non_central_beta_l_handle, long double>;
        out->entropy = nullptr; // not provided in Boost
        out->free = &free_fn<non_central_beta_l_handle>;
        return true;
    } else if (n == "non_central_chi_squared" || n == "noncentral_chi_squared" || n == "noncentralchisquared" || n == "non_central_chi2" || n == "noncentral_chi2" || n == "nc_chi_squared" || n == "ncchisquared") {
        long double v = 0;
        long double lambda = 0;
        const char* dfKeys[] = { "df", "nu", "degreesoffreedom" };
        const char* lambdaKeys[] = { "lambda", "noncentrality", "delta", "nc" };
        if (!find_param(params, count, dfKeys, 3, &v)) return false;
        if (!find_param(params, count, lambdaKeys, 4, &lambda)) return false;
        if (!(v > 0) || !(lambda >= 0)) return false;
        auto* h = new (std::nothrow) non_central_chi_squared_l_handle(v, lambda);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<non_central_chi_squared_l_handle, long double>;
        out->logpdf = &logpdf_fn<non_central_chi_squared_l_handle, long double>;
        out->cdf = &cdf_fn<non_central_chi_squared_l_handle, long double>;
        out->sf = &sf_fn<non_central_chi_squared_l_handle, long double>;
        out->hazard = &hazard_fn<non_central_chi_squared_l_handle, long double>;
        out->chf = &chf_fn<non_central_chi_squared_l_handle, long double>;
        out->quantile = &quantile_fn<non_central_chi_squared_l_handle, long double>;
        out->quantile_complement = &quantile_complement_fn<non_central_chi_squared_l_handle, long double>;
        out->range = &range_fn<non_central_chi_squared_l_handle, bs_range_l, long double>;
        out->mean = &mean_fn<non_central_chi_squared_l_handle, long double>;
        out->variance = &variance_fn<non_central_chi_squared_l_handle, long double>;
        out->skewness = &skewness_fn<non_central_chi_squared_l_handle, long double>;
        out->kurtosis = &kurtosis_fn<non_central_chi_squared_l_handle, long double>;
        out->kurtosis_excess = &kurtosis_excess_fn<non_central_chi_squared_l_handle, long double>;
        out->mode = &mode_fn<non_central_chi_squared_l_handle, long double>;
        out->median = &median_fn<non_central_chi_squared_l_handle, long double>;
        out->entropy = nullptr;
        out->free = &free_fn<non_central_chi_squared_l_handle>;
        return true;
    } else if (n == "non_central_f" || n == "noncentral_f" || n == "noncentralf" || n == "non_central_f_ratio" || n == "noncentral_f_ratio" || n == "nc_f" || n == "ncf") {
        long double df1 = 0;
        long double df2 = 0;
        long double lambda = 0;
        const char* df1Keys[] = { "df1", "d1", "m", "degreesoffreedom1" };
        const char* df2Keys[] = { "df2", "d2", "n", "degreesoffreedom2" };
        const char* lambdaKeys[] = { "lambda", "noncentrality", "delta", "nc" };
        if (!find_param(params, count, df1Keys, 4, &df1)) return false;
        if (!find_param(params, count, df2Keys, 4, &df2)) return false;
        if (!find_param(params, count, lambdaKeys, 4, &lambda)) return false;
        if (!(df1 > 0) || !(df2 > 0) || !(lambda >= 0)) return false;
        auto* h = new (std::nothrow) non_central_f_l_handle(df1, df2, lambda);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<non_central_f_l_handle, long double>;
        out->logpdf = &logpdf_fn<non_central_f_l_handle, long double>;
        out->cdf = &cdf_fn<non_central_f_l_handle, long double>;
        out->sf = &sf_fn<non_central_f_l_handle, long double>;
        out->hazard = &hazard_fn<non_central_f_l_handle, long double>;
        out->chf = &chf_fn<non_central_f_l_handle, long double>;
        out->quantile = &quantile_fn<non_central_f_l_handle, long double>;
        out->quantile_complement = &quantile_complement_fn<non_central_f_l_handle, long double>;
        out->range = &range_fn<non_central_f_l_handle, bs_range_l, long double>;
        out->mean = &mean_fn<non_central_f_l_handle, long double>;
        out->variance = &variance_fn<non_central_f_l_handle, long double>;
        out->skewness = &skewness_fn<non_central_f_l_handle, long double>;
        out->kurtosis = &kurtosis_fn<non_central_f_l_handle, long double>;
        out->kurtosis_excess = &kurtosis_excess_fn<non_central_f_l_handle, long double>;
        out->mode = &mode_fn<non_central_f_l_handle, long double>;
        out->median = &median_fn<non_central_f_l_handle, long double>;
        out->entropy = nullptr;
        out->free = &free_fn<non_central_f_l_handle>;
        return true;
    } else if (n == "non_central_t" || n == "noncentral_t" || n == "noncentralt" || n == "non_central_student_t" || n == "noncentral_student_t" || n == "noncentralstudentt" || n == "nc_t") {
        long double v = 0;
        long double lambda = 0;
        const char* dfKeys[] = { "df", "nu", "degreesoffreedom" };
        const char* lambdaKeys[] = { "lambda", "noncentrality", "delta", "nc" };
        if (!find_param(params, count, dfKeys, 3, &v)) return false;
        if (!find_param(params, count, lambdaKeys, 4, &lambda)) return false;
        if (!(v > 0)) return false;
        if (!std::isfinite(lambda)) return false;
        auto* h = new (std::nothrow) non_central_t_l_handle(v, lambda);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<non_central_t_l_handle, long double>;
        out->logpdf = &logpdf_fn<non_central_t_l_handle, long double>;
        out->cdf = &cdf_fn<non_central_t_l_handle, long double>;
        out->sf = &sf_fn<non_central_t_l_handle, long double>;
        out->hazard = &hazard_fn<non_central_t_l_handle, long double>;
        out->chf = &chf_fn<non_central_t_l_handle, long double>;
        out->quantile = &quantile_fn<non_central_t_l_handle, long double>;
        out->quantile_complement = &quantile_complement_fn<non_central_t_l_handle, long double>;
        out->range = &range_fn<non_central_t_l_handle, bs_range_l, long double>;
        out->mean = &mean_fn<non_central_t_l_handle, long double>;
        out->variance = &variance_fn<non_central_t_l_handle, long double>;
        out->skewness = &skewness_fn<non_central_t_l_handle, long double>;
        out->kurtosis = &kurtosis_fn<non_central_t_l_handle, long double>;
        out->kurtosis_excess = &kurtosis_excess_fn<non_central_t_l_handle, long double>;
        out->mode = &mode_fn<non_central_t_l_handle, long double>;
        out->median = &median_fn<non_central_t_l_handle, long double>;
        out->entropy = nullptr;
        out->free = &free_fn<non_central_t_l_handle>;
        return true;
    } else if (n == "pareto" || n == "pareto_distribution") {
        long double scale = 0;
        long double shape = 0;
        const char* scaleKeys[] = { "scale", "xm", "minimum", "lower", "x0" };
        const char* shapeKeys[] = { "shape", "alpha" };
        if (!find_param(params, count, scaleKeys, 5, &scale)) return false;
        if (!find_param(params, count, shapeKeys, 2, &shape)) return false;
        if (!(scale > 0) || !(shape > 0)) return false;
        auto* h = new (std::nothrow) pareto_l_handle(scale, shape);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<pareto_l_handle, long double>;
        out->logpdf = &logpdf_fn<pareto_l_handle, long double>;
        out->cdf = &cdf_fn<pareto_l_handle, long double>;
        out->sf = &sf_fn<pareto_l_handle, long double>;
        out->hazard = &hazard_fn<pareto_l_handle, long double>;
        out->chf = &chf_fn<pareto_l_handle, long double>;
        out->quantile = &quantile_fn<pareto_l_handle, long double>;
        out->quantile_complement = &quantile_complement_fn<pareto_l_handle, long double>;
        out->range = &range_fn<pareto_l_handle, bs_range_l, long double>;
        out->mean = &mean_fn<pareto_l_handle, long double>;
        out->variance = &variance_fn<pareto_l_handle, long double>;
        out->skewness = &skewness_fn<pareto_l_handle, long double>;
        out->kurtosis = &kurtosis_fn<pareto_l_handle, long double>;
        out->kurtosis_excess = &kurtosis_excess_fn<pareto_l_handle, long double>;
        out->mode = &mode_fn<pareto_l_handle, long double>;
        out->median = &median_fn<pareto_l_handle, long double>;
        out->entropy = nullptr;
        out->free = &free_fn<pareto_l_handle>;
        return true;
    } else if (n == "poisson" || n == "poisson_distribution") {
        long double mean = 0;
        const char* meanKeys[] = { "mean", "lambda", "mu" };
        if (!find_param(params, count, meanKeys, 3, &mean)) return false;
        if (!(mean >= 0)) return false;
        auto* h = new (std::nothrow) poisson_l_handle(mean);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<poisson_l_handle, long double>;
        out->logpdf = &logpdf_fn<poisson_l_handle, long double>;
        out->cdf = &cdf_fn<poisson_l_handle, long double>;
        out->sf = &sf_fn<poisson_l_handle, long double>;
        out->hazard = &hazard_fn<poisson_l_handle, long double>;
        out->chf = &chf_fn<poisson_l_handle, long double>;
        out->quantile = &quantile_fn<poisson_l_handle, long double>;
        out->quantile_complement = &quantile_complement_fn<poisson_l_handle, long double>;
        out->range = &range_fn<poisson_l_handle, bs_range_l, long double>;
        out->mean = &mean_fn<poisson_l_handle, long double>;
        out->variance = &variance_fn<poisson_l_handle, long double>;
        out->skewness = &skewness_fn<poisson_l_handle, long double>;
        out->kurtosis = &kurtosis_fn<poisson_l_handle, long double>;
        out->kurtosis_excess = &kurtosis_excess_fn<poisson_l_handle, long double>;
        out->mode = &mode_fn<poisson_l_handle, long double>;
        out->median = &median_fn<poisson_l_handle, long double>;
        out->entropy = nullptr;
        out->free = &free_fn<poisson_l_handle>;
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
    } else if (n == "inverse_chi_squared" || n == "inversechisquared" || n == "inv_chi_squared" || n == "invchisquared" || n == "inverse_chi2" || n == "inv_chi2") {
        long double df = 0;
        const char* dfKeys[] = { "df", "nu", "degreesoffreedom", "v" };
        if (!find_param(params, count, dfKeys, 4, &df)) return false;
        long double scale = 0;
        const char* scaleKeys[] = { "scale", "sigma2", "xi" };
        bool hasScale = find_param(params, count, scaleKeys, 3, &scale);
        if (!hasScale) {
            if (!(df > 0)) return false;
            scale = 1 / df;
        }
        auto* h = new (std::nothrow) inverse_chi_squared_l_handle(df, scale);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<inverse_chi_squared_l_handle, long double>;
        out->logpdf = &logpdf_fn<inverse_chi_squared_l_handle, long double>;
        out->cdf = &cdf_fn<inverse_chi_squared_l_handle, long double>;
        out->sf = &sf_fn<inverse_chi_squared_l_handle, long double>;
        out->hazard = &hazard_fn<inverse_chi_squared_l_handle, long double>;
        out->chf = &chf_fn<inverse_chi_squared_l_handle, long double>;
        out->quantile = &quantile_fn<inverse_chi_squared_l_handle, long double>;
        out->quantile_complement = &quantile_complement_fn<inverse_chi_squared_l_handle, long double>;
        out->range = &range_fn<inverse_chi_squared_l_handle, bs_range_l, long double>;
        out->mean = &mean_fn<inverse_chi_squared_l_handle, long double>;
        out->variance = &variance_fn<inverse_chi_squared_l_handle, long double>;
        out->skewness = &skewness_fn<inverse_chi_squared_l_handle, long double>;
        out->kurtosis = &kurtosis_fn<inverse_chi_squared_l_handle, long double>;
        out->kurtosis_excess = &kurtosis_excess_fn<inverse_chi_squared_l_handle, long double>;
        out->mode = &mode_fn<inverse_chi_squared_l_handle, long double>;
        out->median = &median_fn<inverse_chi_squared_l_handle, long double>;
        out->entropy = nullptr;
        out->free = &free_fn<inverse_chi_squared_l_handle>;
        return true;
    } else if (n == "non_central_chi_squared" || n == "noncentral_chi_squared" || n == "noncentralchisquared" || n == "non_central_chi2" || n == "noncentral_chi2" || n == "nc_chi_squared" || n == "ncchisquared") {
        long double v = 0;
        long double lambda = 0;
        const char* dfKeys[] = { "df", "nu", "degreesoffreedom" };
        const char* lambdaKeys[] = { "lambda", "noncentrality", "delta", "nc" };
        if (!find_param(params, count, dfKeys, 3, &v)) return false;
        if (!find_param(params, count, lambdaKeys, 4, &lambda)) return false;
        if (!(v > 0)) return false;
        if (!(lambda >= 0)) return false;
        auto* h = new (std::nothrow) non_central_chi_squared_l_handle(v, lambda);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<non_central_chi_squared_l_handle, long double>;
        out->logpdf = &logpdf_fn<non_central_chi_squared_l_handle, long double>;
        out->cdf = &cdf_fn<non_central_chi_squared_l_handle, long double>;
        out->sf = &sf_fn<non_central_chi_squared_l_handle, long double>;
        out->hazard = &hazard_fn<non_central_chi_squared_l_handle, long double>;
        out->chf = &chf_fn<non_central_chi_squared_l_handle, long double>;
        out->quantile = &quantile_fn<non_central_chi_squared_l_handle, long double>;
        out->quantile_complement = &quantile_complement_fn<non_central_chi_squared_l_handle, long double>;
        out->range = &range_fn<non_central_chi_squared_l_handle, bs_range_l, long double>;
        out->mean = &mean_fn<non_central_chi_squared_l_handle, long double>;
        out->variance = &variance_fn<non_central_chi_squared_l_handle, long double>;
        out->skewness = &skewness_fn<non_central_chi_squared_l_handle, long double>;
        out->kurtosis = &kurtosis_fn<non_central_chi_squared_l_handle, long double>;
        out->kurtosis_excess = &kurtosis_excess_fn<non_central_chi_squared_l_handle, long double>;
        out->mode = &mode_fn<non_central_chi_squared_l_handle, long double>;
        out->median = &median_fn<non_central_chi_squared_l_handle, long double>;
        out->entropy = nullptr;
        out->free = &free_fn<non_central_chi_squared_l_handle>;
        return true;
    } else if (n == "inverse_gamma" || n == "inversegamma" || n == "inv_gamma" || n == "invgamma") {
        long double shape = 0;
        const char* shapeKeys[] = { "shape", "alpha", "k" };
        if (!find_param(params, count, shapeKeys, 3, &shape)) return false;
        long double scale = 1;
        const char* scaleKeys[] = { "scale", "theta", "beta" };
        find_param(params, count, scaleKeys, 3, &scale);
        auto* h = new (std::nothrow) inverse_gamma_l_handle(shape, scale);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<inverse_gamma_l_handle, long double>;
        out->logpdf = &logpdf_fn<inverse_gamma_l_handle, long double>;
        out->cdf = &cdf_fn<inverse_gamma_l_handle, long double>;
        out->sf = &sf_fn<inverse_gamma_l_handle, long double>;
        out->hazard = &hazard_fn<inverse_gamma_l_handle, long double>;
        out->chf = &chf_fn<inverse_gamma_l_handle, long double>;
        out->quantile = &quantile_fn<inverse_gamma_l_handle, long double>;
        out->quantile_complement = &quantile_complement_fn<inverse_gamma_l_handle, long double>;
        out->range = &range_fn<inverse_gamma_l_handle, bs_range_l, long double>;
        out->mean = &mean_fn<inverse_gamma_l_handle, long double>;
        out->variance = &variance_fn<inverse_gamma_l_handle, long double>;
        out->skewness = &skewness_fn<inverse_gamma_l_handle, long double>;
        out->kurtosis = &kurtosis_fn<inverse_gamma_l_handle, long double>;
        out->kurtosis_excess = &kurtosis_excess_fn<inverse_gamma_l_handle, long double>;
        out->mode = &mode_fn<inverse_gamma_l_handle, long double>;
        out->median = &median_fn<inverse_gamma_l_handle, long double>;
        out->entropy = nullptr;
        out->free = &free_fn<inverse_gamma_l_handle>;
        return true;
    } else if (n == "inverse_gaussian" || n == "inversegaussian" || n == "inverse_normal" || n == "inversenormal" || n == "wald") {
        long double mean = 0;
        const char* meanKeys[] = { "mean", "mu", "location" };
        if (!find_param(params, count, meanKeys, 3, &mean)) return false;
        long double scale = 1;
        const char* scaleKeys[] = { "scale", "lambda", "shape" };
        find_param(params, count, scaleKeys, 3, &scale);
        auto* h = new (std::nothrow) inverse_gaussian_l_handle(mean, scale);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<inverse_gaussian_l_handle, long double>;
        out->logpdf = &logpdf_fn<inverse_gaussian_l_handle, long double>;
        out->cdf = &cdf_fn<inverse_gaussian_l_handle, long double>;
        out->sf = &sf_fn<inverse_gaussian_l_handle, long double>;
        out->hazard = &hazard_fn<inverse_gaussian_l_handle, long double>;
        out->chf = &chf_fn<inverse_gaussian_l_handle, long double>;
        out->quantile = &quantile_fn<inverse_gaussian_l_handle, long double>;
        out->quantile_complement = &quantile_complement_fn<inverse_gaussian_l_handle, long double>;
        out->range = &range_fn<inverse_gaussian_l_handle, bs_range_l, long double>;
        out->mean = &mean_fn<inverse_gaussian_l_handle, long double>;
        out->variance = &variance_fn<inverse_gaussian_l_handle, long double>;
        out->skewness = &skewness_fn<inverse_gaussian_l_handle, long double>;
        out->kurtosis = &kurtosis_fn<inverse_gaussian_l_handle, long double>;
        out->kurtosis_excess = &kurtosis_excess_fn<inverse_gaussian_l_handle, long double>;
        out->mode = &mode_fn<inverse_gaussian_l_handle, long double>;
        out->median = &median_fn<inverse_gaussian_l_handle, long double>;
        out->entropy = nullptr;
        out->free = &free_fn<inverse_gaussian_l_handle>;
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
    } else if (n == "negative_binomial" || n == "negativebinomial" || n == "neg_binomial" || n == "negative_binomial_distribution" || n == "nbinom") {
        long double successes = 0, p = 0;
        const char* rKeys[] = { "r", "successes", "target", "count" };
        const char* pKeys[] = { "p", "prob", "probability", "success" };
        if (!find_param(params, count, rKeys, 4, &successes)) return false;
        if (!find_param(params, count, pKeys, 4, &p)) return false;
        if (!(successes > 0)) return false;
        if (!(p > 0 && p <= 1)) return false;
        auto* h = new (std::nothrow) negative_binomial_l_handle(successes, p);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<negative_binomial_l_handle, long double>;
        out->logpdf = &logpdf_fn<negative_binomial_l_handle, long double>;
        out->cdf = &cdf_fn<negative_binomial_l_handle, long double>;
        out->sf = &sf_fn<negative_binomial_l_handle, long double>;
        out->hazard = &hazard_fn<negative_binomial_l_handle, long double>;
        out->chf = &chf_fn<negative_binomial_l_handle, long double>;
        out->quantile = &quantile_fn<negative_binomial_l_handle, long double>;
        out->quantile_complement = &quantile_complement_fn<negative_binomial_l_handle, long double>;
        out->range = &range_fn<negative_binomial_l_handle, bs_range_l, long double>;
        out->mean = &mean_fn<negative_binomial_l_handle, long double>;
        out->variance = &variance_fn<negative_binomial_l_handle, long double>;
        out->skewness = &skewness_fn<negative_binomial_l_handle, long double>;
        out->kurtosis = &kurtosis_fn<negative_binomial_l_handle, long double>;
        out->kurtosis_excess = &kurtosis_excess_fn<negative_binomial_l_handle, long double>;
        out->mode = &mode_fn<negative_binomial_l_handle, long double>;
        out->median = &median_via_quantile_fn<negative_binomial_l_handle, long double>;
        out->entropy = nullptr;
        out->free = &free_fn<negative_binomial_l_handle>;
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
    } else if (n == "hyperexponential" || n == "hyper_exponential" || n == "hyperexp" || n == "hyperexponential_distribution") {
        std::vector<long double> probs;
        std::vector<long double> rates;
        if (!collect_hyperex_parameters(params, count, probs, rates)) return false;
        hyperexponential_l_handle* h = nullptr;
        if (probs.empty()) {
            h = new (std::nothrow) hyperexponential_l_handle(rates);
        } else {
            h = new (std::nothrow) hyperexponential_l_handle(probs, rates);
        }
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<hyperexponential_l_handle, long double>;
        out->logpdf = &logpdf_fn<hyperexponential_l_handle, long double>;
        out->cdf = &cdf_fn<hyperexponential_l_handle, long double>;
        out->sf = &sf_fn<hyperexponential_l_handle, long double>;
        out->hazard = &hazard_fn<hyperexponential_l_handle, long double>;
        out->chf = &chf_fn<hyperexponential_l_handle, long double>;
        out->quantile = &quantile_fn<hyperexponential_l_handle, long double>;
        out->quantile_complement = &quantile_complement_fn<hyperexponential_l_handle, long double>;
        out->range = &range_fn<hyperexponential_l_handle, bs_range_l, long double>;
        out->mean = &mean_fn<hyperexponential_l_handle, long double>;
        out->variance = &variance_fn<hyperexponential_l_handle, long double>;
        out->skewness = &skewness_fn<hyperexponential_l_handle, long double>;
        out->kurtosis = &kurtosis_fn<hyperexponential_l_handle, long double>;
        out->kurtosis_excess = &kurtosis_excess_fn<hyperexponential_l_handle, long double>;
        out->mode = &mode_fn<hyperexponential_l_handle, long double>;
        out->median = &median_fn<hyperexponential_l_handle, long double>;
        out->entropy = nullptr;
        out->free = &free_fn<hyperexponential_l_handle>;
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
    } else if (n == "normal" || n == "normal_distribution" || n == "gauss" || n == "gaussian" || n == "gaussian_distribution" || n == "gauss_distribution") {
        double loc = 0, scale = 1;
        const char* locKeys[] = { "location", "loc", "mu", "mean" };
        const char* scaleKeys[] = { "sd", "standard_deviation", "sigma" };
        find_param(params, count, locKeys, 4, &loc);
        find_param(params, count, scaleKeys, 3, &scale);
        auto* h = new (std::nothrow) normal_l_handle(loc, scale);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<normal_l_handle, long double>;
        out->logpdf = &logpdf_fn<normal_l_handle, long double>;
        out->cdf = &cdf_fn<normal_l_handle, long double>;
        out->sf = &sf_fn<normal_l_handle, long double>;
        out->hazard = &hazard_fn<normal_l_handle, long double>;
        out->chf = &chf_fn<normal_l_handle, long double>;
        out->quantile = &quantile_fn<normal_l_handle, long double>;
        out->quantile_complement = &quantile_complement_fn<normal_l_handle, long double>;
        out->range = &range_fn<normal_l_handle, bs_range_l, long double>;
        out->mean = &mean_fn<normal_l_handle, long double>;
        out->variance = &variance_fn<normal_l_handle, long double>;
        out->skewness = &skewness_fn<normal_l_handle, long double>;
        out->kurtosis = &kurtosis_fn<normal_l_handle, long double>;
        out->kurtosis_excess = &kurtosis_excess_fn<normal_l_handle, long double>;;
        out->mode = &mode_fn<normal_l_handle, long double>;
        out->median = &median_fn<normal_l_handle, long double>;
        out->entropy = &entropy_fn<normal_l_handle, long double>;;
        out->free = &free_fn<normal_l_handle>;
        return true;
    } else if (n == "logistic" || n == "logistic_distribution") {
        long double loc = 0, scale = 1;
        const char* locKeys[] = { "location", "loc", "mu", "median" };
        const char* scaleKeys[] = { "scale", "s", "sigma", "diversity" };
        find_param(params, count, locKeys, 4, &loc);
        find_param(params, count, scaleKeys, 4, &scale);
        if (!(scale > 0)) return false;
        auto* h = new (std::nothrow) logistic_l_handle(loc, scale);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<logistic_l_handle, long double>;
        out->logpdf = &logpdf_fn<logistic_l_handle, long double>;
        out->cdf = &cdf_fn<logistic_l_handle, long double>;
        out->sf = &sf_fn<logistic_l_handle, long double>;
        out->hazard = &hazard_fn<logistic_l_handle, long double>;
        out->chf = &chf_fn<logistic_l_handle, long double>;
        out->quantile = &quantile_fn<logistic_l_handle, long double>;
        out->quantile_complement = &quantile_complement_fn<logistic_l_handle, long double>;
        out->range = &range_fn<logistic_l_handle, bs_range_l, long double>;
        out->mean = &mean_fn<logistic_l_handle, long double>;
        out->variance = &variance_fn<logistic_l_handle, long double>;
        out->skewness = &skewness_fn<logistic_l_handle, long double>;
        out->kurtosis = &kurtosis_fn<logistic_l_handle, long double>;
        out->kurtosis_excess = &kurtosis_excess_fn<logistic_l_handle, long double>;
        out->mode = &mode_fn<logistic_l_handle, long double>;
        out->median = &median_fn<logistic_l_handle, long double>;
        out->entropy = &entropy_fn<logistic_l_handle, long double>;
        out->free = &free_fn<logistic_l_handle>;
        return true;
    } else if (n == "lognormal" || n == "log_normal" || n == "lognormal_distribution") {
        long double mu = 0, sigma = 1;
        const char* muKeys[] = { "location", "loc", "mu", "meanlog" };
        const char* sigmaKeys[] = { "scale", "sigma", "sd", "standard_deviation" };
        find_param(params, count, muKeys, 4, &mu);
        find_param(params, count, sigmaKeys, 4, &sigma);
        if (!(sigma > 0)) return false;
        auto* h = new (std::nothrow) lognormal_l_handle(mu, sigma);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<lognormal_l_handle, long double>;
        out->logpdf = &logpdf_fn<lognormal_l_handle, long double>;
        out->cdf = &cdf_fn<lognormal_l_handle, long double>;
        out->sf = &sf_fn<lognormal_l_handle, long double>;
        out->hazard = &hazard_fn<lognormal_l_handle, long double>;
        out->chf = &chf_fn<lognormal_l_handle, long double>;
        out->quantile = &quantile_fn<lognormal_l_handle, long double>;
        out->quantile_complement = &quantile_complement_fn<lognormal_l_handle, long double>;
        out->range = &range_fn<lognormal_l_handle, bs_range_l, long double>;
        out->mean = &mean_fn<lognormal_l_handle, long double>;
        out->variance = &variance_fn<lognormal_l_handle, long double>;
        out->skewness = &skewness_fn<lognormal_l_handle, long double>;
        out->kurtosis = &kurtosis_fn<lognormal_l_handle, long double>;
        out->kurtosis_excess = &kurtosis_excess_fn<lognormal_l_handle, long double>;
        out->mode = &mode_fn<lognormal_l_handle, long double>;
        out->median = &median_fn<lognormal_l_handle, long double>;
        out->entropy = &entropy_fn<lognormal_l_handle, long double>;
        out->free = &free_fn<lognormal_l_handle>;
        return true;
    } else if (n == "mapairy" || n == "map_airy" || n == "mapairy_distribution") {
        long double loc = 0, scale = 1;
        const char* locKeys[] = { "location", "loc", "mu" };
        const char* scaleKeys[] = { "scale", "c", "sigma" };
        find_param(params, count, locKeys, 3, &loc);
        find_param(params, count, scaleKeys, 3, &scale);
        if (!(scale > 0)) return false;
        auto* h = new (std::nothrow) mapairy_l_handle(loc, scale);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<mapairy_l_handle, long double>;
        out->logpdf = &logpdf_fn<mapairy_l_handle, long double>;
        out->cdf = &cdf_fn<mapairy_l_handle, long double>;
        out->sf = &sf_fn<mapairy_l_handle, long double>;
        out->hazard = &hazard_fn<mapairy_l_handle, long double>;
        out->chf = &chf_fn<mapairy_l_handle, long double>;
        out->quantile = &quantile_fn<mapairy_l_handle, long double>;
        out->quantile_complement = &quantile_complement_fn<mapairy_l_handle, long double>;
        out->range = &range_fn<mapairy_l_handle, bs_range_l, long double>;
        out->mean = &mean_fn<mapairy_l_handle, long double>;
        out->variance = &variance_fn<mapairy_l_handle, long double>;
        out->skewness = nullptr;
        out->kurtosis = nullptr;
        out->kurtosis_excess = nullptr;
        out->mode = &mode_fn<mapairy_l_handle, long double>;
        out->median = &median_fn<mapairy_l_handle, long double>;
        out->entropy = &entropy_fn<mapairy_l_handle, long double>;
        out->free = &free_fn<mapairy_l_handle>;
        return true;
    } else if (n == "rayleigh" || n == "rayleigh_distribution") {
        long double scale = 0;
        const char* scaleKeys[] = { "scale", "sigma", "beta" };
        if (!find_param(params, count, scaleKeys, 3, &scale)) return false;
        if (!(scale > 0)) return false;
        auto* h = new (std::nothrow) rayleigh_l_handle(scale);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<rayleigh_l_handle, long double>;
        out->logpdf = &logpdf_fn<rayleigh_l_handle, long double>;
        out->cdf = &cdf_fn<rayleigh_l_handle, long double>;
        out->sf = &sf_fn<rayleigh_l_handle, long double>;
        out->hazard = &hazard_fn<rayleigh_l_handle, long double>;
        out->chf = &chf_fn<rayleigh_l_handle, long double>;
        out->quantile = &quantile_fn<rayleigh_l_handle, long double>;
        out->quantile_complement = &quantile_complement_fn<rayleigh_l_handle, long double>;
        out->range = &range_fn<rayleigh_l_handle, bs_range_l, long double>;
        out->mean = &mean_fn<rayleigh_l_handle, long double>;
        out->variance = &variance_fn<rayleigh_l_handle, long double>;
        out->skewness = &skewness_fn<rayleigh_l_handle, long double>;
        out->kurtosis = &kurtosis_fn<rayleigh_l_handle, long double>;
        out->kurtosis_excess = &kurtosis_excess_fn<rayleigh_l_handle, long double>;
        out->mode = &mode_fn<rayleigh_l_handle, long double>;
        out->median = &median_fn<rayleigh_l_handle, long double>;
        out->entropy = &entropy_fn<rayleigh_l_handle, long double>;
        out->free = &free_fn<rayleigh_l_handle>;
        return true;
    } else if (
        n == "saspoint5" || n == "sas_point5" || n == "sas_point_5" ||
        n == "saspointfive" || n == "sas0.5" || n == "sas_alpha_half" ||
        n == "stable_point5" || n == "stable_alpha_half"
    ) {
        long double location = 0;
        long double scale = 1;
        const char* locKeys[] = { "location", "loc", "mu" };
        const char* scaleKeys[] = { "scale", "sigma", "c" };
        find_param(params, count, locKeys, 3, &location);
        find_param(params, count, scaleKeys, 3, &scale);
        if (!(scale > 0)) return false;
        auto* h = new (std::nothrow) saspoint5_l_handle(location, scale);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<saspoint5_l_handle, long double>;
        out->logpdf = &logpdf_fn<saspoint5_l_handle, long double>;
        out->cdf = &cdf_fn<saspoint5_l_handle, long double>;
        out->sf = &sf_fn<saspoint5_l_handle, long double>;
        out->hazard = &hazard_fn<saspoint5_l_handle, long double>;
        out->chf = &chf_fn<saspoint5_l_handle, long double>;
        out->quantile = &quantile_fn<saspoint5_l_handle, long double>;
        out->quantile_complement = &quantile_complement_fn<saspoint5_l_handle, long double>;
        out->range = &range_fn<saspoint5_l_handle, bs_range_l, long double>;
        out->mean = nullptr;
        out->variance = nullptr;
        out->skewness = nullptr;
        out->kurtosis = nullptr;
        out->kurtosis_excess = nullptr;
        out->mode = &mode_fn<saspoint5_l_handle, long double>;
        out->median = &median_fn<saspoint5_l_handle, long double>;
        out->entropy = &entropy_fn<saspoint5_l_handle, long double>;
        out->free = &free_fn<saspoint5_l_handle>;
        return true;
    } else if (n == "skew_normal" || n == "skewnormal" || n == "skew_normal_distribution") {
        long double location = 0;
        long double scale = 1;
        long double shape = 0;
        const char* locKeys[] = { "location", "loc", "mu" };
        const char* scaleKeys[] = { "scale", "sigma", "omega" };
        const char* shapeKeys[] = { "shape", "alpha", "skew" };
        find_param(params, count, locKeys, 3, &location);
        find_param(params, count, scaleKeys, 3, &scale);
        find_param(params, count, shapeKeys, 3, &shape);
        if (!(scale > 0)) return false;
        auto* h = new (std::nothrow) skew_normal_l_handle(location, scale, shape);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<skew_normal_l_handle, long double>;
        out->logpdf = &logpdf_fn<skew_normal_l_handle, long double>;
        out->cdf = &cdf_fn<skew_normal_l_handle, long double>;
        out->sf = &sf_fn<skew_normal_l_handle, long double>;
        out->hazard = &hazard_fn<skew_normal_l_handle, long double>;
        out->chf = &chf_fn<skew_normal_l_handle, long double>;
        out->quantile = &quantile_fn<skew_normal_l_handle, long double>;
        out->quantile_complement = &quantile_complement_fn<skew_normal_l_handle, long double>;
        out->range = &range_fn<skew_normal_l_handle, bs_range_l, long double>;
        out->mean = &mean_fn<skew_normal_l_handle, long double>;
        out->variance = &variance_fn<skew_normal_l_handle, long double>;
        out->skewness = &skewness_fn<skew_normal_l_handle, long double>;
        out->kurtosis = &kurtosis_fn<skew_normal_l_handle, long double>;
        out->kurtosis_excess = &kurtosis_excess_fn<skew_normal_l_handle, long double>;
        out->mode = &mode_fn<skew_normal_l_handle, long double>;
        out->median = &median_via_quantile_fn<skew_normal_l_handle, long double>;
        out->entropy = nullptr;
        out->free = &free_fn<skew_normal_l_handle>;
        return true;
    } else if (n == "triangular" || n == "triangular_distribution" || n == "triangle") {
        long double lower = 0;
        long double upper = 0;
        long double mode = 0;
        const char* lowerKeys[] = { "lower", "min", "minimum", "a" };
        const char* upperKeys[] = { "upper", "max", "maximum", "b" };
        const char* modeKeys[] = { "mode", "peak", "c" };
        if (!find_param(params, count, lowerKeys, 4, &lower)) return false;
        if (!find_param(params, count, upperKeys, 4, &upper)) return false;
        if (!find_param(params, count, modeKeys, 3, &mode)) return false;
        if (!(upper > lower)) return false;
        if (!(mode >= lower && mode <= upper)) return false;
        auto* h = new (std::nothrow) triangular_l_handle(lower, mode, upper);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<triangular_l_handle, long double>;
        out->logpdf = &logpdf_fn<triangular_l_handle, long double>;
        out->cdf = &cdf_fn<triangular_l_handle, long double>;
        out->sf = &sf_fn<triangular_l_handle, long double>;
        out->hazard = &hazard_fn<triangular_l_handle, long double>;
        out->chf = &chf_fn<triangular_l_handle, long double>;
        out->quantile = &quantile_fn<triangular_l_handle, long double>;
        out->quantile_complement = &quantile_complement_fn<triangular_l_handle, long double>;
        out->range = &range_fn<triangular_l_handle, bs_range_l, long double>;
        out->mean = &mean_fn<triangular_l_handle, long double>;
        out->variance = &variance_fn<triangular_l_handle, long double>;
        out->skewness = &skewness_fn<triangular_l_handle, long double>;
        out->kurtosis = &kurtosis_fn<triangular_l_handle, long double>;
        out->kurtosis_excess = &kurtosis_excess_fn<triangular_l_handle, long double>;
        out->mode = &mode_fn<triangular_l_handle, long double>;
        out->median = &median_fn<triangular_l_handle, long double>;
        out->entropy = &entropy_fn<triangular_l_handle, long double>;
        out->free = &free_fn<triangular_l_handle>;
        return true;
    } else if (
        n == "uniform" || n == "uniform_distribution" || n == "uniform_real" ||
        n == "rectangular" || n == "rectangular_distribution"
    ) {
        long double lower = 0;
        long double upper = 1;
        const char* lowerKeys[] = { "lower", "min", "minimum", "a" };
        const char* upperKeys[] = { "upper", "max", "maximum", "b" };
        find_param(params, count, lowerKeys, 4, &lower);
        find_param(params, count, upperKeys, 4, &upper);
        if (!(upper > lower)) return false;
        auto* h = new (std::nothrow) uniform_l_handle(lower, upper);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<uniform_l_handle, long double>;
        out->logpdf = &logpdf_fn<uniform_l_handle, long double>;
        out->cdf = &cdf_fn<uniform_l_handle, long double>;
        out->sf = &sf_fn<uniform_l_handle, long double>;
        out->hazard = &hazard_fn<uniform_l_handle, long double>;
        out->chf = &chf_fn<uniform_l_handle, long double>;
        out->quantile = &quantile_fn<uniform_l_handle, long double>;
        out->quantile_complement = &quantile_complement_fn<uniform_l_handle, long double>;
        out->range = &range_fn<uniform_l_handle, bs_range_l, long double>;
        out->mean = &mean_fn<uniform_l_handle, long double>;
        out->variance = &variance_fn<uniform_l_handle, long double>;
        out->skewness = &skewness_fn<uniform_l_handle, long double>;
        out->kurtosis = &kurtosis_fn<uniform_l_handle, long double>;
        out->kurtosis_excess = &kurtosis_excess_fn<uniform_l_handle, long double>;
        out->mode = &mode_fn<uniform_l_handle, long double>;
        out->median = &median_fn<uniform_l_handle, long double>;
        out->entropy = &entropy_fn<uniform_l_handle, long double>;
        out->free = &free_fn<uniform_l_handle>;
        return true;
    } else if (n == "weibull" || n == "weibull_distribution") {
        long double shape = 0;
        long double scale = 1;
        const char* shapeKeys[] = { "shape", "k", "alpha" };
        const char* scaleKeys[] = { "scale", "lambda", "beta" };
        if (!find_param(params, count, shapeKeys, 3, &shape)) return false;
        find_param(params, count, scaleKeys, 3, &scale);
        if (!(shape > 0) || !(scale > 0)) return false;
        auto* h = new (std::nothrow) weibull_l_handle(shape, scale);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<weibull_l_handle, long double>;
        out->logpdf = &logpdf_fn<weibull_l_handle, long double>;
        out->cdf = &cdf_fn<weibull_l_handle, long double>;
        out->sf = &sf_fn<weibull_l_handle, long double>;
        out->hazard = &hazard_fn<weibull_l_handle, long double>;
        out->chf = &chf_fn<weibull_l_handle, long double>;
        out->quantile = &quantile_fn<weibull_l_handle, long double>;
        out->quantile_complement = &quantile_complement_fn<weibull_l_handle, long double>;
        out->range = &range_fn<weibull_l_handle, bs_range_l, long double>;
        out->mean = &mean_fn<weibull_l_handle, long double>;
        out->variance = &variance_fn<weibull_l_handle, long double>;
        out->skewness = &skewness_fn<weibull_l_handle, long double>;
        out->kurtosis = &kurtosis_fn<weibull_l_handle, long double>;
        out->kurtosis_excess = &kurtosis_excess_fn<weibull_l_handle, long double>;
        out->mode = &mode_fn<weibull_l_handle, long double>;
        out->median = &median_fn<weibull_l_handle, long double>;
        out->entropy = &entropy_fn<weibull_l_handle, long double>;
        out->free = &free_fn<weibull_l_handle>;
        return true;
    } else if (n == "kolmogorov_smirnov" || n == "kolmogorov-smirnov" || n == "kolmogorovsmirnov" || n == "kolmogorov_smirnov_distribution" || n == "ks" || n == "ks_distribution") {
        long double observations = 0;
        const char* nKeys[] = { "n", "sample_count", "samplecount", "sample_size", "samples", "observations" };
        if (!find_param(params, count, nKeys, 6, &observations)) return false;
        auto* h = new (std::nothrow) kolmogorov_smirnov_l_handle(observations);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<kolmogorov_smirnov_l_handle, long double>;
        out->logpdf = &logpdf_fn<kolmogorov_smirnov_l_handle, long double>;
        out->cdf = &cdf_fn<kolmogorov_smirnov_l_handle, long double>;
        out->sf = &sf_fn<kolmogorov_smirnov_l_handle, long double>;
        out->hazard = &hazard_fn<kolmogorov_smirnov_l_handle, long double>;
        out->chf = &chf_fn<kolmogorov_smirnov_l_handle, long double>;
        out->quantile = &quantile_fn<kolmogorov_smirnov_l_handle, long double>;
        out->quantile_complement = &quantile_complement_fn<kolmogorov_smirnov_l_handle, long double>;
        out->range = &range_fn<kolmogorov_smirnov_l_handle, bs_range_l, long double>;
        out->mean = &mean_fn<kolmogorov_smirnov_l_handle, long double>;
        out->variance = &variance_fn<kolmogorov_smirnov_l_handle, long double>;
        out->skewness = &skewness_fn<kolmogorov_smirnov_l_handle, long double>;
        out->kurtosis = &kurtosis_fn<kolmogorov_smirnov_l_handle, long double>;
        out->kurtosis_excess = &kurtosis_excess_fn<kolmogorov_smirnov_l_handle, long double>;
        out->mode = &mode_fn<kolmogorov_smirnov_l_handle, long double>;
        out->median = &median_via_quantile_fn<kolmogorov_smirnov_l_handle, long double>;
        out->entropy = nullptr;
        out->free = &free_fn<kolmogorov_smirnov_l_handle>;
        return true;
    } else if (n == "landau" || n == "landau_distribution") {
        long double loc = 0, scale = 1;
        const char* locKeys[] = { "location", "loc", "mu" };
        const char* scaleKeys[] = { "scale", "c", "sigma" };
        find_param(params, count, locKeys, 3, &loc);
        find_param(params, count, scaleKeys, 3, &scale);
        auto* h = new (std::nothrow) landau_l_handle(loc, scale);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<landau_l_handle, long double>;
        out->logpdf = &logpdf_fn<landau_l_handle, long double>;
        out->cdf = &cdf_fn<landau_l_handle, long double>;
        out->sf = &sf_fn<landau_l_handle, long double>;
        out->hazard = &hazard_fn<landau_l_handle, long double>;
        out->chf = &chf_fn<landau_l_handle, long double>;
        out->quantile = &quantile_fn<landau_l_handle, long double>;
        out->quantile_complement = &quantile_complement_fn<landau_l_handle, long double>;
        out->range = &range_fn<landau_l_handle, bs_range_l, long double>;
        out->mean = nullptr;
        out->variance = nullptr;
        out->skewness = nullptr;
        out->kurtosis = nullptr;
        out->kurtosis_excess = nullptr;
        out->mode = &mode_fn<landau_l_handle, long double>;
        out->median = &median_fn<landau_l_handle, long double>;
        out->entropy = &entropy_fn<landau_l_handle, long double>;
        out->free = &free_fn<landau_l_handle>;
        return true;
    } else if (n == "laplace" || n == "laplace_distribution" || n == "double_exponential" || n == "doubleexponential") {
        long double loc = 0, scale = 1;
        const char* locKeys[] = { "location", "loc", "mu", "mean" };
        const char* scaleKeys[] = { "scale", "diversity", "b" };
        find_param(params, count, locKeys, 4, &loc);
        find_param(params, count, scaleKeys, 3, &scale);
        auto* h = new (std::nothrow) laplace_l_handle(loc, scale);
        if (!h) return false;
        out->ctx = h;
        out->pdf = &pdf_fn<laplace_l_handle, long double>;
        out->logpdf = &logpdf_fn<laplace_l_handle, long double>;
        out->cdf = &cdf_fn<laplace_l_handle, long double>;
        out->sf = &sf_fn<laplace_l_handle, long double>;
        out->hazard = &hazard_fn<laplace_l_handle, long double>;
        out->chf = &chf_fn<laplace_l_handle, long double>;
        out->quantile = &quantile_fn<laplace_l_handle, long double>;
        out->quantile_complement = &quantile_complement_fn<laplace_l_handle, long double>;
        out->range = &range_fn<laplace_l_handle, bs_range_l, long double>;
        out->mean = &mean_fn<laplace_l_handle, long double>;
        out->variance = &variance_fn<laplace_l_handle, long double>;
        out->skewness = &skewness_fn<laplace_l_handle, long double>;
        out->kurtosis = &kurtosis_fn<laplace_l_handle, long double>;
        out->kurtosis_excess = &kurtosis_excess_fn<laplace_l_handle, long double>;
        out->mode = &mode_fn<laplace_l_handle, long double>;
        out->median = &median_fn<laplace_l_handle, long double>;
        out->entropy = &entropy_fn<laplace_l_handle, long double>;
        out->free = &free_fn<laplace_l_handle>;
        return true;
    }

        return false;
    } catch (...) {
        return false;
    }
}

} // extern "C"
