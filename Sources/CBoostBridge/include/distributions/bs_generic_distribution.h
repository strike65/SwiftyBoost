//
//  Generic distribution vtable (Boost.Math) — C ABI
//  Provides a uniform handle for continuous probability distributions
//  with function pointers for evaluation and summary statistics.
//
#pragma once

#include <stddef.h> // size_t
#include <stdbool.h> // bool

#ifdef __cplusplus
extern "C" {
#endif

// Reuse range PODs defined for other distributions (avoid redefining types)
// C-friendly range PODs for [lower, upper]
typedef struct { float lower, upper; } bs_range_f;
typedef struct { double lower, upper; } bs_range_d;
typedef struct { long double lower, upper; } bs_range_l;

// Parameter key-value pairs
typedef struct { const char* _Nonnull key; float value; } bs_param_f;
typedef struct { const char* _Nonnull key; double value; } bs_param_d;
typedef struct { const char* _Nonnull key; long double value; } bs_param_l;

// Double-precision generic distribution interface
typedef struct {
    void* _Nonnull ctx; // opaque handle
    // Core functions
    double (* _Nullable pdf)(const void* _Nonnull, double);
    double (* _Nullable logpdf)(const void* _Nonnull, double);
    double (* _Nullable cdf)(const void* _Nonnull, double);
    double (* _Nullable sf)(const void* _Nonnull, double);
    double (* _Nullable hazard)(const void* _Nonnull, double);
    double (* _Nullable chf)(const void* _Nonnull, double);
    double (* _Nullable quantile)(const void* _Nonnull, double);
    double (* _Nullable quantile_complement)(const void* _Nonnull, double);
    // Support and stats
    bs_range_d (* _Nullable range)(const void* _Nonnull);
    double (* _Nullable mean)(const void* _Nonnull);
    double (* _Nullable variance)(const void* _Nonnull);
    double (* _Nullable skewness)(const void* _Nonnull);
    double (* _Nullable kurtosis)(const void* _Nonnull);
    double (* _Nullable kurtosis_excess)(const void* _Nonnull);
    double (* _Nullable mode)(const void* _Nonnull);
    double (* _Nullable median)(const void* _Nonnull);
    double (* _Nullable entropy)(const void* _Nonnull);
    // Lifetime
    void (* _Nullable free)(void* _Nonnull);
} bs_dist_d;

// Single-precision generic distribution interface
typedef struct {
    void* _Nonnull ctx; // opaque handle
    float (* _Nullable pdf)(const void* _Nonnull, float);
    float (* _Nullable logpdf)(const void* _Nonnull, float);
    float (* _Nullable cdf)(const void* _Nonnull, float);
    float (* _Nullable sf)(const void* _Nonnull, float);
    float (* _Nullable hazard)(const void* _Nonnull, float);
    float (* _Nullable chf)(const void* _Nonnull, float);
    float (* _Nullable quantile)(const void* _Nonnull, float);
    float (* _Nullable quantile_complement)(const void* _Nonnull, float);
    bs_range_f (* _Nullable range)(const void* _Nonnull);
    float (* _Nullable mean)(const void* _Nonnull);
    float (* _Nullable variance)(const void* _Nonnull);
    float (* _Nullable skewness)(const void* _Nonnull);
    float (* _Nullable kurtosis)(const void* _Nonnull);
    float (* _Nullable kurtosis_excess)(const void* _Nonnull);
    float (* _Nullable mode)(const void* _Nonnull);
    float (* _Nullable median)(const void* _Nonnull);
    float (* _Nullable entropy)(const void* _Nonnull);
    void (* _Nullable free)(void* _Nonnull);
} bs_dist_f;

// Extended-precision (long double) generic distribution interface
typedef struct {
    void* _Nonnull ctx; // opaque handle
    long double (* _Nullable pdf)(const void* _Nonnull, long double);
    long double (* _Nullable logpdf)(const void* _Nonnull, long double);
    long double (* _Nullable cdf)(const void* _Nonnull, long double);
    long double (* _Nullable sf)(const void* _Nonnull, long double);
    long double (* _Nullable hazard)(const void* _Nonnull, long double);
    long double (* _Nullable chf)(const void* _Nonnull, long double);
    long double (* _Nullable quantile)(const void* _Nonnull, long double);
    long double (* _Nullable quantile_complement)(const void* _Nonnull, long double);
    bs_range_l (* _Nullable range)(const void* _Nonnull);
    long double (* _Nullable mean)(const void* _Nonnull);
    long double (* _Nullable variance)(const void* _Nonnull);
    long double (* _Nullable skewness)(const void* _Nonnull);
    long double (* _Nullable kurtosis)(const void* _Nonnull);
    long double (* _Nullable kurtosis_excess)(const void* _Nonnull);
    long double (* _Nullable mode)(const void* _Nonnull);
    long double (* _Nullable median)(const void* _Nonnull);
    long double (* _Nullable entropy)(const void* _Nonnull);
    void (* _Nullable free)(void* _Nonnull);
} bs_dist_l;

// Factories: construct a generic distribution by name with parameters
// Supported names: "gamma", "student_t" (initial set)
// Returns true on success; false on failure (invalid name/parameters)
bool bs_dist_make_d(const char* _Nonnull name, const bs_param_d* _Nonnull params, size_t count, bs_dist_d* _Nonnull out);
bool bs_dist_make_f(const char* _Nonnull name, const bs_param_f* _Nonnull params, size_t count, bs_dist_f* _Nonnull out);
bool bs_dist_make_l(const char* _Nonnull name, const bs_param_l* _Nonnull params, size_t count, bs_dist_l* _Nonnull out);

#ifdef __cplusplus
}
#endif
