# Unified Distribution Factory (VTable) — Design, Rationale, and Usage

This document explains the unified, vtable‑based distribution factory introduced for SwiftyBoost. It covers the design principles, concrete implementation in the C bridge and Swift, supported use cases, extension points, and trade‑offs compared to per‑distribution Swift wrappers.


## Goals and Motivation

- Reduce Swift boilerplate for each new distribution while retaining Boost.Math’s high‑quality numerics.
- Provide one generic runtime constructor: create a distribution by name with a dictionary of parameters.
- Keep the fast path: call directly into Boost with minimal overhead (one function‑pointer dispatch level).
- Preserve the existing strong, typed Swift wrappers for discoverability and clarity; add a generic path in parallel.


## High‑Level Design

The design introduces a single, generic C ABI for “a distribution of reals” using a vtable (set of function pointers) and a Swift wrapper that conforms to your `DistributionProtocol`.

- C layer (CBoostBridge):
  - A factory function builds the appropriate Boost distribution and returns a handle with function pointers for pdf, cdf, quantiles, hazards, stats, etc.
  - The vtable is uniform across precision variants: double (`bs_dist_d`), float (`bs_dist_f`), and long double (`bs_dist_l`).
  - Parameters are passed as key/value arrays; keys are case‑insensitive with aliases per distribution.

- Swift layer:
  - `Distribution.Dynamic<T>` (where `T` is `Double`, `Float`, or `Float80` on x86) wraps the vtable and conforms to `DistributionProtocol`.
  - Callers can construct distributions by name plus a parameter dictionary.
  - All methods forward to the vtable, mapping NaN/∞ to optionals where appropriate.


## What’s Implemented Today

- Supported distributions via the factory (aliases in backticks):
  - **Gamma** — `gamma`, `gamma_distribution`
    - Params: `shape|k` (required), `scale|theta` (optional, defaults to 1).
    - Typed wrapper: ``Distribution/Gamma``.
  - **Beta** — `beta`, `beta_distribution`
    - Params: `alpha|a|p|shape1` (required), `beta|b|q|shape2` (required).
    - Typed wrapper: ``Distribution/Beta``.
  - **Chi-squared** — `chisquared`, `chi_squared`, `chi2`, `chi-squared`, `chisquare`
    - Params: `df|nu|degreesOfFreedom` (required).
    - Typed wrapper: ``Distribution/ChiSquared``.
  - **Student’s t** — `studentt`, `student_t`, `students_t`, `t`, `t_distribution`
    - Params: `df|nu|degreesOfFreedom` (required).
    - Typed wrapper: ``Distribution/StudentT``.
  - **Fisher’s F** — `fisherf`, `f`, `f_distribution`
    - Params: `df1|d1|m|degreesOfFreedom1` (required), `df2|d2|n|degreesOfFreedom2` (required).
    - Typed wrapper: ``Distribution/FisherF``.
  - **Bernoulli** — `bernoulli`, `bernoulli_distribution`
    - Params: `p|prob|probability|success|theta` (required).
    - Typed wrapper: ``Distribution/Bernoulli`` (entropy fallback computed in Swift).
  - **Binomial** — `binomial`, `binomial_distribution`
    - Params: `n|trials` (required, interpreted as trial count), `p|prob|probability|success` (required).
    - Typed wrapper: ``Distribution/Binomial`` plus planning helpers in Swift.
  - **Cauchy** — `cauchy`, `cauchy_distribution`
    - Params: `location|loc|mu|median|x0` (optional, defaults to 0), `scale|gamma|sigma|b` (required, > 0).
    - Typed wrapper: ``Distribution/Cauchy`` (entropy computed in Swift).
  - **Exponential** — `exponential`, `exponential_distribution`, `exp`
    - Params: `lambda|rate` (required, > 0). `scale` is accepted indirectly by typed wrapper.
    - Typed wrapper: ``Distribution/Exponential``.
  - **Extreme value / Gumbel** — `extremevalue`, `extreme_value`, `gumbel`, `extreme_value_distribution`
    - Params: `location|loc|mu` (optional, defaults to 0), `scale|gamma|sigma|b` (required, > 0).
    - Typed wrapper: ``Distribution/ExtremeValueGumpel`` (entropy fallback provided in Swift).
  - **Arcsine** — `arcsine`, `arcsine_distribution`
    - Params: `minX|min|a|lower` (required), `maxX|max|b|upper` (required).
    - Typed wrapper: ``Distribution/Arcsine`` (entropy fallback provided in Swift).
  - **Geometric** — `geometric`, `geometric_distribution`
    - Params: `p|prob|probability|success|theta` (required).
    - Typed wrapper: ``Distribution/Geometric`` (delegates to Dynamic; entropy currently unavailable in Boost).
  - **Holtsmark** — `holtsmark`, `holtsmark_distribution`
    - Params: `location|loc|mu|median|x0` (optional), `scale|gamma|sigma|b` (required).
    - Typed wrapper: ``Distribution/Holtsmark``.

- Swift typed wrappers delegate to the factory internally, so the dynamic and typed paths share the same Boost-backed semantics and alias handling.


## C ABI: Generic Distribution VTable

- Declaration: `Sources/CBoostBridge/include/distributions/bs_generic_distribution.h`
  - Parameter records (key/value pairs):
    - `bs_param_d { const char* key; double value; }`
    - `bs_param_f { const char* key; float value; }`
    - `bs_param_l { const char* key; long double value; }`
  - Vtables:
    - `bs_dist_d`, `bs_dist_f`, `bs_dist_l` each contain:
      - Core: `pdf`, `logpdf`, `cdf`, `sf`, `hazard`, `chf`, `quantile`, `quantile_complement`
      - Support/stats: `range`, `mean`, `variance`, `skewness`, `kurtosis`, `kurtosis_excess`, `mode`, `median`, `entropy`
      - Lifetime: `free(void*)`
  - Nullability annotations:
    - `ctx` is marked `_Nonnull` (always valid for a successfully constructed handle).
    - Function pointers are marked `_Nullable` because some entries (like `entropy`) may not be provided for certain distributions.
    - The first argument of each function pointer is annotated `const void* _Nonnull` to reflect that the context is required.
  - Factories:
    - `bool bs_dist_make_d(const char* name, const bs_param_d* params, size_t count, bs_dist_d* out);`
    - `bool bs_dist_make_f(const char* name, const bs_param_f* params, size_t count, bs_dist_f* out);`
    - `bool bs_dist_make_l(const char* name, const bs_param_l* params, size_t count, bs_dist_l* out);`

- Implementation: `Sources/CBoostBridge/impl/bs_generic_distribution.hxx`
  - Parses `name` case‑insensitively.
  - Pulls params by alias sets (e.g., `shape`/`k`; `scale`/`theta`; `df`/`nu`/`degreesOfFreedom`).
  - Allocates a small handle with the Boost distribution instance as `ctx` and fills vtable function pointers.
  - Uses small inline templates to implement thunks (pdf/cdf/quantile/…); hazards are derived via `pdf` and `cdf(complement)`.
  - Error and exception policy: numeric exceptions from Boost are wrapped by `bs_wrap` to NaN or +∞ (see `Sources/CBoostBridge/internal/bs_internal.hpp:1`).
  - Returns `false` if the name is unknown or a required parameter is missing/invalid.

- Integration:
  - Included by umbrella header `Sources/CBoostBridge/include/CBoostBridge.h:63` to expose the ABI to Swift.
  - Compiled into the bridge via `Sources/CBoostBridge/CBoostBridge.cpp:67` include list.


## Swift Wrapper: `Distribution.Dynamic<T>`

- Declaration: `Sources/SwiftyBoost/Math/Distributions/DynamicDistribution.swift:1`
- Usage:
  - Construct by name with a parameter dictionary:
    - `let g = try Distribution.Dynamic<Double>(distributionName: "gamma", parameters: ["shape": 4.5, "scale": 1.0])`
    - `let t = try Distribution.Dynamic<Float>(distributionName: "student_t", parameters: ["df": 12])`
  - Call as any distribution:
    - `let p = try g.pdf(2.0)`
    - `let x = try t.quantile(0.95)`
- Implementation details:
  - The initializer converts the Swift dictionary to a C array of `bs_param_*` (keys `strdup`’d and freed immediately after the factory call returns).
  - Chooses the precision factory by `T`:
    - `Double` → `bs_dist_make_d`, `Float` → `bs_dist_make_f`; `Float80` on x86 → `bs_dist_make_l` (gated by `#if arch(x86_64) || arch(i386)`).
  - Stores the returned vtable struct in a small `Box` class and cleans up via `free(ctx)` in `deinit`.
  - All protocol methods forward to the corresponding vtable function pointer if present.
  - Swift supplies fallbacks for metrics absent from the vtable (e.g. Bernoulli entropy, Cauchy entropy, Fisher F / Arcsine entropy) so the high-level API remains uniform.
  - For stats (`mean`, `variance`, …) non‑finite values are mapped to `nil`.


## Error Handling and Policies

- Factory errors:
  - Unknown `name` or missing required parameters → `bs_dist_make_*` returns `false`; Swift throws `DistributionError.invalidCombination`.
- Numeric errors during evaluation:
  - Boost throws (domain/overflow) → caught and translated by `bs_wrap` to NaN or +∞; Swift then returns optionals or raw values depending on the API.
- Range and support:
  - `range()` returns the Boost range; may include infinities.


## Performance

- Overhead is a single function‑pointer call per evaluation before entering Boost; negligible relative to the cost of special function and distribution computations.
- Handles are constructed once and reused; no internal allocations per call.


## Concurrency

- The Swift wrapper is `Sendable`; the underlying `ctx` hosts an immutable Boost distribution object once constructed.
- Concurrent reads (pdf/cdf/quantile/…) are safe across threads.


## Extending the Factory (Adding a New Distribution)

To support a new distribution everywhere (Double/Float/Long Double):

1) Include the Boost distribution in the implementation TU

- Edit `Sources/CBoostBridge/impl/bs_generic_distribution.hxx:1` and add the appropriate Boost header, e.g.:
  - `#include <boost/math/distributions/fisher_f.hpp>`

2) Define small handle types per precision

- Example:
  - `struct fisher_f_d_handle { boost::math::fisher_f_distribution<double> dist; fisher_f_d_handle(double m, double n): dist(m,n) {} };`
  - Repeat for `float` and `long double`.

3) In `bs_dist_make_d/f/l`, add a case for the name

- Normalize name (already done); add case labels and parse parameters by alias.
- Validate required parameters (return `false` if missing).
- Allocate the handle with `new (std::nothrow)`; set `out->ctx` and all vtable function pointers using the provided thunk templates:
  - `out->pdf = &pdf_fn<fisher_f_d_handle, double>;` (and similarly for others)
  - Range/stats: `out->range = &range_fn<fisher_f_d_handle, bs_range_d, double>;`, `out->mean = &mean_fn<...>;`, etc.
  - Lifetime: `out->free = &free_fn<fisher_f_d_handle>;`

4) No Swift changes required

- `Distribution.Dynamic<T>` automatically works for the new name and aliases once the C factory supports them.


## Parameter Keys and Aliases

- All parameter keys are matched case‑insensitively.
- Alias sets for each distribution are listed in “What’s Implemented Today”; quick reference:
  - Gamma: `shape|k`, `scale|theta`
  - Beta: `alpha|a|p|shape1`, `beta|b|q|shape2`
  - Chi-squared: `df|nu|degreesOfFreedom`
  - Student’s t: `df|nu|degreesOfFreedom`
  - Fisher’s F: `df1|d1|m|degreesOfFreedom1`, `df2|d2|n|degreesOfFreedom2`
  - Bernoulli & Geometric: `p|prob|probability|success|theta`
  - Binomial: `n|trials`, `p|prob|probability|success`
  - Cauchy & Holtsmark: `location|loc|mu|median|x0`, `scale|gamma|sigma|b`
  - Exponential: `lambda|rate` (scale handled by typed wrapper convenience API)
  - Extreme value / Gumbel: `location|loc|mu`, `scale|gamma|sigma|b`
  - Arcsine: `minX|min|a|lower`, `maxX|max|b|upper`
- Add aliases by extending the key arrays in the factory implementation.


## Comparison: Dynamic vs Typed Swift Wrappers

- Dynamic (this factory):
  - Pros: one Swift type, easy to add distributions (C++ only), flexible runtime construction from names.
  - Cons: parameter validation happens at runtime via strings; less compiler guidance than typed inits.

- Typed wrappers (existing in SwiftyBoost):
  - Pros: strong parameter labeling, better doc and discoverability, direct mapping to domain terms.
  - Cons: one Swift file per distribution.

It’s reasonable to keep both: use typed wrappers for direct, high‑clarity usage; use the dynamic factory for generic tooling, serialization, or configuration‑driven scenarios.


## Examples

Swift (Double precision):

```swift
let g = try Distribution.Dynamic<Double>(
  distributionName: "gamma",
  parameters: ["shape": 4.5, "scale": 1.0]
)
let p = try g.pdf(2.0)
let F = try g.cdf(2.0)
let x95 = try g.quantile(0.95)
```

Swift (Float, Student’s t):

```swift
let t = try Distribution.Dynamic<Float>(
  distributionName: "student_t",
  parameters: ["df": 12]
)
let tail = try t.sf(1.96)
```

Swift (Binomial counts):

```swift
let bin = try Distribution.Dynamic<Double>(
  distributionName: "binomial",
  parameters: ["n": 12, "p": 0.35]
)
let pmf4 = try bin.pdf(4)
let cdf4 = try bin.cdf(4)
```

Error handling:

```swift
do {
  _ = try Distribution.Dynamic<Double>(distributionName: "gamma", parameters: [:])
} catch {
  print("Factory error: \(error)") // missing required parameter
}
```


## Implementation Walkthrough (File by File)

- `Sources/CBoostBridge/include/distributions/bs_generic_distribution.h`
  - Declares the vtable structs `bs_dist_d/f/l`, parameter records, and `bs_dist_make_*` factories.

- `Sources/CBoostBridge/impl/bs_generic_distribution.hxx`
  - Provides helper templates for the vtable thunks and the factory switches for supported names.
  - Uses `bs_wrap` (see `Sources/CBoostBridge/internal/bs_internal.hpp:1`) to translate Boost exceptions to IEEE‑754 results.

- `Sources/CBoostBridge/CBoostBridge.cpp:67`
  - Includes the implementation TU include to ensure linkage for the new symbols.

- `Sources/CBoostBridge/include/CBoostBridge.h:63`
  - Re‑exports the generic distribution headers to Swift as part of the umbrella.

- `Sources/SwiftyBoost/Math/Distributions/DynamicDistribution.swift:1`
  - Swift wrapper that converts a Swift `[String: T]` into `bs_param_*[]`, invokes `bs_dist_make_*`, and stores the returned vtable.
  - Forwards `DistributionProtocol` methods to the vtable; maps non‑finite stats to `nil`.
  - Gates `Float80` support by architecture.


## Memory Management Details

- The factory builds a small C array of parameters. Keys are `strdup`’d and freed immediately after the factory call returns. The vtable/handle does not retain parameter keys.
- The vtable holds an opaque `ctx` with the constructed Boost object. `free(ctx)` is called in the Swift wrapper’s `deinit`.
  - Swift initializes the vtable by allocating an `UnsafeMutablePointer<bs_dist_*>`, invoking `bs_dist_make_*` to populate it, then reading back the struct. This avoids zero-initializing a `_Nonnull` field.


## Known Limitations and Notes

- `Float80` is only available on x86 architectures in Swift; on other platforms, the dynamic wrapper will use `Double` or `Float` variants.
- Some stats (e.g., mode for certain distributions) may be undefined; their thunks return NaN which maps to `nil` in Swift.
- Discrete/lattice distributions are not yet modeled; `latticeStep`/`latticeOrigin` are `nil`.
- The factory currently supports a subset of distributions; adding more is straightforward (see “Extending”).


## Testing Recommendations

- For each distribution supported by the factory:
  - Compare `Distribution.Dynamic<T>` results against the existing typed Swift wrappers for a small grid of inputs.
  - Verify quantile/CDF inverses on a set of probabilities.
  - Exercise parameter alias parsing.


## API Stability

- The vtable layout is considered part of the C ABI surface of CBoostBridge for this project. Add new fields only by appending to the end and provide defaulting behavior to avoid breaking older clients in the same module.


## FAQ

Q: Why a vtable and not a single giant switch in Swift?

A: Centralizing the distribution semantics and parameter parsing near Boost keeps Swift thin and means adding a distribution is a C++‑only change. It also allows consistent exception/numeric error handling within the bridge.

Q: What is the runtime cost?

A: Essentially the cost of one function‑pointer call before entering Boost; negligible relative to distribution evaluations.

Q: Can I keep using the typed wrappers?

A: Yes. They remain supported. The dynamic factory is additive and ideal for configuration‑driven or plugin‑like usage.


---

For questions, open an issue or see the referenced files. The new factory integrates without changing existing APIs, and all tests continue to pass.
