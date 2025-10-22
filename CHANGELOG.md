# Changelog

All notable changes to this project are tracked here, following the principles of [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and adhering to [Semantic Versioning](https://semver.org/).

## Unreleased
- Nothing yet.

## [0.5.0] - 2025-10-22
- Dynamic distribution factory (runtime):
  - Added a unified vtable-based factory in `CBoostBridge` with non-null context pointer (`ctx`) and nullable function pointers for metrics.
  - New Swift wrapper `Distribution.Dynamic<T>` that conforms to `DistributionProtocol` and constructs distributions by name with a parameter dictionary.
  - Supported names and aliases: `gamma|gamma_distribution`, `beta|beta_distribution`, `chi_squared|chisquared|chi2|chi-squared|chisquare`, `student_t|studentt|students_t|t|t_distribution`, `fisherf|f|f_distribution`, `bernoulli|bernoulli_distribution`, `arcsine|arcsine_distribution`.
  - Swift-side fallbacks for metrics not present in vtable: differential entropy for Fisher F and Arcsine.
  - Comprehensive tests comparing Dynamic vs typed wrappers (PDF/CDF/quantile, hazard/CHF, moments), alias handling, error cases, and Float80 smoke on x86_64.
  - New DocC page: “Dynamic Distribution Factory” with usage, aliases, nullability, and initialization notes.
  - Documentation: `DIST-Factory-README.md` details the ABI, design, and extension guide.
  - README updated with a “Dynamic Distribution Factory” section and examples.
  - Typed wrappers refactor: ``Distribution.StudentT``, ``Distribution.FisherF``, and ``Distribution.Arcsine`` now delegate internally to ``Distribution.Dynamic`` instead of holding per-type C handles. This reduces duplication and keeps behavior consistent with the factory. Public APIs and semantics are unchanged.
  - Header nullability tightened: `ctx` marked `_Nonnull`; all vtable function pointers marked `_Nullable` with `const void* _Nonnull` context. Factory prototypes annotated `_Nonnull` for `name`, `params`, and `out`.
  - Robust factory error handling: `bs_dist_make_*` now catch exceptions and return `false` on invalid parameters.
  - Removed the experimental `AnyDistribution` path; use `Distribution.Dynamic<T>` instead. Demo updated accordingly.
- Complex numbers:
  - Replaced the bespoke ``Complex`` struct with the Swift Numerics implementation, re-exported via `SwiftyBoostPrelude`.
  - Retained SwiftyBoost conveniences (`.i`, `fromPolar`, stable magnitude/square root) and Boost-backed elementary functions through extensions.
  - Updated README and DocC guidance to note the Numerics dependency.
- Constants namespace:
  - Introduced `Constants<T>` to expose Boost-backed mathematical constants (π, τ, √π, Euler–Mascheroni, Catalan, ζ(3), φ, logarithms, ratio variants) with dedicated bindings for `Float`, `Double`, and (x86_64) `Float80`.
  - Added inline documentation for every constant accessor plus DocC coverage and README usage examples.
- Mixed-precision promotions:
  - Beta (complete/incomplete/regularized, inverses, derivative, parameter solvers)
  - Gamma (ratios, incomplete/regularized P/Q, inverses, derivative)
  - Bessel (J/Y/I/K + derivatives, integer-order helpers)
  - Elliptic (Legendre F/E/Π; Carlson RC/RF/RD/RJ/RG)
  - Chebyshev Clenshaw ([coefficients] vs `x` mixes)
  - Owen’s T (all pairs)
  - Spherical Harmonics (mixed angle types)
  - Hypergeometric (1F0/0F1 full; 1F1/2F0 selected mixes; pFq arrays with z)
  - Float80 (x86_64) promotions where applicable
- Added Boost-backed helpers:
  - `rsqrt` (reciprocal square root) with domain checks and unit tests for Double/Float/(x86_64) Float80.
  - Sinus cardinalis (π-normalized) functions:
    - `sinc_pi(x) = sin(πx)/(πx)` with removable singularity at 0.
    - `sinhc_pi(x) = sinh(πx)/(πx)` with removable singularity at 0.
    - Complex `sincc_pi(z) = sin(πz)/(πz)` for `ComplexD/F/(x86_64) ComplexL`.
  - Bridge includes and mappings hardened for Boost 1.89.0 (sinc.hpp, sinhc.hpp); π-normalized forms used explicitly to avoid header ambiguities.
- Documentation:
  - Added DocC article “Constants Catalog” documenting the `Constants<T>` namespace, precision choices, and common values.
  - Added distribution references for “Distribution/Beta”, “Distribution/ChiSquared”, and “Distribution/Bernoulli”; refreshed the Quickstart and Dynamic Factory guides with matching examples and alias tables.
  - New DocC page “Result-Type Promotions” (policy, supported APIs, and guidance).
  - README quick reference and usage examples for promotions, `rsqrt`, and Sinus cardinalis.
  - New DocC page “Complex Numbers” documenting `Complex<T>` construction, arithmetic, polar/phase, and Boost-backed elementary functions; linked from the module index.
  - Fixed README typealias to match code: `ComplexL` (Float80) instead of the former custom alias.
  - Added distribution reference pages: “Distribution/FisherF” and “Distribution/Arcsine” with definitions, API overview, usage, and notes.
  - Added special-function references: “Sinus Cardinalis (π-Normalized)” for `sinc_pi`/`sinhc_pi` and complex variants; “Common Helpers and Stable Transforms” for `rsqrt`, `expm1`, `log1p`, `log1pmx`, `powm1`, `cbrt`, `sqrt1pm1`.
  
  Distributions:
  - Added distribution types bridged from Boost:
    - `Distribution.Gamma<T>` (shape/scale) with PDF/CDF/SF, quantiles, hazards, and moments.
    - `Distribution.Beta<T>` (alpha/beta) with PDF/CDF/SF, quantiles, hazards, and moments.
    - `Distribution.ChiSquared<T>` (degrees of freedom) with PDF/CDF/SF, quantiles, hazards, and moments.
    - `Distribution.StudentT<T>` (ν) with PDF/CDF/SF, quantiles, hazards, moments, and power-planning helper.
    - `Distribution.FisherF<T>` (df1, df2) with PDF/CDF/SF, quantiles, hazards, and moments.
    - `Distribution.Bernoulli<T>` (success probability) with PMF/CDF, quantiles, hazards, and entropy fallback.
    - `Distribution.Arcsine<T>` (minX, maxX) with PDF/CDF/SF, quantiles, hazards, and moments.
  - All distributions are generic over `BinaryFloatingPoint` with `Double`/`Float` overloads and `Float80` on x86_64.
- Tests:
  - Added suites for mixed promotions across modules and for new helpers.
  - Real and complex tests for `sinc_pi`, `sinhc_pi`, and `sincc_pi` with identity checks and tolerances.
- Added Swift Numerics-backed complex numbers via typealiasing and extensions
  - Swapped the bespoke `Complex<T>` for `ComplexModule.Complex<T>` and extended it with legacy conveniences (Smith division, stable magnitude, phase/polar helpers).
  - Maintained `ComplexD`, `ComplexF`, and (x86_64) `ComplexL` aliases plus the Boost-powered elementary functions (`bs_cexp`, `bs_clog`, `bs_csin`, `bs_ccos`, `bs_ctan`, `bs_csinh`, `bs_ccosh`, `bs_ctanh`, `bs_catan`).
  - Updated README and DocC to note the Numerics dependency and preserved ExpressibleByFloatLiteral compatibility for tests.
- Corrected README usage examples to match public API naming (`errorFunction`, `besselJ(v:x:)`).
- Aligned supported platforms in README with `Package.swift` (macOS, iOS).
- Updated SPM installation URL to `github.com/strike65/SwiftyBoost`.
- Clarified DocC generation and output paths; added a link to the changelog and versioning policy.
- Added GitHub Actions workflow to build DocC and deploy `Docs/` to GitHub Pages on pushes to `main` (see `.github/workflows/docs.yml`).
- Added README badges: documentation link to `https://strike65.github.io/SwiftyBoost/` and Pages workflow status badge.
- Added direct “Browse online” link under Documentation section in README.
- Docs site root now redirects to the main DocC page for a cleaner entry URL (`https://strike65.github.io/SwiftyBoost` → `/documentation/swiftyboost/`); implemented by writing `Docs/index.html` during Pages deployment.
 - C++ bridge refactor: split the monolithic `CBoostBridge.h` into focused public headers under `Sources/CBoostBridge/include/` and introduced an umbrella header `CBoostBridge.h` that re-exports them. Broke up implementation into maintainable units under `Sources/CBoostBridge/impl/`, composed by `Sources/CBoostBridge/CBoostBridge.cpp`, with shared helpers in `internal/bs_internal.hpp`. No API or symbol changes (all `bs_*` names preserved).

## [0.1.0] - 2025-10-12
- Expanded Boost-backed special function coverage, including additional elliptic, hypergeometric, and Lambert W variants.
- Added the initial Swift Testing suite covering factorial and combinatorics APIs.
- Removed the experimental `ellint3` exposure pending stabilization.

## [0.0.1] - 2025-10-11
- Established the `SwiftyBoost` Swift package with the `CBoostBridge` interop target.
- Vendored Boost.Math headers under `extern/boost`.
- Published core Gamma, Beta, error, Legendre, and Bessel function wrappers with Swift-friendly APIs and error handling.

[0.5.0]: https://github.com/strike65/SwiftyBoost/releases/tag/0.5.0
[0.1.0]: https://github.com/strike65/SwiftyBoost/releases/tag/0.1.0
[0.0.1]: https://github.com/strike65/SwiftyBoost/releases/tag/0.0.1
