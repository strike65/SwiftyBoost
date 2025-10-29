# Changelog

All notable changes to this project are tracked here, following the principles of [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and adhering to [Semantic Versioning](https://semver.org/).

## Unreleased
- _No changes yet._

## [1.0.3] - 2025-10-29
### Removed
- Cleared README and DocC references to result-type promotions now that the feature is no longer supported.

### Changed
- Raised minimum deployment targets to macOS 11, iOS 16, watchOS 9, and tvOS 16 to match the updated package manifest.

## 1.0.2
- romved unnecessary folder
## 1.0.1
- fix: typo in README.md
## [1.0.0] - 2025-10-26
- Features:
  - Added ``Distribution.Empirical`` — a pure Swift empirical distribution with automatic discrete/continuous detection, KNN/KDE estimators, Miller–Madow corrections, bootstrap confidence intervals, and multimodality heuristics.
  - Added Boost-backed Logistic, Log-normal, Map-Airy, and Negative Binomial distributions to the unified factory with typed wrappers ``Distribution.Logistic``, ``Distribution.LogNormal``, ``Distribution.MapAiry``, and ``Distribution.NegativeBinomial`` (parameter validation, runtime aliases, heavy-tail handling).
  - Added Boost-backed Kolmogorov–Smirnov, Laplace (double exponential), and Landau distributions to the unified factory and introduced the typed wrappers ``Distribution.KolmogorovSmirnov``, ``Distribution.Laplace``, and ``Distribution.Landau`` (including parameter validation, moments, and runtime aliases).
  - Added configurable ``Distribution.KLDivergenceOptions`` plus a general-purpose ``Distribution.Dynamic.klDivergence(relativeTo:options:)`` implementation that covers continuous (finite/semi-infinite/infinite) and discrete (finite/unbounded lattice) distributions; typed wrappers now expose ``klDivergence(relativeTo:options:)`` and forward to the dynamic backend whenever a closed form is unavailable.
  - Added Boost-backed quadrature helpers ``SpecialFunctions.Quadrature`` with Gauss–Legendre, Gauss–Kronrod, tanh–sinh, sinh–sinh, and exp–sinh rules plus reusable integrators, metadata, and abscissa extraction.
  - Added Boost-backed Legendre–Stieltjes polynomials to the bridge plus Swift wrappers for evaluation, derivatives, zeros, and norm-squared helpers.
  - Added Hermite polynomial wrappers, including the three-term recurrence helper `hermiteNext` for efficient sequence stepping.
  - Added Jacobi polynomials with first/second/k-th derivative access across precisions.
  - Added the Jacobi Zeta function ``SpecialFunctions.jacobiZeta(_:modulus:)`` across Float/Double/(x86) Float80.
  - Added Airy Ai/Bi and derivative wrappers ``SpecialFunctions.airyAi`` / ``SpecialFunctions.airyBi`` / ``SpecialFunctions.airyAiPrime`` / ``SpecialFunctions.airyBiPrime`` along with zero helpers ``SpecialFunctions.airyAiZero`` / ``SpecialFunctions.airyBiZero`` and bulk retrieval APIs.
  - Added Heuman’s lambda via ``SpecialFunctions.heumanLambda(_:phi:)`` aligned with Boost’s elliptic toolkit.
  - Added Jacobi elliptic functions (`sn`, `cn`, `dn`, quotient variants) through ``SpecialFunctions.jacobiElliptic(_:theta:)`` (returns the triple) plus individual helpers with Float/Float80 support.
  - Added falling factorial ``SpecialFunctions.falling_factorial(_:_:)->T`` with Float/Double/(x86) Float80 overloads, mirroring Boost’s implementation with finiteness checks, zero detection, and overflow/underflow guards.
  - Added Boost-backed hyperexponential distribution support to the factory and a typed wrapper ``Distribution.Hyperexponential`` that accepts arbitrary indexed phase rates plus optional mixing probabilities.
  - Added Boost-backed inverse chi-squared distribution support in the factory and introduced the typed wrapper ``Distribution.InverseChiSquared`` with default and explicit scale initialisers.
  - Added Boost-backed inverse gamma and inverse normal (inverse Gaussian/Wald) distributions to the factory alongside the typed wrappers ``Distribution.InverseGamma`` and ``Distribution.InverseNormal`` (with the ``Distribution.InverseGaussian`` typealias), including parameter validation and runtime alias coverage.
  - Added Boost-backed non-central chi-squared distribution support to the factory and provided the typed wrapper ``Distribution.NonCentralChiSquared`` (analytic moments in Swift, runtime aliases, λ → 0 parity with ``Distribution.ChiSquared``).
  - Added Boost-backed non-central F and non-central Student’s t distributions to the factory with typed wrappers ``Distribution.NonCentralF`` and ``Distribution.NonCentralStudentT`` (parameter validation, runtime aliases, central-case parity).
  - Added Boost-backed Pareto and Poisson distributions to the factory with typed wrappers ``Distribution.Pareto`` and ``Distribution.Poisson`` (analytic moments in Swift, discrete lattice metadata for Poisson).
  - Added Boost-backed Rayleigh, SaS Point5 (alpha = 1/2), Skew-normal, Triangular, Uniform, and Weibull distributions to the unified factory with typed wrappers ``Distribution.Rayleigh``, ``Distribution.SASPoint5``, ``Distribution.SkewNormal``, ``Distribution.Triangular``, ``Distribution.Uniform``, and ``Distribution.Weibull`` (parameter validation, alias coverage, moment/entropy plumbing, and KL divergence interop).
- Documentation:
  - Documented the empirical distribution (DocC article, distribution index entries, README examples, and factory guidance noting it is a pure Swift implementation).
  - Documented the Logistic, Log-normal, Map-Airy, and Negative Binomial distributions (DocC pages, distribution index, README feature list) and updated factory references / design notes.
  - Documented the Rayleigh, SaS Point5, Skew-normal, Triangular, Uniform, and Weibull distributions (DocC coverage, README snippets, factory alias tables, and `DIST-Factory-README.md` updates).
  - Documented the Kolmogorov–Smirnov, Laplace, and Landau distributions (DocC articles, distribution index, README feature list) and refreshed factory alias tables / `DIST-Factory-README.md` coverage for the new names.
  - Documented quadrature helpers across DocC (`<doc:Quadrature>`), README usage, and the new `QUADRATURE-README.md` workflow guide.
  - Extended DocC guidance with an “Orthogonal Polynomials” overview and refreshed README examples to highlight the new wrappers plus Jacobi Zeta.
  - Added documentation references to the Airy helpers.
  - Documented Jacobi theta helpers (q and τ parameterizations) across README and DocC usage guides.
  - Documented falling factorial helpers in DocC and refreshed README combinatorics examples.
  - Documented KL divergence APIs (README usage snippet, `<doc:DistributionsQuickstart>`, and `<doc:DynamicFactory>`) and clarified tuning knobs exposed by ``Distribution.KLDivergenceOptions``.
  - Documented the hyperexponential distribution across DocC (`<doc:DistributionsHyperexponential>`), README usage snippets, and the dynamic factory alias tables.
  - Documented the inverse chi-squared distribution across README, DocC (`<doc:DistributionsInverseChiSquared>`), and the factory reference, including usage examples and alias tables.
  - Documented the inverse gamma and inverse normal (inverse Gaussian) distributions across README samples, DocC articles (`<doc:DistributionsInverseGamma>`, `<doc:DistributionsInverseNormal>`), the dynamic factory reference, and `DIST-Factory-README.md`.
  - Documented the non-central chi-squared distribution (DocC article `<doc:DistributionsNonCentralChiSquared>`, distribution index, README feature list, and factory alias tables).
  - Documented the non-central F and non-central Student’s t distributions (DocC articles `<doc:DistributionsNonCentralF>`, `<doc:DistributionsNonCentralStudentT>`), updated the distribution index, README feature list, and dynamic factory alias tables.
  - Documented the Pareto and Poisson distributions (DocC articles `<doc:DistributionsPareto>`, `<doc:DistributionsPoisson>`), updated the distribution index, README, and factory alias tables.
- Bug Fixes:
  - Corrected ``SpecialFunctions.jacobiTheta2Tau(_:tau:)`` (`Float` overload) to call the θ₂ bridge variant instead of θ₁.
  - Updated tau-parameterized theta helpers to report `tau` as the offending argument when throwing ``SpecialFunctionError.parameterNotFinite``.
  - Fixed the C bridge mappings for θ₂/θ₂τ to call the appropriate Boost `jacobi_theta2{,tau}` kernels.
  - Fixed ``Distribution.Bernoulli.klDivergence(relativeTo:options:)`` and ``Distribution.Exponential.klDivergence(relativeTo:options:)`` to return nats in line with the dynamic backend and corrected the exponential analytic form.
- Tests:
  - Added CSV-backed empirical fixtures (continuous/discrete, n=200) with Swift Testing coverage for entropy, KL divergence estimators, lattice metadata, and multimodality heuristics.
  - Added regression and parity suites for ``Distribution.Logistic``, ``Distribution.LogNormal``, ``Distribution.MapAiry``, and ``Distribution.NegativeBinomial`` (typed vs dynamic behaviour, analytic properties, entropy scaling).
  - Added analytic and dynamic parity tests for ``Distribution.NonCentralChiSquared`` plus regression coverage for the central limit (λ → 0) case.
  - Added parity and regression suites for ``Distribution.NonCentralF`` and ``Distribution.NonCentralStudentT`` (analytic moments, dynamic equality, λ → 0 reductions).
  - Added validation and dynamic parity suites for ``Distribution.Pareto`` and ``Distribution.Poisson`` (analytic moments, discrete behaviour).
  - Added regression and parity tests for ``Distribution.KolmogorovSmirnov``, ``Distribution.Landau``, and ``Distribution.Laplace`` (typed wrappers and dynamic factory) covering analytic moments, entropy scaling, and symmetry properties.
  - Added regression tests covering KL divergence for `Distribution.Dynamic` (exponential/Bernoulli) and the corresponding typed wrappers to ensure they track analytic formulas and exercise both continuous and discrete code paths.
  - Added dedicated suites covering Legendre–Stieltjes, Hermite, Jacobi polynomials, Jacobi Zeta, and Airy functions against the C bridge and invalid-parameter paths.
  - Added regression tests for falling factorial wrappers and Jacobi theta functions (Float/Double/(x86) Float80) including domain-validation paths.
  - Added parity tests for the hyperexponential distribution between ``Distribution.Dynamic`` and the new typed wrapper plus analytical regression coverage (moments, hazards, quantiles).
  - Added regression coverage for the inverse chi-squared distribution (typed vs dynamic parity, analytic moments, and parameter validation).
  - Added regression coverage for the inverse gamma and inverse normal distributions (typed vs dynamic parity, analytic moment checks, parameter validation).

## [0.6.0] - 2025-10-23
- Features:
  - Added ``Distribution.Holtsmark`` as a typed wrapper over Boost.Math’s heavy-tailed Holtsmark distribution (location/scale form) with full PDF/CDF/quantile coverage and entropy support.
- Documentation:
  - Linked all references to Swift Numerics’ ``Real`` protocol directly to the upstream project documentation for easier navigation.
  - Added DocC coverage for ``Distribution/Geometric`` and refreshed README, Quickstart, and factory guides with geometric examples and alias tables.
  - Renamed ``HighPrecisionConstants`` to ``Constants`` throughout the guides to reflect the streamlined API.
- Tests:
  - Added regression coverage for ``Distribution.Geometric`` (closed-form PMF/CDF, quantiles, moments) and ensured the dynamic factory’s geometric entry matches the typed wrapper and alias set.

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
- High-precision constants:
  - Introduced `Constants` (named `HighPrecisionConstants` in the original 0.5.0 release) to expose Boost-backed mathematical constants (π, τ, √π, Euler–Mascheroni, Catalan, ζ(3), φ, logarithms, ratio variants) convertible to any `BinaryFloatingPoint` type at call-site.
  - Added inline documentation for every helper plus DocC coverage and README usage examples.
- Added Boost-backed helpers:
  - `rsqrt` (reciprocal square root) with domain checks and unit tests for Double/Float/(x86_64) Float80.
  - Sinus cardinalis (π-normalized) functions:
    - `sinc_pi(x) = sin(πx)/(πx)` with removable singularity at 0.
    - `sinhc_pi(x) = sinh(πx)/(πx)` with removable singularity at 0.
    - Complex `sincc_pi(z) = sin(πz)/(πz)` for `ComplexD/F/(x86_64) ComplexL`.
  - Bridge includes and mappings hardened for Boost 1.89.0 (sinc.hpp, sinhc.hpp); π-normalized forms used explicitly to avoid header ambiguities.
- Documentation:
  - Added DocC article “Constants Catalog” documenting the constant helpers, precision choices, and common values.
  - Added distribution references for “Distribution/Beta”, “Distribution/ChiSquared”, and “Distribution/Bernoulli”; refreshed the Quickstart and Dynamic Factory guides with matching examples and alias tables.
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

[1.0.0]: https://github.com/strike65/SwiftyBoost/releases/tag/1.0.0
[0.6.0]: https://github.com/strike65/SwiftyBoost/releases/tag/0.6.0
[0.5.0]: https://github.com/strike65/SwiftyBoost/releases/tag/0.5.0
[0.1.0]: https://github.com/strike65/SwiftyBoost/releases/tag/0.1.0
[0.0.1]: https://github.com/strike65/SwiftyBoost/releases/tag/0.0.1
