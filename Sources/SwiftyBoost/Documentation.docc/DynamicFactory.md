# Dynamic Distribution Factory

Construct probability distributions by name at runtime using a small C vtable and a Swift wrapper that conforms to ``DistributionProtocol``.

The dynamic factory complements the typed wrappers (e.g. ``Distribution/Gamma``) and is ideal for configuration‑driven or plugin‑style code where the concrete distribution is chosen at runtime.

## Overview

- Swift entry point: ``Distribution/Dynamic``
- Backed by a C ABI vtable (see `CBoostBridge`) with a non‑null `ctx` pointer and nullable function pointers for metrics.
- Parameters are passed as a dictionary and matched case‑insensitively with sensible aliases per distribution.
- If a metric isn’t provided by the vtable for a distribution (e.g. some entropies), Swift computes a fallback when a closed form is available.

## Usage

```swift
// Gamma(k, theta)
let g = try Distribution.Dynamic<Double>(
  distributionName: "gamma",
  parameters: ["shape": 4.5, "scale": 1.0]
)
let p = try g.pdf(2.0)
let x95 = try g.quantile(0.95)

// Student's t (aliases: student_t, t)
let t = try Distribution.Dynamic<Float>(
  distributionName: "student_t",
  parameters: ["df": 12]
)
let sf = try t.sf(1.96)

// Fisher F (aliases: fisherf, f)
let f = try Distribution.Dynamic<Double>(
  distributionName: "fisherf",
  parameters: ["df1": 10, "df2": 20]
)
let H = try f.chf(2.0)

// Beta
let beta = try Distribution.Dynamic<Double>(
  distributionName: "beta",
  parameters: ["alpha": 2.0, "beta": 5.0]
)
let mean = beta.mean

// Chi-squared
let chi2 = try Distribution.Dynamic<Double>(
  distributionName: "chisquared",
  parameters: ["df": 8.0]
)
let chiSF = try chi2.sf(12.0)

// Bernoulli
let bern = try Distribution.Dynamic<Double>(
  distributionName: "bernoulli",
  parameters: ["p": 0.3]
)
let pmf = try bern.pdf(1)

// Binomial
let bin = try Distribution.Dynamic<Double>(
  distributionName: "binomial",
  parameters: ["n": 12.0, "p": 0.35]
)
let countCdf = try bin.cdf(4.0)

// Arcsine (aliases: arcsine_distribution)
let a = try Distribution.Dynamic<Float>(
  distributionName: "arcsine",
  parameters: ["minX": 0, "maxX": 1]
)
let c = try a.cdf(0.25)

// Continuous location-scale family (Cauchy, includes alias set)
let cauchy = try Distribution.Dynamic<Double>(
  distributionName: "cauchy",
  parameters: ["loc": 0.25, "scale": 1.5]
)
let median = cauchy.median

// Exponential (aliases: exponential_distribution, exp)
let exp = try Distribution.Dynamic<Double>(
  distributionName: "exp",
  parameters: ["rate": 1.2]
)
let tail = try exp.sf(2.0)

// Gumbel / Extreme Value
let gumbel = try Distribution.Dynamic<Double>(
  distributionName: "gumbel",
  parameters: ["loc": 1.0, "scale": 0.75]
)
let peak = gumbel.mode
```

## Supported Names and Aliases

- Gamma: `gamma`, `gamma_distribution`
  - Params: `shape|k` (required), `scale|theta` (default = 1)
- Beta: `beta`, `beta_distribution`
  - Params: `alpha|a|p|shape1` (required), `beta|b|q|shape2` (required)
- Chi-squared: `chisquared`, `chi_squared`, `chi2`, `chi-squared`, `chisquare`
  - Params: `df|nu|degreesOfFreedom` (required)
- StudentT: `studentt`, `student_t`, `students_t`, `t`, `t_distribution`
  - Params: `df|nu|degreesOfFreedom` (required)
- FisherF: `fisherf`, `f`, `f_distribution`
  - Params: `df1|d1|m|degreesOfFreedom1` (required), `df2|d2|n|degreesOfFreedom2` (required)
- Bernoulli: `bernoulli`, `bernoulli_distribution`
  - Params: `p|prob|probability|success|theta` (required)
- Binomial: `binomial`, `binomial_distribution`
  - Params: `n|trials` (required, interpreted as the trial count), `p|prob|probability|success` (required)
- Cauchy: `cauchy`, `cauchy_distribution`
  - Params: `location|loc|mu|median|x0` (optional, defaults to 0), `scale|gamma|sigma|b` (required, > 0)
- Exponential: `exponential`, `exponential_distribution`, `exp`
  - Params: `lambda|rate` (required, > 0)
- Extreme value (Gumbel): `extremevalue`, `extreme_value`, `gumbel`, `extreme_value_distribution`
  - Params: `location|loc|mu` (optional, defaults to 0), `scale|gamma|sigma|b` (required, > 0)
- Geometric: `geometric`, `geometric_distribution`
  - Params: `p|prob|probability|success|theta` (required)
- Holtsmark: `holtsmark`, `holtsmark_distribution`
  - Params: `location|loc|mu|median|x0` (optional, defaults to 0), `scale|gamma|sigma|b` (required, > 0)
- Arcsine: `arcsine`, `arcsine_distribution`
  - Params: `minX|min|a|lower` (required), `maxX|max|b|upper` (required)

### Additional Notes

- ``Distribution/Geometric`` and ``Distribution/Holtsmark`` both reuse the factory under the hood; the dynamic entries remain valuable for configuration-driven scenarios that need to swap distributions at runtime.

## Nullability and Initialization

- The vtable’s `ctx` is `_Nonnull` once constructed; function pointers are `_Nullable` (some metrics may be absent).
- The Swift wrapper allocates an `UnsafeMutablePointer<bs_dist_*>`, calls the C factory to populate it, then reads the struct and releases the pointer.
- Metrics not available in the vtable return `nil` from Swift (e.g. entropy when undefined or not provided). Where a closed form is known, Swift computes it as a fallback (e.g. Fisher F, Arcsine).

## Errors

- Unknown name or missing required parameters → the Swift initializer throws ``DistributionError/invalidCombination(message:value:)`` and, when available, populates ``value`` with the problematic argument.
- Domain violations and numerical errors are mapped to IEEE‑754 values by the bridge; Swift then exposes them as `nil` for optional metrics or returns numeric values for queries.

## Extensibility

Adding new distributions requires only a small extension of the C factory to map names and aliases to Boost objects and wire the vtable pointers. The Swift wrapper does not need changes.

See `DIST-Factory-README.md` for an in‑depth design and an extension guide.
