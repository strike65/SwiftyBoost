#
``SwiftyBoost/Distribution/NonCentralF``

Non-central Fisher–Snedecor F distribution with numerator degrees of freedom `ν₁`, denominator degrees of freedom `ν₂`, and non-centrality parameter `λ`.

## Definition

The non-central F distribution describes ratios of (possibly non-central) chi-squared variates. It is supported on `x ∈ [0, +∞)` with density

```
f(x; ν₁, ν₂, λ) = e^{-(λ + ν₁ x)/(2x)} \sum_{k=0}^{\infty}
  \frac{(λ/2)^k}{k!} \frac{(ν₁ x)^{(ν₁/2 + k) - 1}}{(ν₁/2 + k - 1)!}
  \frac{ν₂^{ν₂/2}}{B(ν₁/2 + k, ν₂/2)} \frac{1}{(ν₂ + ν₁ x)^{(ν₁ + ν₂)/2 + k}},
```

where `B(a, b)` is the beta function. Parameters satisfy `ν₁ > 0`, `ν₂ > 0`, and `λ ≥ 0`.

## Type and Precision

``Distribution/NonCentralF`` is generic over `BinaryFloatingPoint`. Concrete realisations are available for `Float`, `Double`, and (on supported architectures) `Float80`. All evaluations are delegated to Boost.Math via the unified distribution factory; first- and second-order moments are computed analytically in Swift for stability.

## API Overview

- Initialisation: ``Distribution/NonCentralF/init(degreesOfFreedom1:degreesOfFreedom2:nonCentrality:)``
- Density & log-density: ``Distribution/NonCentralF/pdf(_:)``, ``Distribution/NonCentralF/logPdf(_:)``
- Distribution functions: ``Distribution/NonCentralF/cdf(_:)``, ``Distribution/NonCentralF/sf(_:)``
- Inverses: ``Distribution/NonCentralF/quantile(_:)``, ``Distribution/NonCentralF/quantileComplement(_:)``
- Summary statistics: ``Distribution/NonCentralF/mean``, ``Distribution/NonCentralF/variance``,
  ``Distribution/NonCentralF/skewness``, ``Distribution/NonCentralF/kurtosis``,
  ``Distribution/NonCentralF/kurtosisExcess``, ``Distribution/NonCentralF/mode``,
  ``Distribution/NonCentralF/median``
- Hazards: ``Distribution/NonCentralF/hazard(_:)``, ``Distribution/NonCentralF/chf(_:)``
- Divergence: ``Distribution/NonCentralF/klDivergence(relativeTo:options:)``

## Usage

```swift
import SwiftyBoost

let f = try Distribution.NonCentralF<Double>(
    degreesOfFreedom1: 5.0,
    degreesOfFreedom2: 12.0,
    nonCentrality: 3.5
)

let x: Double = 2.0
let density = try f.pdf(x)
let tail = try f.sf(x)
let quantile95 = try f.quantile(0.95)
```

## Mathematical Notes

- Support: `[0, +∞)`
- Mean (ν₂ > 2): `E[X] = ν₂ (ν₁ + λ) / (ν₁ (ν₂ − 2))`
- Variance (ν₂ > 4):
  ```
  Var[X] = 2 (ν₂ / ν₁)^2
           * ((ν₁ + λ)^2 + (ν₁ + 2λ)(ν₂ − 2))
           / ((ν₂ − 2)^2 (ν₂ − 4))
  ```
- Reduces to the central F distribution when `λ = 0`
- No lattice structure: ``Distribution/NonCentralF/latticeStep`` and ``Distribution/NonCentralF/latticeOrigin`` return `nil`
- Higher moments (skewness, kurtosis) are obtained from Boost.Math

## Implementation Details

- Parameter validation is performed in Swift prior to delegating to the dynamic factory.
- The dynamic backend (Boost.Math) supplies PDF/CDF, quantiles, hazards, and higher moments.
- Swift computes mean and variance analytically to match reference formulae and to avoid precision loss when `ν₂` is large.
- Entropy is not currently provided by Boost for the non-central F distribution and surfaces as `nil`.
