#
``SwiftyBoost/Distribution/NonCentralChiSquared``

Non-central chi-squared distribution with degrees of freedom `ν` and non-centrality parameter `λ`.

## Definition

The non-central chi-squared distribution generalises the chi-squared law to allow
non-zero means in the underlying Gaussian components. It is supported on `x ∈ [0, +∞)`
with density

```
f(x; ν, λ) = e^{-(x + λ)/2} \left(\frac{x}{λ}\right)^{ν/4 − 1/2} I_{ν/2 - 1}(\sqrt{λ x}),
```

where `I_k` denotes the modified Bessel function of the first kind. The parameter
constraints are `ν > 0` and `λ ≥ 0`.

## Type and Precision

``Distribution/NonCentralChiSquared`` is generic over `BinaryFloatingPoint`
(`Float`, `Double`, or `Float80` when available). All functionality delegates to
Boost.Math through the dynamic distribution factory, while selected summary
statistics are computed analytically in Swift for improved stability.

## API Overview

- Initialisation: ``Distribution/NonCentralChiSquared/init(degreesOfFreedom:nonCentrality:)``
- Density & log-density: ``Distribution/NonCentralChiSquared/pdf(_:)``,
  ``Distribution/NonCentralChiSquared/logPdf(_:)``
- Distribution functions: ``Distribution/NonCentralChiSquared/cdf(_:)``,
  ``Distribution/NonCentralChiSquared/sf(_:)``
- Inverses: ``Distribution/NonCentralChiSquared/quantile(_:)``,
  ``Distribution/NonCentralChiSquared/quantileComplement(_:)``
- Summary statistics: ``Distribution/NonCentralChiSquared/mean``,
  ``Distribution/NonCentralChiSquared/variance``,
  ``Distribution/NonCentralChiSquared/skewness``,
  ``Distribution/NonCentralChiSquared/kurtosis``,
  ``Distribution/NonCentralChiSquared/kurtosisExcess``,
  ``Distribution/NonCentralChiSquared/mode``,
  ``Distribution/NonCentralChiSquared/median``
- Hazards: ``Distribution/NonCentralChiSquared/hazard(_:)``,
  ``Distribution/NonCentralChiSquared/chf(_:)``
- Divergence: ``Distribution/NonCentralChiSquared/klDivergence(relativeTo:options:)``

## Usage

```swift
import SwiftyBoost

let chi = try Distribution.NonCentralChiSquared<Double>(
    degreesOfFreedom: 6.0,
    nonCentrality: 3.5
)

let x: Double = 4.0
let density = try chi.pdf(x)
let tail = try chi.sf(x)
let quantile90 = try chi.quantile(0.90)

if let skew = chi.skewness {
    print("Skewness:", skew)
}
```

## Mathematical Notes

- Support: `[0, +∞)`
- Mean: `ν + λ`
- Variance: `2(ν + 2λ)`
- Skewness: `2√2 (ν + 3λ) / (ν + 2λ)^{3/2}`
- Excess kurtosis: `12(ν + 4λ) / (ν + 2λ)^2`
- Reduces to the central chi-squared distribution when `λ = 0`
- No lattice structure: ``Distribution/NonCentralChiSquared/latticeStep`` and
  ``Distribution/NonCentralChiSquared/latticeOrigin`` return `nil`

## Implementation Details

- Parameter validation is performed in Swift prior to constructing the dynamic delegate.
- PDF/CDF, quantiles, hazards, and KL divergence are delegated to Boost.Math.
- Analytic expressions for the first four central moments are evaluated in Swift to
  provide high-precision mean, variance, skewness, and kurtosis.
- Entropy is not currently supplied by Boost.Math and surfaces as `nil`.
