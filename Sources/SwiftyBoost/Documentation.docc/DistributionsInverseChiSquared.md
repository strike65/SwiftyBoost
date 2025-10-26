# ``SwiftyBoost/Distribution/InverseChiSquared``

Inverse chi-squared distribution with degrees of freedom `ν` and scale parameter `τ`.

## Definition

The scaled inverse chi-squared distribution is supported on `x ∈ (0, +∞)` with
probability density function

```
f(x; ν, τ) = ((ν τ / 2)^(ν/2) / Γ(ν / 2)) · x^(−ν/2 − 1) · exp(−(ν τ)/(2x)).
```

It is the conjugate prior for a normal variance with known mean and arises as a
special case of the inverse-gamma distribution.

Parameters:

- `degreesOfFreedom` (`ν > 0`)
- `scale` (`τ > 0`, defaults to `1 / ν`)

## Type and Precision

``Distribution/InverseChiSquared`` is generic over `BinaryFloatingPoint`. Concrete
specialisations are available for `Float`, `Double`, and (on x86_64) `Float80`.
All operations delegate to Boost.Math through the dynamic distribution factory.

## API Overview

- Initialisation: ``Distribution/InverseChiSquared/init(degreesOfFreedom:scale:)``
- Density & log-density: ``Distribution/InverseChiSquared/pdf(_:)``,
  ``Distribution/InverseChiSquared/logPdf(_:)``
- Distribution functions: ``Distribution/InverseChiSquared/cdf(_:)``,
  ``Distribution/InverseChiSquared/sf(_:)``
- Inverses: ``Distribution/InverseChiSquared/quantile(_:)``,
  ``Distribution/InverseChiSquared/quantileComplement(_:)``
- Summary statistics: ``Distribution/InverseChiSquared/mean``,
  ``Distribution/InverseChiSquared/variance``,
  ``Distribution/InverseChiSquared/mode``,
  ``Distribution/InverseChiSquared/median``,
  ``Distribution/InverseChiSquared/skewness``,
  ``Distribution/InverseChiSquared/kurtosis``,
  ``Distribution/InverseChiSquared/kurtosisExcess``
- Hazards: ``Distribution/InverseChiSquared/hazard(_:)``,
  ``Distribution/InverseChiSquared/chf(_:)``
- Divergence: ``Distribution/InverseChiSquared/klDivergence(relativeTo:options:)``

## Usage

```swift
import SwiftyBoost

// Degrees of freedom ν = 8, scale τ = 0.5
let invChiSq = try Distribution.InverseChiSquared<Double>(
    degreesOfFreedom: 8,
    scale: 0.5
)

let x: Double = 0.75
let density = try invChiSq.pdf(x)
let tail = try invChiSq.sf(x)
let quantile90 = try invChiSq.quantile(0.90)

// Relative entropy against another configuration
let reference = try Distribution.InverseChiSquared<Double>(
    degreesOfFreedom: 10,
    scale: 0.4
)
let divergence = try invChiSq.klDivergence(relativeTo: reference)
```

## Mathematical Notes

- Support: `(0, +∞)`
- Mean exists for `ν > 2` and equals `τ ν / (ν − 2)`
- Variance exists for `ν > 4` and equals `2 τ² ν² / ((ν − 2)² (ν − 4))`
- The distribution is heavy-tailed; higher-order moments require progressively larger `ν`
- No lattice structure: ``Distribution/InverseChiSquared/latticeStep`` and
  ``Distribution/InverseChiSquared/latticeOrigin`` return `nil`

## Implementation Details

- The Swift wrapper validates input parameters before delegating to
  ``Distribution/Dynamic``.
- Boost.Math provides numerically stable evaluations for PDF/CDF,
  quantiles, and derivatives; entropy is currently unavailable and surfaces as `nil`.
- `klDivergence(relativeTo:options:)` uses adaptive quadrature for continuous
  comparisons with tunable options via ``Distribution/KLDivergenceOptions``.
