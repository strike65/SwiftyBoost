# ``SwiftyBoost/Distribution/InverseGamma``

Inverse gamma distribution Γ⁻¹(α, β) with strictly positive shape `α` and scale `β`.

## Definition

The inverse gamma distribution is supported on `x ∈ (0, +∞)` with probability density

```
f(x; α, β) = β^α / Γ(α) · x^(−α−1) · exp(−β / x).
```

It appears as the conjugate prior for the variance of a normal likelihood and as the continuous analogue of waiting-time models with heavy tails.

Parameters:

- `shape` (`α > 0`)
- `scale` (`β > 0`, defaults to `1`)

## Type and Precision

``Distribution/InverseGamma`` is generic over `BinaryFloatingPoint`. Concrete specialisations are available for `Float`, `Double`, and (on x86_64) `Float80`. All evaluations delegate to Boost.Math through the unified dynamic distribution factory.

## API Overview

- Initialisation: ``Distribution/InverseGamma/init(shape:scale:)``
- Density & log-density: ``Distribution/InverseGamma/pdf(_:)``,
  ``Distribution/InverseGamma/logPdf(_:)``
- Distribution functions: ``Distribution/InverseGamma/cdf(_:)``,
  ``Distribution/InverseGamma/sf(_:)``
- Inverses: ``Distribution/InverseGamma/quantile(_:)``,
  ``Distribution/InverseGamma/quantileComplement(_:)``
- Summary statistics: ``Distribution/InverseGamma/mean``,
  ``Distribution/InverseGamma/variance``,
  ``Distribution/InverseGamma/mode``,
  ``Distribution/InverseGamma/median``,
  ``Distribution/InverseGamma/skewness``,
  ``Distribution/InverseGamma/kurtosis``,
  ``Distribution/InverseGamma/kurtosisExcess``,
  ``Distribution/InverseGamma/entropy``
- Hazards: ``Distribution/InverseGamma/hazard(_:)``,
  ``Distribution/InverseGamma/chf(_:)``
- Divergence: ``Distribution/InverseGamma/klDivergence(relativeTo:options:)``

## Usage

```swift
import SwiftyBoost

// Shape α = 4.5, scale β = 1.2
let invGamma = try Distribution.InverseGamma<Double>(shape: 4.5, scale: 1.2)

let x: Double = 0.75
let density = try invGamma.pdf(x)
let tail = try invGamma.sf(x)
let quantile90 = try invGamma.quantile(0.90)

// Compare against another configuration via KL divergence
let reference = try Distribution.InverseGamma<Double>(shape: 6.0, scale: 0.9)
let divergence = try invGamma.klDivergence(relativeTo: reference)
```

## Mathematical Notes

- Support: `(0, +∞)`
- Mean exists for `α > 1` and equals `β / (α − 1)`
- Variance exists for `α > 2` and equals `β² / [(α − 1)² (α − 2)]`
- Mode is defined for all `α > 0` as `β / (α + 1)`
- Skewness exists for `α > 3`, excess kurtosis for `α > 4`
- Entropy is currently unavailable from Boost and surfaces as `nil`
- No lattice structure: ``Distribution/InverseGamma/latticeStep`` and
  ``Distribution/InverseGamma/latticeOrigin`` return `nil`

## Implementation Details

- The Swift wrapper validates parameters (finiteness and positivity) before delegating to ``Distribution/Dynamic``.
- Boost.Math provides numerically stable evaluations for PDF/CDF, quantiles, and hazard functions; Swift maps non-finite results to `nil` for optional metrics.
- ``Distribution/InverseGamma/klDivergence(relativeTo:options:)`` reuses the dynamic implementation which integrates when no closed form is available.
