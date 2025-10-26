# ``SwiftyBoost/Distribution/MapAiry``

Map-Airy distribution in location–scale form.

## Definition

The Map-Airy distribution appears in optics and diffraction problems. It is parameterised
by a location `μ` and a scale `c > 0`. The distribution is heavy-tailed: only the first
moment is finite; variance and higher moments diverge.

## Type and Precision

``Distribution/MapAiry`` is provided for:

- ``Distribution/MapAiry``<`Double`>
- ``Distribution/MapAiry``<`Float`>
- ``Distribution/MapAiry``<`Float80`> (x86_64 only; falls back to the double backend elsewhere)

## API Overview

- Initialiser: ``Distribution/MapAiry/init(location:scale:)``
- Density & log-density: ``Distribution/MapAiry/pdf(_:)``, ``Distribution/MapAiry/logPdf(_:)``
- Distribution functions: ``Distribution/MapAiry/cdf(_:)``, ``Distribution/MapAiry/sf(_:)``
- Inverse functions: ``Distribution/MapAiry/quantile(_:)``, ``Distribution/MapAiry/quantileComplement(_:)``
- Summary statistics: ``Distribution/MapAiry/mean``, ``Distribution/MapAiry/variance``,
  ``Distribution/MapAiry/mode``, ``Distribution/MapAiry/median``, ``Distribution/MapAiry/skewness``,
  ``Distribution/MapAiry/kurtosis``, ``Distribution/MapAiry/kurtosisExcess``, ``Distribution/MapAiry/entropy``
- Hazard metrics: ``Distribution/MapAiry/hazard(_:)``, ``Distribution/MapAiry/chf(_:)``

## Usage

```swift
import SwiftyBoost

let mapAiry = try Distribution.MapAiry<Double>(location: 0.5, scale: 1.2)

let x: Double = 0.0
let pdf = try mapAiry.pdf(x)
let cdf = try mapAiry.cdf(x)
let sf  = try mapAiry.sf(x)

let quant = try mapAiry.quantile(0.9)
let entropy = mapAiry.entropy
```

## Mathematical Notes

- Support: `x ∈ (-∞, ∞)`
- Mean equals the location parameter `μ`
- Variance, skewness, and kurtosis diverge (returned as `nil`)
- Entropy is finite and scales with `ln c` for `c > 0`

## Implementation Details

- Swift validates `scale > 0` before dispatching to ``Distribution/Dynamic``.
- Boost.Math supplies PDF/CDF/quantile/mode/median/entropy implementations; unavailable
  statistics are surfaced as `nil`.
- KL divergence is evaluated numerically via the dynamic factory.
