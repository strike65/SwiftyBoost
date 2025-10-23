# ``SwiftyBoost/Distribution/Holtsmark``

Holtsmark distribution (stable law with stability index α = 3/2) in location–scale form.

## Definition

The Holtsmark distribution is a symmetric, heavy-tailed stable distribution. In SwiftyBoost it is parameterized by location `μ` and scale `σ > 0`. The probability density function does not admit a simple elementary expression; the implementation delegates to Boost.Math’s numerical routines for PDF, CDF, and quantile evaluations.

## Type and Precision

``Distribution/Holtsmark`` is generic over `BinaryFloatingPoint` with the following specializations:

- ``Distribution/Holtsmark``<`Double`>
- ``Distribution/Holtsmark``<`Float`>
- ``Distribution/Holtsmark``<`Float80`> (x86_64 only; falls back to the double-precision backend elsewhere)

Each instance reuses the unified runtime factory ``Distribution/Dynamic`` to bridge into Boost.Math’s `holtsmark_distribution`.

## API Overview

- Initialization: ``Distribution/Holtsmark/init(location:scale:)``
- Density & log-density: ``Distribution/Holtsmark/pdf(_:)``, ``Distribution/Holtsmark/logPdf(_:)``
- Distribution functions: ``Distribution/Holtsmark/cdf(_:)``, ``Distribution/Holtsmark/sf(_:)``
- Inverse functions: ``Distribution/Holtsmark/quantile(_:)``, ``Distribution/Holtsmark/quantileComplement(_:)``
- Summary statistics: ``Distribution/Holtsmark/mean``, ``Distribution/Holtsmark/variance``, ``Distribution/Holtsmark/mode``, ``Distribution/Holtsmark/median``
- Hazard metrics: ``Distribution/Holtsmark/hazard(_:)``, ``Distribution/Holtsmark/chf(_:)``
- Differential entropy: ``Distribution/Holtsmark/entropy``

## Usage

```swift
import SwiftyBoost

let holts = try Distribution.Holtsmark<Double>(location: 0, scale: 1)

let x: Double = 0.75
let pdf = try holts.pdf(x)
let cdf = try holts.cdf(x)
let sf  = try holts.sf(x)

let q95 = try holts.quantile(0.95)
let mean = holts.mean         // equals location
let entropy = holts.entropy   // differential entropy (nats)
```

## Mathematical Notes

- Support: `x ∈ (−∞, ∞)`
- Mean: equal to the location parameter `μ`
- Variance: infinite (``Distribution/Holtsmark/variance`` returns `nil`)
- Median & mode: equal to `μ`
- Entropy: finite
- Higher-order moments (skewness, kurtosis) diverge; corresponding accessors return `nil`

## Implementation Details

- The Swift wrapper validates `scale > 0` before delegating to ``Distribution/Dynamic``.
- Boost.Math supplies numerically stable implementations for PDF, CDF, quantiles, and entropy; unavailable moments are surfaced as `nil` in Swift.
- As with other location–scale families, hazard and cumulative hazard are computed directly through the bridged vtable.
