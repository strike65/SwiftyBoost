# ``SwiftyBoost/Distribution/Landau``

Landau distribution (strictly stable with heavy tails) in location–scale form.

## Definition

The Landau distribution is a heavy-tailed, strictly stable probability distribution. In SwiftyBoost it is
parameterised by a location `μ` and a strictly positive scale `c`. The distribution has no finite algebraic
moments—mean, variance, skewness, and kurtosis diverge—but its median, mode, and differential entropy are
well defined.

## Type and Precision

``Distribution/Landau`` is generic over `BinaryFloatingPoint` with concrete support for:

- ``Distribution/Landau``<`Double`>
- ``Distribution/Landau``<`Float`>
- ``Distribution/Landau``<`Float80`> (x86_64 only; falls back to the double-precision backend elsewhere)

Every specialization reuses ``Distribution/Dynamic`` with Boost.Math’s `landau_distribution`.

## API Overview

- Initialization: ``Distribution/Landau/init(location:scale:)``
- Location/scale accessors: ``Distribution/Landau/location``, ``Distribution/Landau/scale``
- Density & log-density: ``Distribution/Landau/pdf(_:)``, ``Distribution/Landau/logPdf(_:)``
- Distribution functions: ``Distribution/Landau/cdf(_:)``, ``Distribution/Landau/sf(_:)``
- Inverse functions: ``Distribution/Landau/quantile(_:)``, ``Distribution/Landau/quantileComplement(_:)``
- Summary statistics: ``Distribution/Landau/mode``, ``Distribution/Landau/median``, ``Distribution/Landau/entropy``
- Hazard metrics: ``Distribution/Landau/hazard(_:)``, ``Distribution/Landau/chf(_:)``

## Usage

```swift
import SwiftyBoost

let landau = try Distribution.Landau<Double>(location: 0.0, scale: 1.5)

let x: Double = 2.0
let pdf = try landau.pdf(x)
let cdf = try landau.cdf(x)
let sf  = try landau.sf(x)

let q90 = try landau.quantile(0.90)
let tail = try landau.quantileComplement(1e-4)

let mode = landau.mode
let median = landau.median
let entropy = landau.entropy
```

## Mathematical Notes

- Support: `x ∈ (−∞, ∞)`
- Mean, variance, skewness, and kurtosis diverge; the corresponding accessors return `nil`.
- Median and mode follow location shifts: if ``Distribution/Landau/location`` changes by Δ, both quantities shift by Δ.
- Differential entropy exists and scales as `H(μ, c) = H(0, 1) + ln c`.
- Hazard and cumulative hazard are exposed via the Boost vtable.

## Implementation Details

- The initializer enforces `scale > 0` before delegating to ``Distribution/Dynamic``.
- The C bridge wires Boost’s PDF/CDF/SF/quantile/mode/median/entropy functions; absent metrics are surfaced as `nil`.
- Entropy comes directly from Boost.Math (`landau_entropy_imp`) and respects the expected log-scale shift.
