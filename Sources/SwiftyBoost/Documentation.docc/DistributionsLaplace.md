# ``SwiftyBoost/Distribution/Laplace``

Laplace (double exponential) distribution in location–scale form.

## Definition

The Laplace distribution is a symmetric distribution with exponential tails, parameterised by a location
`μ` and a positive scale `b`. Its density can be written as

```
f(x; μ, b) = (1 / (2b)) · exp(−|x − μ| / b)
```

SwiftyBoost delegates all evaluations to Boost.Math’s `laplace_distribution`, exposing the full set of analytic
moments and entropy.

## Type and Precision

``Distribution/Laplace`` is generic over `BinaryFloatingPoint` with concrete support for:

- ``Distribution/Laplace``<`Double`>
- ``Distribution/Laplace``<`Float`>
- ``Distribution/Laplace``<`Float80`> (x86_64 only; falls back to the double-precision backend elsewhere)

## API Overview

- Initialization: ``Distribution/Laplace/init(location:scale:)``
- Parameters: ``Distribution/Laplace/location``, ``Distribution/Laplace/scale``
- Density & log-density: ``Distribution/Laplace/pdf(_:)``, ``Distribution/Laplace/logPdf(_:)``
- Distribution functions: ``Distribution/Laplace/cdf(_:)``, ``Distribution/Laplace/sf(_:)``
- Inverse functions: ``Distribution/Laplace/quantile(_:)``, ``Distribution/Laplace/quantileComplement(_:)``
- Summary statistics: ``Distribution/Laplace/mean``, ``Distribution/Laplace/variance``,
  ``Distribution/Laplace/mode``, ``Distribution/Laplace/median``, ``Distribution/Laplace/skewness``,
  ``Distribution/Laplace/kurtosis``, ``Distribution/Laplace/kurtosisExcess``, ``Distribution/Laplace/entropy``
- Hazard metrics: ``Distribution/Laplace/hazard(_:)``, ``Distribution/Laplace/chf(_:)``

## Usage

```swift
import SwiftyBoost

let laplace = try Distribution.Laplace<Double>(location: 1.0, scale: 0.75)

let x: Double = 0.6
let pdf = try laplace.pdf(x)
let cdf = try laplace.cdf(x)
let sf  = try laplace.sf(x)

let q95 = try laplace.quantile(0.95)
let upperTail = try laplace.quantileComplement(1e-5)

let mean = laplace.mean         // equals location
let variance = laplace.variance // equals 2 * scale^2
let entropy = laplace.entropy   // ln(2 e b)
```

## Mathematical Notes

- Support: `x ∈ (−∞, ∞)`
- Mean, median, and mode all equal the location parameter `μ`.
- Variance equals `2 b²`.
- Skewness is 0; kurtosis is 6; excess kurtosis is 3.
- Differential entropy is `ln(2 e b)` (returned by ``Distribution/Laplace/entropy``).
- Hazard and cumulative hazard satisfy `h(x) = f(x) / S(x)` and `H(x) = −ln S(x)`, both sourced directly from Boost.

## Implementation Details

- The initializer validates `scale > 0` before forwarding to ``Distribution/Dynamic`` under the name `"laplace"`.
- Boost exposes all metrics directly, so the Swift wrapper simply maps optionals using the unified dynamic backend.
- The distribution is continuous; lattice metadata accessors (`latticeStep`, `latticeOrigin`) return `nil`.
