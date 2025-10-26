# ``SwiftyBoost/Distribution/LogNormal``

Log-normal distribution parameterised by the location and scale of the underlying normal.

## Definition

If `Z ~ N(μ, σ²)` then `X = exp(Z)` follows a log-normal distribution with parameters
`location = μ` and `scale = σ`. The PDF is

```
f(x; μ, σ) = \frac{1}{x σ \sqrt{2π}} \exp\!\left(-\frac{(\ln x - μ)^2}{2σ^2}\right),  x > 0.
```

## Type and Precision

``Distribution/LogNormal`` is available for:

- ``Distribution/LogNormal``<`Double`>
- ``Distribution/LogNormal``<`Float`>
- ``Distribution/LogNormal``<`Float80`> (x86_64 only; falls back to the double backend elsewhere)

## API Overview

- Initialiser: ``Distribution/LogNormal/init(location:scale:)``
- Density & log-density: ``Distribution/LogNormal/pdf(_:)``, ``Distribution/LogNormal/logPdf(_:)``
- Distribution functions: ``Distribution/LogNormal/cdf(_:)``, ``Distribution/LogNormal/sf(_:)``
- Inverse functions: ``Distribution/LogNormal/quantile(_:)``, ``Distribution/LogNormal/quantileComplement(_:)``
- Summary statistics: ``Distribution/LogNormal/mean``, ``Distribution/LogNormal/variance``,
  ``Distribution/LogNormal/mode``, ``Distribution/LogNormal/median``, ``Distribution/LogNormal/skewness``,
  ``Distribution/LogNormal/kurtosis``, ``Distribution/LogNormal/kurtosisExcess``, ``Distribution/LogNormal/entropy``
- Hazard metrics: ``Distribution/LogNormal/hazard(_:)``, ``Distribution/LogNormal/chf(_:)``

## Usage

```swift
import SwiftyBoost

let logNormal = try Distribution.LogNormal<Double>(location: 0.2, scale: 0.6)

let x: Double = 1.5
let pdf = try logNormal.pdf(x)
let cdf = try logNormal.cdf(x)
let quant = try logNormal.quantile(0.95)

let mean = logNormal.mean
let variance = logNormal.variance
let entropy = logNormal.entropy
```

## Mathematical Notes

- Support: `x ∈ (0, ∞)`
- Mean: `exp(μ + σ² / 2)`
- Variance: `(exp(σ²) - 1) · exp(2μ + σ²)`
- Mode: `exp(μ - σ²)`, Median: `exp(μ)`
- Skewness: `(exp(σ²) + 2) √(exp(σ²) - 1)`
- Kurtosis: `exp(4σ²) + 2 exp(3σ²) + 3 exp(2σ²) - 6`
- Differential entropy: `μ + 1/2 + ln(σ √(2π))`

## Implementation Details

- Parameter check ensures `scale > 0` before invoking ``Distribution/Dynamic``.
- Boost.Math supplies all analytic statistics; the Swift wrapper returns them directly.
- KL divergence relies on the shared dynamic quadrature infrastructure.
