# ``SwiftyBoost/Distribution/InverseNormal``

Inverse normal distribution (also known as the inverse Gaussian or Wald distribution) with positive mean `μ` and shape parameter `λ`.

## Definition

The inverse normal distribution is supported on `x ∈ (0, +∞)` with density

```
f(x; μ, λ) = √(λ / (2π x^3)) · exp(−λ (x − μ)^2 / (2 μ^2 x)).
```

It models first-passage times of Brownian motion with drift and appears in reliability analysis and Tweedie compound families.

Parameters:

- `mean` (`μ > 0`)
- `shape` (`λ > 0`)

## Type and Precision

``Distribution/InverseNormal`` is generic over `BinaryFloatingPoint` with concrete specialisations for `Float`, `Double`, and (on x86_64) `Float80`. All computations delegate to Boost.Math via the dynamic distribution factory. ``Distribution/InverseGaussian`` is provided as a convenience typealias.

## API Overview

- Initialisation: ``Distribution/InverseNormal/init(mean:shape:)``
- Density & log-density: ``Distribution/InverseNormal/pdf(_:)``,
  ``Distribution/InverseNormal/logPdf(_:)``
- Distribution functions: ``Distribution/InverseNormal/cdf(_:)``,
  ``Distribution/InverseNormal/sf(_:)``
- Inverses: ``Distribution/InverseNormal/quantile(_:)``,
  ``Distribution/InverseNormal/quantileComplement(_:)``
- Summary statistics: ``Distribution/InverseNormal/mean``,
  ``Distribution/InverseNormal/variance``,
  ``Distribution/InverseNormal/mode``,
  ``Distribution/InverseNormal/median``,
  ``Distribution/InverseNormal/skewness``,
  ``Distribution/InverseNormal/kurtosis``,
  ``Distribution/InverseNormal/kurtosisExcess``,
  ``Distribution/InverseNormal/entropy``
- Hazards: ``Distribution/InverseNormal/hazard(_:)``,
  ``Distribution/InverseNormal/chf(_:)``
- Divergence: ``Distribution/InverseNormal/klDivergence(relativeTo:options:)``

## Usage

```swift
import SwiftyBoost

// Mean μ = 1.25, shape λ = 3.75
let inverseNormal = try Distribution.InverseNormal<Double>(mean: 1.25, shape: 3.75)

let x: Double = 0.9
let density = try inverseNormal.pdf(x)
let tail = try inverseNormal.sf(x)
let median = inverseNormal.median
let q975 = try inverseNormal.quantile(0.975)

// KL divergence against another configuration
let reference = try Distribution.InverseNormal<Double>(mean: 1.2, shape: 2.5)
let divergence = try inverseNormal.klDivergence(relativeTo: reference)
```

## Mathematical Notes

- Support: `(0, +∞)`
- Mean equals `μ` for all valid parameters
- Variance equals `μ³ / λ`
- Mode equals `μ (√(1 + 9 μ² / (4 λ²)) − 3 μ / (2 λ))`
- Skewness equals `3 √(μ / λ)`, excess kurtosis equals `15 μ / λ`
- Entropy is currently unavailable from Boost and surfaces as `nil`
- No lattice structure: ``Distribution/InverseNormal/latticeStep`` and
  ``Distribution/InverseNormal/latticeOrigin`` return `nil`

## Implementation Details

- Parameter validation ensures both `mean` and `shape` are finite and strictly positive before constructing the dynamic handle.
- Boost.Math provides PDF/CDF/quantile/hazard implementations; Swift maps non-finite summaries to `nil`.
- ``Distribution/InverseNormal/klDivergence(relativeTo:options:)`` forwards to the dynamic implementation, which integrates when no analytic form is available.
- ``Distribution/InverseGaussian`` is a typealias of ``Distribution/InverseNormal`` for codebases that prefer the inverse Gaussian terminology.
- Boost.Math reports ``Distribution/InverseNormal/kurtosis`` as `15 μ / λ − 3` (the same value returned by ``Distribution/InverseNormal/kurtosisExcess`` minus 3); the wrapper surfaces the Boost definition unchanged.
