#
``SwiftyBoost/Distribution/Pareto``

Pareto (Type I) distribution with scale parameter `xₘ` and shape parameter `α`.

## Definition

The density is supported on `x ≥ xₘ` and given by

```
f(x; xₘ, α) = α xₘ^{α} / x^{α + 1}.
```

Parameters satisfy `xₘ > 0` and `α > 0`.

## Type and Precision

``Distribution/Pareto`` is generic over `BinaryFloatingPoint`. Concrete instances exist for `Float`, `Double`, and (on supported architectures) `Float80`. All operations delegate to Boost.Math through the unified dynamic distribution factory, while select summary statistics are computed analytically in Swift.

## API Overview

- Initialisation: ``Distribution/Pareto/init(scale:shape:)``
- Density & log-density: ``Distribution/Pareto/pdf(_:)``, ``Distribution/Pareto/logPdf(_:)``
- CDF / SF: ``Distribution/Pareto/cdf(_:)``, ``Distribution/Pareto/sf(_:)``
- Inverses: ``Distribution/Pareto/quantile(_:)``, ``Distribution/Pareto/quantileComplement(_:)``
- Summary statistics: ``Distribution/Pareto/mean``, ``Distribution/Pareto/variance``, ``Distribution/Pareto/mode``, ``Distribution/Pareto/median``, ``Distribution/Pareto/skewness``, ``Distribution/Pareto/kurtosis``
- Hazards: ``Distribution/Pareto/hazard(_:)``, ``Distribution/Pareto/chf(_:)``
- Divergence: ``Distribution/Pareto/klDivergence(relativeTo:options:)``

## Usage

```swift
import SwiftyBoost

let pareto = try Distribution.Pareto<Double>(scale: 2.0, shape: 4.0)
let x: Double = 3.0
let density = try pareto.pdf(x)
let tail = try pareto.sf(x)
let median = pareto.median
```

## Mathematical Notes

- Support: `[xₘ, +∞)`
- Mean exists for `α > 1` with `E[X] = α xₘ / (α − 1)`
- Variance exists for `α > 2` with `Var[X] = α xₘ² / ((α − 1)² (α − 2))`
- Mode equals the scale parameter `xₘ`
- Skewness and kurtosis are supplied by Boost.Math
- No lattice structure (continuous distribution)

## Implementation Details

- Swift validates parameters (positive scale and shape) before delegating to Boost.
- PDF/CDF/quantiles/hazards are handled by Boost.Math via the dynamic factory.
- Mean and variance are computed analytically in Swift for robustness.
- Entropy is currently unavailable from Boost and surfaces as `nil`.
