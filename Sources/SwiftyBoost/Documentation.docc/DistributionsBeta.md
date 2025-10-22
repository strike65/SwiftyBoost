# ``SwiftyBoost/Distribution/Beta``

Beta Distribution B(α, β) with shape parameters α > 0 and β > 0 on the unit interval.

## Definition

Given α > 0 and β > 0, the probability density function is

f(x; α, β) = x^(α−1) (1 − x)^(β−1) / B(α, β) for x ∈ (0, 1),

where B(α, β) is the beta function.

## Type and Precision

``Distribution/Beta`` is generic over `BinaryFloatingPoint` and ships in the usual specializations:

- ``Distribution/Beta``<`Double`>
- ``Distribution/Beta``<`Float`>
- ``Distribution/Beta``<`Float80`> (x86_64 only; falls back to `Double` elsewhere)

Instances allocate a Boost.Math `beta_distribution` once through the dynamic factory and reuse it for every evaluation.

## API Overview

- Initialization: ``Distribution/Beta/init(alpha:beta:)``
- Density & log-density: ``Distribution/Beta/pdf(_:)``, ``Distribution/Beta/logPdf(_:)``
- Distribution functions: ``Distribution/Beta/cdf(_:)``, ``Distribution/Beta/sf(_:)``
- Inverses: ``Distribution/Beta/quantile(_:)``, ``Distribution/Beta/quantileComplement(_:)``
- Moments: ``Distribution/Beta/mean``, ``Distribution/Beta/variance``, ``Distribution/Beta/mode``, ``Distribution/Beta/skewness``
- Hazard metrics: ``Distribution/Beta/hazard(_:)``, ``Distribution/Beta/chf(_:)``

## Usage

```swift
import SwiftyBoost

let beta = try Distribution.Beta<Double>(alpha: 2.0, beta: 5.0)
let pdf = try beta.pdf(0.3)
let cdf = try beta.cdf(0.8)
let q95 = try beta.quantile(0.95)

let mean = beta.mean                  // α / (α + β)
let variance = beta.variance          // αβ / [(α + β)² (α + β + 1)]
let mode = beta.mode                  // (α − 1) / (α + β − 2) for α, β > 1
```

## Mathematical Notes

- Support: x ∈ [0, 1]
- Mean: α / (α + β)
- Variance: α β / [(α + β)² (α + β + 1)]
- Mode: (α − 1) / (α + β − 2) when α, β > 1; undefined otherwise
- Skewness: 2(β − α) √(α + β + 1) / [(α + β + 2) √(α β)]

## Error Handling

- Initialization throws if α ≤ 0 or β ≤ 0, or if parameters are not finite.
- Methods that accept `x` outside [0, 1] throw the underlying ``DistributionError``.
- Quantile helpers validate probability arguments in `[0, 1]`.

## Implementation Details

- Delegates to ``Distribution/Dynamic`` with the `beta` vtable in `CBoostBridge`.
- Entropy is not exposed in Boost’s beta distribution; the Swift wrapper returns `nil`.
- Parameter aliases include `alpha|a|p|shape1` and `beta|b|q|shape2` when using the dynamic factory.
