# ``SwiftyBoost/Distribution/ChiSquared``

Chi-squared distribution χ²(ν) with ν > 0 degrees of freedom.

## Definition

Let ν > 0. The density is

f(x; ν) = 1 / [2^{ν/2} Γ(ν/2)] · x^{ν/2 − 1} e^{−x/2}, for x ≥ 0.

Equivalently, χ²(ν) is the distribution of the sum of squares of ν independent standard normal variables, and is a special case of the gamma distribution Γ(k = ν/2, θ = 2).

## Type and Precision

``Distribution/ChiSquared`` is generic over `BinaryFloatingPoint`:

- ``Distribution/ChiSquared``<`Double`>
- ``Distribution/ChiSquared``<`Float`>
- ``Distribution/ChiSquared``<`Float80`> (x86_64 only; otherwise falls back to `Double`)

Instances allocate a Boost.Math `chi_squared_distribution` one time via the dynamic factory.

## API Overview

- Initialization: ``Distribution/ChiSquared/init(degreesOfFreedom:)``
- Density & log-density: ``Distribution/ChiSquared/pdf(_:)``, ``Distribution/ChiSquared/logPdf(_:)``
- Distribution functions: ``Distribution/ChiSquared/cdf(_:)``, ``Distribution/ChiSquared/sf(_:)``
- Inverses: ``Distribution/ChiSquared/quantile(_:)``, ``Distribution/ChiSquared/quantileComplement(_:)``
- Moments: ``Distribution/ChiSquared/mean``, ``Distribution/ChiSquared/variance``, ``Distribution/ChiSquared/mode``, ``Distribution/ChiSquared/skewness``
- Hazard metrics: ``Distribution/ChiSquared/hazard(_:)``, ``Distribution/ChiSquared/chf(_:)``

## Usage

```swift
import SwiftyBoost

let chi2 = try Distribution.ChiSquared<Double>(degreesOfFreedom: 8)
let pdf = try chi2.pdf(4.0)
let sf  = try chi2.sf(12.0)
let q95 = try chi2.quantile(0.95)

let mean = chi2.mean                // ν
let variance = chi2.variance        // 2ν
let mode = chi2.mode                // ν − 2 (for ν ≥ 2)
```

## Mathematical Notes

- Support: x ∈ [0, ∞)
- Mean: ν
- Variance: 2ν
- Mode: max(ν − 2, 0)
- Skewness: √(8 / ν)
- Kurtosis (excess): 12 / ν

## Error Handling

- Initialization throws if ν ≤ 0 or not finite.
- Density and distribution methods throw when `x < 0`.
- Quantile helpers validate probability arguments in `[0, 1]`.

## Implementation Details

- Delegates to ``Distribution/Dynamic`` with the `chisquared` vtable.
- Entropy is not provided by Boost for χ²; the Swift wrapper returns `nil`.
- The dynamic factory accepts aliases `chisquared`, `chi_squared`, `chi2`, `chi-squared`, and `chisquare`, with parameters under `df|nu|degreesOfFreedom`.
