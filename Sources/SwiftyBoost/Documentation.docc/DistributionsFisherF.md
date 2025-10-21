# ``SwiftyBoost/Distribution/FisherF``

Fisher–Snedecor F distribution with degrees of freedom `df1 > 0`, `df2 > 0`.

## Definition

For `x ≥ 0`, the probability density function (PDF) is

f(x; df1, df2) = d1^(d1/2) d2^(d2/2) x^(d1/2 − 1) / ( B(d1/2, d2/2) (d1 x + d2)^((d1 + d2)/2) ),

where `d1 = df1`, `d2 = df2` and `B(·,·)` is the Beta function. The cumulative distribution function (CDF) can be expressed via the regularized incomplete Beta function as

F(x) = I_{ z }(d1/2, d2/2) with z = d1 x / (d1 x + d2).

## Type and Precision

``Distribution/FisherF`` is generic over `BinaryFloatingPoint` and available in these typical specializations:

- ``Distribution/FisherF``<`Double`>
- ``Distribution/FisherF``<`Float`>
- ``Distribution/FisherF``<`Float80`> (x86_64 only; falls back to `Double` elsewhere)

Each instance delegates to the unified runtime distribution vtable via ``Distribution/Dynamic``. Internally, the dynamic wrapper constructs a Boost.Math `fisher_f_distribution` and reuses that backend for all evaluations.

## API Overview

- Initialization: ``Distribution/FisherF/init(degreesOfFreedom1:degreesOfFreedom2:)``
- Density: ``Distribution/FisherF/pdf(_:)`` and ``Distribution/FisherF/logPdf(_:)``
- Distribution functions: ``Distribution/FisherF/cdf(_:)`` and ``Distribution/FisherF/sf(_:)``
- Inverse functions: ``Distribution/FisherF/quantile(_:)`` and ``Distribution/FisherF/quantileComplement(_:)``
- Moments: ``Distribution/FisherF/mean``, ``Distribution/FisherF/variance``, ``Distribution/FisherF/mode``
- Hazards: ``Distribution/FisherF/hazard(_:)`` and ``Distribution/FisherF/chf(_:)``

## Usage

```swift
import SwiftyBoost

let f = try Distribution.FisherF<Double>(degreesOfFreedom1: 10, degreesOfFreedom2: 20)
let pdf2 = try f.pdf(2.0)
let cdf3 = try f.cdf(3.0)
let q95  = try f.quantile(0.95)
```

## Mathematical Notes

- Support: x ∈ [0, ∞)
- Mean: d2/(d2 − 2) for d2 > 2; undefined otherwise
- Variance: 2 d2² (d1 + d2 − 2) / ( d1 (d2 − 2)² (d2 − 4) ) for d2 > 4; ∞ for 2 < d2 ≤ 4; undefined for d2 ≤ 2
- Mode: ( (d1 − 2)/d1 ) · ( d2 / (d2 + 2) ) for d1 > 2; undefined otherwise
- Skewness exists for d2 > 6; kurtosis exists for d2 > 8 (values provided by the backend)

## Error Handling

- Initialization throws if any degree of freedom is not strictly positive or not finite.
- ``Distribution/FisherF/quantile(_:)`` and ``Distribution/FisherF/quantileComplement(_:)`` throw if their probability arguments are outside `[0, 1]`.

## Implementation Details

All operations are routed through ``Distribution/Dynamic`` which wraps a small C vtable over Boost.Math. The dynamic layer creates a `fisher_f_distribution` in the matching precision and exposes stable function pointers for PDF/CDF/SF, quantiles, hazards, and moments. Entropy is provided by a Swift-side fallback in the dynamic wrapper because Boost does not expose it for F.
