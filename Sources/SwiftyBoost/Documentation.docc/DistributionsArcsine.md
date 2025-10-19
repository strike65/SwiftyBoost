# ``SwiftyBoost/Distribution/Arcsine``

Arcsine distribution on a finite interval `[a, b]` with `a < b`.

## Definition

For `x ∈ (a, b)`, the probability density function (PDF) is

f(x; a, b) = 1 / (π √((x − a)(b − x))).

The cumulative distribution function (CDF) on `[a, b]` can be written in terms of the standard arcsine CDF via the affine map `y = (x − a)/(b − a)`:

F(x) = (2/π) arcsin( √( (x − a)/(b − a) ) ).

## Type and Precision

``Distribution/Arcsine`` is generic over `BinaryFloatingPoint` and available in these typical specializations:

- ``Distribution/Arcsine``<`Double`>
- ``Distribution/Arcsine``<`Float`>
- ``Distribution/Arcsine``<`Float80`> (x86_64 only; falls back to `Double` elsewhere)

Each instance constructs a Boost.Math `arcsine_distribution` once and keeps an internal opaque handle, reused for all evaluations.

## API Overview

- Initialization: ``Distribution/Arcsine/init(minX:maxX:)``
- Density: ``Distribution/Arcsine/pdf(_:)`` and ``Distribution/Arcsine/logPdf(_:)``
- Distribution functions: ``Distribution/Arcsine/cdf(_:)`` and ``Distribution/Arcsine/sf(_:)``
- Inverse functions: ``Distribution/Arcsine/quantile(_:)`` and ``Distribution/Arcsine/quantileComplement(_:)``
- Moments: ``Distribution/Arcsine/mean``, ``Distribution/Arcsine/variance``, ``Distribution/Arcsine/median``
- Hazards: ``Distribution/Arcsine/hazard(_:)`` and ``Distribution/Arcsine/chf(_:)``

## Usage

```swift
import SwiftyBoost

let a = try Distribution.Arcsine<Double>(minX: 0, maxX: 1)
let p2 = try a.pdf(0.2)
let F8 = try a.cdf(0.8)
let q5 = try a.quantile(0.5) // equals median
```

## Mathematical Notes

- Support: x ∈ [a, b]
- Mean: (a + b)/2
- Variance: (b − a)² / 8
- Median: (a + b)/2
- Mode: Density is unbounded at the endpoints; some conventions list `a` and `b` as modes. Backend provides a conventional value.

## Error Handling

- Initialization throws if `a ≥ b` or either bound is not finite.
- ``Distribution/Arcsine/quantile(_:)`` and ``Distribution/Arcsine/quantileComplement(_:)`` throw if their probability arguments are outside `[0, 1]`.

## Implementation Details

All operations delegate to Boost via the C bridge. The Swift type holds an opaque pointer to a `boost::math::arcsine_distribution` in the matching precision. Calls route through stable `bs_` functions (e.g., `bs_arcsine_pdf_h`, `bs_arcsine_cdf_h`, `bs_arcsine_quantile_h`).

