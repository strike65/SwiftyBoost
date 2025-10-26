# ``SwiftyBoost/Distribution/KolmogorovSmirnov``

Kolmogorov–Smirnov distribution for the one-sample statistic with a finite observation count.

## Definition

``Distribution/KolmogorovSmirnov`` models the distribution of the Kolmogorov–Smirnov test statistic
``Dₙ = sup_x |Fₙ(x) - F(x)|`` for a sample of size `n`. The single shape parameter `n > 0` represents the
number of observations in the empirical sample.

Boost.Math supplies numerically stable implementations for the PDF, CDF, quantiles, and analytic
moments up to kurtosis. SwiftyBoost reuses those routines directly through the unified dynamic factory.

## Type and Precision

Generic over any `BinaryFloatingPoint` satisfying the module constraints:

- ``Distribution/KolmogorovSmirnov``<`Double`>
- ``Distribution/KolmogorovSmirnov``<`Float`>
- ``Distribution/KolmogorovSmirnov``<`Float80`> (x86_64 only; falls back to the double-precision backend elsewhere)

## API Overview

- Initialization: ``Distribution/KolmogorovSmirnov/init(numberOfObservations:)``
- Stored parameter: ``Distribution/KolmogorovSmirnov/numberOfObservations``
- Density & log-density: ``Distribution/KolmogorovSmirnov/pdf(_:)``, ``Distribution/KolmogorovSmirnov/logPdf(_:)``
- Distribution functions: ``Distribution/KolmogorovSmirnov/cdf(_:)``, ``Distribution/KolmogorovSmirnov/sf(_:)``
- Inverse functions: ``Distribution/KolmogorovSmirnov/quantile(_:)``, ``Distribution/KolmogorovSmirnov/quantileComplement(_:)``
- Summary statistics: ``Distribution/KolmogorovSmirnov/mean``, ``Distribution/KolmogorovSmirnov/variance``,
  ``Distribution/KolmogorovSmirnov/skewness``, ``Distribution/KolmogorovSmirnov/kurtosis``, ``Distribution/KolmogorovSmirnov/mode``, ``Distribution/KolmogorovSmirnov/median``
- Hazard metrics: ``Distribution/KolmogorovSmirnov/hazard(_:)``, ``Distribution/KolmogorovSmirnov/chf(_:)``

## Usage

```swift
import SwiftyBoost

let ks = try Distribution.KolmogorovSmirnov<Double>(numberOfObservations: 50)

// Evaluate density / distribution functions
let x: Double = 0.25
let pdf = try ks.pdf(x)
let cdf = try ks.cdf(x)
let sf  = try ks.sf(x)

// Invert probabilities
let q95 = try ks.quantile(0.95)
let upper = try ks.quantileComplement(1e-6)

// Analytic statistics
let mean = ks.mean
let variance = ks.variance
let skew = ks.skewness
let kurtosis = ks.kurtosis
```

## Mathematical Notes

- Support: `x ∈ [0, ∞)`
- Mean: `E[Dₙ] = √(π / 2) · ln 2 / √n`
- Variance: `Var[Dₙ] = (π² / 6 − π ln² 2) / (2n)`
- Skewness and kurtosis are finite and sourced directly from Boost.Math’s analytic expressions.
- Median is computed via the 0.5 quantile (Boost does not expose a closed-form median).

## Implementation Details

- Initialization validates `n > 0` and delegates to ``Distribution/Dynamic`` with the canonical parameter key `n`.
- The bridged vtable exposes PDF/CDF/SF/quantiles/moments; hazard and cumulative hazard are derived from those pointers.
- Moments are returned as optionals; SwiftyBoost reports `nil` if Boost indicates divergence or returns a non-finite value.
