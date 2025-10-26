# ``SwiftyBoost/Distribution/NegativeBinomial``

Negative binomial (Pascal) distribution counting failures before a fixed number of successes.

## Definition

Given independent Bernoulli trials with success probability `p`, the negative binomial distribution
models the count of failures observed before `r` successes occur. SwiftyBoost follows the Boost.Math
parameterisation with parameters:

- `successes = r > 0`
- `probabilityOfSuccess = p`, with `0 < p ≤ 1`

## Type and Precision

``Distribution/NegativeBinomial`` is generic over `BinaryFloatingPoint`:

- ``Distribution/NegativeBinomial``<`Double`>
- ``Distribution/NegativeBinomial``<`Float`>
- ``Distribution/NegativeBinomial``<`Float80`> (x86_64 only; falls back to the double backend elsewhere)

## API Overview

- Initialiser: ``Distribution/NegativeBinomial/init(successes:probabilityOfSuccess:)``
- Density & log-density: ``Distribution/NegativeBinomial/pdf(_:)``, ``Distribution/NegativeBinomial/logPdf(_:)``
- Distribution functions: ``Distribution/NegativeBinomial/cdf(_:)``, ``Distribution/NegativeBinomial/sf(_:)``
- Inverse functions: ``Distribution/NegativeBinomial/quantile(_:)``, ``Distribution/NegativeBinomial/quantileComplement(_:)``
- Summary statistics: ``Distribution/NegativeBinomial/mean``, ``Distribution/NegativeBinomial/variance``,
  ``Distribution/NegativeBinomial/mode``, ``Distribution/NegativeBinomial/median``,
  ``Distribution/NegativeBinomial/skewness``, ``Distribution/NegativeBinomial/kurtosis``,
  ``Distribution/NegativeBinomial/kurtosisExcess``, ``Distribution/NegativeBinomial/entropy``
- Hazard metrics: ``Distribution/NegativeBinomial/hazard(_:)``, ``Distribution/NegativeBinomial/chf(_:)``

## Usage

```swift
import SwiftyBoost

let negBin = try Distribution.NegativeBinomial<Double>(successes: 3.5, probabilityOfSuccess: 0.4)

let k: Double = 5
let pmf = try negBin.pdf(k)
let cdf = try negBin.cdf(k)
let q90 = try negBin.quantile(0.90)

let mean = negBin.mean
let variance = negBin.variance
```

## Mathematical Notes

- Support: `k ∈ {0, 1, 2, …}`
- Mean: `r (1 − p) / p`
- Variance: `r (1 − p) / p²`
- Skewness: `(2 − p) / √(r (1 − p))`
- Kurtosis: `3 + 6 / r + (p² / (r (1 − p)))`
- Excess kurtosis: `(6 − p (6 − p)) / (r (1 − p))`
- Entropy is not provided analytically by Boost and is currently unavailable.

## Implementation Details

- The Swift wrapper enforces `successes > 0` and `0 < probabilityOfSuccess ≤ 1` before delegating to
  ``Distribution/Dynamic``.
- Boost.Math exposes all required statistics except entropy; the wrapper forwards optional values
  directly.
- The distribution is discrete with lattice origin 0 and step 1; ``Distribution/Dynamic`` captures the
  lattice metadata automatically.
