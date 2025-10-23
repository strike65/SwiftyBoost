# ``SwiftyBoost/Distribution/Geometric``

Geometric distribution interpreted as the number of failures before the first success in a sequence of independent Bernoulli trials.

## Definition

For success probability `p` with `0 < p ≤ 1`, the probability mass function (PMF) is

```
f(k; p) = (1 - p)^k · p,    for k = 0, 1, 2, …
```

The cumulative distribution function (CDF) can be written as

```
F(k) = 1 - (1 - p)^(k + 1).
```

SwiftyBoost adopts Boost.Math’s “failures before first success” convention. The support is the lattice `{0, 1, 2, …}`.

## Type and Precision

``Distribution/Geometric`` is generic over `BinaryFloatingPoint` with the usual specializations:

- ``Distribution/Geometric``<`Double`>
- ``Distribution/Geometric``<`Float`>
- ``Distribution/Geometric``<`Float80`> (x86_64 only; falls back to the double-precision backend elsewhere)

Each instance delegates to the unified runtime factory ``Distribution/Dynamic`` and reuses a Boost.Math `geometric_distribution` for every evaluation.

## API Overview

- Initialization: ``Distribution/Geometric/init(probabibilityOfSuccess:)``
- Mass & log-mass: ``Distribution/Geometric/pdf(_:)``, ``Distribution/Geometric/logPdf(_:)``
- Distribution functions: ``Distribution/Geometric/cdf(_:)``, ``Distribution/Geometric/sf(_:)``
- Inverse functions: ``Distribution/Geometric/quantile(_:)``, ``Distribution/Geometric/quantileComplement(_:)``
- Moments: ``Distribution/Geometric/mean``, ``Distribution/Geometric/variance``, ``Distribution/Geometric/skewness``, ``Distribution/Geometric/kurtosis``, ``Distribution/Geometric/kurtosisExcess``
- Hazard metrics: ``Distribution/Geometric/hazard(_:)``, ``Distribution/Geometric/chf(_:)``

## Usage

```swift
import SwiftyBoost

// Failures before first success with p = 0.35
let geom = try Distribution.Geometric<Double>(probabibilityOfSuccess: 0.35)

let pmf2 = try geom.pdf(2)          // (1 - p)^2 * p
let cdf3 = try geom.cdf(3)          // 1 - (1 - p)^(3 + 1)
let tail5 = try geom.sf(5)          // (1 - p)^(5 + 1)

let mean = geom.mean                // (1 - p) / p
let variance = geom.variance        // (1 - p) / p^2
```

## Mathematical Notes

- Support: `k ∈ {0, 1, 2, …}`
- Mean: `(1 - p) / p`
- Variance: `(1 - p) / p²`
- Skewness: `(2 - p) / √(1 - p)`
- Excess kurtosis: `6 + p² / (1 - p)`
- Entropy: ``Distribution/Geometric/entropy`` returns `−((1 − p)/p) · log(1 − p) − log p`` (nats).

## Implementation Details

- The Swift wrapper validates probabilities locally before delegating to ``Distribution/Dynamic``.
- All evaluations forward to Boost.Math via the C bridge vtable entries (`bs_geometric_*` thunks).
- Hazard and cumulative hazard follow Boost’s conventions for discrete distributions; for lattice-specific metadata use ``Distribution/Geometric/range`` and ``Distribution/Geometric/supportLowerBound``/``Distribution/Geometric/supportUpperBound``.
