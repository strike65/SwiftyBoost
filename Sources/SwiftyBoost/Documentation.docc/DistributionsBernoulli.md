# ``SwiftyBoost/Distribution/Bernoulli``

Bernoulli distribution modelling a single binary trial with success probability `p`.

## Definition

For `0 ≤ p ≤ 1`, the probability mass function is

- P(X = 1) = p  
- P(X = 0) = 1 − p

The support is the discrete set {0, 1}.

## Type and Precision

``Distribution/Bernoulli`` is generic over `BinaryFloatingPoint`:

- ``Distribution/Bernoulli``<`Double`>
- ``Distribution/Bernoulli``<`Float`>
- ``Distribution/Bernoulli``<`Float80`> (x86_64 only)

Values are bridged to Boost’s `bernoulli_distribution` through the dynamic factory. Although the distribution is discrete, the interface mirrors the continuous distributions for consistency (PDF behaves as the PMF).

## API Overview

- Initialization: ``Distribution/Bernoulli/init(p:)``
- PMF/log-PMF: ``Distribution/Bernoulli/pdf(_:)``, ``Distribution/Bernoulli/logPdf(_:)``
- CDF/SF: ``Distribution/Bernoulli/cdf(_:)``, ``Distribution/Bernoulli/sf(_:)``
- Inverses: ``Distribution/Bernoulli/quantile(_:)``, ``Distribution/Bernoulli/quantileComplement(_:)``
- Summary statistics: ``Distribution/Bernoulli/mean``, ``Distribution/Bernoulli/variance``, ``Distribution/Bernoulli/entropy``
- Hazard metrics: ``Distribution/Bernoulli/hazard(_:)``, ``Distribution/Bernoulli/chf(_:)``

## Usage

```swift
import SwiftyBoost

let bern = try Distribution.Bernoulli<Double>(p: 0.3)
let pmf1 = try bern.pdf(1)              // 0.3
let cdf0 = try bern.cdf(0)              // 0.7
let q80  = try bern.quantile(0.8)       // 1

let mean = bern.mean                    // p
let variance = bern.variance            // p(1 − p)
let entropy = bern.entropy              // −p log(p) − (1 − p) log(1 − p)
```

## Error Handling

- Initialization throws ``DistributionError/invalidCombination(message:)`` when `p` is outside `[0, 1]`.
- PMF/CDF helpers accept `x` values in `{0, 1}`; other inputs throw the underlying error from ``Distribution/Dynamic``.
- Quantile helpers validate probability arguments in `[0, 1]`.

## Implementation Details

- Delegates to ``Distribution/Dynamic`` with the `bernoulli` vtable and aliases `bernoulli`, `bernoulli_distribution`.
- Entropy is not provided by Boost; the Swift wrapper evaluates the closed form in Swift using ``SwiftyBoost/SpecialFunctions/log1p(_:)->T`` for numerical stability.
- The dynamic factory accepts parameter aliases `p`, `prob`, `probability`, `success`, and `theta`.
