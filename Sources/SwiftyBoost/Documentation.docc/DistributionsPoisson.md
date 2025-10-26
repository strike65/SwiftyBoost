#
``SwiftyBoost/Distribution/Poisson``

Poisson distribution with rate (mean) `λ`.

## Definition

The Poisson distribution is discrete on the non-negative integers with mass function

```
P(X = k) = e^{-λ} λ^{k} / k!,  k = 0, 1, 2, …
```

where `λ ≥ 0`.

## Type and Precision

``Distribution/Poisson`` is generic over `BinaryFloatingPoint`. Concrete implementations exist for `Float`, `Double`, and (on supported architectures) `Float80`. All evaluations are delegated to Boost.Math through the dynamic distribution factory.

## API Overview

- Initialisation: ``Distribution/Poisson/init(lambda:)``
- Probability mass/log mass: ``Distribution/Poisson/pdf(_:)``, ``Distribution/Poisson/logPdf(_:)``
- Distribution functions: ``Distribution/Poisson/cdf(_:)``, ``Distribution/Poisson/sf(_:)``
- Inverses: ``Distribution/Poisson/quantile(_:)``, ``Distribution/Poisson/quantileComplement(_:)``
- Summary statistics: ``Distribution/Poisson/mean``, ``Distribution/Poisson/variance``, ``Distribution/Poisson/skewness``, ``Distribution/Poisson/kurtosis``, ``Distribution/Poisson/kurtosisExcess``, ``Distribution/Poisson/mode``, ``Distribution/Poisson/median``
- Hazards: ``Distribution/Poisson/hazard(_:)``, ``Distribution/Poisson/chf(_:)``
- Divergence: ``Distribution/Poisson/klDivergence(relativeTo:options:)``

## Usage

```swift
import SwiftyBoost

let poisson = try Distribution.Poisson<Double>(lambda: 3.5)

let k: Double = 2
let pmf = try poisson.pdf(k)
let cdfValue = try poisson.cdf(k)
let quantile95 = try poisson.quantile(0.95)
```

## Mathematical Notes

- Support: `{0, 1, 2, …}` with lattice step 1 and origin 0
- Mean and variance both equal `λ`
- Skewness = `1 / √λ`, excess kurtosis = `1 / λ`
- Mode is `⌊λ⌋` (when `λ` is an integer, `λ − 1` is also a mode)
- Reduces to a point-mass at 0 for `λ = 0`

## Implementation Details

- Swift validates `λ ≥ 0` and finiteness before delegating to Boost.
- Boost.Math provides PMF/CDF/quantile evaluations and higher-order statistics.
- Lattice metadata is surfaced via ``Distribution/DistributionProtocol/latticeStep`` and ``Distribution/DistributionProtocol/latticeOrigin``.
- Entropy is currently unavailable from Boost for this distribution and surfaces as `nil`.
