# ``SwiftyBoost/Distribution/Logistic``

Logistic distribution in location–scale form.

## Definition

The logistic distribution is a symmetric distribution with heavier tails than the normal.
It is parameterised by a location parameter `μ` and a scale parameter `s > 0`. Its PDF is

```
f(x; μ, s) = \frac{e^{-(x-μ)/s}}{s \bigl(1 + e^{-(x-μ)/s}\bigr)^2}.
```

## Type and Precision

``Distribution/Logistic`` is generic over `BinaryFloatingPoint` and available as:

- ``Distribution/Logistic``<`Double`>
- ``Distribution/Logistic``<`Float`>
- ``Distribution/Logistic``<`Float80`> (x86_64 only; falls back to the double backend elsewhere)

## API Overview

- Initialiser: ``Distribution/Logistic/init(location:scale:)``
- Density & log-density: ``Distribution/Logistic/pdf(_:)``, ``Distribution/Logistic/logPdf(_:)``
- Distribution functions: ``Distribution/Logistic/cdf(_:)``, ``Distribution/Logistic/sf(_:)``
- Inverse functions: ``Distribution/Logistic/quantile(_:)``, ``Distribution/Logistic/quantileComplement(_:)``
- Summary statistics: ``Distribution/Logistic/mean``, ``Distribution/Logistic/variance``,
  ``Distribution/Logistic/mode``, ``Distribution/Logistic/median``, ``Distribution/Logistic/skewness``,
  ``Distribution/Logistic/kurtosis``, ``Distribution/Logistic/kurtosisExcess``, ``Distribution/Logistic/entropy``
- Hazard metrics: ``Distribution/Logistic/hazard(_:)``, ``Distribution/Logistic/chf(_:)``

## Usage

```swift
import SwiftyBoost

let logistic = try Distribution.Logistic<Double>(location: 0.0, scale: 1.2)

let x: Double = 0.75
let pdf = try logistic.pdf(x)
let cdf = try logistic.cdf(x)
let sf  = try logistic.sf(x)

let q95 = try logistic.quantile(0.95)
let entropy = logistic.entropy
```

## Mathematical Notes

- Support: `x ∈ (-∞, ∞)`
- Mean, median, and mode coincide at `μ`.
- Variance: `(π s)^2 / 3`
- Skewness: `0`, Kurtosis: `21/5`, Excess Kurtosis: `6/5`
- Differential entropy: `H = 2 + ln s`

## Implementation Details

- Parameters are validated in Swift (`scale > 0`) before delegating to ``Distribution/Dynamic``.
- All metrics are provided by Boost.Math; the Swift wrapper simply forwards the values.
- KL divergence reuses the shared dynamic integration routine.
