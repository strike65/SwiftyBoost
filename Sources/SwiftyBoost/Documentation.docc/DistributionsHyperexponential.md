# ``SwiftyBoost/Distribution/Hyperexponential``

Hyperexponential distribution: a finite mixture of exponential phases with optional custom mixing probabilities.

## Definition

A hyperexponential distribution represents the waiting time until completion of one of several exponential phases, where the specific phase is chosen at random. In SwiftyBoost you provide:

- `rates`: strictly positive rate parameters `λᵢ` for each exponential phase.
- `probabilities` (optional): non-negative mixing weights. When omitted, phases are mixed uniformly. If supplied, the weights are normalised to sum to one.

Boost.Math performs all heavy lifting (PDF/CDF/log-PDF, survival, quantiles, hazards, and summary statistics), while SwiftyBoost validates parameters and bridges into the Boost implementation.

## Type and Precision

``Distribution/Hyperexponential`` is generic over `BinaryFloatingPoint`, with specialisations for:

- ``Distribution/Hyperexponential``<`Double`>
- ``Distribution/Hyperexponential``<`Float`>
- ``Distribution/Hyperexponential``<`Float80`> (x86_64 only; falls back to the double-precision backend elsewhere)

Each instance reuses the runtime factory ``Distribution/Dynamic`` with the `hyperexponential` entry in the Boost-backed vtable.

## API Overview

- Initialisation: ``Distribution/Hyperexponential/init(probabilities:rates:)``
- Density & log-density: ``Distribution/Hyperexponential/pdf(_:)``, ``Distribution/Hyperexponential/logPdf(_:)``
- Distribution functions: ``Distribution/Hyperexponential/cdf(_:)``, ``Distribution/Hyperexponential/sf(_:)``
- Inverse functions: ``Distribution/Hyperexponential/quantile(_:)``, ``Distribution/Hyperexponential/quantileComplement(_:)``
- Summary statistics: ``Distribution/Hyperexponential/mean``, ``Distribution/Hyperexponential/variance``, ``Distribution/Hyperexponential/skewness``, ``Distribution/Hyperexponential/kurtosis``, ``Distribution/Hyperexponential/kurtosisExcess``, ``Distribution/Hyperexponential/mode``, ``Distribution/Hyperexponential/median``
- Hazard metrics: ``Distribution/Hyperexponential/hazard(_:)``, ``Distribution/Hyperexponential/chf(_:)``
- Phase metadata: ``Distribution/Hyperexponential/probabilities`` (normalised mixture weights), ``Distribution/Hyperexponential/rates``, ``Distribution/Hyperexponential/phaseCount``

## Usage

```swift
import SwiftyBoost

// Mix three exponential phases with custom weights
let hyperexp = try Distribution.Hyperexponential<Double>(
    probabilities: [0.3, 0.2, 0.5],
    rates: [0.6, 1.5, 3.0]
)

let x: Double = 1.25
let pdf = try hyperexp.pdf(x)
let cdf = try hyperexp.cdf(x)
let sf  = try hyperexp.sf(x)

// Quantiles and hazards
let q95 = try hyperexp.quantile(0.95)
let hazard = try hyperexp.hazard(0.5)

// Summary statistics
let mean = hyperexp.mean
let variance = hyperexp.variance
let skewness = hyperexp.skewness
```

## Mathematical Notes

- Support: `x ∈ [0, +∞)`
- PDF: `f(x) = Σᵢ pᵢ λᵢ e^(−λᵢ x)` where `pᵢ` are the normalised mixing probabilities and `λᵢ > 0`.
- CDF: `F(x) = 1 − Σᵢ pᵢ e^(−λᵢ x)`.
- Survival: `S(x) = Σᵢ pᵢ e^(−λᵢ x)`.
- Mean: `E[X] = Σᵢ pᵢ / λᵢ`.
- Variance: `Var[X] = 2 Σᵢ pᵢ / λᵢ² − (E[X])²`.
- Skewness & kurtosis follow closed forms derived from the sums `Σᵢ pᵢ / λᵢ^k` for `k = 1 … 4`; Boost computes them analytically.
- Mode: 0.
- Median: Obtained numerically via the Boost quantile routine.

## Implementation Details

- The Swift initializer validates that all rates are finite and strictly positive, while mixing probabilities (if provided) are finite, non-negative, and contain at least one positive entry. Probabilities are normalised before being forwarded to Boost.
- Parameter dictionaries sent to ``Distribution/Dynamic`` use indexed keys (`rate0`, `rate1`, …) and optional weights (`prob0`, `prob1`, …). When no weights are supplied, the backend uses a uniform mixture.
- Hazard and cumulative hazard are computed via the generic vtable using Boost’s numerically stable survival and PDF evaluations.
