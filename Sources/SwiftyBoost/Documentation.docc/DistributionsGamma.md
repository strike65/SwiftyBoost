# ``SwiftyBoost/Distribution/Gamma``

Gamma Distribution Γ(k, θ) with shape k > 0 and scale θ > 0.

## Definition

Given shape parameter k > 0 and scale parameter θ > 0, the probability density function (PDF) is

f(x; k, θ) = x^(k−1) e^(−x/θ) / (Γ(k) θ^k), for x ≥ 0.

The cumulative distribution function (CDF) equals the regularized lower incomplete gamma function: F(x) = P(k, x/θ). The survival function (upper tail) is S(x) = Q(k, x/θ).

Quantiles are defined by F(x) = p ⇒ x = θ · P⁻¹(k, p).

## Type and Precision

``Distribution/Gamma`` is generic over `BinaryFloatingPoint` and is available in these common specializations:

- ``Distribution/Gamma``<`Double`>
- ``Distribution/Gamma``<`Float`>
- ``Distribution/Gamma``<`Float80`> (x86_64 only; falls back to `Double` on other platforms)

Each instance constructs a Boost.Math `gamma_distribution` once and holds an internal opaque handle. All evaluations reuse that handle for efficiency.

## API Overview

- Initialization: ``Distribution/Gamma/init(shape:scale:)``
- Density: ``Distribution/Gamma/pdf(_:)`` and log-density via ``SwiftyBoost/SpecialFunctions/logGamma(_:)->T`` building blocks
- Distribution functions: ``Distribution/Gamma/cdf(_:)`` and ``Distribution/Gamma/sf(_:)``
- Inverse functions: ``Distribution/Gamma/quantile(_:)`` and ``Distribution/Gamma/quantileComplement(_:)``
- Moments: ``Distribution/Gamma/mean``, ``Distribution/Gamma/variance``, ``Distribution/Gamma/mode``

## Usage

```swift
import SwiftyBoost

// Construct a Gamma distribution once and reuse it.
let g = try Distribution.Gamma<Double>(shape: 2.5, scale: 1.2)

let x: Double = 0.8
let pdf = try g.pdf(x)      // density at x
let cdf = try g.cdf(2.0)    // lower-tail probability
let sf  = try g.sf(2.0)     // upper-tail probability

let q95 = try g.quantile(0.95)
let mean = g.mean           // k θ
let var2 = g.variance       // k θ²
let mode = g.mode           // (k − 1) θ if k ≥ 1; otherwise nil
```

## Mathematical Notes

- Mean: E[X] = k θ
- Variance: Var[X] = k θ²
- Mode: (k − 1) θ if k ≥ 1; undefined for k < 1
- Support: x ∈ [0, ∞)

## Error Handling

- Initialization throws if `shape ≤ 0` or `scale ≤ 0`.
- ``Distribution/Gamma/pdf(_:)``, ``Distribution/Gamma/cdf(_:)``, and ``Distribution/Gamma/sf(_:)`` throw if `x < 0`.
- ``Distribution/Gamma/quantile(_:)`` and ``Distribution/Gamma/quantileComplement(_:)`` throw if their probability arguments are outside `[0, 1]`.

## Implementation Details

All operations delegate to Boost via the C bridge. The Swift type holds an opaque pointer to a `boost::math::gamma_distribution` of the appropriate precision. Calls route through stable `bs_` functions (e.g., `bs_gamma_pdf_h`, `bs_gamma_cdf_h`, `bs_gamma_quantile_h`).

## See Also

- ``SwiftyBoost/SpecialFunctions/regularizedGammaP(_:x:)->T``
- ``SwiftyBoost/SpecialFunctions/regularizedGammaQ(_:x:)->T``
- ``SwiftyBoost/SpecialFunctions/regularizedGammaPInv(_:p:)->T``
- ``SwiftyBoost/SpecialFunctions/regularizedGammaQInv(_:q:)->T``
