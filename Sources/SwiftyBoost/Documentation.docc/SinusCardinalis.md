# Sinus Cardinalis (π-Normalized)

Boost-backed π-normalized cardinal sine and hyperbolic cardinal sine with real and complex variants.

## Definitions

For real `x`:

- ``SwiftyBoost/SpecialFunctions/sinc_pi(_:)->T`` computes `sinc(πx) = sin(πx)/(πx)`, with the removable singularity `sinc(0) = 1`.
- ``SwiftyBoost/SpecialFunctions/sinhc_pi(_:)->T`` computes `sinhc(πx) = sinh(πx)/(πx)`, with the removable singularity `sinhc(0) = 1`.

For complex `z` (double/float/(x86_64) long double):

- ``SwiftyBoost/SpecialFunctions/sincc_pi(_:)`` computes `sinc(πz) = sin(πz)/(πz)`, with the removable singularity `sincc(0) = 1 + 0i`.
- ``SwiftyBoost/SpecialFunctions/sinhcc_pi(_:)`` computes `sinhc(πz) = sinh(πz)/(πz)`, with the removable singularity `sinhcc(0) = 1 + 0i`.

## Type and Precision

- [Real](https://github.com/apple/swift-numerics) overloads are available generically for `T: BinaryFloatingPoint` and specialized for `Double`, `Float`, and (x86_64) `Float80`.
- Complex overloads exist for `ComplexD`, `ComplexF`, and (x86_64) `ComplexL`.

## Usage

This snippet first exercises the [Real](https://github.com/apple/swift-numerics) overloads before demonstrating the complex identity.

```swift
import SwiftyBoost

let s0  = SpecialFunctions.sinc_pi(0.0)     // 1
let s12 = SpecialFunctions.sinc_pi(0.5)     // ≈ 2/π
let sh0 = SpecialFunctions.sinhc_pi(0.0)    // 1

// Complex identity: sincc_pi(i y) = sinhcc_pi(y)
let z = ComplexD(0, 0.7)
let sc = SpecialFunctions.sincc_pi(z)
let sh = SpecialFunctions.sinhc_pi(0.7)
```

## Notes

- Values propagate IEEE-754 NaN and infinities per the backend.
- Implementations delegate to Boost via `CBoostBridge` with stable `bs_` functions: `bs_sinc_pi`, `bs_sinhc_pi`, `bs_sincc_pi`, `bs_sinhcc_pi` (and `_f`, `_l` float/long variants).
