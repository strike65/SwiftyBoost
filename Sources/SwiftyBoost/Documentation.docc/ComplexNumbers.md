# Complex Numbers

SwiftyBoost provides a lightweight, generic complex type with stable arithmetic and Boost-backed elementary functions for common precisions.

## Overview

``Complex`` is generic over `BinaryFloatingPoint` and conforms to `Sendable`, `Equatable`, `Hashable`, and `Codable`. The type stores `real` and `imag` parts and offers:

- Stable arithmetic: addition, subtraction, multiplication, and Smith’s algorithm for division.
- Robust magnitude and principal square root implementations.
- Polar helpers: `fromPolar(radius:phase:)` and the principal `phase` (argument).
- Elementary complex functions for common precisions via the Boost bridge (`bs_*`): `exp`, `log`, `sin`, `cos`, `tan`, `sinh`, `cosh`, `tanh`, and `atan`.

Typealiases for convenience:
- `ComplexD` = `Complex<Double>`
- `ComplexF` = `Complex<Float>`
- `ComplexL` = `Complex<Float80>` (x86_64 only)

## Construction

```swift
let z1 = ComplexD(3, -4)              // 3 - 4i
let z2: ComplexF = .i                  // 0 + 1i
let z3 = Complex<Double>.one           // 1 + 0i
let p  = ComplexD.fromPolar(radius: 2, phase: .pi/6)
```

## Arithmetic

```swift
let a = ComplexD(1, 2)
let b = ComplexD(-3, 0.5)

let sum = a + b
let dif = a - b
let prod = a * b
let quo = a / b                        // Smith’s stable algorithm

let scaled = 2 * a - 0.5               // scalar mixes
let inv = a.reciprocal()
```

## Magnitude, Phase, and Square Root

```swift
let m = a.magnitude                    // robust: max-min scaling
let arg = a.phase                      // principal value in (−π, π]
let r = a.normalized()                 // unit magnitude or unchanged if zero
let s = a.squareRoot                   // principal sqrt
```

## Elementary Complex Functions (Boost-backed)

Available for `ComplexD`, `ComplexF`, and (x86_64) `ComplexL`:

```swift
let e  = a.exp
let l  = a.log
let sn = a.sin
let cs = a.cos
let tn = a.tan
let sh = a.sinh
let ch = a.cosh
let th = a.tanh
let at = a.atan
```

These properties delegate to the C++ bridge (`CBoostBridge`) and use Boost.Math’s implementations. All bridged symbols are prefixed with `bs_` to ensure ABI stability.

## Identities and Helpers

The π-normalized cardinal functions integrate with complex inputs:

```swift
// sincc_pi(0 + iy) = sinhc_pi(y)
let z = ComplexD(0, 0.7)
let sc = SpecialFunctions.sincc_pi(z)   // equals SpecialFunctions.sinhc_pi(0.7)
```

## Availability

- `ComplexL` (Float80) is available on `x86_64` only.
- Polar and elementary functions are specialized for `Double`, `Float`, and (where applicable) `Float80`.

