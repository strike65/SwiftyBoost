# Constants Catalog

High-precision mathematical constants bridged from Boost.Math into Swift.

## Overview

``Constants`` is a generic namespace (`Constants<T>`) that exposes static computed
properties for well-known constants—including π, τ, √π, √2, `e`, Euler–Mascheroni (γ),
Catalan’s constant, Apery’s ζ(3), the golden ratio φ, and a family of logarithms and reciprocal forms.
Each accessor calls through to a dedicated Boost-backed entry point to preserve precision for the selected floating-point type.

## Selecting a Precision

- Use `Constants<Double>` for the default path; it maps to `double` values in Boost.
- Use `Constants<Float>` to shave cycles in GPU-heavy or mobile scenarios; values route to `float`.
- On x86_64 (and i386 simulators) you can access extended precision via `Constants<Float80>`.

All specializations are thread-safe and return pure values with no side effects.

## Pi-Derived Set

- `Constants<T>.pi`, `twoPi`, `halfPi`, `quarterPi`, `thirdPi`, `twoThirdsPi`, `threeQuartersPi`, `sixthPi`
- `Constants<T>.piSqr`, `rootPi`, `oneDivPi`, `oneDivTwoPi`, `oneDivRootPi`, `twoDivPi`, `twoDivRootPi`

## Logarithms and Exponentials

- `Constants<T>.e` (Euler’s number)
- `Constants<T>.lnTwo`, `lnTen`, `lnLnTwo`

## Special Constants

- `Constants<T>.rootTwo`, `rootThree`
- `Constants<T>.euler` (Euler–Mascheroni γ)
- `Constants<T>.catalan`
- `Constants<T>.zetaThree` (Apéry’s constant ζ(3))
- `Constants<T>.phi` (golden ratio)

## Usage

```swift
import SwiftyBoost

let tau = Constants<Double>.twoPi
let sqrtPi = Constants<Float>.rootPi
#if arch(x86_64)
let gamma80 = Constants<Float80>.euler
#endif
```

## Implementation Details

- All constants forward to `bs_const_*` helpers defined in `CBoostBridge`.
- Each helper returns the value computed (or stored) natively in Boost for the requested precision, avoiding intermediate conversions.
- Every accessor carries a DocC symbol reference and inline doc comment for quick lookup in Xcode.
