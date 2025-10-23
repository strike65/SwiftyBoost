# Constants Catalog

High-precision mathematical constants bridged from Boost.Math into Swift.

## Overview

``Constants`` exposes a string-backed catalogue of well-known mathematical
constants—π, τ, √π, √2, `e`, Euler–Mascheroni (γ), Catalan’s constant, Apéry’s ζ(3),
the golden ratio φ, and a family of logarithmic/reciprocal variants. Values are stored
as high-precision decimal literals sourced from Boost.Math’s tables and converted on
demand into any `BinaryFloatingPoint` type.

Use ``Constants/value(_:as:)`` when you need to select constants dynamically at
runtime, or the convenience helpers such as ``Constants/pi(_:)`` and
``Constants/e(_:)`` when the identifier is known at compile time.

## Selecting a Precision

- Call helpers with an explicit static type: `Constants.pi(Double.self)`,
  `Constants.euler(Float.self)`, or—on x86_64—`Constants.zetaThree(Float80.self)`.
- Let Swift infer the type from context: `let tau: Double = Constants.twoPi()`.
- All conversions are thread-safe, side-effect free, and work for any type that conforms to
  both [Real](https://github.com/apple/swift-numerics) and `BinaryFloatingPoint`.

## Pi-Derived Set

- `Constants.pi()`, `twoPi()`, `halfPi()`, `quarterPi()`, `thirdPi()`,
  `twoThirdsPi()`, `threeQuartersPi()`, `sixthPi()`
- `Constants.piSqr()`, `rootPi()`, `oneDivPi()`, `oneDivTwoPi()`,
  `oneDivRootPi()`, `twoDivPi()`, `twoDivRootPi()`

## Logarithms and Exponentials

- `Constants.e()` (Euler’s number)
- `Constants.oneDivE()` (reciprocal of `e`)
- `Constants.lnTwo()`, `lnTen()`, `lnLnTwo()`

## Special Constants

- `Constants.rootTwo()`, `rootThree()`
- `Constants.euler()` (Euler–Mascheroni γ)
- `Constants.catalan()`
- `Constants.zetaThree()` (Apéry’s constant ζ(3))
- `Constants.phi()` (golden ratio)

## Usage

```swift
import SwiftyBoost

let tau: Double = Constants.twoPi()
let sqrtPi: Float = Constants.rootPi()
#if arch(x86_64)
let gamma80: Float80 = Constants.euler(Float80.self)
#endif
```

You can also index into the catalogue dynamically:

```swift
let identifier: Constants.Identifier = .lnTwo
let ln2 = Constants.value(identifier, as: Double.self)
```

## Implementation Details

- The decimal catalogue mirrors Boost.Math’s high-precision constants, trimmed for storage and
  embedded directly in Swift.
- Conversion prefers initializing through `LosslessStringConvertible` when available and falls
  back to a manual parser with binary exponent accumulation to support all `BinaryFloatingPoint` types.
- Utilities such as ``Constants/convert(decimal:as:)`` and ``Constants/powTen(of:)``
  provide deterministic rounding behavior across precisions.
- An exhaustive ``Constants/Identifier`` enum enumerates every supported constant, simplifying
  dynamic lookups and exhaustive switching in client code.
