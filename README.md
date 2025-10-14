# SwiftyBoost

[![Documentation](https://img.shields.io/badge/docs-website-blue?logo=swift)](https://strike65.github.io/SwiftyBoost/) [![Pages](https://github.com/strike65/SwiftyBoost/actions/workflows/docs.yml/badge.svg)](https://github.com/strike65/SwiftyBoost/actions/workflows/docs.yml)

SwiftyBoost gives Swift developers direct access to Boost.Math special functions through a thin, type-safe wrapper. The package keeps Swift ergonomics—generics, errors, and platform availability—while reusing Boost’s battle-tested numerical algorithms via a small C++ bridge target.

## Features
- Gamma, Beta, error, Bessel, Legendre, elliptic (Legendre and Carlson), Lambert W, Owen's T, and other high-precision helpers.
- `CBoostBridge` target that forwards Swift calls into the vendored Boost headers under `extern/boost`.
- Architectural awareness with dedicated `Float`, `Double`, and (x86_64) `Float80` overloads plus generic `BinaryFloatingPoint` entry points.

## Requirements
- Swift 6.2 or later.
- macOS 13 or iOS 16 (per `Package.swift`).
- `Float80` support requires an x86_64 runtime.

## Installation (Swift Package Manager)

Add SwiftyBoost to your project manifest:

```swift
// swift-tools-version: 6.2
import PackageDescription

let package = Package(
    name: "YourApp",
    platforms: [ .iOS(.v16), .macOS(.v13) ],
    dependencies: [
        .package(url: "https://github.com/strike65/SwiftyBoost.git", from: "0.1.0")
    ],
    targets: [
        .target(
            name: "YourApp",
            dependencies: [
                .product(name: "SwiftyBoost", package: "SwiftyBoost")
            ]
        )
    ]
)
```

## Usage

```swift
import SwiftyBoost

// Gamma family
let gamma = try SpecialFunctions.gamma(5.0)              // 24
let logG  = try SpecialFunctions.logGamma(10.0)          // ≈ 12.8018

// Error functions
let erf   = try SpecialFunctions.errorFunction(1.2)
let erfc  = try SpecialFunctions.complementaryErrorFunction(1.2)

// Bessel J_v(x)
let j0    = try SpecialFunctions.besselJ(v: 0.0, x: 2.5)

// Factorials
let nFact: Double = try SpecialFunctions.factorial(10)   // 3,628,800
```

APIs throw `SpecialFunctionError` for invalid inputs or domain violations. See DocC symbol docs for function‑specific bounds mirrored from Boost.Math.

## Development
- Core Swift sources reside in `Sources/SwiftyBoost`, grouped by namespace under `Math/SpecialFunctions`.
- The C++ interop layer lives in `Sources/CBoostBridge` and depends on headers in `extern/boost`. Run `git submodule update --init --recursive` after pulling Boost changes.
- Build with `swift build`. Execute `swift test` (or filter suites with `swift test --filter FactorialAndCombinatoricsTests`) to validate changes in `Tests/SwiftyBoostTests`.
- Use `swift package generate-xcodeproj` if you need an Xcode project for debugging.

## Documentation
- DocC sources live under `Sources/SwiftyBoost/Documentation.docc`.
- Generate a static docs site with: `make documentation`.
- Output is written to `Docs/` with entry `Docs/index.html`.
- To build a `.doccarchive` instead: `make documentation-archive` (outputs `Docs/SwiftyBoost.doccarchive`).
 - Browse online: https://strike65.github.io/SwiftyBoost/

## Contributing

### Project Structure & Module Organization
Core APIs live in `Sources/SwiftyBoost`, with entry points in `SwiftyBoost.swift` and focused implementations under `Math/SpecialFunctions`. The `Sources/CBoostBridge` target hosts the C++ interop layer that forwards into `extern/boost`; run `git submodule update --init --recursive` whenever Boost is refreshed. Swift Testing suites that mirror the public API reside in `Tests/SwiftyBoostTests`, and the optional `SwiftyBoostDemo` workspace exercises the package inside an app target.

### Build, Test, and Development Commands
Use `swift build` for a macOS/iOS-compatible build. Execute `swift test` to run the entire suite, or narrow scope with `swift test --filter FactorialAndCombinatoricsTests`. `swift package generate-xcodeproj` is handy when you need an Xcode project for debugging. When modifying the bridge, rebuild with `swift build -Xcxx -std=c++20` to ensure compiler flags match the manifest.

### Coding Style & Naming Conventions
Follow the Swift API Design Guidelines: types in PascalCase, methods and variables in camelCase, and namespaces grouped with nested structs (see `Sources/SwiftyBoost/Math/Namespaces.swift`). Prefer 4-space indentation and trailing commas in multiline literals. Keep bridge helpers prefixed with `bs_` to align with exported headers and avoid symbol clashes. Run `swift format .` before submitting, and align Boost-facing C++ changes with LLVM-style formatting.

### Testing Guidelines
Add coverage with `@Suite` blocks and targeted `@Test` macros from the `Testing` framework. Mirror Boost edge cases by comparing against bridge helpers such as `bs_rising_factorial`. Use architecture guards like `#if arch(x86_64)` when Float80 paths apply. Pair every bug fix with a regression test, and keep high-precision checks tolerant to floating-point error (e.g. `#expect(abs(a - b) < 1e-12)`).

### Commit & Pull Request Guidelines
Follow the short, prefix-based commit style already in history (`add: …`, `fix: …`, `upd: …`, etc.), keep subject lines under 72 characters, and expand in the body when context matters. Pull requests should explain the motivation, list observable changes, link related issues, and attach benchmarks or screenshots for API updates. Confirm `swift test` passes and call out any platform-specific verification before requesting review.

## Versioning & License
- Versioning follows Semantic Versioning. See `CHANGELOG.md` for released changes.
- License: see `LICENSE.md`. Third-party notices: `THIRD-PARTY-NOTICES.txt`, `LICENSE-BOOST.txt`.
