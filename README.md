# SwiftyBoost

[![Documentation](https://img.shields.io/badge/docs-website-blue?logo=swift)](https://strike65.github.io/SwiftyBoost/documentation/swiftyboost) [![Pages](https://github.com/strike65/SwiftyBoost/actions/workflows/docs.yml/badge.svg)](https://github.com/strike65/SwiftyBoost/actions/workflows/docs.yml)

SwiftyBoost gives Swift developers direct access to Boost.Math special functions through a thin, type-safe wrapper. The package keeps Swift ergonomics—generics, errors, and platform availability—while reusing Boost’s battle-tested numerical algorithms via a small C++ bridge target.

## Features
- Gamma, Beta, error, Bessel, Legendre, elliptic (Legendre and Carlson), Lambert W, Owen's T, and other high-precision helpers.
- Probability distributions with Boost-backed implementations:
  - Gamma, Beta, Chi-squared, Student’s t, Fisher’s F, Bernoulli, and Arcsine (PDF/PMF, CDF/SF, quantiles, hazards, moments).
  - Typed wrappers delegate internally to a unified runtime vtable (`Distribution.Dynamic`) for consistent behavior across precisions.
  - Unified runtime factory to construct distributions by name at runtime (see "Dynamic Distribution Factory").
- High-precision mathematical constants (`π`, `e`, `√2`, Euler–Mascheroni, Catalan, ζ(3), φ, …) exposed via the generic `Constants<T>` namespace with dedicated `Float`, `Double`, and (x86_64) `Float80` entry points.
- `CBoostBridge` target that forwards Swift calls into the vendored Boost headers under `extern/boost`.
- Architectural awareness with dedicated `Float`, `Double`, and (x86_64) `Float80` overloads plus generic `BinaryFloatingPoint` entry points.
- Re-exported Swift Numerics `Complex<T>` type with SwiftyBoost convenience APIs: arithmetic, polar helpers, and Boost-backed elementary functions (`exp`, `log`, `sin`, `cos`, `tan`, `sinh`, `cosh`, `tanh`, `atan`). Double/Float/(x86_64) Float80 specializations call bridge functions (`bs_*`).

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

// Common helpers
let invS = try SpecialFunctions.rsqrt(9.0)               // 1/3
let s0   = try SpecialFunctions.sinc_pi(0.0)             // 1 (limit)
let s12  = try SpecialFunctions.sinc_pi(0.5)             // 2/π
let sh12 = try SpecialFunctions.sinhc_pi(0.5)            // sinh(π/2)/(π/2)
let z    = ComplexD(0, 0.7)
let sc   = SpecialFunctions.sincc_pi(z)                  // equals sinhc_pi(0.7)
```

```swift
// Constants (choose your precision)
let tau   = Constants<Double>.twoPi
let phi   = Constants<Float>.phi
#if arch(x86_64)
let rootPi80 = Constants<Float80>.rootPi
#endif
```

```swift
// Distributions
// Student’s t (ν = 5)
let t = try Distribution.StudentT<Double>(degreesOfFreedom: 5)
let t_q975 = try t.quantile(0.975)

// Fisher’s F (df1 = 10, df2 = 20)
let f = try Distribution.FisherF<Double>(degreesOfFreedom1: 10, degreesOfFreedom2: 20)
let f_cdf = try f.cdf(2.5)

// Beta(α = 2, β = 5)
let beta = try Distribution.Beta<Double>(alpha: 2, beta: 5)
let beta_mean = beta.mean

// Chi-squared (ν = 8)
let chi2 = try Distribution.ChiSquared<Double>(degreesOfFreedom: 8)
let chi2_sf = try chi2.sf(12.0)

// Bernoulli(p = 0.3)
let bern = try Distribution.Bernoulli<Double>(p: 0.3)
let bern_entropy = bern.entropy

// Arcsine on [0, 1]
let arcsine = try Distribution.Arcsine<Double>(minX: 0, maxX: 1)
let arcsine_pdf = try arcsine.pdf(0.2)
```

APIs throw `SpecialFunctionError` for invalid inputs or domain violations. See DocC symbol docs for function‑specific bounds mirrored from Boost.Math.

### Dynamic Distribution Factory

Create distributions by name with a parameter dictionary, without committing to a concrete Swift wrapper type upfront. This is powered by a small C vtable and a Swift wrapper that conforms to `DistributionProtocol`.

```swift
// Gamma(k, theta)
let dynG = try Distribution.Dynamic<Double>(
  distributionName: "gamma",
  parameters: ["shape": 4.5, "scale": 1.0]
)
let p = try dynG.pdf(2.0)
let x = try dynG.quantile(0.95)

// Student's t (aliases: student_t, t)
let dynT = try Distribution.Dynamic<Float>(
  distributionName: "student_t",
  parameters: ["df": 12]
)

// Fisher F (aliases: fisherf, f)
let dynF = try Distribution.Dynamic<Double>(
  distributionName: "fisherf",
  parameters: ["df1": 10, "df2": 20]
)

// Beta (aliases: beta_distribution)
let dynB = try Distribution.Dynamic<Double>(
  distributionName: "beta",
  parameters: ["alpha": 2.0, "beta": 5.0]
)

// Chi-squared (aliases: chi_squared, chisq)
let dynChi = try Distribution.Dynamic<Double>(
  distributionName: "chisquared",
  parameters: ["df": 8.0]
)

// Bernoulli (aliases: bernoulli_distribution)
let dynBern = try Distribution.Dynamic<Double>(
  distributionName: "bernoulli",
  parameters: ["p": 0.3]
)

// Arcsine (aliases: arcsine_distribution)
let dynA = try Distribution.Dynamic<Float>(
  distributionName: "arcsine",
  parameters: ["minX": 0, "maxX": 1]
)
```

Notes:
- The Swift wrapper falls back to Swift-side formulas for certain metrics not provided by Boost (e.g., Fisher F and Arcsine entropy).
- Typed distribution wrappers delegate to `Distribution.Dynamic` internally, so dynamic and typed paths share the same Boost backend and numerical policies.
- Supported names and parameter aliases are documented in `DIST-Factory-README.md`.
- Under the hood, the factory uses a C vtable with a non-null context pointer and nullable function pointers per metric.

See `DIST-Factory-README.md` for the full design, ABI surface, and extension guide.

### Complex Numbers

SwiftyBoost re-exports Swift Numerics’ `Complex<T>` and layers on stable arithmetic helpers plus Boost-backed elementary functions where available.

```swift
import SwiftyBoost

// Construction and constants
let z1 = Complex<Double>(3, -4)
let z2: ComplexD = .i                      // 0 + 1i
let z3 = Complex<Float>.fromPolar(radius: 2, phase: .pi/6)

// Arithmetic
let s = z1 + ComplexD(1, 2)
let p = z1 * z2
let q = z1 / ComplexD(-1, 4)              // Smith’s algorithm for stability

// Properties and transforms
let m  = z1.magnitude                      // 5
let arg = z1.phase                         // atan2(imag, real)
let zr = z1.squareRoot                     // principal sqrt

// Elementary complex functions (Boost-backed for Double/Float/Float80)
let e  = z1.exp
let l  = z1.log
let sn = z1.sin
let cs = z1.cos
```

Typealiases are provided for convenience: `ComplexD` (`Double`), `ComplexF` (`Float`), and `ComplexL` (`Float80`, x86_64 only).

## Development
- Core Swift sources reside in `Sources/SwiftyBoost`, grouped by namespace under `Math/SpecialFunctions`.
- The C++ interop layer lives in `Sources/CBoostBridge` and depends on headers in `extern/boost`. Run `git submodule update --init --recursive` after pulling Boost changes.
  - Public headers are organized by domain under `Sources/CBoostBridge/include/` (e.g. `bs_constants.h`, `bs_gamma.h`, `bs_bessel.h`, …). The umbrella header `CBoostBridge.h` re-exports these and is the module’s public entry point.
  - Implementations are split into focused units under `Sources/CBoostBridge/impl/` and composed by `Sources/CBoostBridge/CBoostBridge.cpp`. Shared helpers live in `Sources/CBoostBridge/internal/bs_internal.hpp`.
  - All exported bridge symbols remain prefixed with `bs_` and simply forward to Boost.Math; no new algorithms are implemented in the bridge.
- Build with `swift build`. Execute `swift test` (or filter suites with `swift test --filter FactorialAndCombinatoricsTests`) to validate changes in `Tests/SwiftyBoostTests`.
- Use `swift package generate-xcodeproj` if you need an Xcode project for debugging.

## Documentation
- DocC sources live under `Sources/SwiftyBoost/Documentation.docc`.
- Generate a static docs site with: `make documentation`.
 - Output is written to `Docs/` and published via GitHub Pages.
- The site root redirects to the main DocC page: https://strike65.github.io/SwiftyBoost/
- Direct landing page (canonical DocC path): https://strike65.github.io/SwiftyBoost/documentation/swiftyboost/
- To build a `.doccarchive` instead: `make documentation-archive` (outputs `Docs/SwiftyBoost.doccarchive`).

### Result-Type Promotions (Quick Reference)

- Float ↔ Double → Double
  - Multi-argument functions promoted from mixed `Float`/`Double` evaluate in `Double` and return `Double`.
- Any presence of Float80 → Float80 (x86_64)
  - On x86_64, if any argument is `Float80`, evaluation uses extended precision and returns `Float80`.
- Single-argument functions
  - Return the same type as the argument; no promotion applies.

Supported mixed promotions (high level):
- Gamma: ratios, incomplete/regularized P/Q (+ inverses, derivative)
- Beta: complete/incomplete/regularized (+ inverses, derivative, parameter solvers)
- Bessel: cylindrical J/Y/I/K and their derivatives (integer-order convenience overloads included)
- Elliptic: Legendre F/E/Π (incomplete + complete Π), Carlson RC/RF/RD/RJ/RG
- Chebyshev series: Clenshaw recursion mixes `[coefficients]` and `x` types
- Owen’s T: all pairs
- Spherical harmonics: mixed angle types (ComplexD/ComplexL return types)
- Hypergeometric: 1F0, 0F1 (full), 1F1/2F0 (selected single‑parameter mixes), pFq (arrays and z)

Advice:
- Prefer `Double` unless you know `Float` is sufficient or `Float80` is required (x86_64 only).
- Use explicit `Float`/`Double`/`Float80` overloads in tight loops to avoid conversions.
- Mixed promotions are for convenience when inputs naturally come in different precisions.
- Promotions don’t change numeric domains; invalid inputs still throw `SpecialFunctionError`.

See DocC: “Result-Type Promotions” for the full list and details.

## Contributing

### Project Structure & Module Organization
Core APIs live in `Sources/SwiftyBoost`, with entry points in `SwiftyBoost.swift` and focused implementations under `Math/SpecialFunctions`. The `Sources/CBoostBridge` target hosts the C++ interop layer that forwards into `extern/boost`; run `git submodule update --init --recursive` whenever Boost is refreshed. Swift Testing suites that mirror the public API reside in `Tests/SwiftyBoostTests`, and the optional `SwiftyBoostDemo` workspace exercises the package inside an app target.

### Build, Test, and Development Commands
Use `swift build` for a macOS/iOS-compatible build. Execute `swift test` to run the entire suite, or narrow scope with `swift test --filter FactorialAndCombinatoricsTests`. `swift package generate-xcodeproj` is handy when you need an Xcode project for debugging. The manifest sets the C++ language standard to match the project; no extra flags are required for normal builds.

### Coding Style & Naming Conventions
Follow the Swift API Design Guidelines: types in PascalCase, methods and variables in camelCase, and namespaces grouped with nested structs (see `Sources/SwiftyBoost/Math/Namespaces.swift`). Prefer 4-space indentation and trailing commas in multiline literals. Keep bridge helpers prefixed with `bs_` to align with exported headers and avoid symbol clashes. Run `swift format .` before submitting, and align Boost-facing C++ changes with LLVM-style formatting.

### Testing Guidelines
Add coverage with `@Suite` blocks and targeted `@Test` macros from the `Testing` framework. Mirror Boost edge cases by comparing against bridge helpers such as `bs_rising_factorial`. Use architecture guards like `#if arch(x86_64)` when Float80 paths apply. Pair every bug fix with a regression test, and keep high-precision checks tolerant to floating-point error (e.g. `#expect(abs(a - b) < 1e-12)`).

### Commit & Pull Request Guidelines
Follow the short, prefix-based commit style already in history (`add: …`, `fix: …`, `upd: …`, etc.), keep subject lines under 72 characters, and expand in the body when context matters. Pull requests should explain the motivation, list observable changes, link related issues, and attach benchmarks or screenshots for API updates. Confirm `swift test` passes and call out any platform-specific verification before requesting review.

## Versioning & License
- Versioning follows Semantic Versioning. See `CHANGELOG.md` for released changes.
- License: see `LICENSE.md`. Third-party notices: `THIRD-PARTY-NOTICES.txt`, `LICENSE-BOOST.txt`.
