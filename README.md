# SwiftyBoost

[![Documentation](https://img.shields.io/badge/docs-website-blue?logo=swift)](https://strike65.github.io/SwiftyBoost/documentation/swiftyboost) [![Pages](https://github.com/strike65/SwiftyBoost/actions/workflows/docs.yml/badge.svg)](https://github.com/strike65/SwiftyBoost/actions/workflows/docs.yml)

SwiftyBoost gives Swift developers direct access to selected functions/methods of Boost.Math through a thin, type-safe wrapper. The package keeps Swift ergonomics—generics, errors, and platform availability—while reusing Boost’s battle-tested numerical algorithms via a small C++ bridge target.

## Disclaimer

Any bugs, crashes, numerical inaccuracies, or documentation mistakes in this Swift wrapper are solely the responsibility of the author of this package. They are not affiliated with, caused by, or endorsed by Boost or Boost.org. The Boost libraries are used here via headers and are respected as upstream dependencies; any wrapper-side behavior or errors should be reported to this project, not to Boost.

## Features
- Gamma, Beta, error, Bessel, Airy, Legendre (including Legendre–Stieltjes quadrature polynomials), Gegenbauer, Jacobi, Jacobi elliptic functions, Jacobi theta functions, Hermite, elliptic (Legendre/Carlson forms, Jacobi Zeta, Heuman’s lambda), Lambert W, Owen's T, and other high-precision helpers.
- Probability distributions with Boost-backed implementations:
  - Gamma, Beta, Chi-squared, Non-central Chi-squared, Student’s t, Non-central Student’s t, Fisher’s F, Non-central F, Pareto, Bernoulli, Geometric, Poisson, Binomial, Negative Binomial, Cauchy, Exponential, Extreme Value (Gumbel), Normal, Logistic, Log-normal, Laplace (double exponential), Landau, Kolmogorov–Smirnov, Map-Airy, Holtsmark, Hyperexponential, Inverse Chi-Squared, Inverse Gamma, Inverse Normal (Inverse Gaussian/Wald), and Arcsine (PDF/PMF, CDF/SF, quantiles, hazards, moments).
  - Typed wrappers delegate internally to a unified runtime vtable (`Distribution.Dynamic`) for consistent behavior across precisions.
  - Continuous and discrete distributions expose `klDivergence(relativeTo:options:)`, powered by configurable quadrature/integration defaults via `Distribution.KLDivergenceOptions`.
  - Unified runtime factory to construct distributions by name at runtime (see "Dynamic Distribution Factory").
- Empirical distribution constructed directly from sample data with automatic discrete/continuous detection, KNN/KDE entropy and KL estimators, and bootstrap confidence intervals.
- High-precision mathematical constants (`π`, `e`, `√2`, Euler–Mascheroni, Catalan, ζ(3), φ, …) available through `Constants` helpers across `Float`, `Double`, and (x86_64) `Float80`.
- `CBoostBridge` target that forwards Swift calls into the vendored Boost headers under `extern/boost`.
- Architectural awareness with dedicated `Float`, `Double`, and (x86_64) `Float80` overloads plus generic `BinaryFloatingPoint` entry points.
- Re-exported Swift Numerics `Complex<T>` type with SwiftyBoost convenience APIs: arithmetic, polar helpers, and Boost-backed elementary functions (`exp`, `log`, `sin`, `cos`, `tan`, `sinh`, `cosh`, `tanh`, `atan`). Double/Float/(x86_64) Float80 specializations call bridge functions (`bs_*`).
- Fixed and adaptive quadrature helpers (`SpecialFunctions.Quadrature`) covering Gauss–Legendre, Gauss–Kronrod, tanh–sinh, sinh–sinh, and exp–sinh rules with reusable integrators, metadata, and abscissa/weight export.
- some features are implemented by the author (e.g. KL-Divergence, entropy in some cases)

## Requirements
- Swift 6.2 or later.
- Minimum platforms: macOS 11, iOS 16, watchOS 9, tvOS 16 (per `Package.swift`).
- `Float80` support requires an x86_64 runtime.

## Installation (Swift Package Manager)

Add SwiftyBoost to your project manifest:

```swift
// swift-tools-version: 6.2
import PackageDescription

let package = Package(
    name: "YourApp",
    platforms: [
        .macOS(.v11),
        .iOS(.v16),
        .watchOS(.v9),
        .tvOS(.v16)
    ],
    dependencies: [
        .package(url: "https://github.com/strike65/SwiftyBoost.git", from: "1.0.0")
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
let fall = try SpecialFunctions.falling_factorial(6.0, 3) // 6 · 5 · 4 = 120
let z    = ComplexD(0, 0.7)
let sc   = SpecialFunctions.sincc_pi(z)                  // equals sinhc_pi(0.7)

// Orthogonal polynomials
let ls = try SpecialFunctions.legendreStieltjes(5, 0.3)  // Legendre–Stieltjes E₅(0.3)
let h2 = try SpecialFunctions.hermite(2, 1.0)            // H₂(1) = 4x² - 2 → 2
let j13 = try SpecialFunctions.jacobi(n: 3, alpha: 0.5, beta: -0.25, x: 0.2)
let jPrime = try SpecialFunctions.jacobiPrime(n: 3, alpha: 0.5, beta: -0.25, x: 0.2)
let zeta = try SpecialFunctions.jacobiZeta(0.5, modulus: 0.8)
let lambda = try SpecialFunctions.heumanLambda(0.8, phi: 0.3)
let theta1 = try SpecialFunctions.jacobiTheta1(0.4, q: 0.2)
let theta3Tau = try SpecialFunctions.jacobiTheta3Tau(0.4, tau: 1.5)
let sn = try SpecialFunctions.jacobiEllipticSn(0.75, theta: 1.1)
let trio = try SpecialFunctions.jacobiElliptic(0.75, theta: 1.1)  // (sn, cn, dn)
let ai = try SpecialFunctions.airyAi(-1.0)
let biPrime = try SpecialFunctions.airyBiPrime(0.75)
let firstAiZero = try SpecialFunctions.airyAiZero(0)    // ≈ −2.3381
let aiZeros = try SpecialFunctions.airyAiZeros(startIndex: 0, count: 3)
```

```swift
// Constants (choose your precision)
let tau: Double = Constants.twoPi()
let phi: Float = Constants.phi()
#if arch(x86_64)
let rootPi80: Float80 = Constants.rootPi()
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

// Geometric(p = 0.35) — failures before first success
let geom = try Distribution.Geometric<Double>(probabibilityOfSuccess: 0.35)
let geom_mean = geom.mean

// Holtsmark(location = 0, scale = 1)
let holtsmark = try Distribution.Holtsmark<Double>(location: 0, scale: 1)
let holtsmark_pdf = try holtsmark.pdf(0.5)

// Relative entropy (KL divergence) between two exponentials
let p = try Distribution.Exponential<Double>(lambda: 1.2)
let q = try Distribution.Exponential<Double>(lambda: 0.75)
let divergence = try p.klDivergence(relativeTo: q)

// Hyperexponential mixture (three phases)
let hyperexp = try Distribution.Hyperexponential<Double>(
    probabilities: [0.25, 0.5, 0.25],
    rates: [0.5, 1.2, 3.0]
)
let hyperexpMean = hyperexp.mean

// Inverse chi-squared (ν = 8, τ = 0.5)
let invChiSq = try Distribution.InverseChiSquared<Double>(
    degreesOfFreedom: 8,
    scale: 0.5
)
let invChiSqMean = invChiSq.mean

// Inverse gamma (α = 4.5, β = 1.2)
let invGamma = try Distribution.InverseGamma<Double>(shape: 4.5, scale: 1.2)
let invGammaVariance = invGamma.variance

// Empirical distribution from raw samples (automatic discrete vs continuous detection)
let samples: [Double] = [1, 2, 2, 4]
let empirical = try Distribution.Empirical(samples: samples)
let empiricalCDF = try empirical.cdf(2)
let empiricalEntropy = empirical.entropy
let maybeEstimate = try empirical.entropyEstimate(bootstrapSamples: 32, confidenceLevel: 0.9)

// KL divergence between two empirical samples via KNN estimator
let other = try Distribution.Empirical(samples: [0, 1, 2, 2])
let kl = try empirical.klDivergence(relativeTo: other, estimator: .knn(k: 3))

// Inverse normal / Gaussian (μ = 1.2, λ = 3.4)
let invNormal = try Distribution.InverseNormal<Double>(mean: 1.2, shape: 3.4)
let invNormalQuantile = try invNormal.quantile(0.95)

// Arcsine on [0, 1]
let arcsine = try Distribution.Arcsine<Double>(minX: 0, maxX: 1)
let arcsine_pdf = try arcsine.pdf(0.2)

// Rayleigh(σ = 0.8)
let rayleigh = try Distribution.Rayleigh<Double>(scale: 0.8)
let rayleighEntropy = rayleigh.entropy

// SaS Point5 (μ = 0, c = 1)
let sasPoint5 = try Distribution.SASPoint5<Double>(location: 0, scale: 1)
let sasMedian = sasPoint5.median

// Skew-normal
let skewNormal = try Distribution.SkewNormal<Double>(location: 0.5, scale: 1.2, shape: -3)
let skewCdf = try skewNormal.cdf(0)

// Triangular on [-1, 2] with mode 0.2
let triangle = try Distribution.Triangular<Double>(lower: -1, mode: 0.2, upper: 2)
let triangleMean = triangle.mean

// Uniform on [-1, 3]
let uniform = try Distribution.Uniform<Double>(lower: -1, upper: 3)
let uniformPdf = try uniform.pdf(1)

// Weibull(k = 1.5, λ = 0.7)
let weibull = try Distribution.Weibull<Double>(shape: 1.5, scale: 0.7)
let weibullQuantile = try weibull.quantile(0.9)

```

```swift
// Quadrature
let kronrod = try SpecialFunctions.Quadrature.Integrator<Double>(rule: .gaussKronrod(points: 61))
let partial = try kronrod.integrate(over: .finite(lower: 0, upper: 1)) { x in 1.0 / (1.0 + x * x) }

let gaussian = try SpecialFunctions.Quadrature.integrate(
    using: .tanhSinh(maxRefinements: 12, tolerance: 1e-12),
    over: .finite(lower: -1, upper: 1)
) { (x: Double) in
    1.0 / (1.0 + x * x)
}
let evaluations = gaussian.functionEvaluations
```

APIs throw `SpecialFunctionError` for invalid inputs or domain violations; each case carries context such as the offending parameter name, and `invalidCombination` echoes the rejected value via its `value` payload. See DocC symbol docs for function‑specific bounds mirrored from Boost.Math.

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

// Geometric (aliases: geometric_distribution)
let dynGeom = try Distribution.Dynamic<Double>(
  distributionName: "geometric",
  parameters: ["success": 0.35]
)

// Holtsmark (aliases: holtsmark_distribution)
let dynHolts = try Distribution.Dynamic<Double>(
  distributionName: "holtsmark",
  parameters: ["loc": 0.0, "scale": 1.0]
)

// Hyperexponential (aliases: hyper_exponential, hyperexp)
let dynHyperexp = try Distribution.Dynamic<Double>(
  distributionName: "hyperexponential",
  parameters: [
    "rate0": 0.5,
    "rate1": 1.5,
    "rate2": 3.0,
    "prob0": 0.2,
    "prob1": 0.5,
    "prob2": 0.3
  ]
)

// Inverse chi-squared (aliases: inv_chi_squared, inverse_chi2)
let dynInvChiSq = try Distribution.Dynamic<Double>(
  distributionName: "inverse_chi_squared",
  parameters: [
    "df": 8.0,
    "scale": 0.5
  ]
)

// Inverse gamma (aliases: inversegamma, inv_gamma)
let dynInvGamma = try Distribution.Dynamic<Double>(
  distributionName: "inverse_gamma",
  parameters: [
    "shape": 4.5,
    "scale": 1.2
  ]
)

// Inverse normal / Gaussian (aliases: inverse_gaussian, inverse_normal, wald)
let dynInvNormal = try Distribution.Dynamic<Double>(
  distributionName: "inverse_gaussian",
  parameters: [
    "mean": 1.2,
    "lambda": 3.4
  ]
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
