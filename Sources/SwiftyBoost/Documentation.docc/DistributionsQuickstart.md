# Distributions Quickstart

Learn how to construct and reuse distribution objects for efficient evaluation.

## Construct once, reuse

```swift
import SwiftyBoost

// Gamma( k = 2.5, θ = 1.2 )
let g = try Distribution.Gamma<Double>(shape: 2.5, scale: 1.2)

let x: Double = 0.8
let pdf = try g.pdf(x)      // f(x)
let cdf = try g.cdf(2.0)    // F(x)
let sf  = try g.sf(2.0)     // 1 − F(x)
let q95 = try g.quantile(0.95)

// Moments
_ = g.mean
_ = g.variance
_ = g.mode
```

```swift
// Student's t with ν = 5
let t = try Distribution.StudentT<Double>(degreesOfFreedom: 5)
let p0  = try t.pdf(0.0)
let cdf = try t.cdf(2.0)
let q97 = try t.quantile(0.975)

// Planning helper (Boost):
let nu = Distribution.StudentT<Double>.findDegreesOfFreedom(
  differenceFromMean: 1.0, alpha: 0.05, beta: 0.2, sd: 1.0)
```

```swift
// Beta(α = 2, β = 5)
let beta = try Distribution.Beta<Double>(alpha: 2, beta: 5)
let betaCDF = try beta.cdf(0.3)
let betaMean = beta.mean

// Chi-squared with ν = 8
let chi2 = try Distribution.ChiSquared<Double>(degreesOfFreedom: 8)
let chi2SF = try chi2.sf(12.0)
let chi2Mode = chi2.mode

// Bernoulli(p = 0.3)
let bern = try Distribution.Bernoulli<Double>(p: 0.3)
let bernPMF1 = try bern.pdf(1)
let bernEntropy = bern.entropy
```

```swift
// Geometric(p = 0.35) — failures before first success
let geom = try Distribution.Geometric<Double>(probabibilityOfSuccess: 0.35)
let geomPMF2 = try geom.pdf(2)
let geomMean = geom.mean
```

```swift
// Fisher’s F (df1 = 10, df2 = 20)
let f = try Distribution.FisherF<Double>(degreesOfFreedom1: 10, degreesOfFreedom2: 20)
let f_pdf = try f.pdf(2.0)
let f_sf  = try f.sf(3.5)
let f_q95 = try f.quantile(0.95)
```

```swift
// Arcsine distribution on [0, 1]
let a = try Distribution.Arcsine<Double>(minX: 0, maxX: 1)
let a_pdf = try a.pdf(0.2)
let a_cdf = try a.cdf(0.8)

// Binomial(n = 12, p = 0.35)
let bin = try Distribution.Binomial<Double>(
  numberOfTrials: 12,
  probabibilityOfSuccess: 0.35
)
let binPMF4 = try bin.pdf(4)
let binVariance = bin.variance

// Exponential(λ = 1.2)
let exp = try Distribution.Exponential<Double>(lambda: 1.2)
let expTail = try exp.sf(2.0)
let expMean = exp.mean

// Cauchy(location = 0.5, scale = 1.5)
let cauchy = try Distribution.Cauchy<Double>(location: 0.5, scale: 1.5)
let cauchyQuantile = try cauchy.quantile(0.95)

// Holtsmark (location = 0, scale = 1)
let holts = try Distribution.Holtsmark<Double>(location: 0, scale: 1)
let holtsPdf = try holts.pdf(0.5)

// Gumbel (extreme value) with loc = 0, scale = 0.75
let gumbel = try Distribution.ExtremeValueGumpel<Double>(location: 0, scale: 0.75)
let gumbelCdf = try gumbel.cdf(1.2)
let gumbelMode = gumbel.mode
```

## Implementation model

Typed distribution wrappers (Gamma, Student’s t, Fisher’s F, Arcsine, Geometric, Binomial, Cauchy, Holtsmark, Exponential, Extreme Value/Gumbel) delegate to ``Distribution/Dynamic``, a unified runtime vtable backed by Boost.Math via the C bridge. Each instance constructs its Boost backend once through the dynamic factory and reuses it for all evaluations to ensure performance and consistent numerical policies.

## Error handling

- Constructors validate parameters and throw when they are not positive.
- Bernoulli, Geometric, and Binomial enforce probabilities in `[0, 1]`; Cauchy, Exponential, and Extreme Value enforce strictly positive scale parameters.
- Methods that accept probabilities (`quantile`, `quantileComplement`) validate inputs in `[0, 1]`.
- Gamma methods throw for negative `x`.

## Precision choices

All distribution types are generic over `BinaryFloatingPoint`. Use `Double` by default for accuracy; use `Float` for performance. On x86_64, `Float80` is available and routes to extended precision backends.
