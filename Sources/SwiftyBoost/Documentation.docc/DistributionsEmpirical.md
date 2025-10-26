# Empirical Distribution

``Distribution/Empirical`` turns a finite sample into a fully fledged distribution object that conforms to ``DistributionProtocol``. The implementation is entirely Swift, adding features that are not provided by Boost.Math:

- Automatic detection of discrete (lattice) versus continuous samples.
- Smoothed probability masses for discrete data (Add-α with Miller–Madow correction).
- K-nearest neighbour (KNN) and Gaussian kernel density (KDE) estimators for entropy and KL divergence when the sample behaves continuously.
- Bootstrap helpers that produce confidence intervals for entropy or KL divergence via percentile or BCa resampling.
- Heuristics for multimodality detection (lattice scans or KDE peak counting) surfaced through ``Distribution/Empirical/isLikelyMultimodal``.

## Creating an empirical distribution

```swift
let samples: [Double] = [-1.2, -0.6, -0.2, 0.0, 0.35, 0.8, 1.2]
let empirical = try Distribution.Empirical(samples: samples)

empirical.isDiscrete          // false → treated as continuous
let density = try empirical.pdf(0.1)
let cdf    = try empirical.cdf(0.5)
let mean   = empirical.mean
let entropy = empirical.entropy   // Uses KNN or KDE depending on the sample
```

For discrete (lattice) samples the API behaves like other discrete distributions:

```swift
let counts = try Distribution.Empirical(samples: [1, 2, 2, 4])
counts.isDiscrete              // true
counts.latticeStep             // Optional(1)
try counts.pdf(2)              // smoothed PMF value
try counts.cdf(2)              // cumulative mass
counts.entropy                 // Miller–Madow corrected entropy
```

## KNN/KDE estimators

The type automatically selects an estimator when possible:

- Discrete data use smoothed frequencies.
- Continuous data fall back to a KNN entropy estimator (default `k = 3`) and KDE densitiy estimates for PDFs.
- When the KNN estimator is not applicable (for example, very small samples), the implementation falls back to a KDE plug-in.

You can choose explicit estimators through ``Distribution/Empirical/DensityEstimator`` when calling
``Distribution/Empirical/klDivergence(relativeTo:estimator:)`` or the bootstrap helpers.

```swift
let p = try Distribution.Empirical(samples: [-1.0, -0.5, 0.0, 0.5])
let q = try Distribution.Empirical(samples: [-0.8, -0.4, 0.2, 0.6])
let knn = try p.klDivergence(relativeTo: q, estimator: .knn(k: 2))
let kde = try p.klDivergence(relativeTo: q, estimator: .kdeGaussian(bandwidth: nil))
```

## Bootstrap intervals

Confidence intervals are available for both entropy and KL divergence:

```swift
let bootstrap = try empirical.entropyEstimate(
    estimator: .automatic,
    bootstrapSamples: 64,
    confidenceLevel: 0.9,
    method: .percentile
)

if let interval = bootstrap.confidenceInterval {
    print("Entropy ≈ \(bootstrap.value) in [\(interval.lower), \(interval.upper)]")
}
```

For KL divergence:

```swift
let klBootstrap = try p.klDivergenceEstimate(
    relativeTo: q,
    estimator: .knn(k: 3),
    bootstrapSamples: 64,
    confidenceLevel: 0.9
)
```

## Multimodality heuristic

``Distribution/Empirical/isLikelyMultimodal`` flags samples that appear to have more than one mode. Lattice samples look for multiple local maxima; continuous samples perform a KDE scan across a grid. This is intended as a coarse diagnostic and does not replace formal tests.

```swift
let unimodal = try Distribution.Empirical(samples: stride(from: -1.0, through: 1.0, by: 0.1))
let multimodal = try Distribution.Empirical(samples: [-3, -2.5, -2.4, -0.2, 0, 0.25, 2.2, 2.7, 3.1])

unimodal.isLikelyMultimodal    // false
multimodal.isLikelyMultimodal  // true
```

## When to use

Use ``Distribution/Empirical`` when you have raw samples and need the convenience of the common distribution protocol without fitting a parametric distribution. The type reuses the same API surface as the Boost-backed wrappers, making it easy to mix empirical and analytic distributions when evaluating PDFs, CDFs, moments, or KL divergence.

