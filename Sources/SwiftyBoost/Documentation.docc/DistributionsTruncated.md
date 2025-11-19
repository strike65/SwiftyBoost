# Truncated Distributions

`Distribution.TruncatedDistribution` clamps any continuous ``Distribution/DistributionProtocol`` conformer to a smaller interval without rewriting its PDF/CDF/quantile logic. The type keeps a reference to the base distribution, precomputes `cdf(lower)` / `cdf(upper)`, and caches the normalization constant so subsequent evaluations simply scale the underlying Boost-backed values.

## Requirements

- The base conformer must be continuous (`isDiscrete == false`).
- Supply at least one of `lower` or `upper`; the initializer intersects those bounds with the base support and throws if the resulting range is empty.
- Bounds can be semi-infinite (e.g. only a lower cutoff) as long as the base distribution supports that region.

If the base CDF reports zero probability between the requested bounds, the initializer throws ``SwiftyBoost/TruncationError/invalidNormalizationConstant`` so you immediately know the truncation interval is unsupported.

## Usage

```swift
import SwiftyBoost

let base = try Distribution.Normal<Double>(mean: 0, sd: 1)
let truncated = try Distribution.TruncatedDistribution(
    base: base,
    lower: -1.0,
    upper:  1.0
)

// PDF/CDF/SF scale the base values by the cached normalization constant.
let pdf = try truncated.pdf(0.2)
let cdf = try truncated.cdf(0.0)
let sf  = try truncated.sf(0.0)

// Quantiles translate the truncated probability back into the base CDF.
let q90 = try truncated.quantile(0.90)
let uq10 = try truncated.quantileComplement(0.10)

// Metadata exposes the truncated support directly.
let support = truncated.range        // (lower: -1.0, upper: 1.0)
let median = truncated.median        // Cached at init via the base quantile.
```

## Notes

- Values outside `[supportLowerBound, supportUpperBound]` always return a PDF of zero, a CDF of zero/one, and a survival function of one/zero respectively.
- ``Distribution/TruncatedDistribution`` conforms to ``Distribution/DistributionProtocol``, so it can participate in KL divergence just like the built-in Boost wrappers and ``Distribution/Empirical``.
- Because the type simply wraps an existing distribution, all Boost-backed precision options remain available (Float, Double, and Float80 on x86_64) and you can reuse the same base instance across multiple truncated views.
