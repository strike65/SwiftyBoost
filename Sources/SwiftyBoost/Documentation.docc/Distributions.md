# Distributions

This section documents probability distributions exposed by SwiftyBoost using Boost.Math backends. Each distribution type is constructed once and reuses an internal opaque handle to a Boost distribution object for efficient evaluation.

## Implemented

- ``Distribution/Gamma``
- ``Distribution/Beta``
- ``Distribution/ChiSquared``
- ``Distribution/StudentT``
- ``Distribution/FisherF``
- ``Distribution/Bernoulli``
- ``Distribution/Geometric``
- ``Distribution/Binomial``
- ``Distribution/Cauchy``
- ``Distribution/Holtsmark``
- ``Distribution/Exponential``
- ``Distribution/ExtremeValueGumpel``
- ``Distribution/Arcsine``

### Tutorials

- <doc:DistributionsQuickstart>

## Notes

- All public APIs are thin wrappers over Boost functions via the C bridge (`bs_` prefixed symbols). No new algorithms are implemented in Swift.
- Generic types support `Double`, `Float`, and, on x86_64, `Float80` (falling back to `Double` where not available).
- Methods that accept probabilities validate input in `[0, 1]`. Domain errors are thrown using ``SwiftyBoost/SpecialFunctionError``.
- ``Distribution/Bernoulli`` and ``Distribution/Geometric`` are discrete but keep the continuous-style interface for consistency; treat `pdf` as the PMF.
- ``Distribution/Holtsmark`` and ``Distribution/Cauchy`` share the location–scale heavy-tail family; both expose the same dynamic aliases documented in <doc:DynamicFactory>.
