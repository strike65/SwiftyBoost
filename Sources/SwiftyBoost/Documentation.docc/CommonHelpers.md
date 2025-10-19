# Common Helpers and Stable Transforms

Selected numerically stable helpers that forward to Boost via the C bridge.

## Reciprocal Square Root

- ``SwiftyBoost/SpecialFunctions/rsqrt(_:)->T`` computes `1/sqrt(x)` with domain checks.
  - Domain: `x ≥ 0` and finite; throws on negative or non-finite input.
  - Backends: `bs_rsqrt` (double); float/long double variants where available.

## Stable Elementaries

- ``SwiftyBoost/SpecialFunctions/expm1(_:)->T``: computes `exp(x) − 1` stably for small `x`; throws if `x` is not finite.
- ``SwiftyBoost/SpecialFunctions/log1p(_:)->T``: computes `log(1 + x)`; domain `x > −1` and finite; throws otherwise.
- ``SwiftyBoost/SpecialFunctions/log1pmx(_:)->T``: computes `log(1 + x) − x`; domain `x > −1` and finite; throws otherwise.
- ``SwiftyBoost/SpecialFunctions/powm1(_:_:)->T``: computes `x^y − 1`; throws on non-finite inputs, or when `x < 0` and `y` is non-integer.
- ``SwiftyBoost/SpecialFunctions/cbrt(_:)->T``: real cube root; throws if input is not finite.
- ``SwiftyBoost/SpecialFunctions/sqrt1pm1(x:)->T``: computes `sqrt(1 + x) − 1`; non-throwing; delegates to Boost.

All helpers provide generic `T: BinaryFloatingPoint` entry points, with `Float`, `Double`, and (x86_64) `Float80` specializations where appropriate. Implementations are thin wrappers over `bs_*` bridge functions and introduce no novel algorithms.

