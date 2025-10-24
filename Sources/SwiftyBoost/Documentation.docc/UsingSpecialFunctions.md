# Using Special Functions

Discover common entry points and best practices when calling into ``SwiftyBoost``.

## Factorials and Combinatorics

- Use ``SpecialFunctions/factorial(_:)->Double`` for `UInt32` inputs that fit within Boost’s factorial table; the API throws when a value exceeds the backend’s supported range.
- Prefer generic overloads such as ``SpecialFunctions/factorial(_:)->Double`` when working with `BinaryFloatingPoint`, and choose the type according to the precision your algorithm requires.
- Rising factorial helpers ``SpecialFunctions/rising_factorial(_:_:)->Double`` / ``SpecialFunctions/pochhammer(_:_:)->Double`` (and Float/Float80/generic overloads) guard against non-finite bases, detect zeros when the sequence crosses 0, and throw ``SpecialFunctionError/invalidCombination`` if the result would overflow the target type.
- Falling factorial helpers ``SpecialFunctions/falling_factorial(_:_:)->Double`` (with matching Float/Float80/generic overloads) offer the same guards: they reject non-finite `x`, short-circuit to 0 when `x` is a non-negative integer smaller than `n`, and pre-check magnitude via `ln Γ(x+1) − ln Γ(x−n+1)` to catch overflow/underflow.

## Orthogonal Polynomials

- Legendre families:
  - ``SpecialFunctions/legendreP(_:_:)->Double`` and ``SpecialFunctions/associatedLegendreP(_:_:_:)->Double`` evaluate Pₙ(x) and Pₙᵐ(x).
  - ``SpecialFunctions/legendrePPrime(_:_:)->Double`` and ``SpecialFunctions/legendrePZeros(degree:)->[Double]`` cover x-derivatives and root finding.
  - ``SpecialFunctions/legendreStieltjes(_:_:)->Double`` and companions (prime, zeros, norm-squared) expose Boost’s Gauss–Kronrod Legendre–Stieltjes polynomials for quadrature construction.
- Hermite polynomials:
  - ``SpecialFunctions/hermite(_:_:)->Double`` evaluates Hₙ(x); ``SpecialFunctions/hermiteNext(_:_:hn:hnm1:)`` advances the three-term recurrence without recomputing the full sequence.
- Jacobi family:
  - ``SpecialFunctions/jacobi(n:alpha:beta:x:)`` evaluates Pₙ^{(α,β)}(x); ``SpecialFunctions/jacobiPrime(n:alpha:beta:x:)``, ``SpecialFunctions/jacobiDoublePrime(n:alpha:beta:x:)``, and ``SpecialFunctions/jacobiDerivative(n:alpha:beta:x:k:)`` expose first/second/k-th derivatives directly.
- Airy functions:
  - ``SpecialFunctions/airyAi(_:)->Double`` / ``SpecialFunctions/airyBi(_:)->Double`` and their derivatives ``SpecialFunctions/airyAiPrime(_:)->Double`` / ``SpecialFunctions/airyBiPrime(_:)->Double`` provide Ai/Bi and slopes with Float/Double/(x86) Float80 overloads.
  - Zero helpers ``SpecialFunctions/airyAiZero(_:)`` / ``SpecialFunctions/airyBiZero(_:)`` return individual roots, while ``SpecialFunctions/airyAiZeros(startIndex:count:)`` and ``SpecialFunctions/airyBiZeros(startIndex:count:)`` fetch consecutive roots in one call.
- Gegenbauer (ultraspherical) polynomials and derivatives are available via ``SpecialFunctions/gegenbauer(n:lambda:x:)`` and related helpers.
- All APIs share validation conventions: integer indices must be non-negative and fit in `UInt32`, while parameters and evaluation points must be finite. Each overload is mirrored for `Float`, `Double`, generic `BinaryFloatingPoint`, and (on x86) `Float80`.

## Elliptic Integrals

- ``SpecialFunctions/bernoulli_b2n_f(_:)`` and companions provide Legendre-form integrals; Carlson symmetric variants such as ``SpecialFunctions/carlsonRF(_:_:_:)->T`` offer improved convergence for certain domains.
- Check availability annotations: some overloads are limited to `Float`/`Double` while others extend to `Float80` on x86_64.
- ``SpecialFunctions/jacobiZeta(_:modulus:)`` exposes the Jacobi Zeta function Z(φ, k), complementing the incomplete Legendre forms with a direct route to Z via Boost’s implementation.
- ``SpecialFunctions/heumanLambda(_:phi:)`` evaluates Heuman’s lambda Λ(φ, k); overloads exist for `Float`, generic `BinaryFloatingPoint`, and (x86) `Float80`.
- Jacobi elliptic functions are available via ``SpecialFunctions/jacobiElliptic(_:theta:)`` (returns `(sn, cn, dn)`) and the individual helpers ``SpecialFunctions/jacobiEllipticSn(_:theta:)`` / `Cn` / `Dn` plus the quotient variants `Sd`, `Sc`, `Ns`, `Nd`, `Nc`, `Ds`, `Dc`, `Cs`, `Cd` with matching Float/Float80 overloads.
- Jacobi theta functions provide θ₁–θ₄ through ``SpecialFunctions/jacobiTheta1(_:q:)`` / ``SpecialFunctions/jacobiTheta2(_:q:)`` / ``SpecialFunctions/jacobiTheta3(_:q:)`` / ``SpecialFunctions/jacobiTheta4(_:q:)`` with Float/Double/generic/(x86) Float80 support and nome restriction `0 < q < 1`. Half-period ratio variants ``SpecialFunctions/jacobiTheta1Tau(_:tau:)`` / ``SpecialFunctions/jacobiTheta2Tau(_:tau:)`` / ``SpecialFunctions/jacobiTheta3Tau(_:tau:)`` / ``SpecialFunctions/jacobiTheta4Tau(_:tau:)`` are available when working directly with τ; all inputs must be finite.

## Error Handling and Validation

- All throwing APIs emit ``SpecialFunctionError`` cases; ``invalidCombination`` includes a `value` payload echoing the offending input when available, while `parameterExceedsMaximumIntegerValue` flags arguments outside documented bounds.
- Use non-throwing helpers (for example `sinhc_pi(_:)`) when you need guarantees about domain coverage and Boost already constrains the evaluation.

## Testing and Verification

- Compare Swift results against bridge helpers such as `bs_rising_factorial` in `CBoostBridge` when you extend functionality or add regression coverage.
- For floating-point comparisons, prefer absolute or relative tolerance checks instead of equality.

## Cardinal Sine (π-Normalized)

- Use ``SwiftyBoost/SpecialFunctions/sinc_pi(_:)->T`` and ``SwiftyBoost/SpecialFunctions/sinhc_pi(_:)->T`` for real inputs. Both handle the removable singularity at 0 by returning 1.
- Complex variants ``SwiftyBoost/SpecialFunctions/sincc_pi(_:)`` and ``SwiftyBoost/SpecialFunctions/sinhcc_pi(_:)`` integrate with ``ComplexD``/``ComplexF``/(x86_64) ``ComplexL``.
- Identity: `sincc_pi(i y) = sinhcc_pi(y)`.

## Reciprocal Square Root and Stable Transforms

- ``SwiftyBoost/SpecialFunctions/rsqrt(_:)->T`` computes `1/√x` with domain checks (`x ≥ 0`).
- Stability helpers: ``SwiftyBoost/SpecialFunctions/expm1(_:)->T``, ``SwiftyBoost/SpecialFunctions/log1p(_:)->T``, ``SwiftyBoost/SpecialFunctions/log1pmx(_:)->T``, ``SwiftyBoost/SpecialFunctions/powm1(_:_:)->T``, ``SwiftyBoost/SpecialFunctions/cbrt(_:)->T``, and ``SwiftyBoost/SpecialFunctions/sqrt1pm1(x:)->T``.
