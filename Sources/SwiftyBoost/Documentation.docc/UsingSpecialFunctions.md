# Using Special Functions

Discover common entry points and best practices when calling into ``SwiftyBoost``.

## Factorials and Combinatorics

- Use ``SpecialFunctions/factorial(_:)->Double`` for `UInt32` inputs that fit within Boost’s factorial table; the API throws when a value exceeds the backend’s supported range.
- Prefer generic overloads such as ``SpecialFunctions/factorial(_:)->Double`` when working with `BinaryFloatingPoint`, and choose the type according to the precision your algorithm requires.

## Elliptic Integrals

- ``SpecialFunctions/bernoulli_b2n_f(_:)`` and companions provide Legendre-form integrals; Carlson symmetric variants such as ``SpecialFunctions/carlsonRF(_:_:_:)->T`` offer improved convergence for certain domains.
- Check availability annotations: some overloads are limited to `Float`/`Double` while others extend to `Float80` on x86_64.

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
