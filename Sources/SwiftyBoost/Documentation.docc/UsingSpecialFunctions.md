# Using Special Functions

Discover common entry points and best practices when calling into ``SwiftyBoost``.

## Factorials and Combinatorics

- Use ``SpecialFunctions/factorial(_:)->Double`` for `UInt32` inputs that fit within Boost’s factorial table; the API throws when a value exceeds the backend’s supported range.
- Prefer generic overloads such as ``SpecialFunctions/factorial(_:)->Double`` when working with `BinaryFloatingPoint`, and choose the type according to the precision your algorithm requires.

## Elliptic Integrals

- ``SpecialFunctions/bernoulli_b2n_f(_:)`` and companions provide Legendre-form integrals; Carlson symmetric variants such as ``SpecialFunctions/ellint_rf(x:y:z:)`` offer improved convergence for certain domains.
- Check availability annotations: some overloads are limited to `Float`/`Double` while others extend to `Float80` on x86_64.

## Error Handling and Validation

- All throwing APIs emit ``SpecialFunctionError`` cases; handle `parameterExceedsMaximumIntegerValue` to detect when inputs fall outside documented bounds.
- Use non-throwing helpers (for example ``SpecialFunctions/sinhc(_:)``) when you need guarantees about domain coverage and Boost already constrains the evaluation.

## Testing and Verification

- Compare Swift results against bridge helpers such as `bs_rising_factorial` in `CBoostBridge` when you extend functionality or add regression coverage.
- For floating-point comparisons, prefer absolute or relative tolerance checks instead of equality.
