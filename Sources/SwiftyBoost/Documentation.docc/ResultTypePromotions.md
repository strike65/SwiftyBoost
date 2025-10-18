# Result-Type Promotions

Learn how SwiftyBoost chooses the result type when function arguments mix `Float`, `Double`, and (x86_64) `Float80`.

## Policy

- Float ↔ Double → Double
  - If a function has more than one floating argument and you pass a mix of `Float` and `Double`, the operation runs in `Double` and returns `Double`.
- Any presence of Float80 → Float80 (x86_64)
  - On x86_64, if any floating argument is `Float80`, the operation runs in extended precision and returns `Float80`.
- Single-argument functions
  - For functions with a single floating argument, the result type matches the argument type. No promotion applies.

All numerical work is delegated to Boost.Math via the `bs_*` bridge; promotions are thin Swift overloads that convert arguments and forward to the appropriate backend.

## APIs With Mixed Promotions

- Gamma family
  - `gammaRatio(a, b)`, `gammaDeltaRatio(a, delta)`
  - Incomplete/regularized and inverses: `incompleteGammaLower/Upper`, `regularizedGammaP/Q`, `regularizedGammaPInv/QInv`, `regularizedGammaPDerivative`
- Beta family
  - `beta(a, b)`
  - Incomplete/regularized and inverses: `incompleteBetaUnnormalized`, `regularizedIncompleteBeta`, `complementaryRegularizedIncompleteBeta`, `inverseRegularizedIncompleteBeta`, `inverseComplementaryRegularizedIncompleteBeta`
  - Solvers/derivative: `solveAForRegularizedIncompleteBeta`, `solveBForRegularizedIncompleteBeta`, `regularizedIncompleteBetaDerivative`
- Bessel (cylindrical + derivatives)
  - `besselJ`, `besselY`, `modifiedBesselI`, `modifiedBesselK`, and their `Prime` derivatives
  - Integer-order convenience overloads follow the same promotion rules
- Elliptic integrals
  - Legendre forms: `incompleteEllipticIntegralF/E`, `incompleteEllipticIntegralPi`, `completeEllipticIntegralPi` (complete K/E are single-argument)
  - Carlson forms: `carlsonRC`, `carlsonRF`, `carlsonRD`, `carlsonRJ`, `carlsonRG`
- Chebyshev series (Clenshaw)
  - `chebyshevClenshawRecurrence(coefficients:x:)` supports mixing coefficient- and `x`-types ([Float] vs Double, etc.)
- Owen’s T
  - `owensT(h:a:)` supports all mixed pairs
- Spherical harmonics
  - `sphericalHarmonic(n:m:theta:phi:)`, `Y_real`, and `Y_imag` support mixed angle types; Float80 mixes return `ComplexL`/`Float80`.
- Hypergeometric
  - `hypergeometric1F0`, `hypergeometric0F1` (full Float/Double and Float80 mixes)
  - `hypergeometric1F1`, `hypergeometric2F0` (select single-parameter mixes to avoid overload explosion)
  - `hypergeometricPFQ` (arrays): `[Float]` + Double z → Double; Float80 variants return `Float80` (x86_64)

## Single-Argument Functions (No Promotion)

- Legendre `Pₙ(x)`, `Pₙᵐ(x)`, `Pₙ′(x)`; Laguerre `Lₙ(x)`, `Lₙ^m(x)`; Lambert W; error functions (`erf`, `erfc`, `erf⁻¹`); exponential integrals `Ei(x)` and `Eₙ(x)`; digamma/polygamma/zeta; cardinal B-splines and derivatives.

Choose the overload for the precision you need: `Float`, `Double`, or (x86_64) `Float80`.

## Guidance

- Pick the narrowest precision that meets your accuracy targets. For many tasks, `Double` is the right default; use `Float80` on x86_64 when you genuinely need more than 53 bits of mantissa.
- Prefer explicit typed overloads in tight loops to avoid incidental conversions.
- Use mixed promotions to simplify call sites when arguments naturally arrive as different types—e.g., a `Float` parameter and a `Double` abscissa—without manual casts.
- Be mindful of functions with restricted domains (e.g., `Y_v(x)`, `K_v(x)`, incomplete/regularized gamma/beta). Domain checks remain in effect regardless of promotion.
- For divergent series such as 2F0, interpret results carefully; promotions don’t change analytic behavior.

## Notes

- Platform: Float80 promotions only compile on x86_64 architectures.
- Backend fidelity: promotions only adapt Swift argument types. All computations are still performed by Boost.Math via `bs_*` bridge functions.

