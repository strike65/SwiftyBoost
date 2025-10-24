# Quadrature Helpers — Usage and Design Guide

This document outlines the Boost-backed quadrature facilities exposed through `SpecialFunctions.Quadrature`. It covers the Swift-facing API, supported rules and parameters, metadata, error semantics, and the underlying C bridge so future contributors can extend the integration stack confidently.

## Overview

- Swift entry point: `Sources/SwiftyBoost/Math/SpecialFunctions/Quadrature.swift`
- C ABI: `Sources/CBoostBridge/include/Quadrature/bs_quadrature.h`
- Implementation: `Sources/CBoostBridge/impl/Quadrature/bs_quadrature.hxx`

The Swift wrapper manages opaque C handles, converts Swift closures into C callbacks, and surfaces result metadata (value, estimated error, L₁ norm, iteration count, function evaluation count, convergence flag). All APIs are generic over `Scalar: Real & BinaryFloatingPoint` with concrete support for `Float`, `Double`, and `Float80` (the latter on x86 targets).

## Goals

- Provide reusable integrators that amortize the cost of constructing Boost quadrature objects.
- Offer a one-shot integration helper for simple use cases.
- Surface rich diagnostics (metadata, function-call counts) without exposing C handles directly.
- Allow callers to export abscissa/weights for fixed Gaussian rules.
- Validate point counts, refinement depths, tolerances, and intervals before hitting the C bridge.

## Swift API Surface

- ``SpecialFunctions/Quadrature`` namespace
  - ``SpecialFunctions/Quadrature/Rule`` — describes the integration rule.
  - ``SpecialFunctions/Quadrature/Interval`` — `.automatic` (native domain) or `.finite(lower:upper:)`.
  - ``SpecialFunctions/Quadrature/Integrator`` — reusable handle with `integrate(over:integrand:)` and `copyAbscissaWeights(into:weights:)`.
  - ``SpecialFunctions/Quadrature/integrate(using:over:integrand:)`` — convenience one-shot entry point.
  - ``SpecialFunctions/Quadrature/Result`` — integral estimate plus diagnostics.
  - ``SpecialFunctions/Quadrature/Metadata`` — describes the instantiated backend (rule, C enum `QuadratureType`, precision, node count, adaptivity, infinite-bound support hint).
  - ``SpecialFunctions/Quadrature/Error`` — Swift error type (`invalidConfiguration`, `invalidInterval`, `backendUnavailable`).

Integrators are `Sendable`, so you can move them across tasks, but each `integrate` call is synchronous and should not be performed concurrently on the same instance.

## Supported Rules and Parameters

| Rule | Swift case | Natural interval | Notes |
| --- | --- | --- | --- |
| Gauss–Legendre | `.gaussLegendre(points:)` | `[-1, 1]` | Supported point counts: 7, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100. Use `.finite(lower:upper:)` for arbitrary `[a, b]` via change of variables. |
| Gauss–Kronrod | `.gaussKronrod(points:)` | `[-1, 1]` | Supported point counts: 15, 21, 31, 41, 51, 61. Exposes embedded Gauss/Kronrod nodes for error estimation workflows. |
| Tanh–Sinh | `.tanhSinh(maxRefinements:tolerance:)` | `[a, b]` (default `[-1, 1]`) | Adaptive double-exponential rule for finite intervals. Defaults: 10 refinements, absolute tolerance `1e-9`. Metadata flag `supportsInfiniteBounds` is `false`. |
| Sinh–Sinh | `.sinhSinh(maxRefinements:tolerance:)` | `(-∞, ∞)` | Adaptive double-exponential rule optimized for even integrands. Metadata flag `supportsInfiniteBounds` is `true`. |
| Exp–Sinh | `.expSinh(maxRefinements:tolerance:)` | `[0, ∞)` | Adaptive double-exponential rule for semi-infinite domains (`[0, ∞)`). Metadata `supportsInfiniteBounds` is `true`. |

Adaptive rules accept optional overrides:

- `maxRefinements`: positive integer (default 10). Must fit in `Int32`.
- `tolerance`: positive absolute tolerance (default `1e-9`). Downcast to the target scalar type.

Invalid configurations throw ``SpecialFunctions/Quadrature/Error/invalidConfiguration(_:)`` before reaching the bridge.

## Result Metadata and Diagnostics

- ``SpecialFunctions/Quadrature/Result/value`` — integral estimate.
- ``SpecialFunctions/Quadrature/Result/estimatedError`` — currently zero-filled; Boost’s default facilities do not expose an error estimate for these rules. Field reserved for future enhancement.
- ``SpecialFunctions/Quadrature/Result/l1Norm`` — zero-filled placeholder mirroring `QuadratureResult*` structs.
- ``SpecialFunctions/Quadrature/Result/iterations`` — 1 for fixed rules, 0 for adaptive ones.
- ``SpecialFunctions/Quadrature/Result/functionEvaluations`` — node count for fixed rules; for adaptive rules the bridge increments this counter on each callback invocation.
- ``SpecialFunctions/Quadrature/Result/converged`` / ``SpecialFunctions/Quadrature/Result/didConverge`` — `true` for fixed rules, `true` when adaptive rules complete without throwing.
- ``SpecialFunctions/Quadrature/Metadata`` — includes the original Swift rule, the bridge enums (`QuadratureType`, `QuadraturePrecision`), reported node count (`points`, 0 for adaptive), `isAdaptive`, and `supportsInfiniteBounds`.

## Abscissa and Weights

Fixed Gaussian rules expose compile-time node/weight tables via `quad_get_abscissa_weights_*`. To fetch them:

```swift
let integrator = try SpecialFunctions.Quadrature.Integrator<Double>(rule: .gaussLegendre(points: 30))
let result = try integrator.integrate { x in 1.0 - x * x }  // Automatic interval [-1, 1]

var nodes = Array(repeating: 0.0, count: result.metadata.points)
var weights = nodes
let success = nodes.withUnsafeMutableBufferPointer { ab in
    weights.withUnsafeMutableBufferPointer { wt in
        integrator.copyAbscissaWeights(into: ab, weights: wt)
    }
}
```

- Buffers must share the same length and match `result.metadata.points`.
- Adaptive rules return `false` because they do not expose a fixed tableau.

## Error Handling

- ``SpecialFunctions/Quadrature/Error/invalidConfiguration(_:)`` — non-positive point counts, unsupported counts, non-positive refinement limits, or non-positive tolerances.
- ``SpecialFunctions/Quadrature/Error/backendUnavailable(_:)`` — rule missing for requested precision (e.g., `Float80` on arm64).
- ``SpecialFunctions/Quadrature/Error/invalidInterval(lower:upper:)`` — bounds are non-finite or `lower >= upper`. Bounds supplied to `.finite` are downcast to the target scalar (`Float`/`Double`/`Float80`); ensure they remain finite after conversion.

Exceptions thrown by Boost during evaluation propagate as Swift errors originating at the bridge boundary.

## Precision and Platforms

- `Double` and `Float` are always available.
- `Float80` support is limited to x86 architectures. Attempting to create an integrator with `Scalar == Float80` on unsupported platforms throws `invalidConfiguration`.
- Metadata includes the `QuadraturePrecision` enum to reflect the actual scalar used by the bridge.

## Extending the Bridge

To add a new rule:

1. Extend `QuadratureType` in `bs_quadrature.h` and update `quad_type_to_string`, `quad_is_adaptive`, and `quad_supports_infinite_bounds`.
2. Implement a new concrete class in `bs_quadrature.hxx` deriving from `QuadratureBase<Real>` for each supported precision.
3. Update the factory switch in `create_quadrature<Real>` (search for `case QUAD_GAUSS_LEGENDRE`) to instantiate the new implementation per point count or parameter set.
4. Expose creation functions in the C header (e.g., `quad_gauss_jacobi_create_d`) and add matching cases in the Swift bridge (`QuadratureBridgeDouble.createHandle`).
5. Document the new rule: update this file, README, DocC (`Quadrature.md`), and the changelog.
6. Add regression tests that compare results against Boost’s direct API or known integrals.

## Testing Tips

- Add unit tests under `Tests/SwiftyBoostTests` targeting representative integrals (polynomials, Gaussians, oscillatory functions) and verifying known values, symmetry, or relationships between rules.
- For adaptive rules, assert monotonic behavior (e.g., increasing `maxRefinements` reduces `estimatedError` when a future bridge surfaces it) and capture function evaluation counts.
- Validate that `copyAbscissaWeights` returns sorted symmetric nodes for Gauss–Legendre and that weights sum to 2 (the length of `[-1, 1]`).

## Related Documentation

- README usage examples (quadrature section).
- DocC page `<doc:Quadrature>` for API-level guidance.
- `Sources/SwiftyBoost/Documentation.docc/UsingSpecialFunctions.md` — quick-start tips referencing quadrature helpers.
