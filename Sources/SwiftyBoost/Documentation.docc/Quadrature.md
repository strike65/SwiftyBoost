# Quadrature

Boost-backed Gaussian and double-exponential quadrature helpers with reusable integrators, metadata, and abscissa accessors.

## Available Rules

- ``SpecialFunctions/Quadrature/Rule/gaussLegendre(points:)`` — fixed Gauss–Legendre nodes on `[-1, 1]`; supported point counts are 7, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, and 100. Use ``SpecialFunctions/Quadrature/Interval/finite(lower:upper:)`` with an affine change of variables for arbitrary finite bounds.
- ``SpecialFunctions/Quadrature/Rule/gaussKronrod(points:)`` — fixed Gauss–Kronrod extensions (15, 21, 31, 41, 51, 61 points). Returns abscissa/weights for embedded Gauss and Kronrod nodes.
- ``SpecialFunctions/Quadrature/Rule/tanhSinh(maxRefinements:tolerance:)`` — adaptive tanh–sinh on finite intervals (default `[-1, 1]`) with optional overrides for refinement depth and absolute tolerance (defaults: 10 refinements, `1e-9` tolerance).
- ``SpecialFunctions/Quadrature/Rule/sinhSinh(maxRefinements:tolerance:)`` — adaptive sinh–sinh, optimized for even integrands on `(-∞, ∞)` with identical parameter semantics.
- ``SpecialFunctions/Quadrature/Rule/expSinh(maxRefinements:tolerance:)`` — adaptive exp–sinh on `[0, ∞)` for semi-infinite integrals; tunable refinements/tolerance.

Fixed rules accept both ``SpecialFunctions/Quadrature/Interval/automatic`` and custom finite intervals. Adaptive rules honor intervals according to the underlying Boost behavior: tanh–sinh evaluates over the supplied `[lower, upper]`, sinh–sinh ignores bounds and always integrates over `(-∞, ∞)`, and exp–sinh allows shifting the lower bound (`[a, ∞)`).

## Integrator Workflow

Use ``SpecialFunctions/Quadrature/Integrator`` to amortize handle creation when evaluating multiple integrals with the same rule.

```swift
import SwiftyBoost

let gl32 = try Quadrature.Integrator<Double>(rule: .gaussLegendre(points: 32))
let result = try gl32.integrate(over: .finite(lower: 0, upper: .pi / 2)) { x in
    sin(x)
}

print(result.value)           // ≈ 1.0 (integral of sin(x) from 0 to π/2)
print(result.estimatedError)  // backend absolute error estimate
print(result.didConverge)     // true for fixed rules
```

``SpecialFunctions/Quadrature/Result`` carries the estimated value, backend error estimate, L₁ norm, iteration count (1 for fixed rules, 0 for adaptive rules), function evaluation tally, and a ``SpecialFunctions/Quadrature/Metadata`` payload describing the instantiated backend (`rule`, `type`, `precision`, `points`, `isAdaptive`, `supportsInfiniteBounds`). The `supportsInfiniteBounds` flag reflects the bridge’s capabilities (currently `true` for sinh–sinh and exp–sinh).

## One-Shot Integration

For single evaluations, call ``SpecialFunctions/Quadrature/integrate(using:over:integrand:)`` which constructs an integrator, evaluates the closure, and disposes the handle automatically.

```swift
let gaussian = try Quadrature.integrate(
    using: .tanhSinh(maxRefinements: 12, tolerance: 1e-12),
    over: .finite(lower: -1, upper: 1)
) { (x: Double) in
    1.0 / (1.0 + x * x)
}
print(gaussian.value)   // ≈ π / 2
let info = gaussian.metadata
```

## Abscissa and Weights

Fixed rules expose their nodes (`xᵢ`) and weights (`wᵢ`). Allocate buffers with capacity equal to ``SpecialFunctions/Quadrature/Metadata/points`` and call ``SpecialFunctions/Quadrature/Integrator/copyAbscissaWeights(into:weights:)``.

```swift
var nodes = Array(repeating: 0.0, count: result.metadata.points)
var weights = nodes
let copied = nodes.withUnsafeMutableBufferPointer { abPtr in
    weights.withUnsafeMutableBufferPointer { wtPtr in
        gl32.copyAbscissaWeights(into: abPtr, weights: wtPtr)
    }
}
```

Adaptive rules return `false` because they do not expose a fixed tableau.

## Error Handling and Intervals

- ``SpecialFunctions/Quadrature/Error/invalidConfiguration(_:)`` flags unsupported point counts or non-positive refinement/tolerance values reported by the bridge.
- ``SpecialFunctions/Quadrature/Error/backendUnavailable(_:)`` appears when a rule is unavailable for the requested precision (for example, Float80 on non-x86 architectures).
- ``SpecialFunctions/Quadrature/Error/invalidInterval(lower:upper:)`` guards against non-finite or reversed bounds when using ``SpecialFunctions/Quadrature/Interval/finite(lower:upper:)``.

When integrating with `Float`, ensure `lower` and `upper` fit in `Float` as well; bounds are downcast before reaching the bridge.

## Precision and Concurrency

- Integrators support `Double`, `Float`, and (on x86_64 or i386) `Float80`. Attempting to construct other scalar types triggers ``SpecialFunctions/Quadrature/Error/invalidConfiguration(_:)``.
- ``SpecialFunctions/Quadrature/Integrator`` is ``Sendable``; each call to ``SpecialFunctions/Quadrature/Integrator/integrate(over:integrand:)`` executes synchronously, so you can create multiple integrators and run them on distinct tasks or threads.
