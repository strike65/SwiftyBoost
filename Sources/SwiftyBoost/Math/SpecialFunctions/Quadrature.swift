//
//  Quadrature.swift
//  Math/SpecialFunctions
//
//  High-level Swift wrappers for the Boost.Math quadrature backends exposed
//  through `CBoostBridge`. These helpers manage opaque C handles, convert Swift
//  closures into C callbacks, and surface integration metadata with support for
//  `Float`, `Double`, and (where available) `Float80`.
//

import SwiftyBoostPrelude

// The C enums are bridged via CBoostBridge and used for metadata reporting only.
// We mark them retroactively Sendable because they are plain C enums.
extension QuadratureType: @retroactive @unchecked Sendable {}
extension QuadraturePrecision: @retroactive @unchecked Sendable {}

public extension SpecialFunctions {
    /// Numerical integration (quadrature) helpers backed by Boost.Math.
    ///
    /// This namespace provides:
    /// - A type-safe description of available quadrature rules.
    /// - A reusable `Integrator` that owns a single C handle and can evaluate multiple integrals.
    /// - A one-shot `integrate` convenience that creates and disposes an integrator for a single call.
    /// - Rich result metadata (estimated error, L1 norm, iterations, function calls, convergence flag).
    ///
    /// Precision
    /// - Generic over `Scalar: Real & BinaryFloatingPoint` with concrete support for:
    ///   - `Double`
    ///   - `Float`
    ///   - `Float80` on x86 (where available)
    ///
    /// Thread-safety
    /// - `Integrator` instances are Sendable; the underlying handle is confined to the
    ///   `Storage` reference and used synchronously during each `integrate` call.
    /// - You can create multiple integrators and use them concurrently from different tasks/threads.
    ///
    /// Intervals
    /// - `.automatic` chooses the natural interval of the rule (e.g. Gauss–Legendre: [-1, 1]).
    /// - `.finite(lower:upper:)` integrates over a user-specified finite interval.
    enum Quadrature {
        /// Describes the quadrature rule to employ.
        ///
        /// Notes per rule:
        /// - `.gaussLegendre(points:)`
        ///   - Fixed, non-adaptive rule on [-1, 1] with positive weights.
        ///   - Use `.finite` with an affine change of variables if you need [a, b].
        /// - `.gaussKronrod(points:)`
        ///   - Fixed Gauss–Kronrod extension (e.g. 15, 21, 31 points depending on availability).
        /// - `.tanhSinh(maxRefinements:tolerance:)`
        ///   - Adaptive double-exponential rule over (-∞, ∞). Natural interval is infinite.
        ///   - Provide parameters to override defaults; otherwise use the backend defaults.
        /// - `.sinhSinh(maxRefinements:tolerance:)`
        ///   - Adaptive double-exponential rule suited for even integrands on (-∞, ∞).
        /// - `.expSinh(maxRefinements:tolerance:)`
        ///   - Adaptive double-exponential rule suited for semi-infinite domains [0, ∞).
        ///
        /// Point count constraints:
        /// - Must be > 0 and fit in Int32; otherwise an error is thrown at construction time.
        public enum Rule: Sendable, Equatable {
            /// Gauss–Legendre rule with a given number of points on [-1, 1].
            case gaussLegendre(points: Int)
            /// Gauss–Kronrod fixed rule on [-1, 1].
            case gaussKronrod(points: Int)
            /// Tanh–Sinh (double-exponential) adaptive rule over (-∞, ∞).
            /// - Parameters override backend defaults if provided.
            case tanhSinh(maxRefinements: Int? = nil, tolerance: Double? = nil)
            /// Sinh–Sinh adaptive rule (often used for even integrands) over (-∞, ∞).
            case sinhSinh(maxRefinements: Int? = nil, tolerance: Double? = nil)
            /// Exp–Sinh adaptive rule over [0, ∞).
            case expSinh(maxRefinements: Int? = nil, tolerance: Double? = nil)
        }

        /// Specifies how to interpret integration bounds.
        ///
        /// - `.automatic`: Use the rule’s natural interval:
        ///   - Gauss–Legendre, Gauss–Jacobi, Gauss–Kronrod: [-1, 1]
        ///   - Gauss–Hermite: (-∞, ∞)
        ///   - Gauss–Laguerre: [0, ∞)
        ///   - Tanh–Sinh, Sinh–Sinh: (-∞, ∞)
        ///   - Exp–Sinh: [0, ∞)
        /// - `.finite(lower:upper:)`: Integrate over a finite, ordered interval [lower, upper].
        ///   - Both bounds must be finite and `lower < upper`. For `Float` integrators the bounds
        ///     must also be representable as finite `Float` values.
        public enum Interval: Sendable, Equatable {
            /// Use the rule’s natural interval (e.g., [-1, 1] for Gauss–Legendre).
            case automatic
            /// Integrate over `[lower, upper]`; bounds must be finite and ordered.
            case finite(lower: Double, upper: Double)
        }

        /// Describes the instantiated quadrature backend and its capabilities.
        ///
        /// Values are queried from the C bridge after creating the handle:
        /// - `type`: The backend rule family (e.g., Gauss–Legendre, Tanh–Sinh).
        /// - `precision`: Underlying scalar precision of the handle.
        /// - `points`: Number of nodes for fixed rules (0 for adaptive rules).
        /// - `isAdaptive`: Whether the backend is adaptive (double-exponential families).
        /// - `supportsInfiniteBounds`: Whether the backend can integrate on infinite domains.
        public struct Metadata: Sendable, Equatable {
            /// The Swift rule used to construct the handle.
            public let rule: Rule
            /// C-side rule family identifier.
            public let type: QuadratureType
            /// C-side precision of the handle (float/double/long double).
            public let precision: QuadraturePrecision
            /// Number of nodes for fixed rules; 0 for adaptive rules.
            public let points: Int
            /// True for adaptive (double-exponential) backends.
            public let isAdaptive: Bool
            /// True if the backend accepts infinite bounds on the C side.
            public let supportsInfiniteBounds: Bool
        }

        /// Encapsulates the result of an integral evaluation.
        ///
        /// Fields mirror the bridge’s `QuadratureResult*` structures with Swift-friendly types.
        /// - value: The estimated integral.
        /// - estimatedError: Backend’s error estimate (absolute).
        /// - l1Norm: Estimated L1 norm of the integrand (if provided by the backend).
        /// - iterations: Iteration count reported by the backend (0 for fixed rules).
        /// - functionEvaluations: Number of integrand calls performed.
        /// - converged: Backend-reported convergence flag (true/false).
        /// - metadata: The `Metadata` describing the instantiated backend.
        public struct Result<Scalar: Real & BinaryFloatingPoint & Sendable>: Sendable, Equatable {
            /// The estimated integral value.
            public let value: Scalar
            /// The backend’s absolute error estimate.
            public let estimatedError: Scalar
            /// Estimated L1 norm of the integrand over the integration domain.
            public let l1Norm: Scalar
            /// Iteration count (adaptive rules) or 0 for fixed rules.
            public let iterations: Int
            /// Number of integrand evaluations performed.
            public let functionEvaluations: Int
            /// Whether the backend reported convergence.
            public let converged: Bool
            /// Metadata describing the backend used for this evaluation.
            public let metadata: Metadata

            /// Convenience alias for `converged`.
            @inlinable public var didConverge: Bool { converged }

            /// Convert the result value to another binary floating-point type.
            @inlinable public func value<T: BinaryFloatingPoint>(as _: T.Type = T.self) -> T {
                T(value)
            }

            fileprivate init(raw: QuadratureBridge.RawResult<Scalar>, metadata: Metadata) {
                self.value = raw.value
                self.estimatedError = raw.estimatedError
                self.l1Norm = raw.l1Norm
                self.iterations = raw.iterations
                self.functionEvaluations = raw.functionEvaluations
                self.converged = raw.converged
                self.metadata = metadata
            }
        }

        /// Errors thrown by quadrature setup or evaluation.
        public enum Error: Swift.Error, Equatable {
            /// Rule-specific validation failed (e.g., non-positive points, α/β ≤ -1).
            case invalidConfiguration(String)
            /// The requested interval was invalid (non-finite or unordered).
            case invalidInterval(lower: Double, upper: Double)
            /// The backend for the requested rule/precision combination is not available.
            case backendUnavailable(String)
        }

        /// Reusable integrator that owns a single quadrature handle.
        ///
        /// Use this type if you plan to evaluate multiple integrals with the same rule
        /// (potentially over different intervals or with different integrands). The
        /// underlying C handle is created once in `init(rule:)` and destroyed when the
        /// integrator is deinitialized.
        ///
        /// Example
        /// - Create a 61-point Gauss–Kronrod integrator and integrate on [0, 1]:
        ///   let gk = try SpecialFunctions.Quadrature.Integrator<Double>(rule: .gaussKronrod(points: 61))
        ///   let result = try gk.integrate(over: .finite(lower: 0, upper: 1)) { x in 1.0 / (1.0 + x*x) }
        ///   print(result.value, result.estimatedError, result.didConverge)
        public struct Integrator<Scalar: Real & BinaryFloatingPoint & Sendable>: Sendable {
            /// Internal reference type that owns the opaque C handle.
            private final class Storage: @unchecked Sendable {
                let handle: QuadratureHandle
                let metadata: Metadata

                /// Create and configure the C handle for the given rule.
                /// - Throws: `Error.invalidConfiguration` or `Error.backendUnavailable`.
                init(rule: Rule) throws {
                    let handle = try QuadratureBridge.createHandle(for: rule, scalar: Scalar.self)
                    let type = quad_get_type(handle)
                    let precision = quad_get_precision(handle)
                    self.handle = handle
                    self.metadata = Metadata(
                        rule: rule,
                        type: type,
                        precision: precision,
                        points: Int(quad_get_points(handle)),
                        isAdaptive: quad_is_adaptive(type) != 0,
                        supportsInfiniteBounds: quad_supports_infinite_bounds(type) != 0
                    )
                }

                deinit {
                    quad_destroy(handle)
                }

                /// Evaluate the integral for the supplied integrand over the given interval.
                /// - Throws: `Error.invalidInterval` for non-finite or unordered bounds.
                func integrate(over interval: Interval,
                               integrand: @escaping @Sendable (Scalar) -> Scalar) throws -> Result<Scalar> {
                    let raw = try QuadratureBridge.integrate(handle: handle,
                                                             over: interval,
                                                             integrand: integrand,
                                                             scalar: Scalar.self)
                    return Result(raw: raw, metadata: metadata)
                }

                /// Copy abscissa and weights for fixed rules into caller-provided buffers.
                /// - Returns: `true` if the backend provided nodes/weights and buffers were large enough.
                /// - Note: Adaptive rules do not expose fixed abscissa/weights and will return `false`.
                func copyAbscissa(into abscissa: UnsafeMutablePointer<Scalar>,
                                  weights: UnsafeMutablePointer<Scalar>,
                                  capacity: Int) -> Bool {
                    QuadratureBridge.copyAbscissa(handle: handle,
                                                  abscissa: abscissa,
                                                  weights: weights,
                                                  capacity: capacity,
                                                  scalar: Scalar.self)
                }
            }

            private let storage: Storage

            /// Metadata describing the underlying backend.
            public var metadata: Metadata { storage.metadata }

            /// Create an integrator for the specified rule.
            /// - Throws:
            ///   - `Error.invalidConfiguration` for invalid parameters (e.g., non-positive points).
            ///   - `Error.backendUnavailable` if the backend is not present for the target precision.
            public init(rule: Rule) throws {
                self.storage = try Storage(rule: rule)
            }

            /// Evaluate the integral for the supplied integrand over the given interval.
            ///
            /// - Parameters:
            ///   - interval: `.automatic` uses the rule’s natural domain; `.finite` integrates over [lower, upper].
            ///   - integrand: A `Sendable` closure from `Scalar` to `Scalar`.
            /// - Returns: A `Result` containing the integral estimate and metadata.
            /// - Throws: `Error.invalidInterval` if bounds are not finite or ordered.
            public func integrate(over interval: Interval = .automatic,
                                  integrand: @escaping @Sendable (Scalar) -> Scalar) throws -> Result<Scalar> {
                try storage.integrate(over: interval, integrand: integrand)
            }

            /// Copy abscissa and weights into caller-provided buffers. Returns false if unavailable.
            ///
            /// - Parameters:
            ///   - abscissa: Mutable buffer to receive nodes (x_i).
            ///   - weights: Mutable buffer to receive weights (w_i).
            /// - Returns:
            ///   - `true` if the backend supplied nodes and weights and both buffers had sufficient capacity.
            ///   - `false` if the backend is adaptive (no fixed nodes), the capacity was insufficient,
            ///     or the underlying C call rejected the request.
            ///
            /// Notes
            /// - For fixed rules (e.g., Gauss–Legendre), the required capacity equals `metadata.points`.
            /// - For adaptive rules, this always returns `false`.
            public func copyAbscissaWeights(into abscissa: UnsafeMutableBufferPointer<Scalar>,
                                            weights: UnsafeMutableBufferPointer<Scalar>) -> Bool {
                guard abscissa.count > 0,
                      abscissa.count == weights.count,
                      let aPtr = abscissa.baseAddress,
                      let wPtr = weights.baseAddress else {
                    return false
                }
                return storage.copyAbscissa(into: aPtr, weights: wPtr, capacity: abscissa.count)
            }
        }

        /// Convenience one-shot entry point that constructs an integrator,
        /// evaluates the integral, and disposes the handle afterwards.
        ///
        /// Prefer this overload for simple, single-use integrations. For repeated
        /// evaluations, construct an `Integrator` to amortize handle creation.
        ///
        /// - Parameters:
        ///   - rule: The quadrature rule to use.
        ///   - interval: Integration interval; defaults to the rule’s natural interval.
        ///   - integrand: A `Sendable` closure mapping `Scalar` to `Scalar`.
        /// - Returns: A `Result` with the integral estimate and metadata.
        /// - Throws: `Error.invalidConfiguration`, `Error.backendUnavailable`, or `Error.invalidInterval`.
        public static func integrate<Scalar: Real & BinaryFloatingPoint & Sendable>(
            using rule: Rule,
            over interval: Interval = .automatic,
            integrand: @escaping @Sendable (Scalar) -> Scalar
        ) throws -> Result<Scalar> {
            try Integrator<Scalar>(rule: rule).integrate(over: interval, integrand: integrand)
        }
    }
}

// MARK: - Rule helpers

private extension SpecialFunctions.Quadrature.Rule {
    /// Validate and downcast a positive point count to Int32 for the C bridge.
    static func checkedCount(_ value: Int, label: String) throws -> Int32 {
        guard value > 0 else {
            throw SpecialFunctions.Quadrature.Error.invalidConfiguration("\(label) requires a positive point count (received \(value)).")
        }
        guard value <= Int(Int32.max) else {
            throw SpecialFunctions.Quadrature.Error.invalidConfiguration("\(label) point count exceeds Int32.max (received \(value)).")
        }
        return Int32(value)
    }

    /// Normalize adaptive parameters (max refinements and tolerance) and validate them.
    ///
    /// - maxRef: Defaults to 10 if nil; must be positive and fit in Int32.
    /// - tolerance: Defaults to 1e-9 if nil; must be positive.
    static func normalizedAdaptiveParameters<Scalar: BinaryFloatingPoint>(
        _ maxRef: Int?,
        _ tolerance: Scalar?,
        label: String
    ) throws -> (Int32, Scalar) {
        let maxRefinements = maxRef ?? 10
        guard maxRefinements > 0 else {
            throw SpecialFunctions.Quadrature.Error.invalidConfiguration("\(label) maxRefinements must be positive (received \(maxRef ?? 0)).")
        }
        guard maxRefinements <= Int(Int32.max) else {
            throw SpecialFunctions.Quadrature.Error.invalidConfiguration("\(label) maxRefinements exceeds Int32.max (received \(maxRefinements)).")
        }
        let tol = tolerance ?? Scalar(1e-9)
        guard tol > 0 else {
            throw SpecialFunctions.Quadrature.Error.invalidConfiguration("\(label) tolerance must be positive (received \(Double(tolerance ?? 0))).")
        }
        return (Int32(maxRefinements), tol)
    }
}

// MARK: - Bridging helpers
//
// Internal adapter that dispatches to precision-specific bridge implementations,
// converts Swift closures to C callbacks, and maps raw results into generic form.

private enum QuadratureBridge {
    /// Precision-agnostic result used to populate `Quadrature.Result`.
    struct RawResult<Scalar: Real & BinaryFloatingPoint & Sendable> {
        let value: Scalar
        let estimatedError: Scalar
        let l1Norm: Scalar
        let iterations: Int
        let functionEvaluations: Int
        let converged: Bool
    }

    /// Create a C-side handle for the given rule and scalar precision.
    static func createHandle<Scalar: Real & BinaryFloatingPoint & Sendable>(
        for rule: SpecialFunctions.Quadrature.Rule,
        scalar _: Scalar.Type
    ) throws -> QuadratureHandle {
        if Scalar.self == Double.self {
            return try QuadratureBridgeDouble.createHandle(for: rule)
        } else if Scalar.self == Float.self {
            return try QuadratureBridgeFloat.createHandle(for: rule)
        } else {
            #if arch(x86_64) || arch(i386)
            if Scalar.self == Float80.self {
                return try QuadratureBridgeFloat80.createHandle(for: rule)
            }
            #endif
            throw SpecialFunctions.Quadrature.Error.invalidConfiguration("Unsupported scalar type: \(Scalar.self)")
        }
    }

    /// Evaluate the integral by forwarding to the precision-specific bridge.
    static func integrate<Scalar: Real & BinaryFloatingPoint & Sendable>(
        handle: QuadratureHandle,
        over interval: SpecialFunctions.Quadrature.Interval,
        integrand: @escaping @Sendable (Scalar) -> Scalar,
        scalar _: Scalar.Type
    ) throws -> RawResult<Scalar> {
        if Scalar.self == Double.self {
            let typed = unsafeBitCast(integrand, to: (@Sendable (Double) -> Double).self)
            let raw = try QuadratureBridgeDouble.integrate(handle: handle, over: interval, integrand: typed)
            return QuadratureBridge.RawResult(
                value: Scalar(raw.value),
                estimatedError: Scalar(raw.estimatedError),
                l1Norm: Scalar(raw.l1Norm),
                iterations: raw.iterations,
                functionEvaluations: raw.functionEvaluations,
                converged: raw.converged
            )
        } else if Scalar.self == Float.self {
            let typed = unsafeBitCast(integrand, to: (@Sendable (Float) -> Float).self)
            let raw = try QuadratureBridgeFloat.integrate(handle: handle, over: interval, integrand: typed)
            return QuadratureBridge.RawResult(
                value: Scalar(raw.value),
                estimatedError: Scalar(raw.estimatedError),
                l1Norm: Scalar(raw.l1Norm),
                iterations: raw.iterations,
                functionEvaluations: raw.functionEvaluations,
                converged: raw.converged
            )
        } else {
            #if arch(x86_64) || arch(i386)
            if Scalar.self == Float80.self {
                let typed = unsafeBitCast(integrand, to: (@Sendable (Float80) -> Float80).self)
                let raw = try QuadratureBridgeFloat80.integrate(handle: handle, over: interval, integrand: typed)
                return QuadratureBridge.RawResult(
                    value: Scalar(raw.value),
                    estimatedError: Scalar(raw.estimatedError),
                    l1Norm: Scalar(raw.l1Norm),
                    iterations: raw.iterations,
                    functionEvaluations: raw.functionEvaluations,
                    converged: raw.converged
                )
            }
            #endif
            throw SpecialFunctions.Quadrature.Error.invalidConfiguration("Unsupported scalar type: \(Scalar.self)")
        }
    }

    /// Copy fixed abscissa/weights for the given precision (if available).
    static func copyAbscissa<Scalar: Real & BinaryFloatingPoint & Sendable>(
        handle: QuadratureHandle,
        abscissa: UnsafeMutablePointer<Scalar>,
        weights: UnsafeMutablePointer<Scalar>,
        capacity: Int,
        scalar _: Scalar.Type
    ) -> Bool {
        if Scalar.self == Double.self {
            return abscissa.withMemoryRebound(to: Double.self, capacity: capacity) { abPtr in
                weights.withMemoryRebound(to: Double.self, capacity: capacity) { wtPtr in
                    QuadratureBridgeDouble.copyAbscissa(handle: handle,
                                                        abscissa: abPtr,
                                                        weights: wtPtr,
                                                        capacity: capacity)
                }
            }
        } else if Scalar.self == Float.self {
            return abscissa.withMemoryRebound(to: Float.self, capacity: capacity) { abPtr in
                weights.withMemoryRebound(to: Float.self, capacity: capacity) { wtPtr in
                    QuadratureBridgeFloat.copyAbscissa(handle: handle,
                                                       abscissa: abPtr,
                                                       weights: wtPtr,
                                                       capacity: capacity)
                }
            }
        } else {
            #if arch(x86_64) || arch(i386)
            if Scalar.self == Float80.self {
                return abscissa.withMemoryRebound(to: Float80.self, capacity: capacity) { abPtr in
                    weights.withMemoryRebound(to: Float80.self, capacity: capacity) { wtPtr in
                        QuadratureBridgeFloat80.copyAbscissa(handle: handle,
                                                             abscissa: abPtr,
                                                             weights: wtPtr,
                                                             capacity: capacity)
                    }
                }
            }
            #endif
            return false
        }
    }
}

// MARK: - Double precision bridge
//
// Precision-specific construction, evaluation, and abscissa copying for Double.

private enum QuadratureBridgeDouble {
    static func createHandle(for rule: SpecialFunctions.Quadrature.Rule) throws -> QuadratureHandle {
        switch rule {
        case let .gaussLegendre(points):
            let n = try SpecialFunctions.Quadrature.Rule.checkedCount(points, label: "Gauss-Legendre")
            guard let handle = quad_gauss_create_d(n) else {
                throw SpecialFunctions.Quadrature.Error.invalidConfiguration("Gauss-Legendre does not support \(points) points.")
            }
            return handle
        case let .gaussKronrod(points):
            let n = try SpecialFunctions.Quadrature.Rule.checkedCount(points, label: "Gauss-Kronrod")
            guard let handle = quad_gauss_kronrod_create_d(n) else {
                throw SpecialFunctions.Quadrature.Error.invalidConfiguration("Gauss-Kronrod does not support \(points) points.")
            }
            return handle
        case let .tanhSinh(maxRefinements, tolerance):
            if maxRefinements == nil && tolerance == nil {
                guard let handle = quad_tanh_sinh_create_d() else {
                    throw SpecialFunctions.Quadrature.Error.backendUnavailable("Tanh-Sinh integrator is unavailable.")
                }
                return handle
            }
            let (maxRef, tol) = try SpecialFunctions.Quadrature.Rule
                .normalizedAdaptiveParameters(maxRefinements, tolerance, label: "Tanh-Sinh")
            guard let handle = quad_tanh_sinh_create_with_params_d(maxRef, tol) else {
                throw SpecialFunctions.Quadrature.Error.invalidConfiguration("Tanh-Sinh integrator rejected maxRef=\(maxRefinements ?? 0), tolerance=\(tolerance ?? 0).")
            }
            return handle
        case let .sinhSinh(maxRefinements, tolerance):
            if maxRefinements == nil && tolerance == nil {
                guard let handle = quad_sinh_sinh_create_d() else {
                    throw SpecialFunctions.Quadrature.Error.backendUnavailable("Sinh-Sinh integrator is unavailable.")
                }
                return handle
            }
            let (maxRef, tol) = try SpecialFunctions.Quadrature.Rule
                .normalizedAdaptiveParameters(maxRefinements, tolerance, label: "Sinh-Sinh")
            guard let handle = quad_sinh_sinh_create_with_params_d(maxRef, tol) else {
                throw SpecialFunctions.Quadrature.Error.invalidConfiguration("Sinh-Sinh integrator rejected maxRef=\(maxRefinements ?? 0), tolerance=\(tolerance ?? 0).")
            }
            return handle
        case let .expSinh(maxRefinements, tolerance):
            if maxRefinements == nil && tolerance == nil {
                guard let handle = quad_exp_sinh_create_d() else {
                    throw SpecialFunctions.Quadrature.Error.backendUnavailable("Exp-Sinh integrator is unavailable.")
                }
                return handle
            }
            let (maxRef, tol) = try SpecialFunctions.Quadrature.Rule
                .normalizedAdaptiveParameters(maxRefinements, tolerance, label: "Exp-Sinh")
            guard let handle = quad_exp_sinh_create_with_params_d(maxRef, tol) else {
                throw SpecialFunctions.Quadrature.Error.invalidConfiguration("Exp-Sinh integrator rejected maxRef=\(maxRefinements ?? 0), tolerance=\(tolerance ?? 0).")
            }
            return handle
        }
    }

    static func integrate(
        handle: QuadratureHandle,
        over interval: SpecialFunctions.Quadrature.Interval,
        integrand: @escaping @Sendable (Double) -> Double
    ) throws -> QuadratureBridge.RawResult<Double> {
        let raw: QuadratureResultD = try withCallbackDouble(integrand) { callback, context in
            switch interval {
            case .automatic:
                return quad_integrate_d(handle, callback, context)
            case let .finite(lower, upper):
                guard lower.isFinite, upper.isFinite, lower < upper else {
                    throw SpecialFunctions.Quadrature.Error.invalidInterval(lower: lower, upper: upper)
                }
                return quad_integrate_interval_d(handle, callback, context, lower, upper)
            }
        }
        return QuadratureBridge.RawResult(
            value: raw.result,
            estimatedError: raw.error,
            l1Norm: raw.l1_norm,
            iterations: Int(raw.iterations),
            functionEvaluations: Int(raw.function_calls),
            converged: raw.converged != 0
        )
    }

    static func copyAbscissa(
        handle: QuadratureHandle,
        abscissa: UnsafeMutablePointer<Double>,
        weights: UnsafeMutablePointer<Double>,
        capacity: Int
    ) -> Bool {
        guard capacity <= Int(Int32.max) else { return false }
        return quad_get_abscissa_weights_d(handle, abscissa, weights, Int32(capacity)) != 0
    }
}

// MARK: - Float precision bridge
//
// Precision-specific construction, evaluation, and abscissa copying for Float.

private enum QuadratureBridgeFloat {
    static func createHandle(for rule: SpecialFunctions.Quadrature.Rule) throws -> QuadratureHandle {
        switch rule {
        case let .gaussLegendre(points):
            let n = try SpecialFunctions.Quadrature.Rule.checkedCount(points, label: "Gauss-Legendre")
            guard let handle = quad_gauss_create_f(n) else {
                throw SpecialFunctions.Quadrature.Error.invalidConfiguration("Gauss-Legendre does not support \(points) points.")
            }
            return handle
        case let .gaussKronrod(points):
            let n = try SpecialFunctions.Quadrature.Rule.checkedCount(points, label: "Gauss-Kronrod")
            guard let handle = quad_gauss_kronrod_create_f(n) else {
                throw SpecialFunctions.Quadrature.Error.invalidConfiguration("Gauss-Kronrod does not support \(points) points.")
            }
            return handle
        case let .tanhSinh(maxRefinements, tolerance):
            if maxRefinements == nil && tolerance == nil {
                guard let handle = quad_tanh_sinh_create_f() else {
                    throw SpecialFunctions.Quadrature.Error.backendUnavailable("Tanh-Sinh integrator is unavailable.")
                }
                return handle
            }
            let (maxRef, tol) = try SpecialFunctions.Quadrature.Rule
                .normalizedAdaptiveParameters(maxRefinements, tolerance.map(Float.init), label: "Tanh-Sinh")
            guard let handle = quad_tanh_sinh_create_with_params_f(maxRef, tol) else {
                throw SpecialFunctions.Quadrature.Error.invalidConfiguration("Tanh-Sinh integrator rejected maxRef=\(maxRefinements ?? 0), tolerance=\(tolerance ?? 0).")
            }
            return handle
        case let .sinhSinh(maxRefinements, tolerance):
            if maxRefinements == nil && tolerance == nil {
                guard let handle = quad_sinh_sinh_create_f() else {
                    throw SpecialFunctions.Quadrature.Error.backendUnavailable("Sinh-Sinh integrator is unavailable.")
                }
                return handle
            }
            let (maxRef, tol) = try SpecialFunctions.Quadrature.Rule
                .normalizedAdaptiveParameters(maxRefinements, tolerance.map(Float.init), label: "Sinh-Sinh")
            guard let handle = quad_sinh_sinh_create_with_params_f(maxRef, tol) else {
                throw SpecialFunctions.Quadrature.Error.invalidConfiguration("Sinh-Sinh integrator rejected maxRef=\(maxRefinements ?? 0), tolerance=\(tolerance ?? 0).")
            }
            return handle
        case let .expSinh(maxRefinements, tolerance):
            if maxRefinements == nil && tolerance == nil {
                guard let handle = quad_exp_sinh_create_f() else {
                    throw SpecialFunctions.Quadrature.Error.backendUnavailable("Exp-Sinh integrator is unavailable.")
                }
                return handle
            }
            let (maxRef, tol) = try SpecialFunctions.Quadrature.Rule
                .normalizedAdaptiveParameters(maxRefinements, tolerance.map(Float.init), label: "Exp-Sinh")
            guard let handle = quad_exp_sinh_create_with_params_f(maxRef, tol) else {
                throw SpecialFunctions.Quadrature.Error.invalidConfiguration("Exp-Sinh integrator rejected maxRef=\(maxRefinements ?? 0), tolerance=\(tolerance ?? 0).")
            }
            return handle
        }
    }

    static func integrate(
        handle: QuadratureHandle,
        over interval: SpecialFunctions.Quadrature.Interval,
        integrand: @escaping @Sendable (Float) -> Float
    ) throws -> QuadratureBridge.RawResult<Float> {
        let raw: QuadratureResultF = try withCallbackFloat(integrand) { callback, context in
            switch interval {
            case .automatic:
                return quad_integrate_f(handle, callback, context)
            case let .finite(lower, upper):
                guard lower.isFinite, upper.isFinite, lower < upper else {
                    throw SpecialFunctions.Quadrature.Error.invalidInterval(lower: lower, upper: upper)
                }
                let lowerScalar = Float(lower)
                let upperScalar = Float(upper)
                guard lowerScalar.isFinite, upperScalar.isFinite else {
                    throw SpecialFunctions.Quadrature.Error.invalidInterval(lower: lower, upper: upper)
                }
                return quad_integrate_interval_f(handle, callback, context, lowerScalar, upperScalar)
            }
        }
        return QuadratureBridge.RawResult(
            value: raw.result,
            estimatedError: raw.error,
            l1Norm: raw.l1_norm,
            iterations: Int(raw.iterations),
            functionEvaluations: Int(raw.function_calls),
            converged: raw.converged != 0
        )
    }

    static func copyAbscissa(
        handle: QuadratureHandle,
        abscissa: UnsafeMutablePointer<Float>,
        weights: UnsafeMutablePointer<Float>,
        capacity: Int
    ) -> Bool {
        guard capacity <= Int(Int32.max) else { return false }
        return quad_get_abscissa_weights_f(handle, abscissa, weights, Int32(capacity)) != 0
    }
}

// MARK: - Float80 bridge (x86 only)
//
// Precision-specific construction, evaluation, and abscissa copying for Float80.

#if arch(x86_64) || arch(i386)
private enum QuadratureBridgeFloat80 {
    static func createHandle(for rule: SpecialFunctions.Quadrature.Rule) throws -> QuadratureHandle {
        switch rule {
        case let .gaussLegendre(points):
            let n = try SpecialFunctions.Quadrature.Rule.checkedCount(points, label: "Gauss-Legendre")
            guard let handle = quad_gauss_create_l(n) else {
                throw SpecialFunctions.Quadrature.Error.invalidConfiguration("Gauss-Legendre does not support \(points) points.")
            }
            return handle
        case let .gaussKronrod(points):
            let n = try SpecialFunctions.Quadrature.Rule.checkedCount(points, label: "Gauss-Kronrod")
            guard let handle = quad_gauss_kronrod_create_l(n) else {
                throw SpecialFunctions.Quadrature.Error.invalidConfiguration("Gauss-Kronrod does not support \(points) points.")
            }
            return handle
        case let .tanhSinh(maxRefinements, tolerance):
            if maxRefinements == nil && tolerance == nil {
                guard let handle = quad_tanh_sinh_create_l() else {
                    throw SpecialFunctions.Quadrature.Error.backendUnavailable("Tanh-Sinh integrator is unavailable.")
                }
                return handle
            }
            let (maxRef, tol) = try SpecialFunctions.Quadrature.Rule
                .normalizedAdaptiveParameters(maxRefinements, tolerance.map(Float80.init), label: "Tanh-Sinh")
            guard let handle = quad_tanh_sinh_create_with_params_l(maxRef, tol) else {
                throw SpecialFunctions.Quadrature.Error.invalidConfiguration("Tanh-Sinh integrator rejected maxRef=\(maxRefinements ?? 0), tolerance=\(tolerance ?? 0).")
            }
            return handle
        case let .sinhSinh(maxRefinements, tolerance):
            if maxRefinements == nil && tolerance == nil {
                guard let handle = quad_sinh_sinh_create_l() else {
                    throw SpecialFunctions.Quadrature.Error.backendUnavailable("Sinh-Sinh integrator is unavailable.")
                }
                return handle
            }
            let (maxRef, tol) = try SpecialFunctions.Quadrature.Rule
                .normalizedAdaptiveParameters(maxRefinements, tolerance.map(Float80.init), label: "Sinh-Sinh")
            guard let handle = quad_sinh_sinh_create_with_params_l(maxRef, tol) else {
                throw SpecialFunctions.Quadrature.Error.invalidConfiguration("Sinh-Sinh integrator rejected maxRef=\(maxRefinements ?? 0), tolerance=\(tolerance ?? 0).")
            }
            return handle
        case let .expSinh(maxRefinements, tolerance):
            if maxRefinements == nil && tolerance == nil {
                guard let handle = quad_exp_sinh_create_l() else {
                    throw SpecialFunctions.Quadrature.Error.backendUnavailable("Exp-Sinh integrator is unavailable.")
                }
                return handle
            }
            let (maxRef, tol) = try SpecialFunctions.Quadrature.Rule
                .normalizedAdaptiveParameters(maxRefinements, tolerance.map(Float80.init), label: "Exp-Sinh")
            guard let handle = quad_exp_sinh_create_with_params_l(maxRef, tol) else {
                throw SpecialFunctions.Quadrature.Error.invalidConfiguration("Exp-Sinh integrator rejected maxRef=\(maxRefinements ?? 0), tolerance=\(tolerance ?? 0).")
            }
            return handle
        }
    }

    static func integrate(
        handle: QuadratureHandle,
        over interval: SpecialFunctions.Quadrature.Interval,
        integrand: @escaping @Sendable (Float80) -> Float80
    ) throws -> QuadratureBridge.RawResult<Float80> {
        let raw: QuadratureResultL = try withCallbackFloat80(integrand) { callback, context in
            switch interval {
            case .automatic:
                return quad_integrate_l(handle, callback, context)
            case let .finite(lower, upper):
                guard lower.isFinite, upper.isFinite, lower < upper else {
                    throw SpecialFunctions.Quadrature.Error.invalidInterval(lower: lower, upper: upper)
                }
                return quad_integrate_interval_l(handle, callback, context, Float80(lower), Float80(upper))
            }
        }
        return QuadratureBridge.RawResult(
            value: raw.result,
            estimatedError: raw.error,
            l1Norm: raw.l1_norm,
            iterations: Int(raw.iterations),
            functionEvaluations: Int(raw.function_calls),
            converged: raw.converged != 0
        )
    }

    static func copyAbscissa(
        handle: QuadratureHandle,
        abscissa: UnsafeMutablePointer<Float80>,
        weights: UnsafeMutablePointer<Float80>,
        capacity: Int
    ) -> Bool {
        guard capacity <= Int(Int32.max) else { return false }
        return quad_get_abscissa_weights_l(handle, abscissa, weights, Int32(capacity)) != 0
    }
}
#endif

// MARK: - Callback plumbing
//
// Convert `@Sendable (Scalar) -> Scalar` closures into C function pointers with an
// opaque context pointer. We retain the callback box for the duration of the C call
// via `withExtendedLifetime` to ensure the closure remains valid.

private typealias CIntegrandDouble = @convention(c) (Double, UnsafeMutableRawPointer?) -> Double
private typealias CIntegrandFloat = @convention(c) (Float, UnsafeMutableRawPointer?) -> Float
#if arch(x86_64) || arch(i386)
private typealias CIntegrandFloat80 = @convention(c) (Float80, UnsafeMutableRawPointer?) -> Float80
#endif

private final class QuadratureCallbackBoxDouble: @unchecked Sendable {
    let closure: @Sendable (Double) -> Double
    init(_ closure: @escaping @Sendable (Double) -> Double) {
        self.closure = closure
    }
}

private final class QuadratureCallbackBoxFloat: @unchecked Sendable {
    let closure: @Sendable (Float) -> Float
    init(_ closure: @escaping @Sendable (Float) -> Float) {
        self.closure = closure
    }
}

#if arch(x86_64) || arch(i386)
private final class QuadratureCallbackBoxFloat80: @unchecked Sendable {
    let closure: @Sendable (Float80) -> Float80
    init(_ closure: @escaping @Sendable (Float80) -> Float80) {
        self.closure = closure
    }
}
#endif

private let quadratureCallbackDouble: CIntegrandDouble = { x, context in
    guard let context else { return .nan }
    let box = Unmanaged<QuadratureCallbackBoxDouble>.fromOpaque(context).takeUnretainedValue()
    return box.closure(x)
}

private let quadratureCallbackFloat: CIntegrandFloat = { x, context in
    guard let context else { return .nan }
    let box = Unmanaged<QuadratureCallbackBoxFloat>.fromOpaque(context).takeUnretainedValue()
    return box.closure(x)
}

#if arch(x86_64) || arch(i386)
private let quadratureCallbackFloat80: CIntegrandFloat80 = { x, context in
    guard let context else { return .nan }
    let box = Unmanaged<QuadratureCallbackBoxFloat80>.fromOpaque(context).takeUnretainedValue()
    return box.closure(x)
}
#endif

private func withCallbackDouble<R>(
    _ integrand: @escaping @Sendable (Double) -> Double,
    _ body: (CIntegrandDouble, UnsafeMutableRawPointer?) throws -> R
) rethrows -> R {
    let box = QuadratureCallbackBoxDouble(integrand)
    let opaque = Unmanaged.passUnretained(box).toOpaque()
    return try withExtendedLifetime(box) {
        try body(quadratureCallbackDouble, opaque)
    }
}

private func withCallbackFloat<R>(
    _ integrand: @escaping @Sendable (Float) -> Float,
    _ body: (CIntegrandFloat, UnsafeMutableRawPointer?) throws -> R
) rethrows -> R {
    let box = QuadratureCallbackBoxFloat(integrand)
    let opaque = Unmanaged.passUnretained(box).toOpaque()
    return try withExtendedLifetime(box) {
        try body(quadratureCallbackFloat, opaque)
    }
}

#if arch(x86_64) || arch(i386)
private func withCallbackFloat80<R>(
    _ integrand: @escaping @Sendable (Float80) -> Float80,
    _ body: (CIntegrandFloat80, UnsafeMutableRawPointer?) throws -> R
) rethrows -> R {
    let box = QuadratureCallbackBoxFloat80(integrand)
    let opaque = Unmanaged.passUnretained(box).toOpaque()
    return try withExtendedLifetime(box) {
        try body(quadratureCallbackFloat80, opaque)
    }
}
#endif
