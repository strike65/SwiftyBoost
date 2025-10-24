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

extension QuadratureType: @retroactive @unchecked Sendable {}
extension QuadraturePrecision: @retroactive @unchecked Sendable {}

public extension SpecialFunctions {
    enum Quadrature {
        /// Describes the quadrature rule to employ.
        public enum Rule: Sendable, Equatable {
            case gaussLegendre(points: Int)
            case gaussHermite(points: Int)
            case gaussLaguerre(points: Int, alpha: Double? = nil)
            case gaussJacobi(points: Int, alpha: Double, beta: Double)
            case gaussKronrod(points: Int)
            case tanhSinh(maxRefinements: Int? = nil, tolerance: Double? = nil)
            case sinhSinh(maxRefinements: Int? = nil, tolerance: Double? = nil)
            case expSinh(maxRefinements: Int? = nil, tolerance: Double? = nil)
        }

        /// Specifies how to interpret integration bounds.
        public enum Interval: Sendable, Equatable {
            /// Use the rule’s natural interval (e.g., [-1, 1] for Gauss–Legendre).
            case automatic
            /// Integrate over `[lower, upper]`; bounds must be finite and ordered.
            case finite(lower: Double, upper: Double)
        }

        /// Describes the instantiated quadrature backend.
        public struct Metadata: Sendable, Equatable {
            public let rule: Rule
            public let type: QuadratureType
            public let precision: QuadraturePrecision
            public let points: Int
            public let isAdaptive: Bool
            public let supportsInfiniteBounds: Bool
        }

        /// Encapsulates the result of an integral evaluation.
        public struct Result<Scalar: Real & BinaryFloatingPoint & Sendable>: Sendable, Equatable {
            public let value: Scalar
            public let estimatedError: Scalar
            public let l1Norm: Scalar
            public let iterations: Int
            public let functionEvaluations: Int
            public let converged: Bool
            public let metadata: Metadata

            @inlinable public var didConverge: Bool { converged }

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

        public enum Error: Swift.Error, Equatable {
            case invalidConfiguration(String)
            case invalidInterval(lower: Double, upper: Double)
            case backendUnavailable(String)
        }

        /// Reusable integrator that owns a single quadrature handle.
        public struct Integrator<Scalar: Real & BinaryFloatingPoint & Sendable>: Sendable {
            private final class Storage: @unchecked Sendable {
                let handle: QuadratureHandle
                let metadata: Metadata

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

                func integrate(over interval: Interval,
                               integrand: @escaping @Sendable (Scalar) -> Scalar) throws -> Result<Scalar> {
                    let raw = try QuadratureBridge.integrate(handle: handle,
                                                             over: interval,
                                                             integrand: integrand,
                                                             scalar: Scalar.self)
                    return Result(raw: raw, metadata: metadata)
                }

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

            public var metadata: Metadata { storage.metadata }

            public init(rule: Rule) throws {
                self.storage = try Storage(rule: rule)
            }

            public func integrate(over interval: Interval = .automatic,
                                  integrand: @escaping @Sendable (Scalar) -> Scalar) throws -> Result<Scalar> {
                try storage.integrate(over: interval, integrand: integrand)
            }

            /// Copy abscissa and weights into caller-provided buffers. Returns false if unavailable.
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
    static func checkedCount(_ value: Int, label: String) throws -> Int32 {
        guard value > 0 else {
            throw SpecialFunctions.Quadrature.Error.invalidConfiguration("\(label) requires a positive point count (received \(value)).")
        }
        guard value <= Int(Int32.max) else {
            throw SpecialFunctions.Quadrature.Error.invalidConfiguration("\(label) point count exceeds Int32.max (received \(value)).")
        }
        return Int32(value)
    }

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

private enum QuadratureBridge {
    struct RawResult<Scalar: Real & BinaryFloatingPoint & Sendable> {
        let value: Scalar
        let estimatedError: Scalar
        let l1Norm: Scalar
        let iterations: Int
        let functionEvaluations: Int
        let converged: Bool
    }

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

    static func integrate<Scalar: Real & BinaryFloatingPoint & Sendable>(
        handle: QuadratureHandle,
        over interval: SpecialFunctions.Quadrature.Interval,
        integrand: @escaping @Sendable (Scalar) -> Scalar,
        scalar _: Scalar.Type
    ) throws -> RawResult<Scalar> {
        if Scalar.self == Double.self {
            let raw = try QuadratureBridgeDouble.integrate(handle: handle,
                                                           over: interval,
                                                           integrand: integrand as! @Sendable (Double) -> Double)
            return QuadratureBridge.RawResult(
                value: Scalar(raw.value),
                estimatedError: Scalar(raw.estimatedError),
                l1Norm: Scalar(raw.l1Norm),
                iterations: raw.iterations,
                functionEvaluations: raw.functionEvaluations,
                converged: raw.converged
            )
        } else if Scalar.self == Float.self {
            let raw = try QuadratureBridgeFloat.integrate(handle: handle,
                                                          over: interval,
                                                          integrand: integrand as! @Sendable (Float) -> Float)
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
                let raw = try QuadratureBridgeFloat80.integrate(handle: handle,
                                                                over: interval,
                                                                integrand: integrand as! @Sendable (Float80) -> Float80)
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

private enum QuadratureBridgeDouble {
    static func createHandle(for rule: SpecialFunctions.Quadrature.Rule) throws -> QuadratureHandle {
        switch rule {
        case let .gaussLegendre(points):
            let n = try SpecialFunctions.Quadrature.Rule.checkedCount(points, label: "Gauss-Legendre")
            guard let handle = quad_gauss_create_d(n) else {
                throw SpecialFunctions.Quadrature.Error.invalidConfiguration("Gauss-Legendre does not support \(points) points.")
            }
            return handle
        case let .gaussHermite(points):
            let n = try SpecialFunctions.Quadrature.Rule.checkedCount(points, label: "Gauss-Hermite")
            guard let handle = quad_gauss_hermite_create_d(n) else {
                throw SpecialFunctions.Quadrature.Error.backendUnavailable("Gauss-Hermite is unavailable for \(points) points (requires Boost ≥ 1.78).")
            }
            return handle
        case let .gaussLaguerre(points, alpha):
            let n = try SpecialFunctions.Quadrature.Rule.checkedCount(points, label: "Gauss-Laguerre")
            if let alpha {
                guard alpha > -1 else {
                    throw SpecialFunctions.Quadrature.Error.invalidConfiguration("Gauss-Laguerre alpha must be greater than -1 (received \(alpha)).")
                }
                guard let handle = quad_gauss_laguerre_create_alpha_d(n, alpha) else {
                    throw SpecialFunctions.Quadrature.Error.invalidConfiguration("Gauss-Laguerre does not support \(points) points with alpha=\(alpha).")
                }
                return handle
            } else {
                guard let handle = quad_gauss_laguerre_create_d(n) else {
                    throw SpecialFunctions.Quadrature.Error.invalidConfiguration("Gauss-Laguerre does not support \(points) points.")
                }
                return handle
            }
        case let .gaussJacobi(points, alpha, beta):
            let n = try SpecialFunctions.Quadrature.Rule.checkedCount(points, label: "Gauss-Jacobi")
            guard alpha > -1 else {
                throw SpecialFunctions.Quadrature.Error.invalidConfiguration("Gauss-Jacobi alpha must be greater than -1 (received \(alpha)).")
            }
            guard beta > -1 else {
                throw SpecialFunctions.Quadrature.Error.invalidConfiguration("Gauss-Jacobi beta must be greater than -1 (received \(beta)).")
            }
            guard let handle = quad_gauss_jacobi_create_d(n, alpha, beta) else {
                throw SpecialFunctions.Quadrature.Error.invalidConfiguration("Gauss-Jacobi does not support \(points) points for α=\(alpha), β=\(beta).")
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

private enum QuadratureBridgeFloat {
    static func createHandle(for rule: SpecialFunctions.Quadrature.Rule) throws -> QuadratureHandle {
        switch rule {
        case let .gaussLegendre(points):
            let n = try SpecialFunctions.Quadrature.Rule.checkedCount(points, label: "Gauss-Legendre")
            guard let handle = quad_gauss_create_f(n) else {
                throw SpecialFunctions.Quadrature.Error.invalidConfiguration("Gauss-Legendre does not support \(points) points.")
            }
            return handle
        case let .gaussHermite(points):
            let n = try SpecialFunctions.Quadrature.Rule.checkedCount(points, label: "Gauss-Hermite")
            guard let handle = quad_gauss_hermite_create_f(n) else {
                throw SpecialFunctions.Quadrature.Error.backendUnavailable("Gauss-Hermite is unavailable for \(points) points (requires Boost ≥ 1.78).")
            }
            return handle
        case let .gaussLaguerre(points, alpha):
            let n = try SpecialFunctions.Quadrature.Rule.checkedCount(points, label: "Gauss-Laguerre")
            if let alpha {
                guard alpha > -1 else {
                    throw SpecialFunctions.Quadrature.Error.invalidConfiguration("Gauss-Laguerre alpha must be greater than -1 (received \(alpha)).")
                }
                guard let handle = quad_gauss_laguerre_create_alpha_f(n, Float(alpha)) else {
                    throw SpecialFunctions.Quadrature.Error.invalidConfiguration("Gauss-Laguerre does not support \(points) points with alpha=\(alpha).")
                }
                return handle
            } else {
                guard let handle = quad_gauss_laguerre_create_f(n) else {
                    throw SpecialFunctions.Quadrature.Error.invalidConfiguration("Gauss-Laguerre does not support \(points) points.")
                }
                return handle
            }
        case let .gaussJacobi(points, alpha, beta):
            let n = try SpecialFunctions.Quadrature.Rule.checkedCount(points, label: "Gauss-Jacobi")
            guard alpha > -1 else {
                throw SpecialFunctions.Quadrature.Error.invalidConfiguration("Gauss-Jacobi alpha must be greater than -1 (received \(alpha)).")
            }
            guard beta > -1 else {
                throw SpecialFunctions.Quadrature.Error.invalidConfiguration("Gauss-Jacobi beta must be greater than -1 (received \(beta)).")
            }
            guard let handle = quad_gauss_jacobi_create_f(n, Float(alpha), Float(beta)) else {
                throw SpecialFunctions.Quadrature.Error.invalidConfiguration("Gauss-Jacobi does not support \(points) points for α=\(alpha), β=\(beta).")
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
        case let .gaussHermite(points):
            let n = try SpecialFunctions.Quadrature.Rule.checkedCount(points, label: "Gauss-Hermite")
            guard let handle = quad_gauss_hermite_create_l(n) else {
                throw SpecialFunctions.Quadrature.Error.backendUnavailable("Gauss-Hermite is unavailable for \(points) points (requires Boost ≥ 1.78).")
            }
            return handle
        case let .gaussLaguerre(points, alpha):
            let n = try SpecialFunctions.Quadrature.Rule.checkedCount(points, label: "Gauss-Laguerre")
            if let alpha {
                guard alpha > -1 else {
                    throw SpecialFunctions.Quadrature.Error.invalidConfiguration("Gauss-Laguerre alpha must be greater than -1 (received \(alpha)).")
                }
                guard let handle = quad_gauss_laguerre_create_alpha_l(n, Float80(alpha)) else {
                    throw SpecialFunctions.Quadrature.Error.invalidConfiguration("Gauss-Laguerre does not support \(points) points with alpha=\(alpha).")
                }
                return handle
            } else {
                guard let handle = quad_gauss_laguerre_create_l(n) else {
                    throw SpecialFunctions.Quadrature.Error.invalidConfiguration("Gauss-Laguerre does not support \(points) points.")
                }
                return handle
            }
        case let .gaussJacobi(points, alpha, beta):
            let n = try SpecialFunctions.Quadrature.Rule.checkedCount(points, label: "Gauss-Jacobi")
            guard alpha > -1 else {
                throw SpecialFunctions.Quadrature.Error.invalidConfiguration("Gauss-Jacobi alpha must be greater than -1 (received \(alpha)).")
            }
            guard beta > -1 else {
                throw SpecialFunctions.Quadrature.Error.invalidConfiguration("Gauss-Jacobi beta must be greater than -1 (received \(beta)).")
            }
            guard let handle = quad_gauss_jacobi_create_l(n, Float80(alpha), Float80(beta)) else {
                throw SpecialFunctions.Quadrature.Error.invalidConfiguration("Gauss-Jacobi does not support \(points) points for α=\(alpha), β=\(beta).")
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
