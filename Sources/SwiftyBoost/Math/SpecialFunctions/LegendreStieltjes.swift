//
//  LegendreStieltjes.swift
//  Math/SpecialFunctions
//
//  Created by Volker Thieme 2025.
//  License: MIT (see project root)
//
//  Swift wrappers for the Legendre–Stieltjes polynomials used in Gauss–Kronrod
//  quadrature extensions. Validation mirrors the domain required by the Boost
//  implementation (order m > 0, finite real x) and delegates evaluation to the
//  C bridge which wraps `boost::math::legendre_stieltjes`.
//

import SwiftyBoostPrelude

public extension SpecialFunctions {

    // MARK: - Generic overloads (Double-backed)

    /// Evaluate the Legendre–Stieltjes polynomial Eₘ(x) of order `m` at `x`.
    ///
    /// The Legendre–Stieltjes polynomials appear in the construction of Gauss–Kronrod
    /// quadrature rules as the “Kronrod extension” factor. This routine forwards to Boost.Math
    /// for the actual computation and returns the result in the requested floating-point type.
    ///
    /// Domain and behavior:
    /// - Requires `m > 0`. An error is thrown otherwise.
    /// - `x` must be finite. An error is thrown if `x` is NaN or ±∞.
    /// - Computation is performed in `Double` precision under the hood for generic `T`, and then
    ///   converted to `T`. If you need native precision, prefer the concrete `Float`/`Float80`
    ///   overloads where available.
    ///
    /// Error handling:
    /// - Throws `SpecialFunctionError<T>.parameterNotPositive(name: "m", value: ...)` if `m <= 0`.
    /// - Throws `SpecialFunctionError<T>.parameterExceedsMaximumIntegerValue(name: "m", max: ...)`
    ///   if `m` cannot be represented as `UInt32`.
    /// - Throws `SpecialFunctionError<T>.parameterNotFinite(name: "x", value: ...)` if `x` is not finite.
    ///
    /// Examples:
    /// ```swift
    /// let y: Double = try SpecialFunctions.legendreStieltjes(3, 0.25)
    /// let z: Float  = try SpecialFunctions.legendreStieltjes(5, Float(0.0)) // includes zero when m is odd
    /// ```
    ///
    /// - Parameters:
    ///   - m: Polynomial order (must satisfy `m > 0`).
    ///   - x: Evaluation point (finite real).
    /// - Returns: `Eₘ(x)` as `T`.
    /// - Throws: A `SpecialFunctionError<T>` describing invalid parameters.
    @inlinable static func legendreStieltjes<T: Real & BinaryFloatingPoint & Sendable>(_ m: Int, _ x: T) throws -> T {
        let mu = try validateLegendreStieltjesOrder(m, as: T.self)
        let dx = D(x)
        guard dx.isFinite else { throw SpecialFunctionError<T>.parameterNotFinite(name: "x", value: x) }
        return T(bs_legendre_stieltjes_d(mu, dx))
    }

    /// Derivative Eₘ′(x) of the Legendre–Stieltjes polynomial of order `m` at `x`.
    ///
    /// Domain and behavior mirror `legendreStieltjes(_:_: )`.
    ///
    /// Error handling:
    /// - Throws `SpecialFunctionError<T>.parameterNotPositive(name: "m", value: ...)` if `m <= 0`.
    /// - Throws `SpecialFunctionError<T>.parameterExceedsMaximumIntegerValue(name: "m", max: ...)`
    ///   if `m` cannot be represented as `UInt32`.
    /// - Throws `SpecialFunctionError<T>.parameterNotFinite(name: "x", value: ...)` if `x` is not finite.
    ///
    /// Example:
    /// ```swift
    /// let d: Double = try SpecialFunctions.legendreStieltjesPrime(4, 0.3)
    /// ```
    ///
    /// - Parameters:
    ///   - m: Polynomial order (must satisfy `m > 0`).
    ///   - x: Evaluation point (finite real).
    /// - Returns: `Eₘ′(x)` as `T`.
    /// - Throws: A `SpecialFunctionError<T>` describing invalid parameters.
    @inlinable static func legendreStieltjesPrime<T: Real & BinaryFloatingPoint & Sendable>(_ m: Int, _ x: T) throws -> T {
        let mu = try validateLegendreStieltjesOrder(m, as: T.self)
        let dx = D(x)
        guard dx.isFinite else { throw SpecialFunctionError<T>.parameterNotFinite(name: "x", value: x) }
        return T(bs_legendre_stieltjes_prime_d(mu, dx))
    }

    /// Squared L² norm ‖Eₘ‖² of the Legendre–Stieltjes polynomial of order `m`.
    ///
    /// Returns the squared norm as defined by Boost.Math for the family used in
    /// Gauss–Kronrod extensions.
    ///
    /// Error handling:
    /// - Throws `SpecialFunctionError<T>.parameterNotPositive(name: "m", value: ...)` if `m <= 0`.
    /// - Throws `SpecialFunctionError<T>.parameterExceedsMaximumIntegerValue(name: "m", max: ...)`
    ///   if `m` cannot be represented as `UInt32`.
    ///
    /// Example:
    /// ```swift
    /// let n2: Double = try SpecialFunctions.legendreStieltjesNormSquared(7)
    /// ```
    ///
    /// - Parameter m: Polynomial order (must satisfy `m > 0`).
    /// - Returns: ‖Eₘ‖² as `T`.
    /// - Throws: A `SpecialFunctionError<T>` describing invalid parameters.
    @inlinable static func legendreStieltjesNormSquared<T: Real & BinaryFloatingPoint & Sendable>(_ m: Int) throws -> T {
        let mu = try validateLegendreStieltjesOrder(m, as: T.self)
        return T(bs_legendre_stieltjes_norm_sq_d(mu))
    }

    /// Zeros of the Legendre–Stieltjes polynomial Eₘ(x) of order `m`.
    ///
    /// Returns the real roots in ascending order as provided by Boost.Math.
    ///
    /// Notes:
    /// - `Eₘ` has `m` zeros when `m` is odd (including 0) and `m − 1` zeros when `m` is even.
    /// - For large `m`, root density increases near the endpoints; numerical ordering is preserved.
    ///
    /// Error handling:
    /// - Throws `SpecialFunctionError<T>.parameterNotPositive(name: "m", value: ...)` if `m <= 0`.
    /// - Throws `SpecialFunctionError<T>.parameterExceedsMaximumIntegerValue(name: "m", max: ...)`
    ///   if `m` cannot be represented as `UInt32`.
    ///
    /// Example:
    /// ```swift
    /// let roots: [Double] = try SpecialFunctions.legendreStieltjesZeros(order: 9)
    /// ```
    ///
    /// - Parameter order: Polynomial order (`m > 0`).
    /// - Returns: The real roots in ascending order as provided by Boost.
    @inlinable static func legendreStieltjesZeros<T: Real & BinaryFloatingPoint & Sendable>(order m: Int) throws -> [T] {
        let mu = try validateLegendreStieltjesOrder(m, as: T.self)
        let count = Int(bs_legendre_stieltjes_zeros_d(mu, nil, 0))
        guard count > 0 else { return [] }
        var tmp = Array<Double>(repeating: .zero, count: count)
        tmp.withUnsafeMutableBufferPointer { buf in
            _ = bs_legendre_stieltjes_zeros_d(mu, buf.baseAddress, count)
        }
        return tmp.map(T.init)
    }

    // MARK: - Float overloads

    /// Float-precision evaluation of Eₘ(x).
    ///
    /// See the generic overload for domain and error details. This variant computes in `Float`
    /// precision end-to-end.
    @inlinable static func legendreStieltjes(_ m: Int, _ x: Float) throws -> Float {
        let mu = try validateLegendreStieltjesOrder(m, as: Float.self)
        guard x.isFinite else { throw SpecialFunctionError<Float>.parameterNotFinite(name: "x", value: x) }
        return bs_legendre_stieltjes_f(mu, x)
    }

    /// Float-precision evaluation of Eₘ′(x).
    ///
    /// See the generic overload for domain and error details. This variant computes in `Float`
    /// precision end-to-end.
    @inlinable static func legendreStieltjesPrime(_ m: Int, _ x: Float) throws -> Float {
        let mu = try validateLegendreStieltjesOrder(m, as: Float.self)
        guard x.isFinite else { throw SpecialFunctionError<Float>.parameterNotFinite(name: "x", value: x) }
        return bs_legendre_stieltjes_prime_f(mu, x)
    }

    /// Float-precision squared norm ‖Eₘ‖².
    @inlinable static func legendreStieltjesNormSquared(_ m: Int) throws -> Float {
        let mu = try validateLegendreStieltjesOrder(m, as: Float.self)
        return bs_legendre_stieltjes_norm_sq_f(mu)
    }

    /// Float-precision zeros of Eₘ in ascending order.
    @inlinable static func legendreStieltjesZeros(order m: Int) throws -> [Float] {
        let mu = try validateLegendreStieltjesOrder(m, as: Float.self)
        let count = Int(bs_legendre_stieltjes_zeros_f(mu, nil, 0))
        guard count > 0 else { return [] }
        var tmp = Array<Float>(repeating: .zero, count: count)
        tmp.withUnsafeMutableBufferPointer { buf in
            _ = bs_legendre_stieltjes_zeros_f(mu, buf.baseAddress, count)
        }
        return tmp
    }

    // MARK: - Float80 overloads (x86 only)

    /// Extended-precision evaluation of Eₘ(x) using `Float80` (x86 only).
    ///
    /// Availability:
    /// - Only on `arch(x86_64)` and `arch(i386)`.
    ///
    /// See the generic overload for domain and error details.
    #if arch(x86_64) || arch(i386)
    @inlinable static func legendreStieltjes(_ m: Int, _ x: Float80) throws -> Float80 {
        let mu = try validateLegendreStieltjesOrder(m, as: Float80.self)
        guard x.isFinite else { throw SpecialFunctionError<Float80>.parameterNotFinite(name: "x", value: x) }
        return bs_legendre_stieltjes_l(mu, x)
    }

    /// Extended-precision evaluation of Eₘ′(x) using `Float80` (x86 only).
    @inlinable static func legendreStieltjesPrime(_ m: Int, _ x: Float80) throws -> Float80 {
        let mu = try validateLegendreStieltjesOrder(m, as: Float80.self)
        guard x.isFinite else { throw SpecialFunctionError<Float80>.parameterNotFinite(name: "x", value: x) }
        return bs_legendre_stieltjes_prime_l(mu, x)
    }

    /// Extended-precision squared norm ‖Eₘ‖² using `Float80` (x86 only).
    @inlinable static func legendreStieltjesNormSquared(_ m: Int) throws -> Float80 {
        let mu = try validateLegendreStieltjesOrder(m, as: Float80.self)
        return bs_legendre_stieltjes_norm_sq_l(mu)
    }

    /// Extended-precision zeros of Eₘ in ascending order using `Float80` (x86 only).
    @inlinable static func legendreStieltjesZeros(order m: Int) throws -> [Float80] {
        let mu = try validateLegendreStieltjesOrder(m, as: Float80.self)
        let count = Int(bs_legendre_stieltjes_zeros_l(mu, nil, 0))
        guard count > 0 else { return [] }
        var tmp = Array<Float80>(repeating: .zero, count: count)
        tmp.withUnsafeMutableBufferPointer { buf in
            _ = bs_legendre_stieltjes_zeros_l(mu, buf.baseAddress, count)
        }
        return tmp
    }
    #endif

    // MARK: - Validation helper

    /// Validate the Legendre–Stieltjes order `m` and convert to `UInt32`.
    ///
    /// - Requires: `m > 0`.
    /// - Throws:
    ///   - `SpecialFunctionError<T>.parameterNotPositive(name: "m", value: ...)` if `m <= 0`.
    ///   - `SpecialFunctionError<T>.parameterExceedsMaximumIntegerValue(name: "m", max: ...)`
    ///     if `m` cannot be represented as `UInt32`.
    /// - Returns: The order as `UInt32` for the C bridge.
    @usableFromInline internal static func validateLegendreStieltjesOrder<T: Real & BinaryFloatingPoint & Sendable>(_ m: Int, as _: T.Type) throws -> UInt32 {
        guard m > 0 else { throw SpecialFunctionError<T>.parameterNotPositive(name: "m", value: T(m)) }
        guard let mu = UInt32(exactly: m) else {
            throw SpecialFunctionError<T>.parameterExceedsMaximumIntegerValue(name: "m", max: Int(UInt32.max))
        }
        return mu
    }
}
