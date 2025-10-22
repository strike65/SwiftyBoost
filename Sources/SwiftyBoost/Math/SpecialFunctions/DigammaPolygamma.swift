//
//  DigammaPolygammaZeta.swift
//  Math/SpecialFunctions
//
//  This file exposes the digamma ψ(x), trigamma ψ₁(x), general polygamma ψ⁽ⁿ⁾(x),
//  and the Riemann zeta ζ(x) functions with type-generic and type-specific
//  overloads. The implementations are backed by the CBoostBridge C shim
//  (wrapping Boost.Math), while performing Swift-side argument validation
//  (domain and finiteness checks).
//
//  All functions throw Swift errors for invalid inputs (e.g., non-finite values)
//  and at poles where the function is not defined (e.g., non-positive integers
//  for digamma/polygamma, and x = 1 for ζ(x)). See each function’s documentation.
//

import SwiftyBoostPrelude
public extension SpecialFunctions {
    
    
    // MARK: - Generic overloads (Double-backed)
    
    /// Compute the digamma function ψ(x) = d/dx ln Γ(x).
    ///
    /// Domain and poles:
    /// - ψ(x) has simple poles at the non-positive integers (x ∈ {0, -1, -2, ...}).
    ///   Calling this function at those points will throw.
    ///
    /// Behavior and precision:
    /// - For finite, non-pole inputs, the result is finite where defined.
    /// - Very large |x| may lose precision depending on `T`.
    /// - This generic overload converts to `Double` for evaluation and returns `T`.
    ///
    /// Thread-safety:
    /// - Pure and thread-safe (no shared mutable state).
    ///
    /// Parameters:
    /// - x: Input value.
    ///
    /// Returns:
    /// - ψ(x) as `T`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotFinite(name: "x")` if `x` is NaN or ±∞.
    /// - `SpecialFunctionError.poleAtNonPositiveInteger(name: "x")` for non-positive integers.
    ///
    /// Example:
    /// ```swift
    /// let d: Double = try digamma(3.0) // ≈ 0.922784...
    /// ```
    @inlinable static func digamma<T: Real & BinaryFloatingPoint>(_ x: T) throws -> T {
        let dx = D(x)
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        if dx <= 0, dx == dx.rounded(.towardZero) { throw SpecialFunctionError.poleAtNonPositiveInteger(name: "x") }
        return T(bs_digamma_d(dx))
    }
    
    /// Compute the trigamma function ψ₁(x) = d/dx ψ(x) = d²/dx² ln Γ(x).
    ///
    /// Domain and poles:
    /// - ψ₁(x) has simple poles at the non-positive integers (x ∈ {0, -1, -2, ...}).
    ///
    /// Behavior and precision:
    /// - For finite, non-pole inputs, the result is positive and finite where defined.
    /// - This generic overload converts to `Double` for evaluation and returns `T`.
    ///
    /// Parameters:
    /// - x: Input value.
    ///
    /// Returns:
    /// - ψ₁(x) as `T`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotFinite(name: "x")` if `x` is NaN or ±∞.
    /// - `SpecialFunctionError.poleAtNonPositiveInteger(name: "x")` for non-positive integers.
    ///
    /// Example:
    /// ```swift
    /// let t: Double = try trigamma(2.0) // ≈ 0.644934...
    /// ```
    @inlinable static func trigamma<T: Real & BinaryFloatingPoint>(_ x: T) throws -> T {
        let dx = D(x)
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        if dx <= 0, dx == dx.rounded(.towardZero) { throw SpecialFunctionError.poleAtNonPositiveInteger(name: "x") }
        return T(bs_trigamma_d(dx))
    }
    
    /// Compute the polygamma function ψ⁽ⁿ⁾(x), the n-th derivative of digamma.
    ///
    /// Domain and poles:
    /// - ψ⁽ⁿ⁾(x) has simple poles at the non-positive integers (x ∈ {0, -1, -2, ...}).
    /// - Order `n` must be ≥ 0.
    ///
    /// Behavior and precision:
    /// - For n = 0 this is digamma; for n = 1 this is trigamma.
    /// - This generic overload converts to `Double` for evaluation and returns `T`.
    ///
    /// Parameters:
    /// - n: The non-negative integer order of the derivative.
    /// - x: Input value.
    ///
    /// Returns:
    /// - ψ⁽ⁿ⁾(x) as `T`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotPositive(name: "order")` if `n < 0`.
    /// - `SpecialFunctionError.parameterNotFinite(name: "x")` if `x` is NaN or ±∞.
    /// - `SpecialFunctionError.poleAtNonPositiveInteger(name: "x")` for non-positive integers.
    ///
    /// Example:
    /// ```swift
    /// let p2: Double = try polygamma(order: 2, 3.0) // second derivative at x=3
    /// ```
    @inlinable static func polygamma<T: Real & BinaryFloatingPoint>(order n: Int, _ x: T) throws -> T {
        guard n >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "order") }
        let dx = D(x)
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        if dx <= 0, dx == dx.rounded(.towardZero) { throw SpecialFunctionError.poleAtNonPositiveInteger(name: "x") }
        return T(bs_polygamma_d(Int32(n), dx))
    }
    
    /// Compute the Riemann zeta function ζ(x).
    ///
    /// Domain and poles:
    /// - ζ(x) has a simple pole at x = 1. Calling this at x = 1 will throw.
    /// - For real x, ζ(x) is defined for all x ≠ 1. (On the reals, the trivial zeros
    ///   are at the negative even integers.)
    ///
    /// Behavior and precision:
    /// - This generic overload converts to `Double` for evaluation and returns `T`.
    ///
    /// Parameters:
    /// - x: Input value.
    ///
    /// Returns:
    /// - ζ(x) as `T`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotFinite(name: "x")` if `x` is NaN or ±∞.
    /// - `SpecialFunctionError.invalidCombination(message: ...)` if `x == 1`.
    ///
    /// Example:
    /// ```swift
    /// let z2: Double = try riemannZeta(2.0) // π^2 / 6 ≈ 1.644934...
    /// ```
    @inlinable static func riemannZeta<T: Real & BinaryFloatingPoint>(_ x: T) throws -> T {
        let dx = D(x)
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        guard dx != 1 else { throw SpecialFunctionError.invalidCombination(message: "riemannZeta has a pole at x = 1") }
        return T(bs_riemann_zeta_d(dx))
    }
    
    // MARK: - Float overloads
    
    /// Digamma ψ(x) for `Float`.
    @inlinable static func digamma(_ x: Float) throws -> Float {
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        if x <= 0, x == x.rounded(.towardZero) { throw SpecialFunctionError.poleAtNonPositiveInteger(name: "x") }
        return bs_digamma_f(x)
    }
    
    /// Trigamma ψ₁(x) for `Float`.
    @inlinable static func trigamma(_ x: Float) throws -> Float {
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        if x <= 0, x == x.rounded(.towardZero) { throw SpecialFunctionError.poleAtNonPositiveInteger(name: "x") }
        return bs_trigamma_f(x)
    }
    
    /// Polygamma ψ⁽ⁿ⁾(x) for `Float`. Order `n` must be ≥ 0.
    @inlinable static func polygamma(order n: Int, _ x: Float) throws -> Float {
        guard n >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "order") }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        if x <= 0, x == x.rounded(.towardZero) { throw SpecialFunctionError.poleAtNonPositiveInteger(name: "x") }
        return bs_polygamma_f(Int32(n), x)
    }
    
    /// Riemann zeta ζ(x) for `Float`. Throws at x = 1.
    @inlinable static func riemannZeta(_ x: Float) throws -> Float {
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        guard x != 1 else { throw SpecialFunctionError.invalidCombination(message: "riemannZeta has a pole at x = 1") }
        return bs_riemann_zeta_f(x)
    }
    
    // MARK: - Float80 overloads (x86_64)
    
#if arch(x86_64)
    /// Digamma ψ(x) for `Float80` (x86_64 only).
    @inlinable static func digamma(_ x: Float80) throws -> Float80 {
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        if x <= 0, x == x.rounded(.towardZero) { throw SpecialFunctionError.poleAtNonPositiveInteger(name: "x") }
        return bs_digamma_l(x)
    }
    
    /// Trigamma ψ₁(x) for `Float80` (x86_64 only).
    @inlinable static func trigamma(_ x: Float80) throws -> Float80 {
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        if x <= 0, x == x.rounded(.towardZero) { throw SpecialFunctionError.poleAtNonPositiveInteger(name: "x") }
        return bs_trigamma_l(x)
    }
    
    /// Polygamma ψ⁽ⁿ⁾(x) for `Float80` (x86_64 only). Order `n` must be ≥ 0.
    @inlinable static func polygamma(order n: Int, _ x: Float80) throws -> Float80 {
        guard n >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "order") }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        if x <= 0, x == x.rounded(.towardZero) { throw SpecialFunctionError.poleAtNonPositiveInteger(name: "x") }
        return bs_polygamma_l(Int32(n), x)
    }
    
    /// Riemann zeta ζ(x) for `Float80` (x86_64 only). Throws at x = 1.
    @inlinable static func riemannZeta(_ x: Float80) throws -> Float80 {
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        guard x != 1 else { throw SpecialFunctionError.invalidCombination(message: "riemannZeta has a pole at x = 1") }
        return bs_riemann_zeta_l(x)
    }
#endif
    
}
