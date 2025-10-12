//
//  Common.swift
//  Math/SpecialFunctions
//
//  Created by VT on 11.10.25.
//  Copyright © 2025 Volker Thieme.
//  License: MIT (see project root)
//
//  This file contains shared types and helpers for the SpecialFunctions module:
//  - Public error enum used by all throwing special functions.
//  - Lightweight conversion helper to `Double` for generic numeric APIs.
//  - A subset of numerically stable wrappers (expm1, log1p, powm1, etc.).
//  - Trigonometric helpers sinPi/cosPi that evaluate sin(πx) and cos(πx).
//
//  All functions validate inputs (e.g., finiteness or range constraints) before
//  delegating to CBoostBridge (Boost.Math) for the numerical implementation.
//

import CBoostBridge
    // MARK: - Errors
    
    /// Errors thrown by special functions when inputs are invalid or outside the domain.
    ///
    /// Use these errors to understand why a computation failed (e.g., argument is
    /// not finite, argument violates domain restrictions, or the function has a
    /// mathematical pole at the requested input).
    public enum SpecialFunctionError: Error & Equatable {
        /// Parameter must be strictly positive (e.g., `n >= 0`).
        case parameterNotPositive(name: String)
        /// Parameter is outside the valid range [min, max].
        case parameterOutOfRange(name: String, min: Double, max: Double)
        /// Parameter is not finite (NaN or ±∞).
        case parameterNotFinite(name: String)
        /// A simple pole occurs at a non-positive integer for this function (e.g., Γ(x)).
        case poleAtNonPositiveInteger(name: String)
        /// Inputs form an invalid combination for the real-valued function.
        case invalidCombination(message: String)
        /// Parameter is greater than max allowed
        case parameterExceedsMaximumIntegerValue(name: String, max: Int)
    }
    
    // MARK: - Helper casts
    
    /// Convert a generic `BinaryFloatingPoint` to `Double`.
    ///
    /// This is a small, inline helper used by generic overloads to funnel
    /// computation through a `Double`-backed C implementation, while preserving
    /// the public API’s generic signature.
    ///
    /// Parameters:
    /// - x: A `BinaryFloatingPoint` value.
    ///
    /// Returns:
    /// - `x` converted to `Double`.
    @usableFromInline internal func D<T: BinaryFloatingPoint>(_ x: T) -> Double { Double(x) }
    
    /// The base of the natural logarithm, e ≈ 2.718281828..., from Boost.Math constants.
    ///
    /// Use this to access a consistent `Double` value sourced from the C bridge.
    @inlinable public var boostE: Double { bs_const_e() }
    
    // MARK: - Numerically stable helpers
    
    /// Compute exp(x) - 1 with improved accuracy for small x.
    ///
    /// Parameters:
    /// - x: The input value `x`.
    ///
    /// Returns:
    /// - `exp(x) - 1` as `T`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotFinite(name: "x")` if `x` is NaN or ±∞.
    ///
    /// Example:
    /// ```swift
    /// let y = try expm1(1e-8 as Double)
    /// ```
    @inlinable public func expm1<T: BinaryFloatingPoint>(_ x: T) throws -> T {
        let dx = D(x)
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return T(bs_expm1(dx))
    }
    
    /// Compute ln(1 + x) with improved accuracy for small x.
    ///
    /// Domain:
    /// - Requires `x > -1` (ln(1+x) undefined for x ≤ -1 in the reals).
    ///
    /// Parameters:
    /// - x: The input value `x`.
    ///
    /// Returns:
    /// - `ln(1 + x)` as `T`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotFinite(name: "x")` if `x` is NaN or ±∞.
    /// - `SpecialFunctionError.parameterOutOfRange(name: "x", ...)` if `x ≤ -1`.
    ///
    /// Example:
    /// ```swift
    /// let y = try log1p(1e-8 as Double)
    /// ```
    @inlinable public func log1p<T: BinaryFloatingPoint>(_ x: T) throws -> T {
        let dx = D(x)
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        guard dx > -1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: -1.0.nextUp, max: Double.infinity) }
        return T(bs_log1p(dx))
    }
    
    /// Compute ln(1 + x) - x with improved accuracy near zero.
    ///
    /// Domain:
    /// - Requires `x > -1`.
    ///
    /// Parameters:
    /// - x: The input value `x`.
    ///
    /// Returns:
    /// - `ln(1 + x) - x` as `T`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotFinite(name: "x")` if `x` is NaN or ±∞.
    /// - `SpecialFunctionError.parameterOutOfRange(name: "x", ...)` if `x ≤ -1`.
    @inlinable public func log1pmx<T: BinaryFloatingPoint>(_ x: T) throws -> T {
        let dx = D(x)
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        guard dx > -1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: -1.0.nextUp, max: Double.infinity) }
        return T(bs_log1pmx(dx))
    }
    
    /// Compute x^y - 1 with improved accuracy, enforcing a real-valued domain.
    ///
    /// Domain:
    /// - Negative base with a non-integer exponent is undefined over the reals.
    ///   This function throws `invalidCombination` in that case.
    ///
    /// Parameters:
    /// - x: The base.
    /// - y: The exponent.
    ///
    /// Returns:
    /// - `x^y - 1` as `T`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotFinite(name: "x"|"y")` if inputs are NaN or ±∞.
    /// - `SpecialFunctionError.invalidCombination(...)` if `x < 0` and `y` is non-integer.
    ///
    /// Example:
    /// ```swift
    /// let v = try powm1(1.000001 as Double, 3.0 as Double)
    /// ```
    @inlinable public func powm1<T: BinaryFloatingPoint>(_ x: T, _ y: T) throws -> T {
        let dx = D(x), dy = D(y)
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        guard dy.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "y") }
        let yIsInteger = dy == dy.rounded(.towardZero)
        guard !(dx < 0 && !yIsInteger) else {
            throw SpecialFunctionError.invalidCombination(message: "powm1 is undefined for negative base with non-integer exponent in the reals")
        }
        return T(bs_powm1(dx, dy))
    }
    
    /// Compute the real cube root of x.
    ///
    /// Parameters:
    /// - x: The input value `x`.
    ///
    /// Returns:
    /// - `cbrt(x)` as `T`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotFinite(name: "x")` if `x` is NaN or ±∞.
    @inlinable public func cbrt<T: BinaryFloatingPoint>(_ x: T) throws -> T {
        let dx = D(x)
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return T(bs_cbrt(dx))
    }
    
    // MARK: - Float overloads
    // These overloads call directly into the Float-precision C implementations for
    // best performance and to avoid unnecessary conversions.
    
    /// Exponential integral E_n(x) for `Float`.
    @inlinable public func exponentialIntegralEn(_ n: Int, _ x: Float) throws -> Float {
        guard n >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "n") }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return bs_expint_En_f(Int32(n), x)
    }
    
    /// Compute exp(x) - 1 for `Float`.
    @inlinable public func expm1(_ x: Float) throws -> Float {
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return bs_expm1_f(x)
    }
    
    /// Compute ln(1 + x) for `Float`. Requires `x > -1`.
    @inlinable public func log1p(_ x: Float) throws -> Float {
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        guard x > -1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: Double(Float(-1).nextUp), max: Double.infinity) }
        return bs_log1p_f(x)
    }
    
    /// Compute ln(1 + x) - x for `Float`. Requires `x > -1`.
    @inlinable public func log1pmx(_ x: Float) throws -> Float {
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        guard x > -1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: Double(Float(-1).nextUp), max: Double.infinity) }
        return bs_log1pmx_f(x)
    }
    
    /// Compute x^y - 1 for `Float`, with real-domain checks.
    @inlinable public func powm1(_ x: Float, _ y: Float) throws -> Float {
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        guard y.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "y") }
        let yIsInteger = y == y.rounded(.towardZero)
        guard !(x < 0 && !yIsInteger) else {
            throw SpecialFunctionError.invalidCombination(message: "powm1 is undefined for negative base with non-integer exponent in the reals")
        }
        return bs_powm1_f(x, y)
    }
    
    /// Compute the real cube root for `Float`.
    @inlinable public func cbrt(_ x: Float) throws -> Float {
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return bs_cbrt_f(x)
    }
    
    // MARK: - Float80 overloads (x86_64)
    // Extended-precision versions for platforms that support Float80.
    
#if arch(x86_64)
    /// Exponential integral E_n(x) for `Float80` (x86_64 only).
    @inlinable public func exponentialIntegralEn(_ n: Int, _ x: Float80) throws -> Float80 {
        guard n >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "n") }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return bs_expint_En_l(Int32(n), x)
    }
    
    /// Compute exp(x) - 1 for `Float80` (x86_64 only).
    @inlinable public func expm1(_ x: Float80) throws -> Float80 {
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return bs_expm1_l(x)
    }
    
    /// Compute ln(1 + x) for `Float80` (x86_64 only). Requires `x > -1`.
    @inlinable public func log1p(_ x: Float80) throws -> Float80 {
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        guard x > -1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: Double(-1).nextUp, max: Double.infinity) }
        return bs_log1p_l(x)
    }
    
    /// Compute ln(1 + x) - x for `Float80` (x86_64 only). Requires `x > -1`.
    @inlinable public func log1pmx(_ x: Float80) throws -> Float80 {
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        guard x > -1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: Double(-1).nextUp, max: Double.infinity) }
        return bs_log1pmx_l(x)
    }
    
    /// Compute x^y - 1 for `Float80` (x86_64 only), with real-domain checks.
    @inlinable public func powm1(_ x: Float80, _ y: Float80) throws -> Float80 {
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        guard y.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "y") }
        let yIsInteger = y == y.rounded(.towardZero)
        guard !(x < 0 && !yIsInteger) else {
            throw SpecialFunctionError.invalidCombination(message: "powm1 is undefined for negative base with non-integer exponent in the reals")
        }
        return bs_powm1_l(x, y)
    }
    
    /// Compute the real cube root for `Float80` (x86_64 only).
    @inlinable public func cbrt(_ x: Float80) throws -> Float80 {
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return bs_cbrt_l(x)
    }
#endif
    
    // MARK: - sinPi / cosPi helpers
    
    /// Compute sin(πx) in a numerically stable way.
    ///
    /// This helper leverages a Boost.Math implementation designed to reduce
    /// catastrophic cancellation near integers.
    ///
    /// Parameters:
    /// - x: The input value `x`.
    ///
    /// Returns:
    /// - `sin(πx)` as `T`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotFinite(name: "x")` if `x` is NaN or ±∞.
    @inlinable public func sinPi<T: BinaryFloatingPoint>(_ x: T) throws -> T {
        let dx = D(x)
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return T(bs_sin_pi(dx))
    }
    
    /// Compute cos(πx) in a numerically stable way.
    ///
    /// Parameters:
    /// - x: The input value `x`.
    ///
    /// Returns:
    /// - `cos(πx)` as `T`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotFinite(name: "x")` if `x` is NaN or ±∞.
    @inlinable public func cosPi<T: BinaryFloatingPoint>(_ x: T) throws -> T {
        let dx = D(x)
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return T(bs_cos_pi(dx))
    }
    
    // MARK: Float overloads for sinPi / cosPi
    
    /// Compute sin(πx) for `Float`.
    @inlinable public func sinPi(_ x: Float) throws -> Float {
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return bs_sin_pi_f(x)
    }
    
    /// Compute cos(πx) for `Float`.
    @inlinable public func cosPi(_ x: Float) throws -> Float {
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return bs_cos_pi_f(x)
    }
    
    // MARK: Float80 overloads for sinPi / cosPi (x86_64)
    
#if arch(x86_64)
    /// Compute sin(πx) for `Float80` (x86_64 only).
    @inlinable public func sinPi(_ x: Float80) throws -> Float80 {
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return bs_sin_pi_l(x)
    }
    
    /// Compute cos(πx) for `Float80` (x86_64 only).
    @inlinable public func cosPi(_ x: Float80) throws -> Float80 {
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return bs_cos_pi_l(x)
    }
#endif
