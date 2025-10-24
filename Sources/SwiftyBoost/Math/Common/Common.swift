//
//  Common.swift
//  Math/SpecialFunctions
//
//  Created by Volker Thieme 2025.
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

import SwiftyBoostPrelude

// MARK: - Errors

/// Errors thrown by special functions when inputs are invalid or outside the domain.
///
/// Use these errors to understand why a computation failed (e.g., an argument
/// is not finite, violates domain restrictions, or the function has a
/// mathematical pole at the requested input).
///
/// You can pattern match on the associated values (e.g., parameter name or
/// min/max bounds) to build user-facing diagnostics.
public enum SpecialFunctionError<T: Real & BinaryFloatingPoint & Sendable>: Error & Equatable {
    /// Parameter must be strictly positive (e.g., `n >= 0`).
    ///
    /// - Parameter name: The name of the offending parameter.
    case parameterNotPositive(name: String, value: T)
    /// Parameter is outside the valid range [min, max].
    ///
    /// - Parameters:
    ///   - name: The name of the offending parameter.
    ///   - min: The inclusive lower bound.
    ///   - max: The inclusive upper bound.
    case parameterOutOfRange(name: String, min: T, max: T)
    /// Parameter is not finite (NaN or ±∞).
    ///
    /// - Parameter name: The name of the offending parameter.
    case parameterNotFinite(name: String, value: T)
    /// A simple pole occurs at a non-positive integer for this function (e.g., Γ(x)).
    ///
    /// - Parameter name: The name of the offending parameter.
    case poleAtNonPositiveInteger(name: String, value: T)
    /// Inputs form an invalid combination for the real-valued function.
    ///
    /// - Parameter message: A human-readable explanation of the invalid combination.
    case invalidCombination(message: String, value: T)
    /// Parameter is greater than the maximum allowed integer value.
    ///
    /// - Parameters:
    ///   - name: The name of the offending parameter.
    ///   - max: The maximum allowed integer value.
    case parameterExceedsMaximumIntegerValue(name: String, max: Int)
    case parameterNotInDomain(name: String, value: T)

}

/// Errors thrown by special functions when inputs are invalid or outside the domain.
///
/// Use these errors to understand why a computation failed (e.g., an argument
/// is not finite, violates domain restrictions, or the function has a
/// mathematical pole at the requested input).
///
/// You can pattern match on the associated values (e.g., parameter name or
/// min/max bounds) to build user-facing diagnostics.
public enum DistributionError<T: Real & BinaryFloatingPoint & Sendable>: Error & Equatable {
    /// Parameter must be strictly positive (e.g., `n >= 0`).
    ///
    /// - Parameter name: The name of the offending parameter.
    case parameterNotPositive(name: String, value: T)
    /// Parameter is outside the valid range [min, max].
    ///
    /// - Parameters:
    ///   - name: The name of the offending parameter.
    ///   - min: The inclusive lower bound.
    ///   - max: The inclusive upper bound.
    case parameterOutOfRange(name: String, min: T, max: T)
    /// Parameter is not finite (NaN or ±∞).
    ///
    /// - Parameter name: The name of the offending parameter.
    case parameterNotFinite(name: String, value: T)
    /// A simple pole occurs at a non-positive integer for this function (e.g., Γ(x)).
    ///
    /// - Parameter name: The name of the offending parameter.
    case poleAtNonPositiveInteger(name: String)
    /// Inputs form an invalid combination for the real-valued function.
    ///
    /// - Parameter message: A human-readable explanation of the invalid combination.
    case invalidCombination(message: String, value: T?)
    /// Parameter is greater than the maximum allowed integer value.
    ///
    /// - Parameters:
    ///   - name: The name of the offending parameter.
    ///   - max: The maximum allowed integer value.
    case parameterExceedsMaximumIntegerValue(name: String, max: Int)
    case generalError(msg: String)
}

// MARK: - Helper casts

/// Converts a generic `BinaryFloatingPoint` to `Double`.
///
/// This small utility enables generic algorithms to call C/Boost-backed
/// implementations that operate on `Double`.
///
/// - Parameter x: A `BinaryFloatingPoint` value.
/// - Returns: `x` converted to `Double`.
@usableFromInline internal func D<T: Real & BinaryFloatingPoint & Sendable>(_ x: T) -> Double { Double(x) }

// MARK: - Numerically stable helpers
extension SpecialFunctions {
    
    /// Computes `exp(x) - 1` with improved numerical stability for small `x`.
    ///
    /// - Parameter x: The exponent.
    /// - Returns: `exp(x) - 1` as `T`.
    /// - Throws: `SpecialFunctionError.parameterNotFinite` if `x` is not finite.
    @inlinable static func expm1<T: Real & BinaryFloatingPoint & Sendable>(_ x: T) throws -> T {
        let dx = D(x)
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        return T(bs_expm1_d(dx))
    }

    /// Computes `log(1 + x)` with improved numerical stability near zero.
    ///
    /// - Parameter x: The input value; must satisfy `x > -1`.
    /// - Returns: `log(1 + x)` as `T`.
    /// - Throws:
    ///   - `SpecialFunctionError.parameterNotFinite` if `x` is not finite.
    ///   - `SpecialFunctionError.parameterOutOfRange` if `x <= -1`.
    @inlinable static func log1p<T: Real & BinaryFloatingPoint & Sendable>(_ x: T) throws -> T {
        let dx = D(x)
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        guard dx > -1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: -1.0.nextUp, max: Double.infinity) }
        return T(bs_log1p_d(dx))
    }
    
    /// Computes `log(1 + x) - x` with improved numerical stability near zero.
    ///
    /// - Parameter x: The input value; must satisfy `x > -1`.
    /// - Returns: `log(1 + x) - x` as `T`.
    /// - Throws:
    ///   - `SpecialFunctionError.parameterNotFinite` if `x` is not finite.
    ///   - `SpecialFunctionError.parameterOutOfRange` if `x <= -1`.
    @inlinable static func log1pmx<T: Real & BinaryFloatingPoint & Sendable>(_ x: T) throws -> T {
        let dx = D(x)
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        guard dx > -1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: -1.0.nextUp, max: Double.infinity) }
        return T(bs_log1pmx_d(dx))
    }
    
    /// Computes `x^y - 1` with improved stability for small results.
    ///
    /// - Parameters:
    ///   - x: The base.
    ///   - y: The exponent.
    /// - Returns: `x^y - 1` as `T`.
    /// - Throws:
    ///   - `SpecialFunctionError.parameterNotFinite` if `x` or `y` is not finite.
    ///   - `SpecialFunctionError.invalidCombination` if `x < 0` and `y` is non-integer
    ///     (undefined in the reals).
    @inlinable static func powm1<T: Real & BinaryFloatingPoint & Sendable>(_ x: T, _ y: T) throws -> T {
        let dx = D(x), dy = D(y)
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        guard dy.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "y", value: y) }
        let yIsInteger = dy == dy.rounded(.towardZero)
        guard !(dx < 0 && !yIsInteger) else {
            throw SpecialFunctionError.invalidCombination(message: "powm1 is undefined for negative base with non-integer exponent in the reals", value: x)
        }
        return T(bs_powm1_d(dx, dy))
    }
    
    /// Computes the real cube root `cbrt(x)`.
    ///
    /// - Parameter x: The input value.
    /// - Returns: The real cube root of `x` as `T`.
    /// - Throws: `SpecialFunctionError.parameterNotFinite` if `x` is not finite.
    @inlinable static func cbrt<T: Real & BinaryFloatingPoint & Sendable>(_ x: T) throws -> T {
        let dx = D(x)
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        return T(bs_cbrt_d(dx))
    }
    
    /// Computes `sqrt(1 + x) - 1` with improved numerical stability near zero.
    ///
    /// - Parameter x: The input value.
    /// - Returns: `sqrt(1 + x) - 1` as `T`.
    /// - Note: This overload does not throw and relies on the underlying Boost implementation.
    @inlinable static func sqrt1pm1<T: Real & BinaryFloatingPoint & Sendable>(x: T) -> T {
        return T(bs_sqrt1pm1_d(Double(x)))
    }
    
    /// Computes the reciprocal square root `1 / sqrt(x)` in a numerically stable way
    /// using Boost.Math’s `rsqrt` backend.
    ///
    /// - Parameter x: Input value (finite). Real-valued branch requires `x ≥ 0`.
    /// - Returns: `1 / sqrt(x)` as `T`.
    /// - Throws:
    ///   - `SpecialFunctionError.parameterNotFinite(name: "x")` if `x` is NaN or ±∞.
    ///   - `SpecialFunctionError.parameterOutOfRange(name: "x", min: 0, max: +∞)` if `x < 0`.
    @inlinable static func rsqrt<T: Real & BinaryFloatingPoint & Sendable>(_ x: T) throws -> T {
        let dx = D(x)
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        guard dx >= 0 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 0.0, max: Double.infinity) }
        return T(bs_rsqrt_d(dx))
    }
    
    //double bs_sqrt1pm1(double x)         { return bs_wrap<double>([&] { return boost::math::sqrt1pm1(x); }); }
    //double bs_hypot_d(double x, double y)  { return bs_wrap<double>([&] { return boost::math::hypot(x, y); }); }
    
    
    // MARK: - Float overloads
    
    /// Computes `exp(x) - 1` for `Float`.
    ///
    /// - Parameter x: The exponent.
    /// - Returns: `exp(x) - 1` as `Float`.
    /// - Throws: `SpecialFunctionError.parameterNotFinite` if `x` is not finite.
    @inlinable static func expm1(_ x: Float) throws -> Float {
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        return bs_expm1_f(x)
    }
    
    /// Computes `log(1 + x)` for `Float`. Requires `x > -1`.
    ///
    /// - Parameter x: The input; must satisfy `x > -1`.
    /// - Returns: `log(1 + x)` as `Float`.
    /// - Throws:
    ///   - `SpecialFunctionError.parameterNotFinite` if `x` is not finite.
    ///   - `SpecialFunctionError.parameterOutOfRange` if `x <= -1`.
    @inlinable static func log1p(_ x: Float) throws -> Float {
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        guard x > -1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: Float(-1).nextUp, max: Float.infinity) }
        return bs_log1p_f(x)
    }
    
    /// Computes `log(1 + x) - x` for `Float`. Requires `x > -1`.
    ///
    /// - Parameter x: The input; must satisfy `x > -1`.
    /// - Returns: `log(1 + x) - x` as `Float`.
    /// - Throws:
    ///   - `SpecialFunctionError.parameterNotFinite` if `x` is not finite.
    ///   - `SpecialFunctionError.parameterOutOfRange` if `x <= -1`.
    @inlinable static func log1pmx(_ x: Float) throws -> Float {
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        guard x > -1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: Float(-1).nextUp, max: Float.infinity) }
        return bs_log1pmx_f(x)
    }
    
    /// Computes `x^y - 1` for `Float`, with real-domain checks.
    ///
    /// - Parameters:
    ///   - x: The base.
    ///   - y: The exponent.
    /// - Returns: `x^y - 1` as `Float`.
    /// - Throws:
    ///   - `SpecialFunctionError.parameterNotFinite` if `x` or `y` is not finite.
    ///   - `SpecialFunctionError.invalidCombination` if `x < 0` and `y` is non-integer.
    @inlinable static func powm1(_ x: Float, _ y: Float) throws -> Float {
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        guard y.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "y", value: y) }
        let yIsInteger = y == y.rounded(.towardZero)
        guard !(x < 0 && !yIsInteger) else {
            throw SpecialFunctionError.invalidCombination(message: "powm1 is undefined for negative base with non-integer exponent in the reals", value: x)
        }
        return bs_powm1_f(x, y)
    }
    
    /// Computes the real cube root for `Float`.
    ///
    /// - Parameter x: The input value.
    /// - Returns: The real cube root of `x` as `Float`.
    /// - Throws: `SpecialFunctionError.parameterNotFinite` if `x` is not finite.
    @inlinable static func cbrt(_ x: Float) throws -> Float {
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        return bs_cbrt_f(x)
    }
    
    /// Computes the reciprocal square root `1 / sqrt(x)` for `Float`.
    ///
    /// - Parameter x: Input value (finite). Real-valued branch requires `x ≥ 0`.
    /// - Returns: `1 / sqrt(x)` as `Float`.
    /// - Throws:
    ///   - `SpecialFunctionError.parameterNotFinite(name: "x")` if `x` is NaN or ±∞.
    ///   - `SpecialFunctionError.parameterOutOfRange(name: "x", min: 0, max: +∞)` if `x < 0`.
    @inlinable static func rsqrt(_ x: Float) throws -> Float {
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        guard x >= 0 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 0.0, max: Float.infinity) }
        return bs_rsqrt_f(x)
    }
    
    // MARK: - Float80 overloads (x86_64)
    
#if arch(x86_64)
    
    /// Computes `exp(x) - 1` for `Float80` (x86_64 only).
    ///
    /// - Parameter x: The exponent.
    /// - Returns: `exp(x) - 1` as `Float80`.
    /// - Throws: `SpecialFunctionError.parameterNotFinite` if `x` is not finite.
    @inlinable static func expm1(_ x: Float80) throws -> Float80 {
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        return bs_expm1_l(x)
    }
    
    /// Computes `log(1 + x)` for `Float80` (x86_64 only). Requires `x > -1`.
    ///
    /// - Parameter x: The input; must satisfy `x > -1`.
    /// - Returns: `log(1 + x)` as `Float80`.
    /// - Throws:
    ///   - `SpecialFunctionError.parameterNotFinite` if `x` is not finite.
    ///   - `SpecialFunctionError.parameterOutOfRange` if `x <= -1`.
    @inlinable static func log1p(_ x: Float80) throws -> Float80 {
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        guard x > -1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: Double(-1).nextUp, max: Double.infinity) }
        return bs_log1p_l(x)
    }
    
    /// Computes `log(1 + x) - x` for `Float80` (x86_64 only). Requires `x > -1`.
    ///
    /// - Parameter x: The input; must satisfy `x > -1`.
    /// - Returns: `log(1 + x) - x` as `Float80`.
    /// - Throws:
    ///   - `SpecialFunctionError.parameterNotFinite` if `x` is not finite.
    ///   - `SpecialFunctionError.parameterOutOfRange` if `x <= -1`.
    @inlinable static func log1pmx(_ x: Float80) throws -> Float80 {
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        guard x > -1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: Double(-1).nextUp, max: Double.infinity) }
        return bs_log1pmx_l(x)
    }
    
    /// Computes `x^y - 1` for `Float80` (x86_64 only), with real-domain checks.
    ///
    /// - Parameters:
    ///   - x: The base.
    ///   - y: The exponent.
    /// - Returns: `x^y - 1` as `Float80`.
    /// - Throws:
    ///   - `SpecialFunctionError.parameterNotFinite` if `x` or `y` is not finite.
    ///   - `SpecialFunctionError.invalidCombination` if `x < 0` and `y` is non-integer.
    @inlinable static func powm1(_ x: Float80, _ y: Float80) throws -> Float80 {
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        guard y.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "y", value: y) }
        let yIsInteger = y == y.rounded(.towardZero)
        guard !(x < 0 && !yIsInteger) else {
            throw SpecialFunctionError.invalidCombination(message: "powm1 is undefined for negative base with non-integer exponent in the reals", value: x)
        }
        return bs_powm1_l(x, y)
    }
    
    /// Computes the real cube root for `Float80` (x86_64 only).
    ///
    /// - Parameter x: The input value.
    /// - Returns: The real cube root of `x` as `Float80`.
    /// - Throws: `SpecialFunctionError.parameterNotFinite` if `x` is not finite.
    @inlinable static func cbrt(_ x: Float80) throws -> Float80 {
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        return bs_cbrt_l(x)
    }
    
    /// Computes the reciprocal square root `1 / sqrt(x)` for `Float80` (x86_64 only).
    ///
    /// - Parameter x: Input value (finite). Real-valued branch requires `x ≥ 0`.
    /// - Returns: `1 / sqrt(x)` as `Float80`.
    /// - Throws:
    ///   - `SpecialFunctionError.parameterNotFinite(name: "x")` if `x` is NaN or ±∞.
    ///   - `SpecialFunctionError.parameterOutOfRange(name: "x", min: 0, max: +∞)` if `x < 0`.
    @inlinable static func rsqrt(_ x: Float80) throws -> Float80 {
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        guard x >= 0 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 0.0, max: Double.infinity) }
        return bs_rsqrt_l(x)
    }
#endif
    
    // MARK: - sinPi / cosPi helpers
    
    /// Computes `sin(πx)` with improved accuracy for half- and integer-aligned arguments.
    ///
    /// - Parameter x: The input value.
    /// - Returns: `sin(πx)` as `T`.
    /// - Throws: `SpecialFunctionError.parameterNotFinite` if `x` is not finite.
    @inlinable static func sinPi<T: Real & BinaryFloatingPoint & Sendable>(_ x: T) throws -> T {
        let dx = D(x)
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        return T(bs_sin_pi_d(dx))
    }
    
    /// Computes `cos(πx)` with improved accuracy for half- and integer-aligned arguments.
    ///
    /// - Parameter x: The input value.
    /// - Returns: `cos(πx)` as `T`.
    /// - Throws: `SpecialFunctionError.parameterNotFinite` if `x` is not finite.
    @inlinable static func cosPi<T: Real & BinaryFloatingPoint & Sendable>(_ x: T) throws -> T {
        let dx = D(x)
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        return T(bs_cos_pi_d(dx))
    }
    
    // MARK: Float overloads for sinPi / cosPi
    
    /// Computes `sin(πx)` for `Float`.
    ///
    /// - Parameter x: The input value.
    /// - Returns: `sin(πx)` as `Float`.
    /// - Throws: `SpecialFunctionError.parameterNotFinite` if `x` is not finite.
    @inlinable static func sinPi(_ x: Float) throws -> Float {
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        return bs_sin_pi_f(x)
    }
    
    /// Computes `cos(πx)` for `Float`.
    ///
    /// - Parameter x: The input value.
    /// - Returns: `cos(πx)` as `Float`.
    /// - Throws: `SpecialFunctionError.parameterNotFinite` if `x` is not finite.
    @inlinable static func cosPi(_ x: Float) throws -> Float {
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        return bs_cos_pi_f(x)
    }
    
    // MARK: Float80 overloads for sinPi / cosPi (x86_64)
    
#if arch(x86_64)
    /// Computes `sin(πx)` for `Float80` (x86_64 only).
    ///
    /// - Parameter x: The input value.
    /// - Returns: `sin(πx)` as `Float80`.
    /// - Throws: `SpecialFunctionError.parameterNotFinite` if `x` is not finite.
    @inlinable static func sinPi(_ x: Float80) throws -> Float80 {
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        return bs_sin_pi_l(x)
    }
    
    /// Computes `cos(πx)` for `Float80` (x86_64 only).
    ///
    /// - Parameter x: The input value.
    /// - Returns: `cos(πx)` as `Float80`.
    /// - Throws: `SpecialFunctionError.parameterNotFinite` if `x` is not finite.
    @inlinable static func cosPi(_ x: Float80) throws -> Float80 {
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        return bs_cos_pi_l(x)
    }
#endif
}

