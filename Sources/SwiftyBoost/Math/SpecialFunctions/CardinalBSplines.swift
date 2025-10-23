//
//  Created by Volker Thieme 2025.
//  Copyright © 2025 Volker Thieme.
//
//  Permission is hereby granted, free of charge, to any person obtaining a copy
//  of this software and associated documentation files (the "Software"), to deal
//  in the Software without restriction, including without limitation the rights
//  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//  copies of the Software, and to permit persons to whom the Software is
//  furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included in
//  all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
//  THE SOFTWARE.
//  
//
import SwiftyBoostPrelude
/// Cardinal B-spline special functions and their derivatives.
///
/// This extension exposes a family of functions for evaluating the cardinal B-spline
/// of order `n` as well as its first and second derivatives (often referred to as prime
/// and double-prime). Implementations are backed by functions provided via `CBoostBridge`
/// (thin wrappers around Boost.Math).
///
/// - Terminology:
///   - The “cardinal B-spline” of order `n` is the compactly-supported, non-negative
///     spline basis function of degree `n - 1`.
///   - The “prime” function computes the first derivative with respect to `x`.
///   - The “double prime” function computes the second derivative with respect to `x`.
///
/// - Domain:
///   - All functions here are defined on the compact support `x ∈ [-1, 1]` in this
///     normalized, centered formulation. Passing values outside this range will throw
///     `SpecialFunctionError.parameterOutOfRange`.
///
/// - Overloads and types:
///   - Generic overloads accept any `BinaryFloatingPoint` and forward to the underlying
///     double-precision implementation.
///   - Specialized overloads are provided for `Float` and, on x86 architectures,
///     `Float80` to avoid intermediate conversions.
///
/// - Throws: `SpecialFunctionError.parameterOutOfRange` if `x` is outside `[-1, 1]`.
///   For derivative functions, an additional range check applies to `n`.
public extension SpecialFunctions {
    
    /// Evaluates the cardinal B-spline of order `n` at position `x`.
    ///
    /// This generic overload accepts any `BinaryFloatingPoint` and internally
    /// evaluates using double precision before converting the result back to `T`.
    ///
    /// - Parameters:
    ///   - n: The spline order (a positive integer). The spline has degree `n - 1`.
    ///   - x: The evaluation point, which must lie in `[-1, 1]` for this normalized form.
    /// - Returns: The cardinal B-spline value `B_n(x)`.
    /// - Throws: `SpecialFunctionError.parameterOutOfRange` if `x ∉ [-1, 1]`.
    ///
    /// - Example:
    ///   ```swift
    ///   let y: Double = try SpecialFunctions.cardinal_B_Spline(4, 0.25)
    ///   ```
    static func cardinal_B_Spline<T: Real & BinaryFloatingPoint & Sendable>(_ n: Int, _ x: T) throws -> T {
        guard x <= 1, x >= -1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 1, max: 1) }
        return T(bs_cardinal_b_spline_d(UInt32(n), D(x)))
    }
    
    /// Evaluates the cardinal B-spline of order `n` at position `x` (single precision).
    ///
    /// - Parameters:
    ///   - n: The spline order (a positive integer). The spline has degree `n - 1`.
    ///   - x: The evaluation point, which must lie in `[-1, 1]`.
    /// - Returns: The cardinal B-spline value `B_n(x)` as `Float`.
    /// - Throws: `SpecialFunctionError.parameterOutOfRange` if `x ∉ [-1, 1]`.
    ///
    /// - Example:
    ///   ```swift
    ///   let y: Float = try SpecialFunctions.cardinal_B_Spline(3, 0.0 as Float)
    ///   ```
    static func cardinal_B_Spline(_ n: Int, _ x: Float) throws -> Float {
        guard x <= 1, x >= -1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 1, max: 1) }
        return bs_cardinal_b_spline_f(UInt32(n), x)
    }
    
    #if arch(i386) || arch(x86_64)
    /// Evaluates the cardinal B-spline of order `n` at position `x` (extended precision).
    ///
    /// Available on x86 architectures that support `Float80`.
    ///
    /// - Parameters:
    ///   - n: The spline order (a positive integer). The spline has degree `n - 1`.
    ///   - x: The evaluation point, which must lie in `[-1, 1]`.
    /// - Returns: The cardinal B-spline value `B_n(x)` as `Float80`.
    /// - Throws: `SpecialFunctionError.parameterOutOfRange` if `x ∉ [-1, 1]`.
    static func cardinal_B_Spline(_ n: Int, _ x: Float80) throws -> Float80 {
        guard x <= 1, x >= -1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 1, max: 1) }
        return bs_cardinal_b_spline_l(UInt32(n), x)
    }
    #endif

    /// Evaluates the first derivative of the cardinal B-spline of order `n` at position `x`.
    ///
    /// This generic overload accepts any `BinaryFloatingPoint` and internally
    /// evaluates using double precision before converting the result back to `T`.
    ///
    /// - Parameters:
    ///   - n: The spline order. Must satisfy `n ≥ 3` for the derivative to be well-defined here.
    ///   - x: The evaluation point, which must lie in `[-1, 1]`.
    /// - Returns: The first derivative `B'_n(x)`.
    /// - Throws:
    ///   - `SpecialFunctionError.parameterOutOfRange` if `x ∉ [-1, 1]`.
    ///   - `SpecialFunctionError.parameterOutOfRange` if `n < 3`.
    ///
    /// - Example:
    ///   ```swift
    ///   let dy: Double = try SpecialFunctions.cardinal_B_Spline_prime(5, 0.1)
    ///   ```
    static func cardinal_B_Spline_prime<T: Real & BinaryFloatingPoint & Sendable>(_ n: Int, _ x: T) throws -> T {
        guard x <= 1, x >= -1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 1, max: 1) }
        guard  n >= 3 else { throw SpecialFunctionError.parameterOutOfRange(name: "n", min: 3, max: .infinity) }
        return T(bs_cardinal_b_spline_prime_d(UInt32(n), D(x)))
    }
    
    /// Evaluates the first derivative of the cardinal B-spline of order `n` at position `x` (single precision).
    ///
    /// - Parameters:
    ///   - n: The spline order. Must satisfy `n ≥ 3`.
    ///   - x: The evaluation point, which must lie in `[-1, 1]`.
    /// - Returns: The first derivative `B'_n(x)` as `Float`.
    /// - Throws:
    ///   - `SpecialFunctionError.parameterOutOfRange` if `x ∉ [-1, 1]`.
    ///   - `SpecialFunctionError.parameterOutOfRange` if `n < 3`.
    static func cardinal_B_Spline_prime(_ n: Int, _ x: Float) throws -> Float {
        guard x <= 1, x >= -1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 1, max: 1) }
        guard  n >= 3 else { throw SpecialFunctionError.parameterOutOfRange(name: "n", min: 3, max: .infinity) }
        return bs_cardinal_b_spline_prime_f(UInt32(n), x)
    }
    
    #if arch(i386) || arch(x86_64)
    /// Evaluates the first derivative of the cardinal B-spline of order `n` at position `x` (extended precision).
    ///
    /// Available on x86 architectures that support `Float80`.
    ///
    /// - Parameters:
    ///   - n: The spline order. Must satisfy `n ≥ 3`.
    ///   - x: The evaluation point, which must lie in `[-1, 1]`.
    /// - Returns: The first derivative `B'_n(x)` as `Float80`.
    /// - Throws:
    ///   - `SpecialFunctionError.parameterOutOfRange` if `x ∉ [-1, 1]`.
    ///   - `SpecialFunctionError.parameterOutOfRange` if `n < 3`.
    static func cardinal_B_Spline_prime(_ n: Int, _ x: Float80) throws -> Float80 {
        guard x <= 1, x >= -1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 1, max: 1) }
        guard  n >= 3 else { throw SpecialFunctionError.parameterOutOfRange(name: "n", min: 3, max: .infinity) }
        return bs_cardinal_b_spline_prime_l(UInt32(n), x)
    }
    #endif

    /// Evaluates the second derivative of the cardinal B-spline of order `n` at position `x`.
    ///
    /// This generic overload accepts any `BinaryFloatingPoint` and internally
    /// evaluates using double precision before converting the result back to `T`.
    ///
    /// - Parameters:
    ///   - n: The spline order. Must satisfy `n ≥ 3`.
    ///   - x: The evaluation point, which must lie in `[-1, 1]`.
    /// - Returns: The second derivative `B''_n(x)`.
    /// - Throws:
    ///   - `SpecialFunctionError.parameterOutOfRange` if `x ∉ [-1, 1]`.
    ///   - `SpecialFunctionError.parameterOutOfRange` if `n < 3`.
    ///
    /// - Example:
    ///   ```swift
    ///   let d2y: Double = try SpecialFunctions.cardinal_B_Spline_double_prime(6, -0.2)
    ///   ```
    static func cardinal_B_Spline_double_prime<T: Real & BinaryFloatingPoint & Sendable>(_ n: Int, _ x: T) throws -> T {
        guard x <= 1, x >= -1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 1, max: 1) }
        guard  n >= 3 else { throw SpecialFunctionError.parameterOutOfRange(name: "n", min: 3, max: .infinity) }
        return T(bs_cardinal_b_spline_double_prime_d(UInt32(n), D(x)))
    }
    
    /// Evaluates the second derivative of the cardinal B-spline of order `n` at position `x` (single precision).
    ///
    /// - Parameters:
    ///   - n: The spline order. Must satisfy `n ≥ 3`.
    ///   - x: The evaluation point, which must lie in `[-1, 1]`.
    /// - Returns: The second derivative `B''_n(x)` as `Float`.
    /// - Throws:
    ///   - `SpecialFunctionError.parameterOutOfRange` if `x ∉ [-1, 1]`.
    ///   - `SpecialFunctionError.parameterOutOfRange` if `n < 3`.
    static func cardinal_B_Spline_double_prime(_ n: Int, _ x: Float) throws -> Float {
        guard x <= 1, x >= -1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 1, max: 1) }
        guard  n >= 3 else { throw SpecialFunctionError.parameterOutOfRange(name: "n", min: 3, max: .infinity) }
        return bs_cardinal_b_spline_double_prime_f(UInt32(n), x)
    }
    
    #if arch(i386) || arch(x86_64)
    /// Evaluates the second derivative of the cardinal B-spline of order `n` at position `x` (extended precision).
    ///
    /// Available on x86 architectures that support `Float80`.
    ///
    /// - Parameters:
    ///   - n: The spline order. Must satisfy `n ≥ 3`.
    ///   - x: The evaluation point, which must lie in `[-1, 1]`.
    /// - Returns: The second derivative `B''_n(x)` as `Float80`.
    /// - Throws:
    ///   - `SpecialFunctionError.parameterOutOfRange` if `x ∉ [-1, 1]`.
    ///   - `SpecialFunctionError.parameterOutOfRange` if `n < 3`.
    static func cardinal_B_Spline_double_prime(_ n: Int, _ x: Float80) throws -> Float80 {
        guard x <= 1, x >= -1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 1, max: 1) }
        guard  n >= 3 else { throw SpecialFunctionError.parameterOutOfRange(name: "n", min: 3, max: .infinity) }
        return bs_cardinal_b_spline_double_prime_l(UInt32(n), x)
    }
    #endif
}
