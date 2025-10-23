//
//  Created by Volker Thieme
//  Copyright © 2025 Volker Thieme. All rights reserved.
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
//  Hyperbolic arc functions backed by Boost.Math via CBoostBridge.
//
//  This file provides overloaded implementations of the inverse hyperbolic
//  functions asinh, acosh, and atanh for standard Swift floating-point types.
//  The generic overloads accept any BinaryFloatingPoint and forward to the
//  Double implementation for evaluation, then convert the result back to the
//  requested type. Dedicated overloads for Double, Float, and (on x86) Float80
//  call through to type-specific C entry points for maximal precision.
//
//  Domain notes:
//  - asinh(x): defined for all real x
//  - acosh(x): defined for x >= 1; returns NaN for x < 1
//  - atanh(x): defined for |x| < 1; returns +/-infinity as x approaches +/-1,
//              and NaN for |x| > 1
//
//  Accuracy and behavior are determined by the underlying Boost implementations.
//  Propagation of NaN, infinities, and domain errors follows the semantics of
//  the C/Boost functions exposed by CBoostBridge.
//
import SwiftyBoostPrelude

public extension SpecialFunctions {
    /// Inverse hyperbolic tangent.
    ///
    /// Computes atanh(x) for any BinaryFloatingPoint by evaluating in Double precision
    /// and converting the result back to the generic type.
    ///
    /// - Parameter x: The input value. The mathematical domain is |x| < 1.
    /// - Returns: The inverse hyperbolic tangent of `x`.
    /// - Note: For |x| > 1 the result is NaN. As x approaches +/-1, the result tends to +/-infinity.
    @inlinable
    static func atanh<T: Real & BinaryFloatingPoint & Sendable>(_ x: T) -> T {
        return T(bs_atanh_d(Double(x)))
    }

    /// Inverse hyperbolic tangent for Double.
    ///
    /// - Parameter x: The input value. The mathematical domain is |x| < 1.
    /// - Returns: The inverse hyperbolic tangent of `x`.
    @inlinable
    static func atanh(_ x: Double) -> Double {
        return bs_atanh_d(x)
    }

    /// Inverse hyperbolic tangent for Float.
    ///
    /// - Parameter x: The input value. The mathematical domain is |x| < 1.
    /// - Returns: The inverse hyperbolic tangent of `x`.
    @inlinable
    static func atanh(_ x: Float) -> Float {
        return bs_atanh_f(x)
    }

    /// Inverse hyperbolic sine.
    ///
    /// Computes asinh(x) for any BinaryFloatingPoint by evaluating in Double precision
    /// and converting the result back to the generic type.
    ///
    /// - Parameter x: The input value (all real numbers are allowed).
    /// - Returns: The inverse hyperbolic sine of `x`.
    @inlinable
    static func asinh<T: Real & BinaryFloatingPoint & Sendable>(_ x: T) -> T {
        return T(bs_asinh_d(Double(x)))
    }

    /// Inverse hyperbolic sine for Double.
    ///
    /// - Parameter x: The input value (all real numbers are allowed).
    /// - Returns: The inverse hyperbolic sine of `x`.
    @inlinable
    static func asinh(_ x: Double) -> Double {
        return bs_asinh_d(x)
    }

    /// Inverse hyperbolic sine for Float.
    ///
    /// - Parameter x: The input value (all real numbers are allowed).
    /// - Returns: The inverse hyperbolic sine of `x`.
    @inlinable
    static func asinh(_ x: Float) -> Float {
        return bs_asinh_f(x)
    }

    /// Inverse hyperbolic cosine.
    ///
    /// Computes acosh(x) for any BinaryFloatingPoint by evaluating in Double precision
    /// and converting the result back to the generic type.
    ///
    /// - Parameter x: The input value. The mathematical domain is x >= 1.
    /// - Returns: The inverse hyperbolic cosine of `x`.
    /// - Note: For x < 1 the result is NaN.
    @inlinable
    static func acosh<T: Real & BinaryFloatingPoint & Sendable>(_ x: T) -> T {
        return T(bs_acosh_d(Double(x)))
    }

    /// Inverse hyperbolic cosine for Double.
    ///
    /// - Parameter x: The input value. The mathematical domain is x >= 1.
    /// - Returns: The inverse hyperbolic cosine of `x`.
    @inlinable
    static func acosh(_ x: Double) -> Double {
        return bs_acosh_d(x)
    }

    /// Inverse hyperbolic cosine for Float.
    ///
    /// - Parameter x: The input value. The mathematical domain is x >= 1.
    /// - Returns: The inverse hyperbolic cosine of `x`.
    @inlinable
    static func acosh(_ x: Float) -> Float {
        return bs_acosh_f(x)
    }

    #if arch(i386) || arch(x86_64)
    /// Inverse hyperbolic tangent for Float80 (x86 only).
    ///
    /// - Parameter x: The input value. The mathematical domain is |x| < 1.
    /// - Returns: The inverse hyperbolic tangent of `x`.
    @inlinable
    static func atanh(_ x: Float80) -> Float80 {
        return bs_atanh_l(x)
    }

    /// Inverse hyperbolic sine for Float80 (x86 only).
    ///
    /// - Parameter x: The input value (all real numbers are allowed).
    /// - Returns: The inverse hyperbolic sine of `x`.
    @inlinable
    static func asinh(_ x: Float80) -> Float80 {
        return bs_asinh_l(x)
    }

    /// Inverse hyperbolic cosine for Float80 (x86 only).
    ///
    /// - Parameter x: The input value. The mathematical domain is x >= 1.
    /// - Returns: The inverse hyperbolic cosine of `x`.
    @inlinable
    static func acosh(_ x: Float80) -> Float80 {
        return bs_acosh_l(x)
    }
    #endif
}
