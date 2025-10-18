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
//  Swift wrappers for the (generalized) Laguerre polynomials L_n(x) and L_n^m(x).
//  These APIs validate common real-domain constraints and delegate numerical work
//  to CBoostBridge (wrapping Boost.Math).
//
//  Functions (real-valued):
//  - laguerre(n, x):        L_n(x)          — Laguerre polynomial of degree n.
//  - assocLaguerre(n, m, x): L_n^m(x)       — Associated (generalized) Laguerre polynomial.
//
//  Domain (real-valued use):
//  - Indices n (and m for the associated form) must be integers with n ≥ 0, m ≥ 0.
//  - The argument x must be a finite real.
//
//  Overloads:
//  - Generic BinaryFloatingPoint overloads convert to Double using D(_:) and
//    return results as the same generic type T.
//  - Type-specific overloads for Float and (on x86_64) Float80 call directly
//    into the matching C functions for performance.
//
//  Error model:
//  - n < 0 (or m < 0) throws SpecialFunctionError.parameterNotPositive.
//  - Non-finite x throws SpecialFunctionError.parameterNotFinite.
//
//  References:
//  - NIST DLMF §18 (Orthogonal Polynomials): https://dlmf.nist.gov/18
//  - Boost.Math Laguerre: https://www.boost.org/doc/libs/release/libs/math/doc/html/math_toolkit/sf_laguerre/laguerre.html
//

import CBoostBridge

public extension SpecialFunctions {
    // MARK: - Generic overloads (Double-backed)

    /// Compute the Laguerre polynomial L_n(x).
    ///
    /// Definition (one of several equivalent forms):
    /// - L_n(x) = Σ_{k=0}^n (n choose k) [(-x)^k / k!]
    ///
    /// Domain:
    /// - Requires an integer degree n ≥ 0 and finite real x.
    ///
    /// Parameters:
    /// - n: Degree of the polynomial (n ≥ 0).
    /// - x: Evaluation point (finite real).
    ///
    /// Returns:
    /// - L_n(x) as `T`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotPositive(name: "n")` if `n < 0`.
    /// - `SpecialFunctionError.parameterNotFinite(name: "x")` if `x` is NaN or ±∞.
    @inlinable static func laguerre<T: BinaryFloatingPoint>(_ n: Int, _ x: T) throws -> T {
        guard n >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "n") }
        let dx = D(x)
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return T(bs_laguerre(UInt32(n), dx))
    }

    /// Compute the associated (generalized) Laguerre polynomial L_n^m(x).
    ///
    /// Common definition:
    /// - L_n^m(x) = Σ_{k=0}^n (n + m choose n − k) [(-x)^k / k!], for integers n ≥ 0, m ≥ 0.
    ///
    /// Domain:
    /// - Requires integer indices n ≥ 0, m ≥ 0 and finite real x.
    ///
    /// Parameters:
    /// - n: Degree (n ≥ 0).
    /// - m: Order (m ≥ 0).
    /// - x: Evaluation point (finite real).
    ///
    /// Returns:
    /// - L_n^m(x) as `T`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotPositive(name: "n")` if `n < 0`.
    /// - `SpecialFunctionError.parameterNotPositive(name: "m")` if `m < 0`.
    /// - `SpecialFunctionError.parameterNotFinite(name: "x")` if `x` is NaN or ±∞.
    @inlinable static func assocLaguerre<T: BinaryFloatingPoint>(_ n: Int, _ m: Int, _ x: T) throws -> T {
        guard n >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "n") }
        guard m >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "m") }
        let dx = D(x)
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return T(bs_assoc_laguerre(UInt32(n), UInt32(m), dx))
    }

    // MARK: - Float overloads

    /// L_n(x) for `Float`.
    ///
    /// Throws if `n < 0` or `x` is not finite.
    @inlinable static func laguerre(_ n: Int, _ x: Float) throws -> Float {
        guard n >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "n") }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return bs_laguerre_f(UInt32(n), x)
    }

    /// L_n^m(x) for `Float`.
    ///
    /// Throws if `n < 0`, `m < 0`, or `x` is not finite.
    @inlinable static func assocLaguerre(_ n: Int, _ m: Int, _ x: Float) throws -> Float {
        guard n >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "n") }
        guard m >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "m") }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return bs_assoc_laguerre_f(UInt32(n), UInt32(m), x)
    }

    // MARK: - Float80 overloads (x86_64 only)

    #if arch(x86_64)
    /// L_n(x) for `Float80` (x86_64 only).
    ///
    /// Throws if `n < 0` or `x` is not finite.
    @inlinable static func laguerre(_ n: Int, _ x: Float80) throws -> Float80 {
        guard n >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "n") }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return bs_laguerre_l(UInt32(n), x)
    }

    /// L_n^m(x) for `Float80` (x86_64 only).
    ///
    /// Throws if `n < 0`, `m < 0`, or `x` is not finite.
    @inlinable static func assocLaguerre(_ n: Int, _ m: Int, _ x: Float80) throws -> Float80 {
        guard n >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "n") }
        guard m >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "m") }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return bs_assoc_laguerre_l(UInt32(n), UInt32(m), x)
    }
    #endif
}
