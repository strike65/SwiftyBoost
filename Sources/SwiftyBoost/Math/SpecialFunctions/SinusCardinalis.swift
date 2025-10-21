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
//  SinusCardinalis.swift
//  ------------------------------------
//  Normalized cardinal sine and hyperbolic cardinal sine utilities.
//
//  This file provides real and complex variants of the normalized cardinal
//  sine and its hyperbolic counterpart in π-normalized form:
//
//    - sinc(πx)  = sin(πx)  / (πx) with removable singularity sinc(0) = 1
//    - sinhc(πx) = sinh(πx) / (πx) with removable singularity sinhc(0) = 1
//
//  Overloads are provided for Float, Double and, on Intel architectures,
//  Float80. Complex overloads operate on ComplexF, ComplexD, and ComplexL.
//
//  The implementations delegate to highly-accurate routines bridged from
//  CBoostBridge (Boost special functions), while preserving Swift-friendly
//  APIs and semantics.
//

import CBoostBridge

public extension SpecialFunctions {
    
    /// Normalized cardinal sine, π-scaled, for real arguments.
    ///
    /// Computes
    ///   sinc(πx) = sin(πx) / (πx),
    /// with the removable singularity at `x = 0` defined by the limit `1`.
    ///
    /// - Parameter x: The real input value.
    /// - Returns: The value of `sin(πx) / (πx)`, or `1` when `x == 0`.
    /// - Notes:
    ///   - This generic overload supports `Double`, `Float`, and, on Intel architectures,
    ///     `Float80`. The computation is performed in `Double` precision and converted
    ///     back to `T`.
    ///   - Propagates IEEE-754 `NaN` and `±∞` according to the underlying implementation.
    /// - SeeAlso: ``SpecialFunctions/sinc_pi(_:)->T`` (generic), ``SpecialFunctions/sinc_pi(_:)->Double`` (Double), ``SpecialFunctions/sinc_pi(_:)->Float`` (Float)
    @inlinable
    static func sinc_pi<T: BinaryFloatingPoint>(_ x: T) -> T {
            return T(bs_sinc_pi_d(Double(x)))
    }

    /// Normalized cardinal sine, π-scaled, specialized for `Double`.
    ///
    /// - Parameter x: The real input value.
    /// - Returns: `sin(πx) / (πx)`, or `1` when `x == 0`.
    @inlinable
    static func sinc_pi(_ x: Double) -> Double {
            return bs_sinc_pi_d(x)
    }

    /// Normalized cardinal sine, π-scaled, specialized for `Float`.
    ///
    /// - Parameter x: The real input value.
    /// - Returns: `sin(πx) / (πx)`, or `1` when `x == 0`.
    @inlinable
    static func sinc_pi(_ x: Float) -> Float {
            return bs_sinc_pi_f(x)
    }

    /// Normalized hyperbolic cardinal sine, π-scaled, for real arguments.
    ///
    /// Computes
    ///   sinhc(πx) = sinh(πx) / (πx),
    /// with the removable singularity at `x = 0` defined by the limit `1`.
    ///
    /// - Parameter x: The real input value.
    /// - Returns: The value of `sinh(πx) / (πx)`, or `1` when `x == 0`.
    /// - Notes:
    ///   - This generic overload supports `Double`, `Float`, and, on Intel architectures,
    ///     `Float80`. The computation is performed in `Double` precision and converted
    ///     back to `T`.
    ///   - Propagates IEEE-754 `NaN` and `±∞` according to the underlying implementation.
    /// - SeeAlso: ``SpecialFunctions/sinhc_pi(_:)->T`` (generic), ``SpecialFunctions/sinhc_pi(_:)->Double`` (Double), ``SpecialFunctions/sinhc_pi(_:)->Float`` (Float)
    @inlinable
    static func sinhc_pi<T: BinaryFloatingPoint>(_ x: T) -> T {
            return T(bs_sinhc_pi_d(Double(x)))
    }

    /// Normalized hyperbolic cardinal sine, π-scaled, specialized for `Double`.
    ///
    /// - Parameter x: The real input value.
    /// - Returns: `sinh(πx) / (πx)`, or `1` when `x == 0`.
    @inlinable
    static func sinhc_pi(_ x: Double) -> Double {
            return bs_sinhc_pi_d(x)
    }

    /// Normalized hyperbolic cardinal sine, π-scaled, specialized for `Float`.
    ///
    /// - Parameter x: The real input value.
    /// - Returns: `sinh(πx) / (πx)`, or `1` when `x == 0`.
    @inlinable
    static func sinhc_pi(_ x: Float) -> Float {
            return bs_sinhc_pi_f(x)
    }

    // MARK: - Complex cardinal sine sinc(πz)

    /// Complex normalized cardinal sine: `sinc(πz) = sin(πz) / (πz)` for double-precision complex values.
    ///
    /// The removable singularity at `z = 0` is defined by the limit `1 + 0i`.
    ///
    /// - Parameter z: Complex input as `ComplexD` (double-precision).
    /// - Returns: The complex value of `sin(πz) / (πz)`, or `1 + 0i` when `z == 0`.
    /// - Note: Delegates to a Boost-based implementation via `CBoostBridge`.
    @inlinable
    static func sincc_pi(_ z: ComplexD) -> ComplexD {
        let r = bs_sincc_pi_d(complex_d(re: z.real, im: z.imag))
        return ComplexD(r.re, r.im)
    }

    /// Complex normalized cardinal sine: `sinc(πz) = sin(πz) / (πz)` for single-precision complex values.
    ///
    /// The removable singularity at `z = 0` is defined by the limit `1 + 0i`.
    ///
    /// - Parameter z: Complex input as `ComplexF` (single-precision).
    /// - Returns: The complex value of `sin(πz) / (πz)`, or `1 + 0i` when `z == 0`.
    @inlinable
    static func sincc_pi(_ z: ComplexF) -> ComplexF {
        let r = bs_sincc_pi_f(complex_f(re: z.real, im: z.imag))
        return ComplexF(r.re, r.im)
    }

    /// Complex normalized hyperbolic cardinal sine: `sinhc(πz) = sinh(πz) / (πz)` for double-precision complex values.
    ///
    /// The removable singularity at `z = 0` is defined by the limit `1 + 0i`.
    ///
    /// - Parameter z: Complex input as `ComplexD` (double-precision).
    /// - Returns: The complex value of `sinh(πz) / (πz)`, or `1 + 0i` when `z == 0`.
    @inlinable
    static func sinhcc_pi(_ z: ComplexD) -> ComplexD {
        let r = bs_sinhcc_pi_d(complex_d(re: z.real, im: z.imag))
        return ComplexD(r.re, r.im)
    }

    /// Complex normalized hyperbolic cardinal sine: `sinhc(πz) = sinh(πz) / (πz)` for single-precision complex values.
    ///
    /// The removable singularity at `z = 0` is defined by the limit `1 + 0i`.
    ///
    /// - Parameter z: Complex input as `ComplexF` (single-precision).
    /// - Returns: The complex value of `sinh(πz) / (πz)`, or `1 + 0i` when `z == 0`.
    @inlinable
    static func sinhcc_pi(_ z: ComplexF) -> ComplexF {
        let r = bs_sinhcc_pi_f(complex_f(re: z.real, im: z.imag))
        return ComplexF(r.re, r.im)
    }


    #if arch(i386) || arch(x86_64)
    /// Normalized cardinal sine, π-scaled, specialized for `Float80` (Intel only).
    ///
    /// - Parameter x: The real input value.
    /// - Returns: `sin(πx) / (πx)`, or `1` when `x == 0`.
    @inlinable
    static func sinc_pi(_ x: Float80) -> Float80 {
            return bs_sinc_pi_l(x)
    }

    /// Normalized hyperbolic cardinal sine, π-scaled, specialized for `Float80` (Intel only).
    ///
    /// - Parameter x: The real input value.
    /// - Returns: `sinh(πx) / (πx)`, or `1` when `x == 0`.
    @inlinable
    static func sinhc_pi(_ x: Float80) -> Float80 {
            return bs_sinhc_pi_l(x)
    }

    /// Complex normalized cardinal sine: `sinc(πz)` for extended-precision complex values (Intel only).
    ///
    /// - Parameter z: Complex input as `ComplexL` (extended precision).
    /// - Returns: The complex value of `sin(πz) / (πz)`, or `1 + 0i` when `z == 0`.
    /// - Availability: Only on `i386` and `x86_64` architectures.
    @inlinable
    static func sincc_pi(_ z: ComplexL) -> ComplexL {
        let r = bs_sincc_pi_l(complex_l(re: z.real, im: z.imag))
        return ComplexL(r.re, r.im)
    }

    /// Complex normalized hyperbolic cardinal sine: `sinhc(πz)` for extended-precision complex values (Intel only).
    ///
    /// - Parameter z: Complex input as `ComplexL` (extended precision).
    /// - Returns: The complex value of `sinh(πz) / (πz)`, or `1 + 0i` when `z == 0`.
    /// - Availability: Only on `i386` and `x86_64` architectures.
    @inlinable
    static func sinhcc_pi(_ z: ComplexL) -> ComplexL {
        let r = bs_sinhcc_pi_l(complex_l(re: z.real, im: z.imag))
        return ComplexL(r.re, r.im)
    }
    #endif
}
