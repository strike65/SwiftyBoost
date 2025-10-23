//
//  Created by Volker Thieme 2025.
//  Copyright © 2025 Volker Thieme.
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

import SwiftyBoostPrelude

public extension SpecialFunctions {
    
    
    // MARK: - 1F0(a; ; z)
    
    /// Gauss-type hypergeometric function 1F0(a; ; z).
    ///
    /// Definition:
    /// - 1F0(a; ; z) = Σ_{k=0}^∞ (a)_k z^k / k! = (1 − z)^{−a} for |z| < 1, with analytic continuation elsewhere.
    ///
    /// Domain and behavior:
    /// - The series diverges for |z| ≥ 1, but the function is defined via analytic continuation.
    /// - This wrapper checks only for finiteness of inputs. Poles/branch behavior are handled by the backend.
    ///
    /// Parameters:
    /// - a: Upper parameter (finite real).
    /// - z: Argument (finite real).
    ///
    /// Returns:
    /// - 1F0(a; ; z) as `T`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotFinite(name: "a")` if `a` is NaN or ±∞.
    /// - `SpecialFunctionError.parameterNotFinite(name: "z")` if `z` is NaN or ±∞.
    ///
    /// References:
    /// - NIST DLMF §16.2 (Generalized Hypergeometric Function)
    /// - Boost.Math hypergeometric_1F0
    @inlinable static func hypergeometric1F0<T: Real & BinaryFloatingPoint & Sendable>(a: T, z: T) throws -> T {
        let da = D(a), dz = D(z)
        guard da.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a", value: a) }
        guard dz.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "z", value: z) }
        return T(bs_hypergeometric_1F0(da, dz))
    }

    /// 1F0(a; ; z) with mixed `Float`/`Double` arguments; returns `Double`.
    ///
    /// See ``SpecialFunctions/hypergeometric1F0(a:z:)->T`` for definition and domain notes.
    @inlinable static func hypergeometric1F0(a: Float, z: Double) throws -> Double { try hypergeometric1F0(a: Double(a), z: z) }
    /// 1F0(a; ; z) with mixed `Double`/`Float` arguments; returns `Double`.
    ///
    /// See ``SpecialFunctions/hypergeometric1F0(a:z:)->T`` for definition and domain notes.
    @inlinable static func hypergeometric1F0(a: Double, z: Float) throws -> Double { try hypergeometric1F0(a: a, z: Double(z)) }
    
    /// 1F0(a; ; z) for `Float`.
    @inlinable static func hypergeometric1F0(a: Float, z: Float) throws -> Float {
        guard a.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a", value: a) }
        guard z.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "z", value: z) }
        return bs_hypergeometric_1F0_f(a, z)
    }
    
#if arch(x86_64)
    /// 1F0(a; ; z) for `Float80` (x86_64 only).
    @inlinable static func hypergeometric1F0(a: Float80, z: Float80) throws -> Float80 {
        guard a.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a", value: a) }
        guard z.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "z", value: z) }
        return bs_hypergeometric_1F0_l(a, z)
    }

    /// 1F0(a; ; z) with one `Float80` mixed argument; returns `Float80` (x86_64).
    @inlinable static func hypergeometric1F0(a: Float80, z: Double) throws -> Float80 { try hypergeometric1F0(a: a, z: Float80(z)) }
    /// 1F0(a; ; z) with one `Float80` mixed argument; returns `Float80` (x86_64).
    @inlinable static func hypergeometric1F0(a: Double, z: Float80) throws -> Float80 { try hypergeometric1F0(a: Float80(a), z: z) }
    /// 1F0(a; ; z) with one `Float80` mixed argument; returns `Float80` (x86_64).
    @inlinable static func hypergeometric1F0(a: Float80, z: Float) throws -> Float80 { try hypergeometric1F0(a: a, z: Float80(z)) }
    /// 1F0(a; ; z) with one `Float80` mixed argument; returns `Float80` (x86_64).
    @inlinable static func hypergeometric1F0(a: Float, z: Float80) throws -> Float80 { try hypergeometric1F0(a: Float80(a), z: z) }
#endif
    
    
    // MARK: - 0F1(; b; z)
    
    /// Confluent hypergeometric limit function 0F1(; b; z).
    ///
    /// Definition:
    /// - 0F1(; b; z) = Σ_{k=0}^∞ z^k / [(b)_k k!], entire in z for b ∉ {0, −1, −2, …}.
    ///
    /// Domain and behavior:
    /// - Poles occur at b ∈ {0, −1, −2, …}. This wrapper does not pre-validate those discrete parameter poles.
    /// - Inputs must be finite. Numerical stability depends on parameter ranges; Boost.Math handles evaluation strategy.
    ///
    /// Parameters:
    /// - b: Lower parameter (finite real; avoid non‑positive integers to prevent poles).
    /// - z: Argument (finite real).
    ///
    /// Returns:
    /// - 0F1(; b; z) as `T`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotFinite(name: "b")` if `b` is NaN or ±∞.
    /// - `SpecialFunctionError.parameterNotFinite(name: "z")` if `z` is NaN or ±∞.
    ///
    /// References:
    /// - NIST DLMF §13 (Bessel connections), §16.2
    /// - Boost.Math hypergeometric_0F1
    @inlinable static func hypergeometric0F1<T: Real & BinaryFloatingPoint & Sendable>(b: T, z: T) throws -> T {
        let db = D(b), dz = D(z)
        guard db.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "b", value: b) }
        guard dz.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "z", value: z) }
        return T(bs_hypergeometric_0F1(db, dz))
    }

    /// 0F1(; b; z) with mixed `Float`/`Double` arguments; returns `Double`.
    ///
    /// See ``SpecialFunctions/hypergeometric0F1(b:z:)->T`` for definition and domain notes.
    @inlinable static func hypergeometric0F1(b: Float, z: Double) throws -> Double { try hypergeometric0F1(b: Double(b), z: z) }
    /// 0F1(; b; z) with mixed `Double`/`Float` arguments; returns `Double`.
    ///
    /// See ``SpecialFunctions/hypergeometric0F1(b:z:)->T`` for definition and domain notes.
    @inlinable static func hypergeometric0F1(b: Double, z: Float) throws -> Double { try hypergeometric0F1(b: b, z: Double(z)) }
    
    /// 0F1(; b; z) for `Float`.
    @inlinable static func hypergeometric0F1(b: Float, z: Float) throws -> Float {
        guard b.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "b", value: b) }
        guard z.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "z", value: z) }
        return bs_hypergeometric_0F1_f(b, z)
    }
    
#if arch(x86_64)
    /// 0F1(; b; z) for `Float80` (x86_64 only).
    @inlinable static func hypergeometric0F1(b: Float80, z: Float80) throws -> Float80 {
        guard b.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "b", value: b) }
        guard z.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "z", value: z) }
        return bs_hypergeometric_0F1_l(b, z)
    }

    /// 0F1(; b; z) with one `Float80` mixed argument; returns `Float80` (x86_64).
    @inlinable static func hypergeometric0F1(b: Float80, z: Double) throws -> Float80 { try hypergeometric0F1(b: b, z: Float80(z)) }
    /// 0F1(; b; z) with one `Float80` mixed argument; returns `Float80` (x86_64).
    @inlinable static func hypergeometric0F1(b: Double, z: Float80) throws -> Float80 { try hypergeometric0F1(b: Float80(b), z: z) }
    /// 0F1(; b; z) with one `Float80` mixed argument; returns `Float80` (x86_64).
    @inlinable static func hypergeometric0F1(b: Float80, z: Float) throws -> Float80 { try hypergeometric0F1(b: b, z: Float80(z)) }
    /// 0F1(; b; z) with one `Float80` mixed argument; returns `Float80` (x86_64).
    @inlinable static func hypergeometric0F1(b: Float, z: Float80) throws -> Float80 { try hypergeometric0F1(b: Float80(b), z: z) }
#endif
    
    
    // MARK: - 2F0(a, b; ; z)
    
    /// Hypergeometric 2F0(a, b; ; z).
    ///
    /// Series form and interpretation:
    /// - 2F0(a, b; ; z) = Σ_{k=0}^∞ (a)_k (b)_k z^k / k!.
    /// - The series diverges for any z ≠ 0. In practice, 2F0 is interpreted as an asymptotic expansion or via analytic continuation.
    ///
    /// Parameters:
    /// - a: Upper parameter (finite real).
    /// - b: Upper parameter (finite real).
    /// - z: Argument (finite real).
    ///
    /// Returns:
    /// - 2F0(a, b; ; z) as `T`, when defined by the backend.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotFinite(name: ...)` if any input is NaN or ±∞.
    ///
    /// Notes:
    /// - Boost.Math provides an evaluation consistent with its documentation; outside supported regions results may be NaN/Inf.
    /// - No attempt is made here to pre-screen convergence regions.
    ///
    /// References:
    /// - NIST DLMF §16.11
    /// - Boost.Math hypergeometric_2F0
    @inlinable static func hypergeometric2F0<T: Real & BinaryFloatingPoint & Sendable>(a: T, b: T, z: T) throws -> T {
        let da = D(a), db = D(b), dz = D(z)
        guard da.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a", value: a) }
        guard db.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "b", value: b) }
        guard dz.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "z", value: z) }
        return T(bs_hypergeometric_2F0(da, db, dz))
    }

    /// 2F0(a, b; z) with a single `Float` mixed argument; returns `Double`.
    ///
    /// See ``SpecialFunctions/hypergeometric2F0(a:b:z:)->T`` for behavior and domain notes.
    @inlinable static func hypergeometric2F0(a: Float, b: Double, z: Double) throws -> Double { try hypergeometric2F0(a: Double(a), b: b, z: z) }
    /// 2F0(a, b; z) with a single `Float` mixed argument; returns `Double`.
    @inlinable static func hypergeometric2F0(a: Double, b: Float, z: Double) throws -> Double { try hypergeometric2F0(a: a, b: Double(b), z: z) }
    /// 2F0(a, b; z) with a single `Float` mixed argument; returns `Double`.
    @inlinable static func hypergeometric2F0(a: Double, b: Double, z: Float) throws -> Double { try hypergeometric2F0(a: a, b: b, z: Double(z)) }
    
    /// 2F0(a, b; ; z) for `Float`.
    @inlinable static func hypergeometric2F0(a: Float, b: Float, z: Float) throws -> Float {
        guard a.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a", value: a) }
        guard b.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "b", value: b) }
        guard z.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "z", value: z) }
        return bs_hypergeometric_2F0_f(a, b, z)
    }
    
#if arch(x86_64)
    /// 2F0(a, b; ; z) for `Float80` (x86_64 only).
    @inlinable static func hypergeometric2F0(a: Float80, b: Float80, z: Float80) throws -> Float80 {
        guard a.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a", value: a) }
        guard b.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "b", value: b) }
        guard z.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "z", value: z) }
        return bs_hypergeometric_2F0_l(a, b, z)
    }

    /// 2F0(a, b; z) with one `Float80` mixed argument; returns `Float80` (x86_64).
    @inlinable static func hypergeometric2F0(a: Float80, b: Double, z: Double) throws -> Float80 { try hypergeometric2F0(a: a, b: Float80(b), z: Float80(z)) }
    /// 2F0(a, b; z) with one `Float80` mixed argument; returns `Float80` (x86_64).
    @inlinable static func hypergeometric2F0(a: Double, b: Float80, z: Double) throws -> Float80 { try hypergeometric2F0(a: Float80(a), b: b, z: Float80(z)) }
    /// 2F0(a, b; z) with one `Float80` mixed argument; returns `Float80` (x86_64).
    @inlinable static func hypergeometric2F0(a: Double, b: Double, z: Float80) throws -> Float80 { try hypergeometric2F0(a: Float80(a), b: Float80(b), z: z) }
    /// 2F0(a, b; z) with one `Float80` mixed argument; returns `Float80` (x86_64).
    @inlinable static func hypergeometric2F0(a: Float80, b: Float, z: Float) throws -> Float80 { try hypergeometric2F0(a: a, b: Float80(b), z: Float80(z)) }
    /// 2F0(a, b; z) with one `Float80` mixed argument; returns `Float80` (x86_64).
    @inlinable static func hypergeometric2F0(a: Float, b: Float80, z: Float) throws -> Float80 { try hypergeometric2F0(a: Float80(a), b: b, z: Float80(z)) }
    /// 2F0(a, b; z) with one `Float80` mixed argument; returns `Float80` (x86_64).
    @inlinable static func hypergeometric2F0(a: Float, b: Float, z: Float80) throws -> Float80 { try hypergeometric2F0(a: Float80(a), b: Float80(b), z: z) }
#endif
    
    
    // MARK: - 1F1(a; b; z) (Kummer’s M)
    
    /// Confluent hypergeometric function of the first kind 1F1(a; b; z) (Kummer’s M).
    ///
    /// Definition:
    /// - 1F1(a; b; z) = Σ_{k=0}^∞ (a)_k z^k / [(b)_k k!].
    ///
    /// Domain and behavior:
    /// - Poles occur at b ∈ {0, −1, −2, …}. This wrapper does not pre-validate those discrete parameter poles.
    /// - Entire in z for admissible b. Evaluation strategy depends on parameter ranges; Boost.Math selects stable forms.
    ///
    /// Parameters:
    /// - a: Upper parameter (finite real).
    /// - b: Lower parameter (finite real; avoid non‑positive integers).
    /// - z: Argument (finite real).
    ///
    /// Returns:
    /// - 1F1(a; b; z) as `T`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotFinite(name: ...)` if any input is NaN or ±∞.
    ///
    /// References:
    /// - NIST DLMF §13, §16.2
    /// - Boost.Math hypergeometric_1F1
    @inlinable static func hypergeometric1F1<T: Real & BinaryFloatingPoint & Sendable>(a: T, b: T, z: T) throws -> T {
        let da = D(a), db = D(b), dz = D(z)
        guard da.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a", value: a) }
        guard db.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "b", value: b) }
        guard dz.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "z", value: z) }
        return T(bs_hypergeometric_1F1(da, db, dz))
    }

    /// 1F1(a; b; z) with a single `Float` mixed argument; returns `Double`.
    ///
    /// See ``SpecialFunctions/hypergeometric1F1(a:b:z:)->T`` for definition and domain notes.
    @inlinable static func hypergeometric1F1(a: Float, b: Double, z: Double) throws -> Double { try hypergeometric1F1(a: Double(a), b: b, z: z) }
    /// 1F1(a; b; z) with a single `Float` mixed argument; returns `Double`.
    @inlinable static func hypergeometric1F1(a: Double, b: Float, z: Double) throws -> Double { try hypergeometric1F1(a: a, b: Double(b), z: z) }
    /// 1F1(a; b; z) with a single `Float` mixed argument; returns `Double`.
    @inlinable static func hypergeometric1F1(a: Double, b: Double, z: Float) throws -> Double { try hypergeometric1F1(a: a, b: b, z: Double(z)) }
    
    /// 1F1(a; b; z) for `Float`.
    @inlinable static func hypergeometric1F1(a: Float, b: Float, z: Float) throws -> Float {
        guard a.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a", value: a) }
        guard b.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "b", value: b) }
        guard z.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "z", value: z) }
        return bs_hypergeometric_1F1_f(a, b, z)
    }
    
#if arch(x86_64)
    /// 1F1(a; b; z) for `Float80` (x86_64 only).
    @inlinable static func hypergeometric1F1(a: Float80, b: Float80, z: Float80) throws -> Float80 {
        guard a.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a", value: a) }
        guard b.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "b", value: b) }
        guard z.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "z", value: z) }
        return bs_hypergeometric_1F1_l(a, b, z)
    }

    /// 1F1(a; b; z) with one `Float80` mixed argument; returns `Float80` (x86_64).
    @inlinable static func hypergeometric1F1(a: Float80, b: Double, z: Double) throws -> Float80 { try hypergeometric1F1(a: a, b: Float80(b), z: Float80(z)) }
    /// 1F1(a; b; z) with one `Float80` mixed argument; returns `Float80` (x86_64).
    @inlinable static func hypergeometric1F1(a: Double, b: Float80, z: Double) throws -> Float80 { try hypergeometric1F1(a: Float80(a), b: b, z: Float80(z)) }
    /// 1F1(a; b; z) with one `Float80` mixed argument; returns `Float80` (x86_64).
    @inlinable static func hypergeometric1F1(a: Double, b: Double, z: Float80) throws -> Float80 { try hypergeometric1F1(a: Float80(a), b: Float80(b), z: z) }
    /// 1F1(a; b; z) with one `Float80` mixed argument; returns `Float80` (x86_64).
    @inlinable static func hypergeometric1F1(a: Float80, b: Float, z: Float) throws -> Float80 { try hypergeometric1F1(a: a, b: Float80(b), z: Float80(z)) }
    /// 1F1(a; b; z) with one `Float80` mixed argument; returns `Float80` (x86_64).
    @inlinable static func hypergeometric1F1(a: Float, b: Float80, z: Float) throws -> Float80 { try hypergeometric1F1(a: Float80(a), b: b, z: Float80(z)) }
    /// 1F1(a; b; z) with one `Float80` mixed argument; returns `Float80` (x86_64).
    @inlinable static func hypergeometric1F1(a: Float, b: Float, z: Float80) throws -> Float80 { try hypergeometric1F1(a: Float80(a), b: Float80(b), z: z) }
#endif
    
    
    // MARK: - Generalized hypergeometric pFq(a; b; z)
    
    /// Generalized hypergeometric function pFq(a; b; z) with Double parameters.
    ///
    /// Definition:
    /// - pFq(a; b; z) = Σ_{k=0}^∞ [(a₁)_k … (a_p)_k] / [(b₁)_k … (b_q)_k] · z^k / k!
    ///
    /// Convergence:
    /// - If p ≤ q, the series converges for all finite z (entire when p < q).
    /// - If p = q + 1, converges for |z| < 1 and conditionally on |z| = 1 depending on parameters.
    /// - If p > q + 1, the series diverges for all z ≠ 0; interpreted as an asymptotic expansion/analytic continuation.
    ///
    /// Domain and notes:
    /// - Poles occur if any b_j ∈ {0, −1, −2, …}. This wrapper does not pre-screen discrete poles.
    /// - Inputs must be finite. Very large parameter magnitudes may be challenging; Boost.Math handles selection of methods.
    ///
    /// Parameters:
    /// - a: Array of upper parameters a₁…a_p.
    /// - b: Array of lower parameters b₁…b_q.
    /// - z: Argument.
    ///
    /// Returns:
    /// - pFq(a; b; z) as `Double`.
    @inlinable static func hypergeometricPFQ(a: [Double], b: [Double], z: Double) -> Double {
        a.withUnsafeBufferPointer { ap in
            b.withUnsafeBufferPointer { bp in
                bs_hypergeometric_pFq_d(ap.baseAddress, ap.count, bp.baseAddress, bp.count, z)
            }
        }
    }

    /// Generalized hypergeometric function pFq(a; b; z) with Float parameters.
    @inlinable static func hypergeometricPFQ(a: [Float], b: [Float], z: Float) -> Float {
        a.withUnsafeBufferPointer { ap in
            b.withUnsafeBufferPointer { bp in
                bs_hypergeometric_pFq_f(ap.baseAddress, ap.count, bp.baseAddress, bp.count, z)
            }
        }
    }

    /// pFq(a; b; z) with `[Float]` parameters and `Double` z; returns `Double`.
    @inlinable static func hypergeometricPFQ(a: [Float], b: [Float], z: Double) -> Double {
        let ad = a.map(Double.init), bd = b.map(Double.init)
        return hypergeometricPFQ(a: ad, b: bd, z: z)
    }
    /// pFq(a; b; z) with `[Double]` parameters and `Float` z; returns `Double`.
    @inlinable static func hypergeometricPFQ(a: [Double], b: [Double], z: Float) -> Double {
        return hypergeometricPFQ(a: a, b: b, z: Double(z))
    }
    
#if arch(x86_64)
    /// Generalized hypergeometric function pFq(a; b; z) with Float80 parameters (x86_64 only).
    @inlinable static func hypergeometricPFQ(a: [Float80], b: [Float80], z: Float80) -> Float80 {
        a.withUnsafeBufferPointer { ap in
            b.withUnsafeBufferPointer { bp in
                bs_hypergeometric_pFq_l(ap.baseAddress, ap.count, bp.baseAddress, bp.count, z)
            }
        }
    }

    /// pFq(a; b; z) with `[Float80]` parameters and `Double` z; returns `Float80` (x86_64).
    @inlinable static func hypergeometricPFQ(a: [Float80], b: [Float80], z: Double) -> Float80 {
        a.withUnsafeBufferPointer { ap in
            b.withUnsafeBufferPointer { bp in
                bs_hypergeometric_pFq_l(ap.baseAddress, ap.count, bp.baseAddress, bp.count, Float80(z))
            }
        }
    }
    /// pFq(a; b; z) with `[Double]` parameters and `Float80` z; returns `Float80` (x86_64).
    @inlinable static func hypergeometricPFQ(a: [Double], b: [Double], z: Float80) -> Float80 {
        let al = a.map(Float80.init), bl = b.map(Float80.init)
        return hypergeometricPFQ(a: al, b: bl, z: z)
    }
    /// pFq(a; b; z) with `[Float]` parameters and `Float80` z; returns `Float80` (x86_64).
    @inlinable static func hypergeometricPFQ(a: [Float], b: [Float], z: Float80) -> Float80 {
        let al = a.map(Float80.init), bl = b.map(Float80.init)
        return hypergeometricPFQ(a: al, b: bl, z: z)
    }
#endif
}
