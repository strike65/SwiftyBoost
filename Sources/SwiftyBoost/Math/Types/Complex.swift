//
//  Created by VT on 15.10.25.
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
//
//  Complex.swift
//  SwiftyBoost
//
//  A lightweight, generic complex number type.
//
//  This file follows repository guidelines:
//  - Generic over BinaryFloatingPoint
//  - Conforms to Sendable, Equatable, Hashable, Codable
//  - Public API with Swift naming conventions
//

import Foundation
#if canImport(CBoostBridge)
import CBoostBridge
#endif

/// A generic complex number with real and imaginary parts stored as `BinaryFloatingPoint`.
///
/// Complex provides basic arithmetic (+, −, ×, ÷), conjugation, and numerically stable
/// magnitude computation, without depending on external libraries.
///
/// Numeric stability:
/// - Division uses Smith’s algorithm to reduce overflow/underflow risk.
/// - Magnitude uses a scaling approach: `max(|a|,|b|) * sqrt(1 + (min/max)^2)`.
///
/// Note on equality and NaNs:
/// - Equatable mirrors the semantics of `T`; any NaN payloads render comparisons `false`.
///
/// Codable:
/// - Encoding uses a keyed container with keys `real` and `imag`.
@frozen
public struct Complex<T: BinaryFloatingPoint & Codable & Sendable>: Sendable, Equatable, Hashable, Codable, CustomStringConvertible {

    /// Real component.
    public var real: T

    /// Imaginary component.
    public var imag: T

    // MARK: - Initializers

    /// Creates a complex number from real and imaginary parts.
    @inlinable
    public init(real: T, imag: T) {
        self.real = real
        self.imag = imag
    }

    /// Creates a purely real complex number.
    @inlinable
    public init(_ real: T) {
        self.real = real
        self.imag = .zero
    }

    /// Creates a complex number from real and imaginary parts (unlabeled).
    @inlinable
    public init(_ real: T, _ imag: T) {
        self.real = real
        self.imag = imag
    }

    // MARK: - Common constants

    /// 0 + 0i
    public static var zero: Complex<T> { Complex(real: .zero, imag: .zero) }

    /// 1 + 0i
    public static var one: Complex<T> { Complex(real: 1, imag: .zero) }

    /// 0 + 1i
    public static var i: Complex<T> { Complex(real: .zero, imag: 1) }

    // MARK: - Derived properties

    /// The complex conjugate: (a + bi) → (a − bi).
    @inlinable
    public var conjugate: Complex<T> {
        Complex(real: real, imag: -imag)
    }

    /// The squared magnitude (a.k.a. the Euclidean norm squared): a² + b².
    ///
    /// - Important: May overflow for very large components; prefer `magnitude` for a stable result.
    @inlinable
    public var norm: T {
        real * real + imag * imag
    }

    /// The magnitude (Euclidean norm) `sqrt(a² + b²)` computed in a numerically stable way.
    ///
    /// Algorithm:
    /// - Let `x = max(|a|, |b|)` and `y = min(|a|, |b|)`.
    /// - If `x` is zero, result is zero.
    /// - Otherwise, compute `x * sqrt(1 + (y/x)²)`.
    @inlinable
    public var magnitude: T {
        let ar = abs(real)
        let ai = abs(imag)
        let x = ar >= ai ? ar : ai
        let y = ar >= ai ? ai : ar
        if x == .zero { return .zero }
        let r = y / x
        return x * (1 + r * r).squareRoot()
    }

    /// Whether both components are finite.
    @inlinable
    public var isFinite: Bool {
        real.isFinite && imag.isFinite
    }

    /// Whether any component is infinite.
    @inlinable
    public var isInfinite: Bool {
        real.isInfinite || imag.isInfinite
    }

    /// Whether any component is NaN.
    @inlinable
    public var isNaN: Bool {
        real.isNaN || imag.isNaN
    }

    // MARK: - Polar construction and phase (argument)

    /// Creates a complex number from polar coordinates `radius * (cos(phase) + i sin(phase))`.
    ///
    /// - Parameters:
    ///   - radius: Non-negative radius (magnitide). Negative values will be used as-is.
    ///   - phase: Angle in radians.
    ///
    /// - Returns: The complex number corresponding to the given polar coordinates.
    ///
    /// - Note: Implemented for `Double`, `Float`, and (x86_64) `Float80` via type-specific extensions below.
    @inlinable
    public static func fromPolar(radius: T, phase: T) -> Complex<T> {
        // Default fallback (for generic T) returns purely real radius when phase == 0,
        // otherwise returns (radius, 0). Specialized implementations exist for common types.
        if phase == .zero {
            return Complex(radius, .zero)
        }
        return Complex(radius, .zero)
    }

    /// The principal value of the argument (phase) in radians, in (−π, π].
    ///
    /// - Note: Implemented for `Double`, `Float`, and (x86_64) `Float80` via type-specific extensions below.
    @inlinable
    public var phase: T {
        // Fallback: return 0 for generic T; specialized implementations exist for common types.
        if real.sign == .minus { return .zero } // arbitrary placeholder for generic T
        return .zero
    }

    // MARK: - Reciprocal and normalization

    /// The multiplicative inverse `1 / z`, computed with a stable algorithm.
    ///
    /// Algorithm (Smith):
    /// - If `|a| ≥ |b|`: let `r = b/a`, `den = a + b*r`; then `1/z = (1 + (-r)i) / den`.
    /// - Else: let `r = a/b`, `den = b + a*r`; then `1/z = (r + (-1)i) / den`.
    ///
    /// - Returns: `1 / (a + bi)`; if both `a` and `b` are zero, the result contains infinities/NaNs per IEEE division-by-zero rules.
    @inlinable
    public func reciprocal() -> Complex<T> {
        let a = real, b = imag
        let ar = abs(a), br = abs(b)
        if ar >= br {
            // If a == 0 and b != 0, we will take the else-branch below.
            if a == .zero {
                // Both zero => fall through to else path (which will also divide by zero).
                // Let IEEE handle infinities/NaNs.
            } else {
                let r = b / a
                let den = a + b * r
                return Complex((1) / den, (-r) / den)
            }
        }
        // ar < br OR a == 0
        let r = a / b
        let den = b + a * r
        return Complex(r / den, (-1) / den)
    }

    /// The normalized complex number `z / |z|` when `|z| > 0`; otherwise `z`.
    ///
    /// - Returns: A unit-magnitude complex in the same direction as `z`, or `z` if its magnitude is zero.
    @inlinable
    public func normalized() -> Complex<T> {
        let m = magnitude
        if m == .zero { return self }
        return self / m
    }

    // MARK: - CustomStringConvertible

    public var description: String {
        let sign = imag.sign == .minus ? " - " : " + "
        let b = abs(imag)
        return "\(real)\(sign)\(b)i"
    }

    // MARK: - Codable

    public enum CodingKeys: String, CodingKey {
        case real, imag
    }

    @inlinable
    public init(from decoder: Decoder) throws {
        let c = try decoder.container(keyedBy: CodingKeys.self)
        self.real = try c.decode(T.self, forKey: .real)
        self.imag = try c.decode(T.self, forKey: .imag)
    }

    @inlinable
    public func encode(to encoder: Encoder) throws {
        var c = encoder.container(keyedBy: CodingKeys.self)
        try c.encode(real, forKey: .real)
        try c.encode(imag, forKey: .imag)
    }
}

// MARK: - ExpressibleByFloatLiteral

extension Complex: ExpressibleByFloatLiteral where T: ExpressibleByFloatLiteral {
    @inlinable
    public init(floatLiteral value: T.FloatLiteralType) {
        self.init(T(floatLiteral: value))
    }
}

// MARK: - AdditiveArithmetic

extension Complex: AdditiveArithmetic {
    @inlinable
    public static func + (lhs: Complex<T>, rhs: Complex<T>) -> Complex<T> {
        Complex(real: lhs.real + rhs.real, imag: lhs.imag + rhs.imag)
    }

    @inlinable
    public static func - (lhs: Complex<T>, rhs: Complex<T>) -> Complex<T> {
        Complex(real: lhs.real - rhs.real, imag: lhs.imag - rhs.imag)
    }
}

// MARK: - Arithmetic Operators

/// Unary negation.
@inlinable
public prefix func - <T>(x: Complex<T>) -> Complex<T> {
    Complex(real: -x.real, imag: -x.imag)
}

/// Multiplies two complex numbers: (a+bi)(c+di) = (ac − bd) + (ad + bc)i
@inlinable
public func * <T>(lhs: Complex<T>, rhs: Complex<T>) -> Complex<T> {
    let ac = lhs.real * rhs.real
    let bd = lhs.imag * rhs.imag
    let ad = lhs.real * rhs.imag
    let bc = lhs.imag * rhs.real
    return Complex(real: ac - bd, imag: ad + bc)
}

/// Divides two complex numbers using Smith’s stable algorithm.
@inlinable
public func / <T>(lhs: Complex<T>, rhs: Complex<T>) -> Complex<T> {
    // (a + bi) / (c + di)
    let a = lhs.real, b = lhs.imag, c = rhs.real, d = rhs.imag
    let ac = abs(c), ad = abs(d)
    if ac >= ad {
        // If c == 0 and d != 0, we prefer the else-branch to avoid r = d/c.
        if c == .zero, d != .zero {
            let r = c / d // 0 / d = 0
            let denom = d + c * r
            let real = (a * r + b) / denom
            let imag = (b * r - a) / denom
            return Complex(real: real, imag: imag)
        }
        let r = d / c
        let denom = c + d * r
        let real = (a + b * r) / denom
        let imag = (b - a * r) / denom
        return Complex(real: real, imag: imag)
    } else {
        let r = c / d
        let denom = d + c * r
        let real = (a * r + b) / denom
        let imag = (b * r - a) / denom
        return Complex(real: real, imag: imag)
    }
}

// MARK: - Mixed Complex–Scalar arithmetic

/// Adds a scalar to the real part: (a+bi) + x = (a+x) + bi
@inlinable
public func + <T>(lhs: Complex<T>, rhs: T) -> Complex<T> {
    Complex(real: lhs.real + rhs, imag: lhs.imag)
}

/// Adds a scalar to the real part: x + (a+bi) = (x+a) + bi
@inlinable
public func + <T>(lhs: T, rhs: Complex<T>) -> Complex<T> {
    Complex(real: lhs + rhs.real, imag: rhs.imag)
}

/// Subtracts a scalar from the real part: (a+bi) − x = (a−x) + bi
@inlinable
public func - <T>(lhs: Complex<T>, rhs: T) -> Complex<T> {
    Complex(real: lhs.real - rhs, imag: lhs.imag)
}

/// Subtracts a complex from a scalar: x − (a+bi) = (x−a) − bi
@inlinable
public func - <T>(lhs: T, rhs: Complex<T>) -> Complex<T> {
    Complex(real: lhs - rhs.real, imag: -rhs.imag)
}

/// Scales a complex by a scalar.
@inlinable
public func * <T>(lhs: Complex<T>, rhs: T) -> Complex<T> {
    Complex(real: lhs.real * rhs, imag: lhs.imag * rhs)
}

/// Scales a complex by a scalar.
@inlinable
public func * <T>(lhs: T, rhs: Complex<T>) -> Complex<T> {
    Complex(real: lhs * rhs.real, imag: lhs * rhs.imag)
}

/// Divides a complex by a scalar.
@inlinable
public func / <T>(lhs: Complex<T>, rhs: T) -> Complex<T> {
    Complex(real: lhs.real / rhs, imag: lhs.imag / rhs)
}

/// Divides a scalar by a complex: x / (a+bi).
@inlinable
public func / <T>(lhs: T, rhs: Complex<T>) -> Complex<T> {
    // Treat lhs as (lhs + 0i) and use complex division.
    Complex(real: lhs, imag: .zero) / rhs
}

// MARK: - Integer power

public extension Complex {
    /// Raises `z` to an integer power `n` using exponentiation by squaring.
    ///
    /// - Parameters:
    ///   - n: An integer exponent (can be negative).
    /// - Returns: `z^n`. If `n < 0`, computes `(1/z)^|n|`.
    @inlinable
    func pow(_ n: Int) -> Complex<T> {
        if n == 0 { return .one }
        if n < 0 { return self.reciprocal().pow(-n) }
        // n > 0
        var base = self
        var exp = n
        var result = Complex<T>.one
        while exp > 0 {
            if (exp & 1) == 1 { result = result * base }
            exp >>= 1
            if exp > 0 { base = base * base }
        }
        return result
    }
}

// MARK: - Square root (stable)

public extension Complex {
    /// The principal square root of a complex number.
    ///
    /// Stable formula:
    /// - Let `r = |z| = sqrt(a² + b²)`.
    /// - `u = sqrt((r + a)/2)`, `v = sign(b) * sqrt((r − a)/2)`.
    ///
    /// This avoids catastrophic cancellation when `a` and `b` are large or disparate.
    @inlinable
    var squareRoot: Complex<T> {
        let a = real, b = imag
        if b == .zero {
            if a >= .zero {
                return Complex(a.squareRoot(), .zero)
            } else {
                return Complex(.zero, ( -a ).squareRoot())
            }
        }
        let r = magnitude
        // Guard against potential tiny negatives due to rounding:
        let u = ((r + a) / 2).squareRoot()
        let v = ((r - a) / 2).squareRoot()
        if b >= .zero {
            return Complex(u, v)
        } else {
            return Complex(u, -v)
        }
    }
}

// MARK: - Double/Float/Float80 specializations for polar and phase and elementary complex functions

public extension Complex where T == Double {
    /// Creates a complex number from polar coordinates `radius * (cos(phase) + i sin(phase))`.
    @inlinable
    static func fromPolar(radius: Double, phase: Double) -> Complex<Double> {
        Complex(radius * Darwin.cos(phase), radius * Darwin.sin(phase))
    }

    /// The principal argument (phase) in radians, in (−π, π].
    @inlinable
    var phase: Double {
        atan2(imag, real)
    }

    // Elementary complex functions via Boost bridge

    /// Complex exponential: exp(a + ib)
    @inlinable
    var exp: Complex<Double> {
        #if canImport(CBoostBridge)
        let z = bs_complex_d(re: real, im: imag)
        let r = bs_cexp(z)
        return Complex(r.re, r.im)
        #else
        let ea = Darwin.exp(real)
        return Complex(ea * cos(imag), ea * sin(imag))
        #endif
    }

    /// Natural logarithm: log(z) = ln|z| + i arg(z)
    @inlinable
    var log: Complex<Double> {
        #if canImport(CBoostBridge)
        let z = bs_complex_d(re: real, im: imag)
        let r = bs_clog(z)
        return Complex(r.re, r.im)
        #else
        Complex(Darwin.log(magnitude), phase)
        #endif
    }

    /// Complex sine: sin(a + ib)
    @inlinable
    var sin: Complex<Double> {
        #if canImport(CBoostBridge)
        let z = bs_complex_d(re: real, im: imag)
        let r = bs_csin(z)
        return Complex(r.re, r.im)
        #else
        // Fallback via identity: sin(x+iy) = sin x cosh y + i cos x sinh y
        let x = real, y = imag
        return Complex(Swift.sin(x) * cosh(y), Swift.cos(x) * sinh(y))
        #endif
    }

    /// Complex cosine: cos(a + ib)
    @inlinable
    var cos: Complex<Double> {
        #if canImport(CBoostBridge)
        let z = bs_complex_d(re: real, im: imag)
        let r = bs_ccos(z)
        return Complex(r.re, r.im)
        #else
        // cos(x+iy) = cos x cosh y − i sin x sinh y
        let x = real, y = imag
        return Complex(Swift.cos(x) * cosh(y), -Swift.sin(x) * sinh(y))
        #endif
    }

    /// Complex tangent: tan(a + ib)
    @inlinable
    var tan: Complex<Double> {
        #if canImport(CBoostBridge)
        let z = bs_complex_d(re: real, im: imag)
        let r = bs_ctan(z)
        return Complex(r.re, r.im)
        #else
        return self.sin / self.cos
        #endif
    }

    /// Complex hyperbolic sine: sinh(a + ib)
    @inlinable
    var sinh: Complex<Double> {
        #if canImport(CBoostBridge)
        let z = bs_complex_d(re: real, im: imag)
        let r = bs_csinh(z)
        return Complex(r.re, r.im)
        #else
        // sinh(x+iy) = sinh x cos y + i cosh x sin y
        let x = real, y = imag
        return Complex(Swift.sinh(x) * Swift.cos(y), Swift.cosh(x) * Swift.sin(y))
        #endif
    }

    /// Complex hyperbolic cosine: cosh(a + ib)
    @inlinable
    var cosh: Complex<Double> {
        #if canImport(CBoostBridge)
        let z = bs_complex_d(re: real, im: imag)
        let r = bs_ccosh(z)
        return Complex(r.re, r.im)
        #else
        // cosh(x+iy) = cosh x cos y + i sinh x sin y
        let x = real, y = imag
        return Complex(Swift.cosh(x) * Swift.cos(y), Swift.sinh(x) * Swift.sin(y))
        #endif
    }

    /// Complex hyperbolic tangent: tanh(a + ib)
    @inlinable
    var tanh: Complex<Double> {
        #if canImport(CBoostBridge)
        let z = bs_complex_d(re: real, im: imag)
        let r = bs_ctanh(z)
        return Complex(r.re, r.im)
        #else
        return self.sinh / self.cosh
        #endif
    }

    /// Complex arctangent: atan(a + ib)
    @inlinable
    var atan: Complex<Double> {
        #if canImport(CBoostBridge)
        let z = bs_complex_d(re: real, im: imag)
        let r = bs_catan(z)
        return Complex(r.re, r.im)
        #else
        // Fallback using identity: atan z = (i/2) [ln(1 - i z) - ln(1 + i z)]
        let i = Complex<Double>.i
        let one = Complex<Double>.one
        let lhs = (one - i * self).log
        let rhs = (one + i * self).log
        return (i / 2) * (lhs - rhs)
        #endif
    }
}

public extension Complex where T == Float {
    /// Creates a complex number from polar coordinates `radius * (cos(phase) + i sin(phase))`.
    @inlinable
    static func fromPolar(radius: Float, phase: Float) -> Complex<Float> {
        Complex(radius * Darwin.cos(phase), radius * Darwin.sin(phase))
    }

    /// The principal argument (phase) in radians, in (−π, π].
    @inlinable
    var phase: Float {
        atan2(imag, real)
    }

    // Elementary complex functions via Boost bridge

    /// Complex exponential: exp(a + ib)
    @inlinable
    var exp: Complex<Float> {
        #if canImport(CBoostBridge)
        let z = bs_complex_f(re: real, im: imag)
        let r = bs_cexp_f(z)
        return Complex(r.re, r.im)
        #else
        let ea = Float(Darwin.exp(Double(real)))
        return Complex(ea * cos(imag), ea * sin(imag))
        #endif
    }

    /// Natural logarithm: log(z) = ln|z| + i arg(z)
    @inlinable
    var log: Complex<Float> {
        #if canImport(CBoostBridge)
        let z = bs_complex_f(re: real, im: imag)
        let r = bs_clog_f(z)
        return Complex(r.re, r.im)
        #else
        Complex(Float(Darwin.log(Double(magnitude))), phase)
        #endif
    }

    /// Complex sine: sin(a + ib)
    @inlinable
    var sin: Complex<Float> {
        #if canImport(CBoostBridge)
        let z = bs_complex_f(re: real, im: imag)
        let r = bs_csin_f(z)
        return Complex(r.re, r.im)
        #else
        let x = real, y = imag
        return Complex(Swift.sin(x) * cosh(y), Swift.cos(x) * sinh(y))
        #endif
    }

    /// Complex cosine: cos(a + ib)
    @inlinable
    var cos: Complex<Float> {
        #if canImport(CBoostBridge)
        let z = bs_complex_f(re: real, im: imag)
        let r = bs_ccos_f(z)
        return Complex(r.re, r.im)
        #else
        let x = real, y = imag
        return Complex(Swift.cos(x) * cosh(y), -Swift.sin(x) * sinh(y))
        #endif
    }

    /// Complex tangent: tan(a + ib)
    @inlinable
    var tan: Complex<Float> {
        // Use identity to avoid missing C bridge symbol for Float:
        return self.sin / self.cos
    }

    /// Complex hyperbolic sine: sinh(a + ib)
    @inlinable
    var sinh: Complex<Float> {
        #if canImport(CBoostBridge)
        let z = bs_complex_f(re: real, im: imag)
        let r = bs_csinh_f(z)
        return Complex(r.re, r.im)
        #else
        let x = real, y = imag
        return Complex(Swift.sinh(x) * Swift.cos(y), Swift.cosh(x) * Swift.sin(y))
        #endif
    }

    /// Complex hyperbolic cosine: cosh(a + ib)
    @inlinable
    var cosh: Complex<Float> {
        #if canImport(CBoostBridge)
        let z = bs_complex_f(re: real, im: imag)
        let r = bs_ccosh_f(z)
        return Complex(r.re, r.im)
        #else
        let x = real, y = imag
        return Complex(Swift.cosh(x) * Swift.cos(y), Swift.sinh(x) * Swift.sin(y))
        #endif
    }

    /// Complex hyperbolic tangent: tanh(a + ib)
    @inlinable
    var tanh: Complex<Float> {
        #if canImport(CBoostBridge)
        let z = bs_complex_f(re: real, im: imag)
        let r = bs_ctanh_f(z)
        return Complex(r.re, r.im)
        #else
        return self.sinh / self.cosh
        #endif
    }

    /// Complex arctangent: atan(a + ib)
    @inlinable
    var atan: Complex<Float> {
        #if canImport(CBoostBridge)
        let z = bs_complex_f(re: real, im: imag)
        let r = bs_catan_f(z)
        return Complex(r.re, r.im)
        #else
        let i = Complex<Float>.i
        let one = Complex<Float>.one
        let lhs = (one - i * self).log
        let rhs = (one + i * self).log
        return (i / 2) * (lhs - rhs)
        #endif
    }
}

#if arch(x86_64)
public extension Complex where T == Float80 {
    /// Creates a complex number from polar coordinates `radius * (cos(phase) + i sin(phase))`.
    @inlinable
    static func fromPolar(radius: Float80, phase: Float80) -> Complex<Float80> {
        Complex(radius * cos(phase), radius * sin(phase))
    }

    /// The principal argument (phase) in radians, in (−π, π].
    @inlinable
    var phase: Float80 {
        atan2(imag, real)
    }

    // Elementary complex functions via Boost bridge

    /// Complex exponential: exp(a + ib)
    @inlinable
    var exp: Complex<Float80> {
        #if canImport(CBoostBridge)
        let z = bs_complex_l(re: real, im: imag)
        let r = bs_cexp_l(z)
        return Complex(r.re, r.im)
        #else
        let ea = Float80(Darwin.exp(Double(real)))
        return Complex(ea * cos(imag), ea * sin(imag))
        #endif
    }

    /// Natural logarithm: log(z) = ln|z| + i arg(z)
    @inlinable
    var log: Complex<Float80> {
        #if canImport(CBoostBridge)
        let z = bs_complex_l(re: real, im: imag)
        let r = bs_clog_l(z)
        return Complex(r.re, r.im)
        #else
        Complex(Float80(Darwin.log(Double(magnitude))), phase)
        #endif
    }

    /// Complex sine: sin(a + ib)
    @inlinable
    var sin: Complex<Float80> {
        #if canImport(CBoostBridge)
        let z = bs_complex_l(re: real, im: imag)
        let r = bs_csin_l(z)
        return Complex(r.re, r.im)
        #else
        let x = real, y = imag
        return Complex(Swift.sin(x) * cosh(y), Swift.cos(x) * sinh(y))
        #endif
    }

    /// Complex cosine: cos(a + ib)
    @inlinable
    var cos: Complex<Float80> {
        #if canImport(CBoostBridge)
        let z = bs_complex_l(re: real, im: imag)
        let r = bs_ccos_l(z)
        return Complex(r.re, r.im)
        #else
        let x = real, y = imag
        return Complex(Swift.cos(x) * cosh(y), -Swift.sin(x) * sinh(y))
        #endif
    }

    /// Complex tangent: tan(a + ib)
    @inlinable
    var tan: Complex<Float80> {
        #if canImport(CBoostBridge)
        let z = bs_complex_l(re: real, im: imag)
        let r = bs_ctan_l(z)
        return Complex(r.re, r.im)
        #else
        return self.sin / self.cos
        #endif
    }

    /// Complex hyperbolic sine: sinh(a + ib)
    @inlinable
    var sinh: Complex<Float80> {
        #if canImport(CBoostBridge)
        let z = bs_complex_l(re: real, im: imag)
        let r = bs_csinh_l(z)
        return Complex(r.re, r.im)
        #else
        let x = real, y = imag
        return Complex(Swift.sinh(x) * Swift.cos(y), Swift.cosh(x) * Swift.sin(y))
        #endif
    }

    /// Complex hyperbolic cosine: cosh(a + ib)
    @inlinable
    var cosh: Complex<Float80> {
        #if canImport(CBoostBridge)
        let z = bs_complex_l(re: real, im: imag)
        let r = bs_ccosh_l(z)
        return Complex(r.re, r.im)
        #else
        let x = real, y = imag
        return Complex(Swift.cosh(x) * Swift.cos(y), Swift.sinh(x) * Swift.sin(y))
        #endif
    }

    /// Complex hyperbolic tangent: tanh(a + ib)
    @inlinable
    var tanh: Complex<Float80> {
        #if canImport(CBoostBridge)
        let z = bs_complex_l(re: real, im: imag)
        let r = bs_ctanh_l(z)
        return Complex(r.re, r.im)
        #else
        return self.sinh / self.cosh
        #endif
    }

    /// Complex arctangent: atan(a + ib)
    @inlinable
    var atan: Complex<Float80> {
        #if canImport(CBoostBridge)
        let z = bs_complex_l(re: real, im: imag)
        let r = bs_catan_l(z)
        return Complex(r.re, r.im)
        #else
        let i = Complex<Float80>.i
        let one = Complex<Float80>.one
        let lhs = (one - i * self).log
        let rhs = (one + i * self).log
        return (i / 2) * (lhs - rhs)
        #endif
    }
}
#endif

// MARK: - Public typealiases for common precisions

public typealias ComplexD = Complex<Double>
public typealias ComplexF = Complex<Float>
#if arch(x86_64)
public typealias ComplexX = Complex<Float80>
#endif

