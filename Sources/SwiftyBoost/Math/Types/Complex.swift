//
//  Complex.swift
//  SwiftyBoost
//
//  Re-exports Swift Numerics' `Complex` type and layers SwiftyBoost-specific
//  conveniences plus Boost-backed elementary functions.
//

import SwiftyBoostPrelude

#if canImport(Darwin)
import Darwin
#endif

/// Alias Swift Numerics' ``ComplexModule/Complex`` so existing call sites can keep
/// referring to ``SwiftyBoost/Complex``.
public typealias Complex<T: Real & BinaryFloatingPoint> = ComplexModule.Complex<T>

/// Double-precision complex number.
public typealias ComplexD = Complex<Double>
/// Single-precision complex number.
public typealias ComplexF = Complex<Float>
#if arch(x86_64)
/// Extended-precision complex number (x86 only).
public typealias ComplexL = Complex<Float80>
#endif

// MARK: - General conveniences for Real floating-point scalars

public extension Complex where RealType: Real & BinaryFloatingPoint {

    /// Imaginary component (alias that mirrors the legacy implementation).
    @inlinable
    var imag: RealType {
        @inline(__always) get { rawStorage.y }
        @inline(__always) set {
            var storage = rawStorage
            storage.y = newValue
            rawStorage = storage
        }
    }

    /// 0 + 0i
    @inlinable
    static var zero: Self { Self(.zero, .zero) }

    /// 1 + 0i
    @inlinable
    static var one: Self { Self(1, .zero) }

    /// 0 + 1i
    @inlinable
    static var i: Self { Self(.zero, 1) }

    /// Complex conjugate `(a + bi) → (a − bi)`.
    @inlinable
    var conjugate: Self {
        let (a, b) = rawStorage
        return Self(a, -b)
    }

    /// Squared magnitude `a² + b²`.
    @inlinable
    var norm: RealType {
        let (a, b) = rawStorage
        return a * a + b * b
    }

    /// Euclidean magnitude computed via a numerically stable scaling trick.
    @inlinable
    var magnitude: RealType {
        let (a, b) = rawStorage
        let ar = abs(a)
        let ai = abs(b)
        let x = max(ar, ai)
        let y = min(ar, ai)
        if x == .zero { return .zero }
        let r = y / x
        return x * (1 + r * r).squareRoot()
    }

    /// Whether both components are finite.
    @inlinable
    var isFinite: Bool {
        let (a, b) = rawStorage
        return a.isFinite && b.isFinite
    }

    /// Whether any component is infinite.
    @inlinable
    var isInfinite: Bool {
        let (a, b) = rawStorage
        return a.isInfinite || b.isInfinite
    }

    /// Whether any component is NaN.
    @inlinable
    var isNaN: Bool {
        let (a, b) = rawStorage
        return a.isNaN || b.isNaN
    }

    /// Multiplicative inverse computed with Smith’s algorithm for stability.
    @inlinable
    func reciprocal() -> Self {
        let (a, b) = rawStorage
        let ar = abs(a)
        let br = abs(b)
        if ar >= br {
            if a == .zero {
                // Fall through to the |b| branch and let IEEE rules handle infinities.
            } else {
                let r = b / a
                let den = a + b * r
                return Self(1 / den, -r / den)
            }
        }
        let r = a / b
        let den = b + a * r
        return Self(r / den, -1 / den)
    }

    /// Normalised complex number `z / |z|` (returns `z` when `|z| == 0`).
    @inlinable
    func normalized() -> Self {
        let m = magnitude
        if m == .zero { return self }
        let (a, b) = rawStorage
        return Self(a / m, b / m)
    }

    /// Raise `z` to an integer power using exponentiation by squaring.
    @inlinable
    func pow(_ n: Int) -> Self {
        if n == 0 { return .one }
        if n < 0 { return reciprocal().pow(-n) }
        var base = self
        var exp = n
        var result = Self.one
        while exp > 0 {
            if (exp & 1) == 1 {
                result = result * base
            }
            exp >>= 1
            if exp > 0 {
                base = base * base
            }
        }
        return result
    }

    /// Principal square root of the complex number.
    @inlinable
    var squareRoot: Self {
        let (a, b) = rawStorage
        if b == .zero {
            if a >= .zero {
                return Self(a.squareRoot(), .zero)
            } else {
                return Self(.zero, (-a).squareRoot())
            }
        }
        let r = magnitude
        let u = ((r + a) / 2).squareRoot()
        let v = ((r - a) / 2).squareRoot()
        if b >= .zero {
            return Self(u, v)
        } else {
            return Self(u, -v)
        }
    }

    /// String formatting `a ± bi`, matching the legacy implementation for tests/docs.
    @inlinable
    var description: String {
        let (a, b) = rawStorage
        let sign = b.sign == .minus ? " - " : " + "
        let bi = abs(b)
        return "\(a)\(sign)\(bi)i"
    }
}

// MARK: - Mixed complex/scalar operators

@inlinable
public func + <T>(lhs: Complex<T>, rhs: T) -> Complex<T> where T: Real & BinaryFloatingPoint {
    let (a, b) = lhs.rawStorage
    return Complex(a + rhs, b)
}

@inlinable
public func + <T>(lhs: T, rhs: Complex<T>) -> Complex<T> where T: Real & BinaryFloatingPoint {
    let (a, b) = rhs.rawStorage
    return Complex(lhs + a, b)
}

@inlinable
public func - <T>(lhs: Complex<T>, rhs: T) -> Complex<T> where T: Real & BinaryFloatingPoint {
    let (a, b) = lhs.rawStorage
    return Complex(a - rhs, b)
}

@inlinable
public func - <T>(lhs: T, rhs: Complex<T>) -> Complex<T> where T: Real & BinaryFloatingPoint {
    let (a, b) = rhs.rawStorage
    return Complex(lhs - a, -b)
}

@inlinable
public func * <T>(lhs: Complex<T>, rhs: T) -> Complex<T> where T: Real & BinaryFloatingPoint {
    let (a, b) = lhs.rawStorage
    return Complex(a * rhs, b * rhs)
}

@inlinable
public func * <T>(lhs: T, rhs: Complex<T>) -> Complex<T> where T: Real & BinaryFloatingPoint {
    let (a, b) = rhs.rawStorage
    return Complex(lhs * a, lhs * b)
}

@inlinable
public func / <T>(lhs: Complex<T>, rhs: T) -> Complex<T> where T: Real & BinaryFloatingPoint {
    let (a, b) = lhs.rawStorage
    return Complex(a / rhs, b / rhs)
}

@inlinable
public func / <T>(lhs: T, rhs: Complex<T>) -> Complex<T> where T: Real & BinaryFloatingPoint {
    Complex(lhs, .zero) / rhs
}

/// Division of two complex numbers using Smith’s stable algorithm.
@inlinable
public func / <T>(lhs: Complex<T>, rhs: Complex<T>) -> Complex<T> where T: Real & BinaryFloatingPoint {
    let (a, b) = lhs.rawStorage
    let (c, d) = rhs.rawStorage
    let ac = abs(c)
    let ad = abs(d)
    if ac >= ad {
        if c == .zero, d != .zero {
            let r = c / d
            let denom = d + c * r
            let real = (a * r + b) / denom
            let imag = (b * r - a) / denom
            return Complex(real, imag)
        }
        let r = d / c
        let denom = c + d * r
        let real = (a + b * r) / denom
        let imag = (b - a * r) / denom
        return Complex(real, imag)
    } else {
        let r = c / d
        let denom = d + c * r
        let real = (a * r + b) / denom
        let imag = (b * r - a) / denom
        return Complex(real, imag)
    }
}

// MARK: - Double specialisation

public extension Complex where RealType == Double {

    /// Creates a complex number from polar coordinates.
    @inlinable
    static func fromPolar(radius: Double, phase: Double) -> Complex<Double> {
        Complex(radius * Darwin.cos(phase), radius * Darwin.sin(phase))
    }

    /// Principal argument (phase) in radians, in (−π, π].
    @inlinable
    var phase: Double { atan2(imag, real) }

    /// Complex exponential via Boost.Math.
    @inlinable
    var exp: Complex<Double> {
        let r = bs_cexp_d(complex_d(re: real, im: imag))
        return Complex(r.re, r.im)
    }

    /// Natural logarithm via Boost.Math.
    @inlinable
    var log: Complex<Double> {
        let r = bs_clog_d(complex_d(re: real, im: imag))
        return Complex(r.re, r.im)
    }

    /// Complex sine via Boost.Math.
    @inlinable
    var sin: Complex<Double> {
        let r = bs_csin_d(complex_d(re: real, im: imag))
        return Complex(r.re, r.im)
    }

    /// Complex cosine via Boost.Math.
    @inlinable
    var cos: Complex<Double> {
        let r = bs_ccos_d(complex_d(re: real, im: imag))
        return Complex(r.re, r.im)
    }

    /// Complex tangent via Boost.Math.
    @inlinable
    var tan: Complex<Double> {
        let r = bs_ctan_d(complex_d(re: real, im: imag))
        return Complex(r.re, r.im)
    }

    /// Complex hyperbolic sine via Boost.Math.
    @inlinable
    var sinh: Complex<Double> {
        let r = bs_csinh_d(complex_d(re: real, im: imag))
        return Complex(r.re, r.im)
    }

    /// Complex hyperbolic cosine via Boost.Math.
    @inlinable
    var cosh: Complex<Double> {
        let r = bs_ccosh_d(complex_d(re: real, im: imag))
        return Complex(r.re, r.im)
    }

    /// Complex hyperbolic tangent via Boost.Math.
    @inlinable
    var tanh: Complex<Double> {
        let r = bs_ctanh_d(complex_d(re: real, im: imag))
        return Complex(r.re, r.im)
    }

    /// Complex arctangent via Boost.Math.
    @inlinable
    var atan: Complex<Double> {
        let r = bs_catan_d(complex_d(re: real, im: imag))
        return Complex(r.re, r.im)
    }
}

// MARK: - Float specialisation

public extension Complex where RealType == Float {

    /// Creates a complex number from polar coordinates.
    @inlinable
    static func fromPolar(radius: Float, phase: Float) -> Complex<Float> {
        Complex(radius * Darwin.cosf(phase), radius * Darwin.sinf(phase))
    }

    /// Principal argument (phase) in radians, in (−π, π].
    @inlinable
    var phase: Float { atan2f(imag, real) }

    @inlinable
    var exp: Complex<Float> {
        let r = bs_cexp_f(complex_f(re: real, im: imag))
        return Complex(r.re, r.im)
    }

    @inlinable
    var log: Complex<Float> {
        let r = bs_clog_f(complex_f(re: real, im: imag))
        return Complex(r.re, r.im)
    }

    @inlinable
    var sin: Complex<Float> {
        let r = bs_csin_f(complex_f(re: real, im: imag))
        return Complex(r.re, r.im)
    }

    @inlinable
    var cos: Complex<Float> {
        let r = bs_ccos_f(complex_f(re: real, im: imag))
        return Complex(r.re, r.im)
    }

    @inlinable
    var tan: Complex<Float> {
        let r = bs_ctan_f(complex_f(re: real, im: imag))
        return Complex(r.re, r.im)
    }

    @inlinable
    var sinh: Complex<Float> {
        let r = bs_csinh_f(complex_f(re: real, im: imag))
        return Complex(r.re, r.im)
    }

    @inlinable
    var cosh: Complex<Float> {
        let r = bs_ccosh_f(complex_f(re: real, im: imag))
        return Complex(r.re, r.im)
    }

    @inlinable
    var tanh: Complex<Float> {
        let r = bs_ctanh_f(complex_f(re: real, im: imag))
        return Complex(r.re, r.im)
    }

    @inlinable
    var atan: Complex<Float> {
        let r = bs_catan_f(complex_f(re: real, im: imag))
        return Complex(r.re, r.im)
    }
}

#if arch(x86_64)
// MARK: - Float80 specialisation

public extension Complex where RealType == Float80 {

    /// Creates a complex number from polar coordinates.
    @inlinable
    static func fromPolar(radius: Float80, phase: Float80) -> Complex<Float80> {
        Complex(radius * Darwin.cosl(phase), radius * Darwin.sinl(phase))
    }

    /// Principal argument (phase) in radians, in (−π, π].
    @inlinable
    var phase: Float80 { atan2(imag, real) }

    @inlinable
    var exp: Complex<Float80> {
        let r = bs_cexp_l(complex_l(re: real, im: imag))
        return Complex(r.re, r.im)
    }

    @inlinable
    var log: Complex<Float80> {
        let r = bs_clog_l(complex_l(re: real, im: imag))
        return Complex(r.re, r.im)
    }

    @inlinable
    var sin: Complex<Float80> {
        let r = bs_csin_l(complex_l(re: real, im: imag))
        return Complex(r.re, r.im)
    }

    @inlinable
    var cos: Complex<Float80> {
        let r = bs_ccos_l(complex_l(re: real, im: imag))
        return Complex(r.re, r.im)
    }

    @inlinable
    var tan: Complex<Float80> {
        let r = bs_ctan_l(complex_l(re: real, im: imag))
        return Complex(r.re, r.im)
    }

    @inlinable
    var sinh: Complex<Float80> {
        let r = bs_csinh_l(complex_l(re: real, im: imag))
        return Complex(r.re, r.im)
    }

    @inlinable
    var cosh: Complex<Float80> {
        let r = bs_ccosh_l(complex_l(re: real, im: imag))
        return Complex(r.re, r.im)
    }

    @inlinable
    var tanh: Complex<Float80> {
        let r = bs_ctanh_l(complex_l(re: real, im: imag))
        return Complex(r.re, r.im)
    }

    @inlinable
    var atan: Complex<Float80> {
        let r = bs_catan_l(complex_l(re: real, im: imag))
        return Complex(r.re, r.im)
    }
}
#endif
