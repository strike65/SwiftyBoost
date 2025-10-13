//
//  Legendre.swift
//  Math/SpecialFunctions
//
//  Legendre polynomials and associated Legendre functions (first kind).
//
//  This file provides Swift wrappers around Boost.Math’s real-valued Legendre
//  functions via the CBoostBridge C shim. The API surface mirrors the common
//  mathematical notation and offers:
//  - P_n(x): Legendre polynomial of non‑negative integer degree n.
//  - P_n^m(x): Associated Legendre function of the first kind with integer
//    degree n ≥ 0 and order m with |m| ≤ n.
//
//  Design goals:
//  - Type-generic entry points for any BinaryFloatingPoint (Double-backed).
//  - Type-specific overloads for Float and (on x86_64) Float80, avoiding
//    intermediate conversions for performance.
//  - Swift-side validation of common constraints (n ≥ 0, |m| ≤ n, finiteness).
//  - Consistent error model using SpecialFunctionError.
//  - Thorough DocC comments describing domains, behavior, and references.
//
//  References:
//  - NIST DLMF §14: https://dlmf.nist.gov/14
//  - Boost.Math Legendre documentation:
//    https://www.boost.org/doc/libs/release/libs/math/doc/html/math_toolkit/sf_legendre.html
//  - Orthogonal Polynomials (general background): G. Szegő, “Orthogonal Polynomials”.
//

import CBoostBridge
public extension SpecialFunctions {
    
    
    // MARK: - Generic BinaryFloatingPoint overloads (Double-backed)
    
    /// Legendre polynomial Pₙ(x) for integer n ≥ 0.
    ///
    /// Overview:
    /// - Computes the (unassociated) Legendre polynomial of degree `n`.
    /// - This generic overload accepts any `BinaryFloatingPoint` input `x` and returns
    ///   a result of the same generic type `T`.
    /// - Internally, the computation is delegated to Boost.Math via a `Double`-backed
    ///   C shim and then converted back to `T`.
    ///
    /// Domain and behavior:
    /// - Mathematically, Pₙ(x) is defined for all real `x`. In many applications,
    ///   `x ∈ [-1, 1]` due to orthogonality on that interval, but this wrapper does
    ///   not enforce that restriction. Very large |x| or large `n` can lead to large
    ///   magnitudes and loss of precision depending on `T`.
    /// - `n` must be a non‑negative integer (n ≥ 0). This is enforced by the wrapper.
    /// - `x` must be finite (not NaN or ±∞).
    ///
    /// Normalization and conventions:
    /// - This is the standard Legendre polynomial Pₙ(x) with P₀(x) = 1, P₁(x) = x,
    ///   and the standard recurrence Pₙ₊₁(x) = ((2n+1)x Pₙ(x) − n Pₙ₋₁(x)) / (n+1).
    ///
    /// Thread-safety:
    /// - Pure and thread-safe; no shared mutable state.
    ///
    /// Parameters:
    /// - n: The degree (non‑negative integer).
    /// - x: The evaluation point (finite real).
    ///
    /// Returns:
    /// - The value of Pₙ(x) as `T`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotPositive(name: "n")` if `n < 0`.
    /// - `SpecialFunctionError.parameterNotFinite(name: "x")` if `x` is NaN or ±∞.
    ///
    /// Example:
    /// ```swift
    /// let p3: Double = try legendreP(3, 0.5) // P3(0.5) = (5x^3 − 3x)/2 at x=0.5
    /// ```
    @inlinable static func legendreP<T: BinaryFloatingPoint>(_ n: Int, _ x: T) throws -> T {
        // Degree must be non-negative.
        guard n >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "n") }
        // Convert to Double for the C backend. Keep a single conversion point to minimize rounding steps.
        let dx = D(x)
        // Validate finiteness prior to calling the C layer.
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        // Delegate to Boost-backed implementation and convert back to T.
        return T(bs_legendre_p(Int32(n), dx))
    }
    
    /// Associated Legendre function of the first kind Pₙᵐ(x) for integer n ≥ 0 and |m| ≤ n.
    ///
    /// Overview:
    /// - Computes the associated Legendre function Pₙᵐ(x) with integer degree `n` and
    ///   order `m`. This generic overload accepts any `BinaryFloatingPoint` input `x`
    ///   and returns a result of the same generic type `T`.
    /// - Internally, the computation is delegated to Boost.Math via a `Double`-backed
    ///   C shim and then converted back to `T`.
    ///
    /// Domain and behavior:
    /// - Mathematically, Pₙᵐ(x) is defined for `n ≥ 0`, integer `m` with |m| ≤ n,
    ///   and real `x`. In many applications, `x ∈ [-1, 1]` due to orthogonality on
    ///   that interval, but this wrapper does not enforce that restriction.
    /// - This wrapper enforces:
    ///   - `n ≥ 0`
    ///   - `|m| ≤ n`
    ///   - `x` must be finite (not NaN or ±∞)
    ///
    /// Conventions:
    /// - Boost.Math’s associated Legendre uses the Condon–Shortley phase, i.e.,
    ///   Pₙᵐ(x) includes the factor (−1)ᵐ.
    /// - For negative `m`, the standard relation holds:
    ///   Pₙ^{−m}(x) = (−1)ᵐ (n−m)! / (n+m)! · Pₙᵐ(x)
    ///
    /// Thread-safety:
    /// - Pure and thread-safe; no shared mutable state.
    ///
    /// Parameters:
    /// - n: Degree (non‑negative integer).
    /// - m: Order (integer with |m| ≤ n).
    /// - x: Evaluation point (finite real).
    ///
    /// Returns:
    /// - The value of Pₙᵐ(x) as `T`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotPositive(name: "n")` if `n < 0`.
    /// - `SpecialFunctionError.parameterOutOfRange(name: "m", min: -n, max: n)` if `|m| > n`.
    /// - `SpecialFunctionError.parameterNotFinite(name: "x")` if `x` is NaN or ±∞.
    ///
    /// Example:
    /// ```swift
    /// let p32: Double = try associatedLegendreP(3, 2, 0.3) // P3^2(0.3)
    /// ```
    @inlinable static func associatedLegendreP<T: BinaryFloatingPoint>(_ n: Int, _ m: Int, _ x: T) throws -> T {
        // Degree/order constraints.
        guard n >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "n") }
        guard abs(m) <= n else { throw SpecialFunctionError.parameterOutOfRange(name: "m", min: Double(-n), max: Double(n)) }
        // Convert to Double for the C backend. Keep a single conversion point to minimize rounding steps.
        let dx = D(x)
        // Validate finiteness prior to calling the C layer.
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        // Delegate to Boost-backed implementation and convert back to T.
        return T(bs_assoc_legendre_p(Int32(n), Int32(m), dx))
    }
    
    /// Derivative of the Legendre polynomial, Pₙ′(x), for integer n ≥ 0.
    ///
    /// Parameters:
    /// - n: Degree (non‑negative integer).
    /// - x: Evaluation point (finite real).
    ///
    /// Returns:
    /// - Pₙ′(x) as `T`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotPositive(name: "n")` if `n < 0`.
    /// - `SpecialFunctionError.parameterNotFinite(name: "x")` if `x` is NaN or ±∞.
    @inlinable static func legendrePPrime<T: BinaryFloatingPoint>(_ n: Int, _ x: T) throws -> T {
        guard n >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "n") }
        let dx = D(x)
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return T(bs_legendre_p_prime(Int32(n), dx))
    }
    
    /// Zeros of the Legendre polynomial Pₗ(x), i.e. the l simple roots in (−1, 1) for l ≥ 1.
    ///
    /// Behavior:
    /// - For l ≤ 0, returns an empty array.
    /// - For l ≥ 1, returns an array of length l with the roots in ascending order as provided by Boost.
    ///
    /// Parameters:
    /// - degree: l, the degree (non‑negative integer).
    ///
    /// Returns:
    /// - An array [x_i] of length l containing the zeros of Pₗ(x), converted to `T`.
    @inlinable static func legendrePZeros<T: BinaryFloatingPoint>(degree l: Int) throws -> [T] {
        guard l >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "degree") }
        guard l > 0 else { return [] }
        var tmp = Array<Double>(repeating: .zero, count: l)
        tmp.withUnsafeMutableBufferPointer { buf in
            bs_legendre_p_zeros(Int32(l), buf.baseAddress!)
        }
        return tmp.map(T.init)
    }
    
    // MARK: - Float overloads
    // These overloads call directly into the Float-precision C implementations for
    // performance and to avoid intermediate conversions.
    
    /// Pₙ(x) for `Float`. Requires `n ≥ 0` and finite `x`.
    ///
    /// Parameters:
    /// - n: Degree (non‑negative integer).
    /// - x: Evaluation point (finite `Float`).
    ///
    /// Returns:
    /// - `Pₙ(x)` as `Float`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotPositive(name: "n")` if `n < 0`.
    /// - `SpecialFunctionError.parameterNotFinite(name: "x")` if `x` is NaN or ±∞.
    ///
    /// Example:
    /// ```swift
    /// let pf = try legendreP(2, 0.1 as Float)
    /// ```
    @inlinable static func legendreP(_ n: Int, _ x: Float) throws -> Float {
        guard n >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "n") }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return bs_legendre_p_f(Int32(n), x)
    }
    
    /// Pₙᵐ(x) for `Float`. Requires `n ≥ 0`, `|m| ≤ n`, and finite `x`.
    ///
    /// Parameters:
    /// - n: Degree (non‑negative integer).
    /// - m: Order (integer with |m| ≤ n).
    /// - x: Evaluation point (finite `Float`).
    ///
    /// Returns:
    /// - `Pₙᵐ(x)` as `Float`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotPositive(name: "n")` if `n < 0`.
    /// - `SpecialFunctionError.parameterOutOfRange(name: "m", min: -n, max: n)` if `|m| > n`.
    /// - `SpecialFunctionError.parameterNotFinite(name: "x")` if `x` is NaN or ±∞.
    ///
    /// Example:
    /// ```swift
    /// let paf = try associatedLegendreP(4, 1, 0.25 as Float)
    /// ```
    @inlinable static func associatedLegendreP(_ n: Int, _ m: Int, _ x: Float) throws -> Float {
        guard n >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "n") }
        guard abs(m) <= n else { throw SpecialFunctionError.parameterOutOfRange(name: "m", min: Double(-n), max: Double(n)) }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return bs_assoc_legendre_p_f(Int32(n), Int32(m), x)
    }
    
    /// Pₙ′(x) for `Float`. Requires `n ≥ 0` and finite `x`.
    @inlinable static func legendrePPrime(_ n: Int, _ x: Float) throws -> Float {
        guard n >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "n") }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return bs_legendre_p_prime_f(Int32(n), x)
    }
    
    /// Zeros of Pₗ(x) for `Float`. For l ≤ 0 returns [].
    @inlinable static func legendrePZeros(degree l: Int) throws -> [Float] {
        guard l >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "degree") }
        guard l > 0 else { return [] }
        var out = Array<Float>(repeating: .zero, count: l)
        out.withUnsafeMutableBufferPointer { buf in
            bs_legendre_p_zeros_f(Int32(l), buf.baseAddress!)
        }
        return out
    }
    
    // MARK: - Float80 overloads (x86_64)
    // Extended-precision entry points for platforms that support Float80.
    
#if arch(x86_64)
    /// Pₙ(x) for `Float80` (x86_64 only). Requires `n ≥ 0` and finite `x`.
    ///
    /// Prefer this overload when you need more precision than `Double` on x86_64.
    ///
    /// Parameters:
    /// - n: Degree (non‑negative integer).
    /// - x: Evaluation point (finite `Float80`).
    ///
    /// Returns:
    /// - `Pₙ(x)` as `Float80`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotPositive(name: "n")` if `n < 0`.
    /// - `SpecialFunctionError.parameterNotFinite(name: "x")` if `x` is NaN or ±∞.
    ///
    /// Availability:
    /// - Only on x86_64 architectures where `Float80` is available.
    @inlinable static func legendreP(_ n: Int, _ x: Float80) throws -> Float80 {
        guard n >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "n") }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return bs_legendre_p_l(Int32(n), x)
    }
    
    /// Pₙᵐ(x) for `Float80` (x86_64 only). Requires `n ≥ 0`, `|m| ≤ n`, and finite `x`.
    ///
    /// Prefer this overload when you need more precision than `Double` on x86_64.
    ///
    /// Parameters:
    /// - n: Degree (non‑negative integer).
    /// - m: Order (integer with |m| ≤ n).
    /// - x: Evaluation point (finite `Float80`).
    ///
    /// Returns:
    /// - `Pₙᵐ(x)` as `Float80`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotPositive(name: "n")` if `n < 0`.
    /// - `SpecialFunctionError.parameterOutOfRange(name: "m", min: -n, max: n)` if `|m| > n`.
    /// - `SpecialFunctionError.parameterNotFinite(name: "x")` if `x` is NaN or ±∞.
    ///
    /// Availability:
    /// - Only on x86_64 architectures where `Float80` is available.
    @inlinable static func associatedLegendreP(_ n: Int, _ m: Int, _ x: Float80) throws -> Float80 {
        guard n >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "n") }
        guard abs(m) <= n else { throw SpecialFunctionError.parameterOutOfRange(name: "m", min: Double(-n), max: Double(n)) }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return bs_assoc_legendre_p_l(Int32(n), Int32(m), x)
    }
    
    /// Pₙ′(x) for `Float80` (x86_64 only). Requires `n ≥ 0` and finite `x`.
    @inlinable static func legendrePPrime(_ n: Int, _ x: Float80) throws -> Float80 {
        guard n >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "n") }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return bs_legendre_p_prime_l(Int32(n), x)
    }
    
    /// Zeros of Pₗ(x) for `Float80` (x86_64 only). For l ≤ 0 returns [].
    @inlinable static func legendrePZeros(degree l: Int) throws -> [Float80] {
        guard l >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "degree") }
        guard l > 0 else { return [] }
        var out = Array<Float80>(repeating: .zero, count: l)
        out.withUnsafeMutableBufferPointer { buf in
            bs_legendre_p_zeros_l(Int32(l), buf.baseAddress!)
        }
        return out
    }
#endif
}

