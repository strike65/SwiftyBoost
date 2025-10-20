//
//  Created by Volker Thieme 2025.
//  Copyright © 2025 Volker Thieme.
//  License: MIT (see project root)
//
//  Gegenbauer (ultraspherical) polynomials Cₙ^(λ)(x) and their derivatives.
//
//  This file provides Swift wrappers around Boost.Math’s real-valued Gegenbauer
//  polynomials via the CBoostBridge C shim. The API surface includes:
//  - Cₙ^(λ)(x): Gegenbauer polynomial of degree n with parameter λ.
//  - d/dx Cₙ^(λ)(x): first derivative (prime).
//  - dᵏ/dxᵏ Cₙ^(λ)(x): k-th derivative.
//
//  Design goals:
//  - Type-generic entry points for any BinaryFloatingPoint (Double-backed).
//  - Type-specific overloads for Float and (on x86_64) Float80.
//  - Swift-side validation of common constraints (n, k ≥ 0; finiteness).
//  - Consistent error model using SpecialFunctionError.
//  - Thorough DocC comments describing domains, behavior, and references.
//
//  Notes on domain:
//  - Gegenbauer polynomials are classically associated with λ > −1/2 for
//    orthogonality on [−1, 1]. Boost’s implementation is defined more broadly.
//  - We validate only finiteness of λ and x and non-negativity/range of indices,
//    mirroring the style used in other wrappers. The recommended domain λ > −1/2
//    is documented but not enforced here.
//
//  References:
//  - NIST DLMF §18: https://dlmf.nist.gov/18
//  - Boost.Math Gegenbauer documentation
//

import CBoostBridge

public extension SpecialFunctions {
    
    // MARK: - Generic BinaryFloatingPoint overloads (Double-backed)
    
    /// Gegenbauer (ultraspherical) polynomial Cₙ^(λ)(x).
    ///
    /// Overview:
    /// - Computes Cₙ^(λ)(x) for degree `n ≥ 0`, parameter `λ`, and real `x`.
    /// - This generic overload accepts any `BinaryFloatingPoint` inputs and returns
    ///   a result of the same generic type `T`.
    /// - Internally, the computation is delegated to Boost.Math via a `Double`-backed
    ///   C shim and then converted back to `T`.
    ///
    /// Domain and behavior:
    /// - Mathematically defined for real `x` and many `λ`; commonly λ > −1/2 for
    ///   orthogonality on [−1, 1]. This wrapper does not enforce λ > −1/2, but
    ///   documents it as a recommended range.
    /// - `n` must be a non‑negative integer and fit into `UInt32` for the C shim.
    /// - `λ` and `x` must be finite (not NaN or ±∞).
    ///
    /// Parameters:
    /// - n: Degree (non‑negative integer; must be ≤ UInt32.max).
    /// - lambda: Gegenbauer parameter λ (finite real).
    /// - x: Evaluation point (finite real).
    ///
    /// Returns:
    /// - The value of Cₙ^(λ)(x) as `T`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotPositive(name: "n")` if `n < 0`.
    /// - `SpecialFunctionError.parameterExceedsMaximumIntegerValue(name: "n", max: ...)` if `n > UInt32.max`.
    /// - `SpecialFunctionError.parameterNotFinite(name: "lambda"|"x")` if inputs are NaN or ±∞.
    @inlinable static func gegenbauer<T: BinaryFloatingPoint>(n: Int, lambda: T, x: T) throws -> T {
        guard n >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "n") }
        guard n <= Int(UInt32.max) else { throw SpecialFunctionError.parameterExceedsMaximumIntegerValue(name: "n", max: Int(UInt32.max)) }
        let dl = D(lambda), dx = D(x)
        guard dl.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "lambda") }
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return T(bs_gegenbauer_d(UInt32(n), dl, dx))
    }

    // Mixed-precision promotions (Float ↔ Double) → Double
    /// Gegenbauer C_n^(λ)(x) with mixed `Float`/`Double` inputs; returns `Double`.
    @inlinable static func gegenbauer(n: Int, lambda: Float, x: Double) throws -> Double { try gegenbauer(n: n, lambda: Double(lambda), x: x) }
    /// Gegenbauer C_n^(λ)(x) with mixed `Double`/`Float` inputs; returns `Double`.
    @inlinable static func gegenbauer(n: Int, lambda: Double, x: Float) throws -> Double { try gegenbauer(n: n, lambda: lambda, x: Double(x)) }
    
    /// First derivative (with respect to x) of the Gegenbauer polynomial, d/dx Cₙ^(λ)(x).
    ///
    /// Parameters:
    /// - n: Degree (non‑negative integer; must be ≤ UInt32.max).
    /// - lambda: Gegenbauer parameter λ (finite real).
    /// - x: Evaluation point (finite real).
    ///
    /// Returns:
    /// - d/dx Cₙ^(λ)(x) as `T`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotPositive(name: "n")` if `n < 0`.
    /// - `SpecialFunctionError.parameterExceedsMaximumIntegerValue(name: "n", max: ...)` if `n > UInt32.max`.
    /// - `SpecialFunctionError.parameterNotFinite(name: "lambda"|"x")` if inputs are NaN or ±∞.
    @inlinable static func gegenbauerPrime<T: BinaryFloatingPoint>(n: Int, lambda: T, x: T) throws -> T {
        guard n >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "n") }
        guard n <= Int(UInt32.max) else { throw SpecialFunctionError.parameterExceedsMaximumIntegerValue(name: "n", max: Int(UInt32.max)) }
        let dl = D(lambda), dx = D(x)
        guard dl.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "lambda") }
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return T(bs_gegenbauer_prime_d(UInt32(n), dl, dx))
    }

    // Mixed-precision promotions (Float ↔ Double) → Double
    /// d/dx C_n^(λ)(x) with mixed `Float`/`Double` inputs; returns `Double`.
    @inlinable static func gegenbauerPrime(n: Int, lambda: Float, x: Double) throws -> Double { try gegenbauerPrime(n: n, lambda: Double(lambda), x: x) }
    /// d/dx C_n^(λ)(x) with mixed `Double`/`Float` inputs; returns `Double`.
    @inlinable static func gegenbauerPrime(n: Int, lambda: Double, x: Float) throws -> Double { try gegenbauerPrime(n: n, lambda: lambda, x: Double(x)) }
    
    /// k-th derivative (with respect to x) of the Gegenbauer polynomial, dᵏ/dxᵏ Cₙ^(λ)(x).
    ///
    /// Parameters:
    /// - n: Degree (non‑negative integer; must be ≤ UInt32.max).
    /// - lambda: Gegenbauer parameter λ (finite real).
    /// - x: Evaluation point (finite real).
    /// - k: Derivative order (non‑negative integer; must be ≤ UInt32.max).
    ///
    /// Returns:
    /// - dᵏ/dxᵏ Cₙ^(λ)(x) as `T`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotPositive(name: "n")` if `n < 0`.
    /// - `SpecialFunctionError.parameterExceedsMaximumIntegerValue(name: "n", max: ...)` if `n > UInt32.max`.
    /// - `SpecialFunctionError.parameterNotPositive(name: "k")` if `k < 0`.
    /// - `SpecialFunctionError.parameterExceedsMaximumIntegerValue(name: "k", max: ...)` if `k > UInt32.max`.
    /// - `SpecialFunctionError.parameterNotFinite(name: "lambda"|"x")` if inputs are NaN or ±∞.
    @inlinable static func gegenbauerDerivative<T: BinaryFloatingPoint>(n: Int, lambda: T, x: T, k: Int) throws -> T {
        guard n >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "n") }
        guard k >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "k") }
        guard n <= Int(UInt32.max) else { throw SpecialFunctionError.parameterExceedsMaximumIntegerValue(name: "n", max: Int(UInt32.max)) }
        guard k <= Int(UInt32.max) else { throw SpecialFunctionError.parameterExceedsMaximumIntegerValue(name: "k", max: Int(UInt32.max)) }
        let dl = D(lambda), dx = D(x)
        guard dl.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "lambda") }
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return T(bs_gegenbauer_derivative_d(UInt32(n), dl, dx, UInt32(k)))
    }

    // Mixed-precision promotions (Float ↔ Double) → Double
    /// d^k/dx^k C_n^(λ)(x) with mixed `Float`/`Double` inputs; returns `Double`.
    @inlinable static func gegenbauerDerivative(n: Int, lambda: Float, x: Double, k: Int) throws -> Double { try gegenbauerDerivative(n: n, lambda: Double(lambda), x: x, k: k) }
    /// d^k/dx^k C_n^(λ)(x) with mixed `Double`/`Float` inputs; returns `Double`.
    @inlinable static func gegenbauerDerivative(n: Int, lambda: Double, x: Float, k: Int) throws -> Double { try gegenbauerDerivative(n: n, lambda: lambda, x: Double(x), k: k) }
    
    
    // MARK: - Float overloads
    
    /// Cₙ^(λ)(x) for `Float`. Requires `n ≥ 0`, finite `lambda` and `x`.
    @inlinable static func gegenbauer(n: Int, lambda: Float, x: Float) throws -> Float {
        guard n >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "n") }
        guard n <= Int(UInt32.max) else { throw SpecialFunctionError.parameterExceedsMaximumIntegerValue(name: "n", max: Int(UInt32.max)) }
        guard lambda.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "lambda") }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return bs_gegenbauer_f(UInt32(n), lambda, x)
    }
    
    /// d/dx Cₙ^(λ)(x) for `Float`. Requires `n ≥ 0`, finite `lambda` and `x`.
    @inlinable static func gegenbauerPrime(n: Int, lambda: Float, x: Float) throws -> Float {
        guard n >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "n") }
        guard n <= Int(UInt32.max) else { throw SpecialFunctionError.parameterExceedsMaximumIntegerValue(name: "n", max: Int(UInt32.max)) }
        guard lambda.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "lambda") }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return bs_gegenbauer_prime_f(UInt32(n), lambda, x)
    }
    
    /// dᵏ/dxᵏ Cₙ^(λ)(x) for `Float`. Requires `n, k ≥ 0`, finite `lambda` and `x`.
    @inlinable static func gegenbauerDerivative(n: Int, lambda: Float, x: Float, k: Int) throws -> Float {
        guard n >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "n") }
        guard k >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "k") }
        guard n <= Int(UInt32.max) else { throw SpecialFunctionError.parameterExceedsMaximumIntegerValue(name: "n", max: Int(UInt32.max)) }
        guard k <= Int(UInt32.max) else { throw SpecialFunctionError.parameterExceedsMaximumIntegerValue(name: "k", max: Int(UInt32.max)) }
        guard lambda.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "lambda") }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return bs_gegenbauer_derivative_f(UInt32(n), lambda, x, UInt32(k))
    }
    
    
    // MARK: - Float80 overloads (x86_64 only)
    
#if arch(x86_64)
    /// Cₙ^(λ)(x) for `Float80` (x86_64 only). Requires `n ≥ 0`, finite `lambda` and `x`.
    @inlinable static func gegenbauer(n: Int, lambda: Float80, x: Float80) throws -> Float80 {
        guard n >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "n") }
        guard n <= Int(UInt32.max) else { throw SpecialFunctionError.parameterExceedsMaximumIntegerValue(name: "n", max: Int(UInt32.max)) }
        guard lambda.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "lambda") }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return bs_gegenbauer_l(UInt32(n), lambda, x)
    }

    // Mixed promotions with Float80 → Float80
    @inlinable static func gegenbauer(n: Int, lambda: Float80, x: Double) throws -> Float80 { try gegenbauer(n: n, lambda: lambda, x: Float80(x)) }
    @inlinable static func gegenbauer(n: Int, lambda: Double, x: Float80) throws -> Float80 { try gegenbauer(n: n, lambda: Float80(lambda), x: x) }
    @inlinable static func gegenbauer(n: Int, lambda: Float80, x: Float) throws -> Float80 { try gegenbauer(n: n, lambda: lambda, x: Float80(x)) }
    @inlinable static func gegenbauer(n: Int, lambda: Float, x: Float80) throws -> Float80 { try gegenbauer(n: n, lambda: Float80(lambda), x: x) }
    
    /// d/dx Cₙ^(λ)(x) for `Float80` (x86_64 only). Requires `n ≥ 0`, finite `lambda` and `x`.
    @inlinable static func gegenbauerPrime(n: Int, lambda: Float80, x: Float80) throws -> Float80 {
        guard n >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "n") }
        guard n <= Int(UInt32.max) else { throw SpecialFunctionError.parameterExceedsMaximumIntegerValue(name: "n", max: Int(UInt32.max)) }
        guard lambda.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "lambda") }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return bs_gegenbauer_prime_l(UInt32(n), lambda, x)
    }

    // Mixed promotions with Float80 → Float80
    @inlinable static func gegenbauerPrime(n: Int, lambda: Float80, x: Double) throws -> Float80 { try gegenbauerPrime(n: n, lambda: lambda, x: Float80(x)) }
    @inlinable static func gegenbauerPrime(n: Int, lambda: Double, x: Float80) throws -> Float80 { try gegenbauerPrime(n: n, lambda: Float80(lambda), x: x) }
    @inlinable static func gegenbauerPrime(n: Int, lambda: Float80, x: Float) throws -> Float80 { try gegenbauerPrime(n: n, lambda: lambda, x: Float80(x)) }
    @inlinable static func gegenbauerPrime(n: Int, lambda: Float, x: Float80) throws -> Float80 { try gegenbauerPrime(n: n, lambda: Float80(lambda), x: x) }
    
    /// dᵏ/dxᵏ Cₙ^(λ)(x) for `Float80` (x86_64 only). Requires `n, k ≥ 0`, finite `lambda` and `x`.
    @inlinable static func gegenbauerDerivative(n: Int, lambda: Float80, x: Float80, k: Int) throws -> Float80 {
        guard n >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "n") }
        guard k >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "k") }
        guard n <= Int(UInt32.max) else { throw SpecialFunctionError.parameterExceedsMaximumIntegerValue(name: "n", max: Int(UInt32.max)) }
        guard k <= Int(UInt32.max) else { throw SpecialFunctionError.parameterExceedsMaximumIntegerValue(name: "k", max: Int(UInt32.max)) }
        guard lambda.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "lambda") }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return bs_gegenbauer_derivative_l(UInt32(n), lambda, x, UInt32(k))
    }

    // Mixed promotions with Float80 → Float80
    @inlinable static func gegenbauerDerivative(n: Int, lambda: Float80, x: Double, k: Int) throws -> Float80 { try gegenbauerDerivative(n: n, lambda: lambda, x: Float80(x), k: k) }
    @inlinable static func gegenbauerDerivative(n: Int, lambda: Double, x: Float80, k: Int) throws -> Float80 { try gegenbauerDerivative(n: n, lambda: Float80(lambda), x: x, k: k) }
    @inlinable static func gegenbauerDerivative(n: Int, lambda: Float80, x: Float, k: Int) throws -> Float80 { try gegenbauerDerivative(n: n, lambda: lambda, x: Float80(x), k: k) }
    @inlinable static func gegenbauerDerivative(n: Int, lambda: Float, x: Float80, k: Int) throws -> Float80 { try gegenbauerDerivative(n: n, lambda: Float80(lambda), x: x, k: k) }
#endif
}
