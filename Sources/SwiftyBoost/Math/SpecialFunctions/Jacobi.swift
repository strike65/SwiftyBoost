//
//  Jacobi.swift
//  Math/SpecialFunctions
//
//  Created by Volker Thieme 2025.
//  License: MIT (see project root)
//
//  Swift wrappers for the Jacobi polynomials Pₙ^{(α,β)}(x) and their derivatives.
//  These wrappers validate indices and finiteness before delegating to Boost.Math
//  through the CBoostBridge layer.
//
//  Overview
//  --------
//  The Jacobi polynomials Pₙ^{(α,β)}(x) form a family of classical orthogonal
//  polynomials on x ∈ [−1, 1] with respect to the weight
//      w(x) = (1 − x)^α (1 + x)^β,
//  where α > −1 and β > −1 ensure orthogonality and finite norms.
//  They generalize several other families (e.g., Legendre for α = β = 0,
//  Gegenbauer for α = β = λ − 1/2).
//
//  What this file provides
//  -----------------------
//  - Pₙ^{(α,β)}(x) evaluation
//  - First derivative d/dx Pₙ^{(α,β)}(x)
//  - Second derivative d²/dx² Pₙ^{(α,β)}(x)
//  - k-th derivative dᵏ/dxᵏ Pₙ^{(α,β)}(x)
//  Each is exposed generically over BinaryFloatingPoint (computed via Double in the bridge),
//  plus dedicated Float overloads and Float80 overloads on x86.
//
//  Validation & Errors
//  -------------------
//  - n (degree) must satisfy n ≥ 0 and fit into UInt32.
//  - k (derivative order) must satisfy k ≥ 0 and fit into UInt32.
//  - alpha, beta, and x must be finite reals.
//  Violations throw SpecialFunctionError<T> with parameter names and offending values.
//
//  Numeric notes
//  -------------
//  - Orthogonality is guaranteed for α, β > −1; outside this region,
//    Pₙ^{(α,β)} is still defined by analytic continuation except at singularities.
//  - Large degrees or extreme parameter values may amplify rounding error;
//    Boost.Math uses stable recurrences and internal scaling to mitigate this,
//    but callers should still consider conditioning and tolerance in downstream use.
//  - For k = 1 or k = 2, prefer the dedicated jacobiPrime / jacobiDoublePrime
//    for slightly less argument handling overhead.
//
//  References
//  ----------
//  - NIST DLMF §18 (Orthogonal Polynomials): https://dlmf.nist.gov/18
//  - Szegő, “Orthogonal Polynomials” (AMS).
//  - Boost.Math documentation for Jacobi polynomials.
//
//  Examples
//  --------
//    // Evaluate P₃^{(1,2)}(0.25)
//    let p: Double = try SpecialFunctions.jacobi(n: 3, alpha: 1.0, beta: 2.0, x: 0.25)
//
//    // First derivative at x = 0
//    let dp: Double = try SpecialFunctions.jacobiPrime(n: 5, alpha: 0.5, beta: 0.5, x: 0.0)
//
//    // Second derivative in Float
//    let d2p: Float = try SpecialFunctions.jacobiDoublePrime(n: 2, alpha: 0.0, beta: 0.0, x: 0.1)
//
//    // k-th derivative (k = 4)
//    let dk: Double = try SpecialFunctions.jacobiDerivative(n: 12, alpha: 1.5, beta: 0.25, x: 0.9, k: 4)
//

import SwiftyBoostPrelude

public extension SpecialFunctions {

    // MARK: - Generic overloads (Double-backed)

    /// Jacobi polynomial Pₙ^{(α,β)}(x).
    ///
    /// Computes the degree-`n` Jacobi polynomial at `x` for real parameters `alpha` (α) and `beta` (β).
    ///
    /// - Mathematical background:
    ///   Pₙ^{(α,β)}(x) are orthogonal on [−1, 1] with weight (1 − x)^α (1 + x)^β for α, β > −1.
    ///   Outside this range the function is defined by analytic continuation where applicable.
    ///
    /// - Parameters:
    ///   - n: Degree (requires `n ≥ 0` and representable as `UInt32`).
    ///   - alpha: Parameter α (finite real; typical orthogonality domain α > −1).
    ///   - beta: Parameter β (finite real; typical orthogonality domain β > −1).
    ///   - x: Evaluation point (finite real).
    ///
    /// - Returns: `Pₙ^{(α,β)}(x)` as `T`.
    ///
    /// - Throws:
    ///   - `SpecialFunctionError.parameterNotPositive("n", …)` if `n < 0`.
    ///   - `SpecialFunctionError.parameterExceedsMaximumIntegerValue("n", …)` if `n > UInt32.max`.
    ///   - `SpecialFunctionError.parameterNotFinite` if any of `alpha`, `beta`, or `x` is NaN/±∞.
    ///
    /// - Complexity: O(1) per call (delegates to Boost.Math).
    ///
    /// - SeeAlso: `jacobiPrime(n:alpha:beta:x:)`, `jacobiDoublePrime(n:alpha:beta:x:)`, `jacobiDerivative(n:alpha:beta:x:k:)`
    @inlinable static func jacobi<T: Real & BinaryFloatingPoint & Sendable>(n: Int, alpha: T, beta: T, x: T) throws -> T {
        let nu = try validateJacobiIndex(n, as: T.self)
        let da = D(alpha), db = D(beta), dx = D(x)
        guard da.isFinite else { throw SpecialFunctionError<T>.parameterNotFinite(name: "alpha", value: alpha) }
        guard db.isFinite else { throw SpecialFunctionError<T>.parameterNotFinite(name: "beta", value: beta) }
        guard dx.isFinite else { throw SpecialFunctionError<T>.parameterNotFinite(name: "x", value: x) }
        return T(bs_jacobi_d(nu, da, db, dx))
    }

    /// First derivative d/dx Pₙ^{(α,β)}(x).
    ///
    /// - Parameters:
    ///   - n: Degree (requires `n ≥ 0`).
    ///   - alpha: Parameter α (finite real).
    ///   - beta: Parameter β (finite real).
    ///   - x: Evaluation point (finite real).
    ///
    /// - Returns: `d/dx Pₙ^{(α,β)}(x)` as `T`.
    ///
    /// - Throws: `SpecialFunctionError<T>` for invalid `n` or any non-finite parameter.
    ///
    /// - Discussion:
    ///   Consistent with identities linking Pₙ^{(α,β)} and Pₙ₋₁^{(α+1,β+1)}; Boost.Math evaluates the derivative directly with stable kernels.
    @inlinable static func jacobiPrime<T: Real & BinaryFloatingPoint & Sendable>(n: Int, alpha: T, beta: T, x: T) throws -> T {
        let nu = try validateJacobiIndex(n, as: T.self)
        let da = D(alpha), db = D(beta), dx = D(x)
        guard da.isFinite else { throw SpecialFunctionError<T>.parameterNotFinite(name: "alpha", value: alpha) }
        guard db.isFinite else { throw SpecialFunctionError<T>.parameterNotFinite(name: "beta", value: beta) }
        guard dx.isFinite else { throw SpecialFunctionError<T>.parameterNotFinite(name: "x", value: x) }
        return T(bs_jacobi_prime_d(nu, da, db, dx))
    }

    /// Second derivative d²/dx² Pₙ^{(α,β)}(x).
    ///
    /// - Parameters:
    ///   - n: Degree (requires `n ≥ 0`).
    ///   - alpha: Parameter α (finite real).
    ///   - beta: Parameter β (finite real).
    ///   - x: Evaluation point (finite real).
    ///
    /// - Returns: `d²/dx² Pₙ^{(α,β)}(x)` as `T`.
    ///
    /// - Throws: `SpecialFunctionError<T>` for invalid `n` or any non-finite parameter.
    @inlinable static func jacobiDoublePrime<T: Real & BinaryFloatingPoint & Sendable>(n: Int, alpha: T, beta: T, x: T) throws -> T {
        let nu = try validateJacobiIndex(n, as: T.self)
        let da = D(alpha), db = D(beta), dx = D(x)
        guard da.isFinite else { throw SpecialFunctionError<T>.parameterNotFinite(name: "alpha", value: alpha) }
        guard db.isFinite else { throw SpecialFunctionError<T>.parameterNotFinite(name: "beta", value: beta) }
        guard dx.isFinite else { throw SpecialFunctionError<T>.parameterNotFinite(name: "x", value: x) }
        return T(bs_jacobi_double_prime_d(nu, da, db, dx))
    }

    /// k-th derivative dᵏ/dxᵏ Pₙ^{(α,β)}(x).
    ///
    /// - Parameters:
    ///   - n: Degree (requires `n ≥ 0`).
    ///   - alpha: Parameter α (finite real).
    ///   - beta: Parameter β (finite real).
    ///   - x: Evaluation point (finite real).
    ///   - k: Derivative order (requires `k ≥ 0` and representable as `UInt32`).
    ///
    /// - Returns: `dᵏ/dxᵏ Pₙ^{(α,β)}(x)` as `T`.
    ///
    /// - Throws:
    ///   - `SpecialFunctionError.parameterNotPositive("n" | "k", …)` if `n < 0` or `k < 0`.
    ///   - `SpecialFunctionError.parameterExceedsMaximumIntegerValue("n" | "k", …)` if the value exceeds `UInt32.max`.
    ///   - `SpecialFunctionError.parameterNotFinite` if any of `alpha`, `beta`, or `x` is NaN/±∞.
    ///
    /// - Discussion:
    ///   For `k == 0`, this reduces to `jacobi`. For `k == 1` or `k == 2`, prefer the dedicated helpers for slightly lower overhead.
    @inlinable static func jacobiDerivative<T: Real & BinaryFloatingPoint & Sendable>(n: Int, alpha: T, beta: T, x: T, k: Int) throws -> T {
        let nu = try validateJacobiIndex(n, as: T.self)
        let ku = try validateJacobiDerivativeOrder(k, as: T.self)
        let da = D(alpha), db = D(beta), dx = D(x)
        guard da.isFinite else { throw SpecialFunctionError<T>.parameterNotFinite(name: "alpha", value: alpha) }
        guard db.isFinite else { throw SpecialFunctionError<T>.parameterNotFinite(name: "beta", value: beta) }
        guard dx.isFinite else { throw SpecialFunctionError<T>.parameterNotFinite(name: "x", value: x) }
        return T(bs_jacobi_derivative_d(nu, da, db, dx, ku))
    }

    // MARK: - Float overloads

    /// Jacobi polynomial Pₙ^{(α,β)}(x) for `Float`.
    ///
    /// See the generic overload for mathematical background and error conditions.
    @inlinable static func jacobi(n: Int, alpha: Float, beta: Float, x: Float) throws -> Float {
        let nu = try validateJacobiIndex(n, as: Float.self)
        guard alpha.isFinite else { throw SpecialFunctionError<Float>.parameterNotFinite(name: "alpha", value: alpha) }
        guard beta.isFinite else { throw SpecialFunctionError<Float>.parameterNotFinite(name: "beta", value: beta) }
        guard x.isFinite else { throw SpecialFunctionError<Float>.parameterNotFinite(name: "x", value: x) }
        return bs_jacobi_f(nu, alpha, beta, x)
    }

    /// First derivative d/dx Pₙ^{(α,β)}(x) for `Float`.
    @inlinable static func jacobiPrime(n: Int, alpha: Float, beta: Float, x: Float) throws -> Float {
        let nu = try validateJacobiIndex(n, as: Float.self)
        guard alpha.isFinite else { throw SpecialFunctionError<Float>.parameterNotFinite(name: "alpha", value: alpha) }
        guard beta.isFinite else { throw SpecialFunctionError<Float>.parameterNotFinite(name: "beta", value: beta) }
        guard x.isFinite else { throw SpecialFunctionError<Float>.parameterNotFinite(name: "x", value: x) }
        return bs_jacobi_prime_f(nu, alpha, beta, x)
    }

    /// Second derivative d²/dx² Pₙ^{(α,β)}(x) for `Float`.
    @inlinable static func jacobiDoublePrime(n: Int, alpha: Float, beta: Float, x: Float) throws -> Float {
        let nu = try validateJacobiIndex(n, as: Float.self)
        guard alpha.isFinite else { throw SpecialFunctionError<Float>.parameterNotFinite(name: "alpha", value: alpha) }
        guard beta.isFinite else { throw SpecialFunctionError<Float>.parameterNotFinite(name: "beta", value: beta) }
        guard x.isFinite else { throw SpecialFunctionError<Float>.parameterNotFinite(name: "x", value: x) }
        return bs_jacobi_double_prime_f(nu, alpha, beta, x)
    }

    /// k-th derivative dᵏ/dxᵏ Pₙ^{(α,β)}(x) for `Float`.
    @inlinable static func jacobiDerivative(n: Int, alpha: Float, beta: Float, x: Float, k: Int) throws -> Float {
        let nu = try validateJacobiIndex(n, as: Float.self)
        let ku = try validateJacobiDerivativeOrder(k, as: Float.self)
        guard alpha.isFinite else { throw SpecialFunctionError<Float>.parameterNotFinite(name: "alpha", value: alpha) }
        guard beta.isFinite else { throw SpecialFunctionError<Float>.parameterNotFinite(name: "beta", value: beta) }
        guard x.isFinite else { throw SpecialFunctionError<Float>.parameterNotFinite(name: "x", value: x) }
        return bs_jacobi_derivative_f(nu, alpha, beta, x, ku)
    }

    // MARK: - Float80 overloads (x86 only)

    #if arch(x86_64) || arch(i386)
    /// Jacobi polynomial Pₙ^{(α,β)}(x) for `Float80` (x86 only).
    @inlinable static func jacobi(n: Int, alpha: Float80, beta: Float80, x: Float80) throws -> Float80 {
        let nu = try validateJacobiIndex(n, as: Float80.self)
        guard alpha.isFinite else { throw SpecialFunctionError<Float80>.parameterNotFinite(name: "alpha", value: alpha) }
        guard beta.isFinite else { throw SpecialFunctionError<Float80>.parameterNotFinite(name: "beta", value: beta) }
        guard x.isFinite else { throw SpecialFunctionError<Float80>.parameterNotFinite(name: "x", value: x) }
        return bs_jacobi_l(nu, alpha, beta, x)
    }

    /// First derivative d/dx Pₙ^{(α,β)}(x) for `Float80` (x86 only).
    @inlinable static func jacobiPrime(n: Int, alpha: Float80, beta: Float80, x: Float80) throws -> Float80 {
        let nu = try validateJacobiIndex(n, as: Float80.self)
        guard alpha.isFinite else { throw SpecialFunctionError<Float80>.parameterNotFinite(name: "alpha", value: alpha) }
        guard beta.isFinite else { throw SpecialFunctionError<Float80>.parameterNotFinite(name: "beta", value: beta) }
        guard x.isFinite else { throw SpecialFunctionError<Float80>.parameterNotFinite(name: "x", value: x) }
        return bs_jacobi_prime_l(nu, alpha, beta, x)
    }

    /// Second derivative d²/dx² Pₙ^{(α,β)}(x) for `Float80` (x86 only).
    @inlinable static func jacobiDoublePrime(n: Int, alpha: Float80, beta: Float80, x: Float80) throws -> Float80 {
        let nu = try validateJacobiIndex(n, as: Float80.self)
        guard alpha.isFinite else { throw SpecialFunctionError<Float80>.parameterNotFinite(name: "alpha", value: alpha) }
        guard beta.isFinite else { throw SpecialFunctionError<Float80>.parameterNotFinite(name: "beta", value: beta) }
        guard x.isFinite else { throw SpecialFunctionError<Float80>.parameterNotFinite(name: "x", value: x) }
        return bs_jacobi_double_prime_l(nu, alpha, beta, x)
    }

    /// k-th derivative dᵏ/dxᵏ Pₙ^{(α,β)}(x) for `Float80` (x86 only).
    @inlinable static func jacobiDerivative(n: Int, alpha: Float80, beta: Float80, x: Float80, k: Int) throws -> Float80 {
        let nu = try validateJacobiIndex(n, as: Float80.self)
        let ku = try validateJacobiDerivativeOrder(k, as: Float80.self)
        guard alpha.isFinite else { throw SpecialFunctionError<Float80>.parameterNotFinite(name: "alpha", value: alpha) }
        guard beta.isFinite else { throw SpecialFunctionError<Float80>.parameterNotFinite(name: "beta", value: beta) }
        guard x.isFinite else { throw SpecialFunctionError<Float80>.parameterNotFinite(name: "x", value: x) }
        return bs_jacobi_derivative_l(nu, alpha, beta, x, ku)
    }
    #endif

    // MARK: - Validation helpers

    /// Validates and converts the Jacobi degree `n` to `UInt32`.
    ///
    /// - Throws:
    ///   - `SpecialFunctionError.parameterNotPositive("n", …)` if `n < 0`.
    ///   - `SpecialFunctionError.parameterExceedsMaximumIntegerValue("n", …)` if `n > UInt32.max`.
    @usableFromInline internal static func validateJacobiIndex<T: Real & BinaryFloatingPoint & Sendable>(_ n: Int, as _: T.Type) throws -> UInt32 {
        guard n >= 0 else { throw SpecialFunctionError<T>.parameterNotPositive(name: "n", value: T(n)) }
        guard let nu = UInt32(exactly: n) else {
            throw SpecialFunctionError<T>.parameterExceedsMaximumIntegerValue(name: "n", max: Int(UInt32.max))
        }
        return nu
    }

    /// Validates and converts the derivative order `k` to `UInt32`.
    ///
    /// - Throws:
    ///   - `SpecialFunctionError.parameterNotPositive("k", …)` if `k < 0`.
    ///   - `SpecialFunctionError.parameterExceedsMaximumIntegerValue("k", …)` if `k > UInt32.max`.
    @usableFromInline internal static func validateJacobiDerivativeOrder<T: Real & BinaryFloatingPoint & Sendable>(_ k: Int, as _: T.Type) throws -> UInt32 {
        guard k >= 0 else { throw SpecialFunctionError<T>.parameterNotPositive(name: "k", value: T(k)) }
        guard let ku = UInt32(exactly: k) else {
            throw SpecialFunctionError<T>.parameterExceedsMaximumIntegerValue(name: "k", max: Int(UInt32.max))
        }
        return ku
    }
}
