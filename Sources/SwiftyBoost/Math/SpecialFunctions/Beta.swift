//
//  Beta.swift
//  Math/SpecialFunctions
//
//  Swift wrappers around the (complete) Beta function B(a, b), the incomplete Beta,
//  its regularized forms, and inverse utilities. These APIs validate arguments and
//  surface domain issues via a Swifty error model, while delegating numerical work
//  to CBoostBridge (Boost.Math).
//
//  Functions included (real-valued):
//  - beta(a, b) = B(a, b) = Γ(a)Γ(b) / Γ(a + b), for a > 0, b > 0.
//  - incompleteBetaUnnormalized(a, b, x) = B_x(a, b) = ∫₀ˣ t^{a-1} (1−t)^{b-1} dt, for a > 0, b > 0, x ∈ [0, 1].
//  - regularizedIncompleteBeta(a, b, x) = I_x(a, b) = B_x(a, b) / B(a, b), for a > 0, b > 0, x ∈ [0, 1].
//  - complementaryRegularizedIncompleteBeta(a, b, x) = 1 − I_x(a, b), same domain.
//  - inverseRegularizedIncompleteBeta(a, b, p): solve x in I_x(a, b) = p, for a > 0, b > 0, p ∈ [0, 1].
//  - inverseComplementaryRegularizedIncompleteBeta(a, b, p): solve x in 1 − I_x(a, b) = p, for a > 0, b > 0, p ∈ [0, 1].
//  - solveAForRegularizedIncompleteBeta(b, x, p): solve a in I_x(a, b) = p (no throws; see notes).
//  - solveBForRegularizedIncompleteBeta(a, x, p): solve b in I_x(a, b) = p (no throws; see notes).
//  - regularizedIncompleteBetaDerivative(a, b, x): d/dx I_x(a, b), for a > 0, b > 0, x ∈ [0, 1].
//
//  Overloads:
//  - Generic BinaryFloatingPoint overloads funnel through a Double-backed implementation
//    using the internal helper D(_:), and convert results back to T.
//  - Fast paths are provided for Float and (on x86_64) Float80, directly invoking
//    the corresponding C functions to avoid conversions.
//
//  Domain notes (real-valued):
//  - beta(a, b) requires a > 0 and b > 0.
//  - Incomplete/regularized forms require a > 0, b > 0, and x ∈ [0, 1].
//  - Inverse forms require a > 0, b > 0, and p ∈ [0, 1].
//  - The “solveA/solveB” utilities are root-finders for parameters; they do not throw.
//    Provide valid inputs (b > 0 or a > 0, x ∈ [0,1], p ∈ [0,1]). Invalid inputs may
//    result in NaN or undefined behavior from the backend.
//
//  Error model:
//  - Generic and Float/Float80 overloads validate inputs and throw SpecialFunctionError
//    for domain violations: parameterNotPositive, parameterOutOfRange, parameterNotFinite.
//  - The solveA/solveB helpers do not validate and do not throw; callers should ensure
//    valid inputs.
//
//  References:
//  - NIST DLMF §8 (Beta Function): https://dlmf.nist.gov/8
//  - Boost.Math Beta functions:
//    https://www.boost.org/doc/libs/release/libs/math/doc/html/math_toolkit/sf_beta.html
//

import SwiftyBoostPrelude

public extension SpecialFunctions {
    
    
    // MARK: - Beta family (generic Double-backed)
    
    /// Compute the complete Beta function B(a, b) = Γ(a)Γ(b) / Γ(a + b).
    ///
    /// Throws:
    /// - SpecialFunctionError.parameterNotFinite(name: ...) if inputs are NaN or ±∞.
    /// - SpecialFunctionError.parameterNotPositive(name: "a", value: a) if a ≤ 0.
    /// - SpecialFunctionError.parameterNotPositive(name: "b", value: b) if b ≤ 0.
    @inlinable static func beta<T: Real & BinaryFloatingPoint & Sendable>(_ a: T, _ b: T) throws -> T {
        guard a.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a", value: a) }
        guard b.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "b", value: b) }
        guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a", value: a) }
        guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b", value: b) }
        return T(bs_beta_d(D(a), D(b)))
    }

    /// Unnormalized incomplete Beta B_x(a, b) = ∫₀ˣ t^{a−1} (1−t)^{b−1} dt.
    ///
    /// Throws:
    /// - SpecialFunctionError.parameterNotFinite(name: ...) if inputs are NaN or ±∞.
    /// - SpecialFunctionError.parameterNotPositive(name: "a", value: a) if a ≤ 0.
    /// - SpecialFunctionError.parameterNotPositive(name: "b", value: b) if b ≤ 0.
    /// - SpecialFunctionError.parameterOutOfRange(name: "x", min: 0, max: 1) if x ∉ [0, 1].
    @inlinable static func incompleteBetaUnnormalized<T: Real & BinaryFloatingPoint & Sendable>(_ a: T, _ b: T, x: T) throws -> T {
        guard a.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a", value: a) }
        guard b.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "b", value: b) }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a", value: a) }
        guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b", value: b) }
        guard x >= 0 && x <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 0.0, max: 1.0) }
        return T(bs_fullBeta_d(D(a), D(b), D(x)))
    }

    /// Regularized incomplete Beta I_x(a, b) = B_x(a, b) / B(a, b).
    ///
    /// Throws:
    /// - SpecialFunctionError.parameterNotFinite(name: ...) if inputs are NaN or ±∞.
    /// - SpecialFunctionError.parameterNotPositive(name: "a", value: a) if a ≤ 0.
    /// - SpecialFunctionError.parameterNotPositive(name: "b", value: b) if b ≤ 0.
    /// - SpecialFunctionError.parameterOutOfRange(name: "x", min: 0, max: 1) if x ∉ [0, 1].
    @inlinable static func regularizedIncompleteBeta<T: Real & BinaryFloatingPoint & Sendable>(_ a: T, _ b: T, x: T) throws -> T {
        guard a.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a", value: a) }
        guard b.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "b", value: b) }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a", value: a) }
        guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b", value: b) }
        guard x >= 0 && x <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 0.0, max: 1.0) }
        return T(bs_ibeta_d(D(a), D(b), D(x)))
    }

    /// Complementary regularized incomplete Beta 1 − I_x(a, b).
    ///
    /// Throws:
    /// - SpecialFunctionError.parameterNotFinite(name: ...) if inputs are NaN or ±∞.
    /// - SpecialFunctionError.parameterNotPositive(name: "a", value: a) if a ≤ 0.
    /// - SpecialFunctionError.parameterNotPositive(name: "b", value: b) if b ≤ 0.
    /// - SpecialFunctionError.parameterOutOfRange(name: "x", min: 0, max: 1) if x ∉ [0, 1].
    @inlinable static func complementaryRegularizedIncompleteBeta<T: Real & BinaryFloatingPoint & Sendable>(_ a: T, _ b: T, x: T) throws -> T {
        guard a.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a", value: a) }
        guard b.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "b", value: b) }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a", value: a) }
        guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b", value: b) }
        guard x >= 0 && x <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 0.0, max: 1.0) }
        return T(bs_ibetac_d(D(a), D(b), D(x)))
    }

    /// Invert the regularized incomplete Beta: find x such that I_x(a, b) = p.
    ///
    /// Throws:
    /// - SpecialFunctionError.parameterNotFinite(name: ...) if inputs are NaN or ±∞.
    /// - SpecialFunctionError.parameterNotPositive(name: "a", value: a) if a ≤ 0.
    /// - SpecialFunctionError.parameterNotPositive(name: "b", value: b) if b ≤ 0.
    /// - SpecialFunctionError.parameterOutOfRange(name: "p", min: 0, max: 1) if p ∉ [0, 1].
    @inlinable static func inverseRegularizedIncompleteBeta<T: Real & BinaryFloatingPoint & Sendable>(_ a: T, _ b: T, p: T) throws -> T {
        guard a.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a", value: a) }
        guard b.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "b", value: b) }
        guard p.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "p", value: p) }
        guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a", value: a) }
        guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b", value: b) }
        guard p >= 0 && p <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "p", min: 0.0, max: 1.0) }
        return T(bs_ibeta_inv_d(D(a), D(b), D(p)))
    }

    /// Invert the complementary regularized incomplete Beta: find x such that 1 − I_x(a, b) = p.
    ///
    /// Throws:
    /// - SpecialFunctionError.parameterNotFinite(name: ...) if inputs are NaN or ±∞.
    /// - SpecialFunctionError.parameterNotPositive(name: "a", value: a) if a ≤ 0.
    /// - SpecialFunctionError.parameterNotPositive(name: "b", value: b) if b ≤ 0.
    /// - SpecialFunctionError.parameterOutOfRange(name: "p", min: 0, max: 1) if p ∉ [0, 1].
    @inlinable static func inverseComplementaryRegularizedIncompleteBeta<T: Real & BinaryFloatingPoint & Sendable>(_ a: T, _ b: T, p: T) throws -> T {
        guard a.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a", value: a) }
        guard b.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "b", value: b) }
        guard p.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "p", value: p) }
        guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a", value: a) }
        guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b", value: b) }
        guard p >= 0 && p <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "p", min: 0.0, max: 1.0) }
        return T(bs_ibetac_inv_d(D(a), D(b), D(p)))
    }

   
    /// Solve for a in I_x(a, b) = p, given b, x, p. No throws; see notes.
    @inlinable static func solveAForRegularizedIncompleteBeta<T: Real & BinaryFloatingPoint & Sendable>(b: T, x: T, p: T) -> T {
        T(bs_ibeta_inva_d(D(b), D(x), D(p)))
    }

    /// Solve for b in I_x(a, b) = p, given a, x, p. No throws; see notes.
    @inlinable static func solveBForRegularizedIncompleteBeta<T: Real & BinaryFloatingPoint & Sendable>(a: T, x: T, p: T) -> T {
        T(bs_ibeta_invb_d(D(a), D(x), D(p)))
    }
   
    // MARK: - Derivative of the regularized incomplete Beta (generic)
    
    /// Derivative with respect to x of the regularized incomplete Beta I_x(a, b).
    ///
    /// Throws:
    /// - SpecialFunctionError.parameterNotFinite(name: ...) if inputs are NaN or ±∞.
    /// - SpecialFunctionError.parameterNotPositive(name: "a", value: a) if a ≤ 0.
    /// - SpecialFunctionError.parameterNotPositive(name: "b", value: b) if b ≤ 0.
    /// - SpecialFunctionError.parameterOutOfRange(name: "x", min: 0, max: 1) if x ∉ [0, 1].
    @inlinable static func regularizedIncompleteBetaDerivative<T: Real & BinaryFloatingPoint & Sendable>(_ a: T, _ b: T, x: T) throws -> T {
        guard a.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a", value: a) }
        guard b.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "b", value: b) }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a", value: a) }
        guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b", value: b) }
        guard x >= 0 && x <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 0.0, max: 1.0) }
        return T(bs_ibeta_derivative_d(D(a), D(b), D(x)))
    }
    
    // MARK: - Float overloads (Beta)
    
    @inlinable static func beta(_ a: Float, _ b: Float) throws -> Float {
        guard a.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a", value: a) }
        guard b.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "b", value: b) }
        guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a", value: a) }
        guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b", value: b) }
        return bs_beta_f(a, b)
    }
    
    @inlinable static func incompleteBetaUnnormalized(_ a: Float, _ b: Float, x: Float) throws -> Float {
        guard a.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a", value: a) }
        guard b.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "b", value: b) }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a", value: a) }
        guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b", value: b) }
        guard x >= 0 && x <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 0.0, max: 1.0) }
        return bs_fullBeta_f(a, b, x)
    }
    
    @inlinable static func regularizedIncompleteBeta(_ a: Float, _ b: Float, x: Float) throws -> Float {
        guard a.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a", value: a) }
        guard b.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "b", value: b) }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a", value: a) }
        guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b", value: b) }
        guard x >= 0 && x <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 0.0, max: 1.0) }
        return bs_ibeta_f(a, b, x)
    }
    
    @inlinable static func complementaryRegularizedIncompleteBeta(_ a: Float, _ b: Float, x: Float) throws -> Float {
        guard a.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a", value: a) }
        guard b.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "b", value: b) }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a", value: a) }
        guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b", value: b) }
        guard x >= 0 && x <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 0.0, max: 1.0) }
        return bs_ibetac_f(a, b, x)
    }
    
    @inlinable static func inverseRegularizedIncompleteBeta(_ a: Float, _ b: Float, p: Float) throws -> Float {
        guard a.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a", value: a) }
        guard b.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "b", value: b) }
        guard p.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "p", value: p) }
        guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a", value: a) }
        guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b", value: b) }
        guard p >= 0 && p <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "p", min: 0.0, max: 1.0) }
        return bs_ibeta_inv_f(a, b, p)
    }
    
    @inlinable static func inverseComplementaryRegularizedIncompleteBeta(_ a: Float, _ b: Float, p: Float) throws -> Float {
        guard a.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a", value: a) }
        guard b.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "b", value: b) }
        guard p.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "p", value: p) }
        guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a", value: a) }
        guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b", value: b) }
        guard p >= 0 && p <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "p", min: 0.0, max: 1.0) }
        return bs_ibetac_inv_f(a, b, p)
    }
    
    @inlinable static func solveAForRegularizedIncompleteBeta(b: Float, x: Float, p: Float) -> Float {
        bs_ibeta_inva_f(b, x, p)
    }
    
    @inlinable static func solveBForRegularizedIncompleteBeta(a: Float, x: Float, p: Float) -> Float {
        bs_ibeta_invb_f(a, x, p)
    }
    
    @inlinable static func regularizedIncompleteBetaDerivative(_ a: Float, _ b: Float, x: Float) throws -> Float {
        guard a.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a", value: a) }
        guard b.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "b", value: b) }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a", value: a) }
        guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b", value: b) }
        guard x >= 0 && x <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 0.0, max: 1.0) }
        return bs_ibeta_derivative_f(a, b, x)
    }
    
    // MARK: - Float80 overloads (x86_64)
    
#if arch(x86_64)
    @inlinable static func beta(_ a: Float80, _ b: Float80) throws -> Float80 {
        guard a.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a", value: a) }
        guard b.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "b", value: b) }
        guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a", value: a) }
        guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b", value: b) }
        return bs_beta_l(a, b)
    }
    
    @inlinable static func incompleteBetaUnnormalized(_ a: Float80, _ b: Float80, x: Float80) throws -> Float80 {
        guard a.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a", value: a) }
        guard b.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "b", value: b) }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a", value: a) }
        guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b", value: b) }
        guard x >= 0 && x <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 0.0, max: 1.0) }
        return bs_fullBeta_l(a, b, x)
    }

    @inlinable static func regularizedIncompleteBeta(_ a: Float80, _ b: Float80, x: Float80) throws -> Float80 {
        guard a.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a", value: a) }
        guard b.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "b", value: b) }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a", value: a) }
        guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b", value: b) }
        guard x >= 0 && x <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 0.0, max: 1.0) }
        return bs_ibeta_l(a, b, x)
    }

    @inlinable static func complementaryRegularizedIncompleteBeta(_ a: Float80, _ b: Float80, x: Float80) throws -> Float80 {
        guard a.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a", value: a) }
        guard b.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "b", value: b) }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a", value: a) }
        guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b", value: b) }
        guard x >= 0 && x <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 0.0, max: 1.0) }
        return bs_ibetac_l(a, b, x)
    }

    @inlinable static func inverseRegularizedIncompleteBeta(_ a: Float80, _ b: Float80, p: Float80) throws -> Float80 {
        guard a.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a", value: a) }
        guard b.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "b", value: b) }
        guard p.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "p", value: p) }
        guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a", value: a) }
        guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b", value: b) }
        guard p >= 0 && p <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "p", min: 0.0, max: 1.0) }
        return bs_ibeta_inv_l(a, b, p)
    }
    
    @inlinable static func inverseComplementaryRegularizedIncompleteBeta(_ a: Float80, _ b: Float80, p: Float80) throws -> Float80 {
        guard a.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a", value: a) }
        guard b.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "b", value: b) }
        guard p.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "p", value: p) }
        guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a", value: a) }
        guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b", value: b) }
        guard p >= 0 && p <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "p", min: 0.0, max: 1.0) }
        return bs_ibetac_inv_l(a, b, p)
    }

    @inlinable static func solveAForRegularizedIncompleteBeta(b: Float80, x: Float80, p: Float80) -> Float80 {
        bs_ibeta_inva_l(b, x, p)
    }

    @inlinable static func solveBForRegularizedIncompleteBeta(a: Float80, x: Float80, p: Float80) -> Float80 {
        bs_ibeta_invb_l(a, x, p)
    }

    @inlinable static func regularizedIncompleteBetaDerivative(_ a: Float80, _ b: Float80, x: Float80) throws -> Float80 {
        guard a.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a", value: a) }
        guard b.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "b", value: b) }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a", value: a) }
        guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b", value: b) }
        guard x >= 0 && x <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 0.0, max: 1.0) }
        return bs_ibeta_derivative_l(a, b, x)
    }
#endif
}
