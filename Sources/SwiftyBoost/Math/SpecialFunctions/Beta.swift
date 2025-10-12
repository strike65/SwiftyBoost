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

import CBoostBridge

public extension SpecialFunctions {
    
    
    // MARK: - Beta family (generic Double-backed)
    
    /// Compute the complete Beta function B(a, b) = Γ(a)Γ(b) / Γ(a + b).
    ///
    /// Domain:
    /// - Requires a > 0 and b > 0.
    ///
    /// Behavior:
    /// - Converts inputs to Double, evaluates via Boost.Math, and converts back to T.
    ///
    /// Parameters:
    /// - a: First shape parameter (must be > 0).
    /// - b: Second shape parameter (must be > 0).
    ///
    /// Returns:
    /// - B(a, b) as T.
    ///
    /// Throws:
    /// - SpecialFunctionError.parameterNotPositive(name: "a") if a ≤ 0.
    /// - SpecialFunctionError.parameterNotPositive(name: "b") if b ≤ 0.
    ///
    /// See also:
    /// - NIST DLMF §8.1–8.4.
    @inlinable static func beta<T: BinaryFloatingPoint>(_ a: T, _ b: T) throws -> T {
        guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
        guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b") }
        return T(bs_beta(D(a), D(b)))
    }
    
    /// Compute the unnormalized incomplete Beta B_x(a, b) = ∫₀ˣ t^{a−1} (1−t)^{b−1} dt.
    ///
    /// Domain:
    /// - Requires a > 0, b > 0, and x ∈ [0, 1].
    ///
    /// Behavior:
    /// - Converts to Double, evaluates via Boost.Math, converts back to T.
    ///
    /// Parameters:
    /// - a: First shape parameter (must be > 0).
    /// - b: Second shape parameter (must be > 0).
    /// - x: Upper integration limit (must be within [0, 1]).
    ///
    /// Returns:
    /// - B_x(a, b) as T.
    ///
    /// Throws:
    /// - SpecialFunctionError.parameterNotPositive(name: "a") if a ≤ 0.
    /// - SpecialFunctionError.parameterNotPositive(name: "b") if b ≤ 0.
    /// - SpecialFunctionError.parameterOutOfRange(name: "x", min: 0, max: 1) if x ∉ [0, 1].
    ///
    /// See also:
    /// - NIST DLMF §8.17(i).
    @inlinable static func incompleteBetaUnnormalized<T: BinaryFloatingPoint>(_ a: T, _ b: T, x: T) throws -> T {
        guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
        guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b") }
        guard x >= 0 && x <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 0.0, max: 1.0) }
        return T(bs_fullBeta(D(a), D(b), D(x)))
    }
    
    /// Compute the regularized incomplete Beta I_x(a, b) = B_x(a, b) / B(a, b).
    ///
    /// Domain:
    /// - Requires a > 0, b > 0, and x ∈ [0, 1].
    ///
    /// Behavior:
    /// - Converts to Double, evaluates via Boost.Math, converts back to T.
    ///
    /// Parameters:
    /// - a: First shape parameter (must be > 0).
    /// - b: Second shape parameter (must be > 0).
    /// - x: Upper integration limit (must be within [0, 1]).
    ///
    /// Returns:
    /// - I_x(a, b) as T within [0, 1].
    ///
    /// Throws:
    /// - SpecialFunctionError.parameterNotPositive(name: "a") if a ≤ 0.
    /// - SpecialFunctionError.parameterNotPositive(name: "b") if b ≤ 0.
    /// - SpecialFunctionError.parameterOutOfRange(name: "x", min: 0, max: 1) if x ∉ [0, 1].
    ///
    /// See also:
    /// - NIST DLMF §8.17(ii).
    @inlinable static func regularizedIncompleteBeta<T: BinaryFloatingPoint>(_ a: T, _ b: T, x: T) throws -> T {
        guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
        guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b") }
        guard x >= 0 && x <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 0.0, max: 1.0) }
        return T(bs_ibeta(D(a), D(b), D(x)))
    }
    
    /// Compute the complementary regularized incomplete Beta 1 − I_x(a, b).
    ///
    /// Domain:
    /// - Requires a > 0, b > 0, and x ∈ [0, 1].
    ///
    /// Parameters:
    /// - a: First shape parameter (must be > 0).
    /// - b: Second shape parameter (must be > 0).
    /// - x: Upper integration limit (must be within [0, 1]).
    ///
    /// Returns:
    /// - 1 − I_x(a, b) as T within [0, 1].
    ///
    /// Throws:
    /// - SpecialFunctionError.parameterNotPositive(name: "a") if a ≤ 0.
    /// - SpecialFunctionError.parameterNotPositive(name: "b") if b ≤ 0.
    /// - SpecialFunctionError.parameterOutOfRange(name: "x", min: 0, max: 1) if x ∉ [0, 1].
    @inlinable static func complementaryRegularizedIncompleteBeta<T: BinaryFloatingPoint>(_ a: T, _ b: T, x: T) throws -> T {
        guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
        guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b") }
        guard x >= 0 && x <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 0.0, max: 1.0) }
        return T(bs_ibetac(D(a), D(b), D(x)))
    }
    
    /// Invert the regularized incomplete Beta: find x such that I_x(a, b) = p.
    ///
    /// Domain:
    /// - Requires a > 0, b > 0, and p ∈ [0, 1].
    ///
    /// Parameters:
    /// - a: First shape parameter (must be > 0).
    /// - b: Second shape parameter (must be > 0).
    /// - p: Probability in [0, 1].
    ///
    /// Returns:
    /// - x ∈ [0, 1] such that I_x(a, b) = p.
    ///
    /// Throws:
    /// - SpecialFunctionError.parameterNotPositive(name: "a") if a ≤ 0.
    /// - SpecialFunctionError.parameterNotPositive(name: "b") if b ≤ 0.
    /// - SpecialFunctionError.parameterOutOfRange(name: "p", min: 0, max: 1) if p ∉ [0, 1].
    @inlinable static func inverseRegularizedIncompleteBeta<T: BinaryFloatingPoint>(_ a: T, _ b: T, p: T) throws -> T {
        guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
        guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b") }
        guard p >= 0 && p <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "p", min: 0.0, max: 1.0) }
        return T(bs_ibeta_inv(D(a), D(b), D(p)))
    }
    
    /// Invert the complementary regularized incomplete Beta: find x such that 1 − I_x(a, b) = p.
    ///
    /// Domain:
    /// - Requires a > 0, b > 0, and p ∈ [0, 1].
    ///
    /// Parameters:
    /// - a: First shape parameter (must be > 0).
    /// - b: Second shape parameter (must be > 0).
    /// - p: Probability in [0, 1].
    ///
    /// Returns:
    /// - x ∈ [0, 1] such that 1 − I_x(a, b) = p.
    ///
    /// Throws:
    /// - SpecialFunctionError.parameterNotPositive(name: "a") if a ≤ 0.
    /// - SpecialFunctionError.parameterNotPositive(name: "b") if b ≤ 0.
    /// - SpecialFunctionError.parameterOutOfRange(name: "p", min: 0, max: 1) if p ∉ [0, 1].
    @inlinable static func inverseComplementaryRegularizedIncompleteBeta<T: BinaryFloatingPoint>(_ a: T, _ b: T, p: T) throws -> T {
        guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
        guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b") }
        guard p >= 0 && p <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "p", min: 0.0, max: 1.0) }
        return T(bs_ibetac_inv(D(a), D(b), D(p)))
    }
    
    /// Solve for a in I_x(a, b) = p, given b, x, p.
    ///
    /// Notes:
    /// - This utility does not throw and does not perform Swift-side validation.
    ///   Callers should ensure b > 0, x ∈ [0, 1], p ∈ [0, 1].
    /// - On invalid inputs or if no solution exists, the backend may return NaN.
    ///
    /// Parameters:
    /// - b: Second shape parameter (expected > 0).
    /// - x: Upper integration limit (expected in [0, 1]).
    /// - p: Probability (expected in [0, 1]).
    ///
    /// Returns:
    /// - The solved a as T (may be NaN if inputs are invalid or the solver fails).
    @inlinable static func solveAForRegularizedIncompleteBeta<T: BinaryFloatingPoint>(b: T, x: T, p: T) -> T {
        T(bs_ibeta_inva(D(b), D(x), D(p)))
    }
    
    /// Solve for b in I_x(a, b) = p, given a, x, p.
    ///
    /// Notes:
    /// - This utility does not throw and does not perform Swift-side validation.
    ///   Callers should ensure a > 0, x ∈ [0, 1], p ∈ [0, 1].
    /// - On invalid inputs or if no solution exists, the backend may return NaN.
    ///
    /// Parameters:
    /// - a: First shape parameter (expected > 0).
    /// - x: Upper integration limit (expected in [0, 1]).
    /// - p: Probability (expected in [0, 1]).
    ///
    /// Returns:
    /// - The solved b as T (may be NaN if inputs are invalid or the solver fails).
    @inlinable static func solveBForRegularizedIncompleteBeta<T: BinaryFloatingPoint>(a: T, x: T, p: T) -> T {
        T(bs_ibeta_invb(D(a), D(x), D(p)))
    }
    
    // MARK: - Derivative of the regularized incomplete Beta (generic)
    
    /// Derivative with respect to x of the regularized incomplete Beta I_x(a, b).
    ///
    /// Definition:
    /// - For a > 0, b > 0 and x ∈ [0, 1], the derivative is:
    ///   d/dx I_x(a, b) = x^{a−1} (1 − x)^{b−1} / B(a, b).
    ///
    /// Properties:
    /// - Nonnegative on [0, 1] for all valid a, b.
    /// - Integrates to I_x(a, b) up to a constant of integration.
    /// - This derivative equals the probability density function of the Beta(a, b)
    ///   distribution evaluated at x.
    ///
    /// Domain:
    /// - Requires a > 0, b > 0, and x ∈ [0, 1].
    ///
    /// Numeric behavior:
    /// - Delegates to Boost.Math’s ibeta_derivative for robust evaluation across
    ///   the parameter space, including extreme a, b and x near {0, 1}.
    ///
    /// Parameters:
    /// - a: First shape parameter (strictly positive).
    /// - b: Second shape parameter (strictly positive).
    /// - x: Evaluation point (in [0, 1]).
    ///
    /// Returns:
    /// - The value of d/dx I_x(a, b) as T.
    ///
    /// Throws:
    /// - SpecialFunctionError.parameterNotPositive(name: "a") if a ≤ 0.
    /// - SpecialFunctionError.parameterNotPositive(name: "b") if b ≤ 0.
    /// - SpecialFunctionError.parameterOutOfRange(name: "x", min: 0, max: 1) if x ∉ [0, 1].
    ///
    /// SeeAlso:
    /// - regularizedIncompleteBeta(_:_:x:) for I_x(a, b).
    /// - beta(_:_) for B(a, b).
    @inlinable static func regularizedIncompleteBetaDerivative<T: BinaryFloatingPoint>(_ a: T, _ b: T, x: T) throws -> T {
        guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
        guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b") }
        guard x >= 0 && x <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 0.0, max: 1.0) }
        return T(bs_ibeta_derivative(D(a), D(b), D(x)))
    }
    
    // MARK: - Float overloads (Beta)
    // Direct Float-precision entry points that avoid generic conversions and call
    // the corresponding C implementations for speed.
    
    /// Complete Beta B(a, b) for Float. Requires a > 0 and b > 0.
    @inlinable static func beta(_ a: Float, _ b: Float) throws -> Float {
        guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
        guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b") }
        return bs_beta_f(a, b)
    }
    
    /// Unnormalized incomplete Beta B_x(a, b) for Float. Requires a > 0, b > 0, x ∈ [0, 1].
    @inlinable static func incompleteBetaUnnormalized(_ a: Float, _ b: Float, x: Float) throws -> Float {
        guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
        guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b") }
        guard x >= 0 && x <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 0.0, max: 1.0) }
        return bs_fullBeta_f(a, b, x)
    }
    
    /// Regularized incomplete Beta I_x(a, b) for Float. Requires a > 0, b > 0, x ∈ [0, 1].
    @inlinable static func regularizedIncompleteBeta(_ a: Float, _ b: Float, x: Float) throws -> Float {
        guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
        guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b") }
        guard x >= 0 && x <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 0.0, max: 1.0) }
        return bs_ibeta_f(a, b, x)
    }
    
    /// Complementary regularized incomplete Beta 1 − I_x(a, b) for Float. Requires a > 0, b > 0, x ∈ [0, 1].
    @inlinable static func complementaryRegularizedIncompleteBeta(_ a: Float, _ b: Float, x: Float) throws -> Float {
        guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
        guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b") }
        guard x >= 0 && x <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 0.0, max: 1.0) }
        return bs_ibetac_f(a, b, x)
    }
    
    /// Inverse regularized incomplete Beta for Float: solve x in I_x(a, b) = p. Requires a > 0, b > 0, p ∈ [0, 1].
    @inlinable static func inverseRegularizedIncompleteBeta(_ a: Float, _ b: Float, p: Float) throws -> Float {
        guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
        guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b") }
        guard p >= 0 && p <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "p", min: 0.0, max: 1.0) }
        return bs_ibeta_inv_f(a, b, p)
    }
    
    /// Inverse complementary regularized incomplete Beta for Float: solve x in 1 − I_x(a, b) = p. Requires a > 0, b > 0, p ∈ [0, 1].
    @inlinable static func inverseComplementaryRegularizedIncompleteBeta(_ a: Float, _ b: Float, p: Float) throws -> Float {
        guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
        guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b") }
        guard p >= 0 && p <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "p", min: 0.0, max: 1.0) }
        return bs_ibetac_inv_f(a, b, p)
    }
    
    /// Solve for a in I_x(a, b) = p for Float. No throws; see notes above for expected domains.
    @inlinable static func solveAForRegularizedIncompleteBeta(b: Float, x: Float, p: Float) -> Float {
        bs_ibeta_inva_f(b, x, p)
    }
    
    /// Solve for b in I_x(a, b) = p for Float. No throws; see notes above for expected domains.
    @inlinable static func solveBForRegularizedIncompleteBeta(a: Float, x: Float, p: Float) -> Float {
        bs_ibeta_invb_f(a, x, p)
    }
    
    /// Derivative d/dx I_x(a, b) for Float.
    ///
    /// See regularizedIncompleteBetaDerivative(_:_:x:) for definition, domain, and discussion.
    @inlinable static func regularizedIncompleteBetaDerivative(_ a: Float, _ b: Float, x: Float) throws -> Float {
        guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
        guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b") }
        guard x >= 0 && x <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 0.0, max: 1.0) }
        return bs_ibeta_derivative_f(a, b, x)
    }
    
    // MARK: - Float80 overloads (x86_64)
    // Extended-precision versions for platforms that support Float80.
    
#if arch(x86_64)
    /// Complete Beta B(a, b) for Float80 (x86_64 only). Requires a > 0 and b > 0.
    @inlinable static func beta(_ a: Float80, _ b: Float80) throws -> Float80 {
        guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
        guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b") }
        return bs_beta_l(a, b)
    }
    
    /// Unnormalized incomplete Beta B_x(a, b) for Float80 (x86_64 only). Requires a > 0, b > 0, x ∈ [0, 1].
    @inlinable static func incompleteBetaUnnormalized(_ a: Float80, _ b: Float80, x: Float80) throws -> Float80 {
        guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
        guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b") }
        guard x >= 0 && x <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 0.0, max: 1.0) }
        return bs_fullBeta_l(a, b, x)
    }
    
    /// Regularized incomplete Beta I_x(a, b) for Float80 (x86_64 only). Requires a > 0, b > 0, x ∈ [0, 1].
    @inlinable static func regularizedIncompleteBeta(_ a: Float80, _ b: Float80, x: Float80) throws -> Float80 {
        guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
        guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b") }
        guard x >= 0 && x <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 0.0, max: 1.0) }
        return bs_ibeta_l(a, b, x)
    }
    
    /// Complementary regularized incomplete Beta 1 − I_x(a, b) for Float80 (x86_64 only). Requires a > 0, b > 0, x ∈ [0, 1].
    @inlinable static func complementaryRegularizedIncompleteBeta(_ a: Float80, _ b: Float80, x: Float80) throws -> Float80 {
        guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
        guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b") }
        guard x >= 0 && x <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 0.0, max: 1.0) }
        return bs_ibetac_l(a, b, x)
    }
    
    /// Inverse regularized incomplete Beta for Float80 (x86_64 only): solve x in I_x(a, b) = p. Requires a > 0, b > 0, p ∈ [0, 1].
    @inlinable static func inverseRegularizedIncompleteBeta(_ a: Float80, _ b: Float80, p: Float80) throws -> Float80 {
        guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
        guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b") }
        guard p >= 0 && p <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "p", min: 0.0, max: 1.0) }
        return bs_ibeta_inv_l(a, b, p)
    }
    
    /// Inverse complementary regularized incomplete Beta for Float80 (x86_64 only): solve x in 1 − I_x(a, b) = p. Requires a > 0, b > 0, p ∈ [0, 1].
    @inlinable static func inverseComplementaryRegularizedIncompleteBeta(_ a: Float80, _ b: Float80, p: Float80) throws -> Float80 {
        guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
        guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b") }
        guard p >= 0 && p <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "p", min: 0.0, max: 1.0) }
        return bs_ibetac_inv_l(a, b, p)
    }
    
    /// Solve for a in I_x(a, b) = p for Float80 (x86_64 only). No throws; see notes above for expected domains.
    @inlinable static func solveAForRegularizedIncompleteBeta(b: Float80, x: Float80, p: Float80) -> Float80 {
        bs_ibeta_inva_l(b, x, p)
    }
    
    /// Solve for b in I_x(a, b) = p for Float80 (x86_64 only). No throws; see notes above for expected domains.
    @inlinable static func solveBForRegularizedIncompleteBeta(a: Float80, x: Float80, p: Float80) -> Float80 {
        bs_ibeta_invb_l(a, x, p)
    }
    
    /// Derivative d/dx I_x(a, b) for Float80 (x86_64 only).
    ///
    /// See regularizedIncompleteBetaDerivative(_:_:x:) for definition, domain, and discussion.
    @inlinable static func regularizedIncompleteBetaDerivative(_ a: Float80, _ b: Float80, x: Float80) throws -> Float80 {
        guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
        guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b") }
        guard x >= 0 && x <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 0.0, max: 1.0) }
        return bs_ibeta_derivative_l(a, b, x)
    }
#endif
}

