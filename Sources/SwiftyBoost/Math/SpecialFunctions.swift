//
//  Created by VT on 11.10.25.
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

import CBoostBridge

// Errors for invalid inputs to special functions
public enum SpecialFunctionError: Error {
    case parameterNotPositive(name: String)
    case parameterOutOfRange(name: String, min: Double, max: Double)
    case parameterNotFinite(name: String)
    case poleAtNonPositiveInteger(name: String)
    case invalidCombination(message: String)
}

// Helper casts
@usableFromInline internal func D<T: BinaryFloatingPoint>(_ x: T) -> Double { Double(x) }

/// Base of the natural logarithm e ≈ 2.71828, provided by Boost.Math constants.
///
/// This is returned as `Double`. See also the Float and Float80 constant wrappers exposed via the C interface if you need those types.
///
/// - Returns: The constant e as `Double`.
@inlinable public var boostE: Double { bs_const_e() }

// MARK: - Gamma / Error (generic Double-backed)

/// Gamma function Γ(x).
///
/// Uses Boost.Math’s `tgamma`. This is the analytic continuation of the factorial to real values.
/// Poles occur at non-positive integers (0, −1, −2, ...).
///
/// - Parameter x: Real input.
/// - Returns: Γ(x).
/// - Throws: `SpecialFunctionError.poleAtNonPositiveInteger` if `x` is a non-positive integer; `SpecialFunctionError.parameterNotFinite` if `x` is not finite.
/// - SeeAlso: `logGamma(_:)`
@inlinable public func gamma<T: BinaryFloatingPoint>(_ x: T) throws -> T {
    let dx = D(x)
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    if dx <= 0, dx == dx.rounded(.towardZero) {
        throw SpecialFunctionError.poleAtNonPositiveInteger(name: "x")
    }
    return T(bs_tgamma(dx))
}

/// Logarithm of the gamma function ln Γ(x).
///
/// Uses Boost.Math’s `lgamma`. Poles occur at non-positive integers (0, −1, −2, ...).
///
/// - Parameter x: Real input.
/// - Returns: ln Γ(x).
/// - Throws: `SpecialFunctionError.poleAtNonPositiveInteger` if `x` is a non-positive integer; `SpecialFunctionError.parameterNotFinite` if `x` is not finite.
/// - SeeAlso: `gamma(_:)`
@inlinable public func logGamma<T: BinaryFloatingPoint>(_ x: T) throws -> T {
    let dx = D(x)
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    if dx <= 0, dx == dx.rounded(.towardZero) {
        throw SpecialFunctionError.poleAtNonPositiveInteger(name: "x")
    }
    return T(bs_lgamma(dx))
}

/// Error function erf(x).
///
/// Uses Boost.Math’s `erf`.
///
/// - Parameter x: Real input.
/// - Returns: erf(x).
/// - Throws: `SpecialFunctionError.parameterNotFinite` if `x` is not finite.
/// - SeeAlso: `complementaryErrorFunction(_:)`
@inlinable public func errorFunction<T: BinaryFloatingPoint>(_ x: T) throws -> T {
    let dx = D(x)
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    return T(bs_erf(dx))
}

/// Complementary error function erfc(x) = 1 − erf(x).
///
/// Uses Boost.Math’s `erfc`.
///
/// - Parameter x: Real input.
/// - Returns: erfc(x).
/// - Throws: `SpecialFunctionError.parameterNotFinite` if `x` is not finite.
/// - SeeAlso: `errorFunction(_:)`
@inlinable public func complementaryErrorFunction<T: BinaryFloatingPoint>(_ x: T) throws -> T {
    let dx = D(x)
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    return T(bs_erfc(dx))
}

// MARK: - Float overloads

/// Float overload of `gamma(_:)`.
@inlinable public func gamma(_ x: Float) throws -> Float {
    guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    if x <= 0, x == x.rounded(.towardZero) {
        throw SpecialFunctionError.poleAtNonPositiveInteger(name: "x")
    }
    return bs_tgamma_f(x)
}

/// Float overload of `logGamma(_:)`.
@inlinable public func logGamma(_ x: Float) throws -> Float {
    guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    if x <= 0, x == x.rounded(.towardZero) {
        throw SpecialFunctionError.poleAtNonPositiveInteger(name: "x")
    }
    return bs_lgamma_f(x)
}

/// Float overload of `errorFunction(_:)`.
@inlinable public func errorFunction(_ x: Float) throws -> Float {
    guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    return bs_erf_f(x)
}

/// Float overload of `complementaryErrorFunction(_:)`.
@inlinable public func complementaryErrorFunction(_ x: Float) throws -> Float {
    guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    return bs_erfc_f(x)
}

// MARK: - Float80 overloads (x86_64 only)

/// Float80 overload of `gamma(_:)`.
#if arch(x86_64)
@inlinable public func gamma(_ x: Float80) throws -> Float80 {
    guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    if x <= 0, x == x.rounded(.towardZero) {
        throw SpecialFunctionError.poleAtNonPositiveInteger(name: "x")
    }
    return bs_tgamma_l(x)
}
#endif

/// Float80 overload of `logGamma(_:)`.
#if arch(x86_64)
@inlinable public func logGamma(_ x: Float80) throws -> Float80 {
    guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    if x <= 0, x == x.rounded(.towardZero) {
        throw SpecialFunctionError.poleAtNonPositiveInteger(name: "x")
    }
    return bs_lgamma_l(x)
}
#endif

/// Float80 overload of `errorFunction(_:)`.
#if arch(x86_64)
@inlinable public func errorFunction(_ x: Float80) throws -> Float80 {
    guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    return bs_erf_l(x)
}
#endif

/// Float80 overload of `complementaryErrorFunction(_:)`.
#if arch(x86_64)
@inlinable public func complementaryErrorFunction(_ x: Float80) throws -> Float80 {
    guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    return bs_erfc_l(x)
}
#endif

// MARK: - Beta family (generic Double-backed)

/// Euler’s beta function B(a, b) = Γ(a)Γ(b)/Γ(a+b).
///
/// Uses Boost.Math’s `beta`.
///
/// - Parameters:
///   - a: First shape parameter, must be > 0.
///   - b: Second shape parameter, must be > 0.
/// - Returns: B(a, b).
/// - Throws: `SpecialFunctionError.parameterNotPositive` if `a` or `b` ≤ 0.
@inlinable public func beta<T: BinaryFloatingPoint>(_ a: T, _ b: T) throws -> T {
    guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
    guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b") }
    return T(bs_beta(D(a), D(b)))
}

/// Unnormalized incomplete beta B(a, b; x) = ∫₀ˣ t^{a−1}(1−t)^{b−1} dt.
///
/// Uses Boost.Math’s `beta(a,b,x)`.
///
/// - Parameters:
///   - a: First shape parameter, must be > 0.
///   - b: Second shape parameter, must be > 0.
///   - x: Integration limit in [0, 1].
/// - Returns: B(a, b; x).
/// - Throws: `SpecialFunctionError.parameterNotPositive` if `a` or `b` ≤ 0; `SpecialFunctionError.parameterOutOfRange` if `x ∉ [0,1]`.
@inlinable public func incompleteBetaUnnormalized<T: BinaryFloatingPoint>(_ a: T, _ b: T, x: T) throws -> T {
    guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
    guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b") }
    guard x >= 0 && x <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 0.0, max: 1.0) }
    return T(bs_fullBeta(D(a), D(b), D(x)))
}

/// Regularized incomplete beta Iₓ(a, b).
///
/// Uses Boost.Math’s `ibeta`.
///
/// - Parameters:
///   - a: First shape parameter, must be > 0.
///   - b: Second shape parameter, must be > 0.
///   - x: Evaluation point in [0, 1].
/// - Returns: Iₓ(a, b).
/// - Throws: `SpecialFunctionError.parameterNotPositive` if `a` or `b` ≤ 0; `SpecialFunctionError.parameterOutOfRange` if `x ∉ [0,1]`.
@inlinable public func regularizedIncompleteBeta<T: BinaryFloatingPoint>(_ a: T, _ b: T, x: T) throws -> T {
    guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
    guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b") }
    guard x >= 0 && x <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 0.0, max: 1.0) }
    return T(bs_ibeta(D(a), D(b), D(x)))
}

/// Complement of regularized incomplete beta I_{1−x}(b, a).
///
/// Uses Boost.Math’s `ibetac`.
///
/// - Parameters:
///   - a: First shape parameter, must be > 0.
///   - b: Second shape parameter, must be > 0.
///   - x: Evaluation point in [0, 1].
/// - Returns: I_{1−x}(b, a).
/// - Throws: `SpecialFunctionError.parameterNotPositive` if `a` or `b` ≤ 0; `SpecialFunctionError.parameterOutOfRange` if `x ∉ [0,1]`.
@inlinable public func complementaryRegularizedIncompleteBeta<T: BinaryFloatingPoint>(_ a: T, _ b: T, x: T) throws -> T {
    guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
    guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b") }
    guard x >= 0 && x <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 0.0, max: 1.0) }
    return T(bs_ibetac(D(a), D(b), D(x)))
}

/// Inverse of regularized incomplete beta: solve for x in Iₓ(a, b) = p.
///
/// Uses Boost.Math’s `ibeta_inv`.
///
/// - Parameters:
///   - a: First shape parameter, must be > 0.
///   - b: Second shape parameter, must be > 0.
///   - p: Probability in [0, 1].
/// - Returns: The solution x ∈ [0, 1].
/// - Throws: `SpecialFunctionError.parameterNotPositive` if `a` or `b` ≤ 0; `SpecialFunctionError.parameterOutOfRange` if `p ∉ [0,1]`.
@inlinable public func inverseRegularizedIncompleteBeta<T: BinaryFloatingPoint>(_ a: T, _ b: T, p: T) throws -> T {
    guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
    guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b") }
    guard p >= 0 && p <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "p", min: 0.0, max: 1.0) }
    return T(bs_ibeta_inv(D(a), D(b), D(p)))
}

/// Inverse of the complement: solve for x in I_{1−x}(a, b) = p.
///
/// Uses Boost.Math’s `ibetac_inv`.
///
/// - Parameters:
///   - a: First shape parameter, must be > 0.
///   - b: Second shape parameter, must be > 0.
///   - p: Probability in [0, 1].
/// - Returns: The solution x ∈ [0, 1].
/// - Throws: `SpecialFunctionError.parameterNotPositive` if `a` or `b` ≤ 0; `SpecialFunctionError.parameterOutOfRange` if `p ∉ [0,1]`.
@inlinable public func inverseComplementaryRegularizedIncompleteBeta<T: BinaryFloatingPoint>(_ a: T, _ b: T, p: T) throws -> T {
    guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
    guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b") }
    guard p >= 0 && p <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "p", min: 0.0, max: 1.0) }
    return T(bs_ibetac_inv(D(a), D(b), D(p)))
}

/// Solve for parameter a in Iₓ(a, b) = p, given b, x, p.
///
/// Uses Boost.Math’s `ibeta_inva`.
///
/// - Parameters:
///   - b: Second shape parameter (> 0).
///   - x: Evaluation point in [0, 1].
///   - p: Probability in [0, 1].
/// - Returns: The solution a.
/// - Note: No explicit validation beyond domain checks in the C/Boost layer is performed for numerical feasibility when multiple solutions exist.
@inlinable public func solveAForRegularizedIncompleteBeta<T: BinaryFloatingPoint>(b: T, x: T, p: T) -> T {
    T(bs_ibeta_inva(D(b), D(x), D(p)))
}

/// Solve for parameter b in Iₓ(a, b) = p, given a, x, p.
///
/// Uses Boost.Math’s `ibeta_invb`.
///
/// - Parameters:
///   - a: First shape parameter (> 0).
///   - x: Evaluation point in [0, 1].
///   - p: Probability in [0, 1].
/// - Returns: The solution b.
/// - Note: No explicit validation beyond domain checks in the C/Boost layer is performed for numerical feasibility when multiple solutions exist.
@inlinable public func solveBForRegularizedIncompleteBeta<T: BinaryFloatingPoint>(a: T, x: T, p: T) -> T {
    T(bs_ibeta_invb(D(a), D(x), D(p)))
}

// MARK: - Float overloads (Beta)

/// Float overload of `beta(_: _:)`.
@inlinable public func beta(_ a: Float, _ b: Float) throws -> Float {
    guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
    guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b") }
    return bs_beta_f(a, b)
}

/// Float overload of `incompleteBetaUnnormalized(_:_:x:)`.
@inlinable public func incompleteBetaUnnormalized(_ a: Float, _ b: Float, x: Float) throws -> Float {
    guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
    guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b") }
    guard x >= 0 && x <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 0.0, max: 1.0) }
    return bs_fullBeta_f(a, b, x)
}

/// Float overload of `regularizedIncompleteBeta(_:_:x:)`.
@inlinable public func regularizedIncompleteBeta(_ a: Float, _ b: Float, x: Float) throws -> Float {
    guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
    guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b") }
    guard x >= 0 && x <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 0.0, max: 1.0) }
    return bs_ibeta_f(a, b, x)
}

/// Float overload of `complementaryRegularizedIncompleteBeta(_:_:x:)`.
@inlinable public func complementaryRegularizedIncompleteBeta(_ a: Float, _ b: Float, x: Float) throws -> Float {
    guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
    guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b") }
    guard x >= 0 && x <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 0.0, max: 1.0) }
    return bs_ibetac_f(a, b, x)
}

/// Float overload of `inverseRegularizedIncompleteBeta(_:_:p:)`.
@inlinable public func inverseRegularizedIncompleteBeta(_ a: Float, _ b: Float, p: Float) throws -> Float {
    guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
    guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b") }
    guard p >= 0 && p <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "p", min: 0.0, max: 1.0) }
    return bs_ibeta_inv_f(a, b, p)
}

/// Float overload of `inverseComplementaryRegularizedIncompleteBeta(_:_:p:)`.
@inlinable public func inverseComplementaryRegularizedIncompleteBeta(_ a: Float, _ b: Float, p: Float) throws -> Float {
    guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
    guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b") }
    guard p >= 0 && p <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "p", min: 0.0, max: 1.0) }
    return bs_ibetac_inv_f(a, b, p)
}

/// Float overload of `solveAForRegularizedIncompleteBeta(b:x:p:)`.
@inlinable public func solveAForRegularizedIncompleteBeta(b: Float, x: Float, p: Float) -> Float {
    bs_ibeta_inva_f(b, x, p)
}

/// Float overload of `solveBForRegularizedIncompleteBeta(a:x:p:)`.
@inlinable public func solveBForRegularizedIncompleteBeta(a: Float, x: Float, p: Float) -> Float {
    bs_ibeta_invb_f(a, x, p)
}

// MARK: - Float80 overloads (Beta) — x86_64 only

/// Float80 overload of `beta(_: _:)`.
#if arch(x86_64)
@inlinable public func beta(_ a: Float80, _ b: Float80) throws -> Float80 {
    guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
    guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b") }
    return bs_beta_l(a, b)
}
#endif

/// Float80 overload of `incompleteBetaUnnormalized(_:_:x:)`.
#if arch(x86_64)
@inlinable public func incompleteBetaUnnormalized(_ a: Float80, _ b: Float80, x: Float80) throws -> Float80 {
    guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
    guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b") }
    guard x >= 0 && x <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 0.0, max: 1.0) }
    return bs_fullBeta_l(a, b, x)
}
#endif

/// Float80 overload of `regularizedIncompleteBeta(_:_:x:)`.
#if arch(x86_64)
@inlinable public func regularizedIncompleteBeta(_ a: Float80, _ b: Float80, x: Float80) throws -> Float80 {
    guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
    guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b") }
    guard x >= 0 && x <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 0.0, max: 1.0) }
    return bs_ibeta_l(a, b, x)
}
#endif

/// Float80 overload of `complementaryRegularizedIncompleteBeta(_:_:x:)`.
#if arch(x86_64)
@inlinable public func complementaryRegularizedIncompleteBeta(_ a: Float80, _ b: Float80, x: Float80) throws -> Float80 {
    guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
    guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b") }
    guard x >= 0 && x <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 0.0, max: 1.0) }
    return bs_ibetac_l(a, b, x)
}
#endif

/// Float80 overload of `inverseRegularizedIncompleteBeta(_:_:p:)`.
#if arch(x86_64)
@inlinable public func inverseRegularizedIncompleteBeta(_ a: Float80, _ b: Float80, p: Float80) throws -> Float80 {
    guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
    guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b") }
    guard p >= 0 && p <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "p", min: 0.0, max: 1.0) }
    return bs_ibeta_inv_l(a, b, p)
}
#endif

/// Float80 overload of `inverseComplementaryRegularizedIncompleteBeta(_:_:p:)`.
#if arch(x86_64)
@inlinable public func inverseComplementaryRegularizedIncompleteBeta(_ a: Float80, _ b: Float80, p: Float80) throws -> Float80 {
    guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
    guard b > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "b") }
    guard p >= 0 && p <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "p", min: 0.0, max: 1.0) }
    return bs_ibetac_inv_l(a, b, p)
}
#endif

/// Float80 overload of `solveAForRegularizedIncompleteBeta(b:x:p:)`.
#if arch(x86_64)
@inlinable public func solveAForRegularizedIncompleteBeta(b: Float80, x: Float80, p: Float80) -> Float80 {
    bs_ibeta_inva_l(b, x, p)
}
#endif

/// Float80 overload of `solveBForRegularizedIncompleteBeta(a:x:p:)`.
#if arch(x86_64)
@inlinable public func solveBForRegularizedIncompleteBeta(a: Float80, x: Float80, p: Float80) -> Float80 {
    bs_ibeta_invb_l(a, x, p)
}
#endif

// MARK: - Digamma / Polygamma / Zeta (generic Double-backed)

/// Digamma ψ(x) = d/dx ln Γ(x).
///
/// Uses Boost.Math’s `digamma`. Poles at non-positive integers.
///
/// - Parameter x: Real input.
/// - Returns: ψ(x).
/// - Throws: `SpecialFunctionError.poleAtNonPositiveInteger` if `x` is a non-positive integer; `SpecialFunctionError.parameterNotFinite` if `x` is not finite.
@inlinable public func digamma<T: BinaryFloatingPoint>(_ x: T) throws -> T {
    let dx = D(x)
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    if dx <= 0, dx == dx.rounded(.towardZero) {
        throw SpecialFunctionError.poleAtNonPositiveInteger(name: "x")
    }
    return T(bs_digamma(dx))
}

/// Trigamma ψ₁(x) = d/dx ψ(x).
///
/// Uses Boost.Math’s `trigamma`. Poles at non-positive integers.
///
/// - Parameter x: Real input.
/// - Returns: ψ₁(x).
/// - Throws: `SpecialFunctionError.poleAtNonPositiveInteger` if `x` is a non-positive integer; `SpecialFunctionError.parameterNotFinite` if `x` is not finite.
@inlinable public func trigamma<T: BinaryFloatingPoint>(_ x: T) throws -> T {
    let dx = D(x)
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    if dx <= 0, dx == dx.rounded(.towardZero) {
        throw SpecialFunctionError.poleAtNonPositiveInteger(name: "x")
    }
    return T(bs_trigamma(dx))
}

/// Polygamma ψ⁽ⁿ⁾(x), the n-th derivative of digamma.
///
/// Uses Boost.Math’s `polygamma`.
///
/// - Parameters:
///   - n: Order n ≥ 0. `n = 0` corresponds to digamma.
///   - x: Real input.
/// - Returns: ψ⁽ⁿ⁾(x).
/// - Throws: `SpecialFunctionError.parameterNotPositive` if `n < 0`; `SpecialFunctionError.poleAtNonPositiveInteger` if `x` is a non-positive integer; `SpecialFunctionError.parameterNotFinite` if inputs are not finite.
@inlinable public func polygamma<T: BinaryFloatingPoint>(order n: Int, _ x: T) throws -> T {
    guard n >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "order") }
    let dx = D(x)
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    if dx <= 0, dx == dx.rounded(.towardZero) {
        throw SpecialFunctionError.poleAtNonPositiveInteger(name: "x")
    }
    return T(bs_polygamma(Int32(n), dx))
}

/// Riemann zeta ζ(x).
///
/// Uses Boost.Math’s `zeta`. Has a simple pole at x = 1.
///
/// - Parameter x: Real input.
/// - Returns: ζ(x).
/// - Throws: `SpecialFunctionError.invalidCombination` if `x == 1`; `SpecialFunctionError.parameterNotFinite` if `x` is not finite.
@inlinable public func riemannZeta<T: BinaryFloatingPoint>(_ x: T) throws -> T {
    let dx = D(x)
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    guard dx != 1 else { throw SpecialFunctionError.invalidCombination(message: "riemannZeta has a pole at x = 1") }
    return T(bs_riemann_zeta(dx))
}

// Float overloads (digamma/polygamma/zeta)

/// Float overload of `digamma(_:)`.
@inlinable public func digamma(_ x: Float) throws -> Float {
    guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    if x <= 0, x == x.rounded(.towardZero) { throw SpecialFunctionError.poleAtNonPositiveInteger(name: "x") }
    return bs_digamma_f(x)
}

/// Float overload of `trigamma(_:)`.
@inlinable public func trigamma(_ x: Float) throws -> Float {
    guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    if x <= 0, x == x.rounded(.towardZero) { throw SpecialFunctionError.poleAtNonPositiveInteger(name: "x") }
    return bs_trigamma_f(x)
}

/// Float overload of `polygamma(order:_:)`.
@inlinable public func polygamma(order n: Int, _ x: Float) throws -> Float {
    guard n >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "order") }
    guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    if x <= 0, x == x.rounded(.towardZero) { throw SpecialFunctionError.poleAtNonPositiveInteger(name: "x") }
    return bs_polygamma_f(Int32(n), x)
}

/// Float overload of `riemannZeta(_:)`.
@inlinable public func riemannZeta(_ x: Float) throws -> Float {
    guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    guard x != 1 else { throw SpecialFunctionError.invalidCombination(message: "riemannZeta has a pole at x = 1") }
    return bs_riemann_zeta_f(x)
}

// Float80 overloads (x86_64)

/// Float80 overload of `digamma(_:)`.
#if arch(x86_64)
@inlinable public func digamma(_ x: Float80) throws -> Float80 {
    guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    if x <= 0, x == x.rounded(.towardZero) { throw SpecialFunctionError.poleAtNonPositiveInteger(name: "x") }
    return bs_digamma_l(x)
}
#endif

/// Float80 overload of `trigamma(_:)`.
#if arch(x86_64)
@inlinable public func trigamma(_ x: Float80) throws -> Float80 {
    guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    if x <= 0, x == x.rounded(.towardZero) { throw SpecialFunctionError.poleAtNonPositiveInteger(name: "x") }
    return bs_trigamma_l(x)
}
#endif

/// Float80 overload of `polygamma(order:_:)`.
#if arch(x86_64)
@inlinable public func polygamma(order n: Int, _ x: Float80) throws -> Float80 {
    guard n >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "order") }
    guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    if x <= 0, x == x.rounded(.towardZero) { throw SpecialFunctionError.poleAtNonPositiveInteger(name: "x") }
    return bs_polygamma_l(Int32(n), x)
}
#endif

/// Float80 overload of `riemannZeta(_:)`.
#if arch(x86_64)
@inlinable public func riemannZeta(_ x: Float80) throws -> Float80 {
    guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    guard x != 1 else { throw SpecialFunctionError.invalidCombination(message: "riemannZeta has a pole at x = 1") }
    return bs_riemann_zeta_l(x)
}
#endif

// MARK: - Owen's T

/// Owen’s T function T(h, a).
///
/// Uses Boost.Math’s `owens_t`.
///
/// - Parameters:
///   - h: Real parameter.
///   - a: Real parameter.
/// - Returns: T(h, a).
/// - Throws: `SpecialFunctionError.parameterNotFinite` if inputs are not finite.
@inlinable public func owensT<T: BinaryFloatingPoint>(h: T, a: T) throws -> T {
    let dh = D(h), da = D(a)
    guard dh.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "h") }
    guard da.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a") }
    return T(bs_owens_t(dh, da))
}

/// Float overload of `owensT(h:a:)`.
@inlinable public func owensT(h: Float, a: Float) throws -> Float {
    guard h.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "h") }
    guard a.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a") }
    return bs_owens_t_f(h, a)
}

/// Float80 overload of `owensT(h:a:)`.
#if arch(x86_64)
@inlinable public func owensT(h: Float80, a: Float80) throws -> Float80 {
    guard h.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "h") }
    guard a.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a") }
    return bs_owens_t_l(h, a)
}
#endif

// MARK: - Exponential integrals and related

/// Generalized exponential integral Eₙ(n, x).
///
/// Uses Boost.Math’s `expint`.
///
/// - Parameters:
///   - n: Order n ≥ 0.
///   - x: Real input.
/// - Returns: Eₙ(n, x).
/// - Throws: `SpecialFunctionError.parameterNotPositive` if `n < 0`; `SpecialFunctionError.parameterNotFinite` if `x` is not finite.
@inlinable public func exponentialIntegralEn<T: BinaryFloatingPoint>(_ n: Int, _ x: T) throws -> T {
    guard n >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "n") }
    let dx = D(x)
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    return T(bs_expint_En(Int32(n), dx))
}

/// exp(x) − 1 computed with improved accuracy for small x.
///
/// Uses Boost.Math’s `expm1`.
///
/// - Parameter x: Real input.
/// - Returns: exp(x) − 1.
/// - Throws: `SpecialFunctionError.parameterNotFinite` if `x` is not finite.
@inlinable public func expm1<T: BinaryFloatingPoint>(_ x: T) throws -> T {
    let dx = D(x)
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    return T(bs_expm1(dx))
}

/// log(1 + x) computed with improved accuracy for small x.
///
/// Uses Boost.Math’s `log1p`.
///
/// - Parameter x: Real input (must be > −1).
/// - Returns: log(1 + x).
/// - Throws: `SpecialFunctionError.parameterNotFinite` if `x` is not finite; `SpecialFunctionError.parameterOutOfRange` if `x ≤ −1`.
@inlinable public func log1p<T: BinaryFloatingPoint>(_ x: T) throws -> T {
    let dx = D(x)
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    guard dx > -1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: -1.0.nextUp, max: Double.infinity) }
    return T(bs_log1p(dx))
}

/// log(1 + x) − x computed with improved accuracy for small x.
///
/// Uses Boost.Math’s `log1pmx`.
///
/// - Parameter x: Real input (must be > −1).
/// - Returns: log(1 + x) − x.
/// - Throws: `SpecialFunctionError.parameterNotFinite` if `x` is not finite; `SpecialFunctionError.parameterOutOfRange` if `x ≤ −1`.
@inlinable public func log1pmx<T: BinaryFloatingPoint>(_ x: T) throws -> T {
    let dx = D(x)
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    guard dx > -1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: -1.0.nextUp, max: Double.infinity) }
    return T(bs_log1pmx(dx))
}

/// xʸ − 1 computed with improved accuracy near y = 0.
///
/// Uses Boost.Math’s `powm1`.
///
/// - Parameters:
///   - x: Base (real). For negative bases, `y` must be an integer to remain real.
///   - y: Exponent (real).
/// - Returns: xʸ − 1.
/// - Throws: `SpecialFunctionError.invalidCombination` for negative base with non-integer exponent; `SpecialFunctionError.parameterNotFinite` if inputs are not finite.
@inlinable public func powm1<T: BinaryFloatingPoint>(_ x: T, _ y: T) throws -> T {
    let dx = D(x), dy = D(y)
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    guard dy.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "y") }
    let yIsInteger = dy == dy.rounded(.towardZero)
    guard !(dx < 0 && !yIsInteger) else {
        throw SpecialFunctionError.invalidCombination(message: "powm1 is undefined for negative base with non-integer exponent in the reals")
    }
    return T(bs_powm1(dx, dy))
}

/// Real cube root cbrt(x), accurate for all finite x.
///
/// Uses Boost.Math’s `cbrt`.
///
/// - Parameter x: Real input.
/// - Returns: The real cube root of `x`.
/// - Throws: `SpecialFunctionError.parameterNotFinite` if `x` is not finite.
@inlinable public func cbrt<T: BinaryFloatingPoint>(_ x: T) throws -> T {
    let dx = D(x)
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    return T(bs_cbrt(dx))
}

// Float overloads (exp-related)

/// Float overload of `exponentialIntegralEn(_:_: )`.
@inlinable public func exponentialIntegralEn(_ n: Int, _ x: Float) throws -> Float {
    guard n >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "n") }
    guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    return bs_expint_En_f(Int32(n), x)
}

/// Float overload of `expm1(_:)`.
@inlinable public func expm1(_ x: Float) throws -> Float {
    guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    return bs_expm1_f(x)
}

/// Float overload of `log1p(_:)`.
@inlinable public func log1p(_ x: Float) throws -> Float {
    guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    guard x > -1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: Double(Float(-1).nextUp), max: Double.infinity) }
    return bs_log1p_f(x)
}

/// Float overload of `log1pmx(_:)`.
@inlinable public func log1pmx(_ x: Float) throws -> Float {
    guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    guard x > -1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: Double(Float(-1).nextUp), max: Double.infinity) }
    return bs_log1pmx_f(x)
}

/// Float overload of `powm1(_: _:)`.
@inlinable public func powm1(_ x: Float, _ y: Float) throws -> Float {
    guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    guard y.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "y") }
    let yIsInteger = y == y.rounded(.towardZero)
    guard !(x < 0 && !yIsInteger) else {
        throw SpecialFunctionError.invalidCombination(message: "powm1 is undefined for negative base with non-integer exponent in the reals")
    }
    return bs_powm1_f(x, y)
}

/// Float overload of `cbrt(_:)`.
@inlinable public func cbrt(_ x: Float) throws -> Float {
    guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    return bs_cbrt_f(x)
}

// Float80 overloads (x86_64)

/// Float80 overload of `exponentialIntegralEn(_:_: )`.
#if arch(x86_64)
@inlinable public func exponentialIntegralEn(_ n: Int, _ x: Float80) throws -> Float80 {
    guard n >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "n") }
    guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    return bs_expint_En_l(Int32(n), x)
}
#endif

/// Float80 overload of `expm1(_:)`.
#if arch(x86_64)
@inlinable public func expm1(_ x: Float80) throws -> Float80 {
    guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    return bs_expm1_l(x)
}
#endif

/// Float80 overload of `log1p(_:)`.
#if arch(x86_64)
@inlinable public func log1p(_ x: Float80) throws -> Float80 {
    guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    guard x > -1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: Double(-1).nextUp, max: Double.infinity) }
    return bs_log1p_l(x)
}
#endif

/// Float80 overload of `log1pmx(_:)`.
#if arch(x86_64)
@inlinable public func log1pmx(_ x: Float80) throws -> Float80 {
    guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    guard x > -1 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: Double(-1).nextUp, max: Double.infinity) }
    return bs_log1pmx_l(x)
}
#endif

/// Float80 overload of `powm1(_: _:)`.
#if arch(x86_64)
@inlinable public func powm1(_ x: Float80, _ y: Float80) throws -> Float80 {
    guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    guard y.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "y") }
    let yIsInteger = y == y.rounded(.towardZero)
    guard !(x < 0 && !yIsInteger) else {
        throw SpecialFunctionError.invalidCombination(message: "powm1 is undefined for negative base with non-integer exponent in the reals")
    }
    return bs_powm1_l(x, y)
}
#endif

/// Float80 overload of `cbrt(_:)`.
#if arch(x86_64)
@inlinable public func cbrt(_ x: Float80) throws -> Float80 {
    guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    return bs_cbrt_l(x)
}
#endif

// MARK: - Trig helpers

/// sin(πx) computed with argument reduction for improved accuracy near integers.
///
/// Uses Boost.Math’s `sin_pi`.
///
/// - Parameter x: Real input.
/// - Returns: sin(πx).
/// - Throws: `SpecialFunctionError.parameterNotFinite` if `x` is not finite.
@inlinable public func sinPi<T: BinaryFloatingPoint>(_ x: T) throws -> T {
    let dx = D(x)
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    return T(bs_sin_pi(dx))
}

/// cos(πx) computed with argument reduction for improved accuracy near half-integers.
///
/// Uses Boost.Math’s `cos_pi`.
///
/// - Parameter x: Real input.
/// - Returns: cos(πx).
/// - Throws: `SpecialFunctionError.parameterNotFinite` if `x` is not finite.
@inlinable public func cosPi<T: BinaryFloatingPoint>(_ x: T) throws -> T {
    let dx = D(x)
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    return T(bs_cos_pi(dx))
}

/// Float overload of `sinPi(_:)`.
@inlinable public func sinPi(_ x: Float) throws -> Float {
    guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    return bs_sin_pi_f(x)
}

/// Float overload of `cosPi(_:)`.
@inlinable public func cosPi(_ x: Float) throws -> Float {
    guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    return bs_cos_pi_f(x)
}

/// Float80 overload of `sinPi(_:)`.
#if arch(x86_64)
@inlinable public func sinPi(_ x: Float80) throws -> Float80 {
    guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    return bs_sin_pi_l(x)
}
#endif

/// Float80 overload of `cosPi(_:)`.
#if arch(x86_64)
@inlinable public func cosPi(_ x: Float80) throws -> Float80 {
    guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    return bs_cos_pi_l(x)
}
#endif

// MARK: - Airy

/// Airy function Ai(x).
///
/// Uses Boost.Math’s `airy_ai`.
///
/// - Parameter x: Real input.
/// - Returns: Ai(x).
/// - Throws: `SpecialFunctionError.parameterNotFinite` if `x` is not finite.
@inlinable public func airyAi<T: BinaryFloatingPoint>(_ x: T) throws -> T {
    let dx = D(x)
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    return T(bs_airy_ai(dx))
}

/// Airy function Bi(x).
///
/// Uses Boost.Math’s `airy_bi`.
///
/// - Parameter x: Real input.
/// - Returns: Bi(x).
/// - Throws: `SpecialFunctionError.parameterNotFinite` if `x` is not finite.
@inlinable public func airyBi<T: BinaryFloatingPoint>(_ x: T) throws -> T {
    let dx = D(x)
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    return T(bs_airy_bi(dx))
}

/// Derivative Ai′(x) of the Airy function Ai.
///
/// Uses Boost.Math’s `airy_ai_prime`.
///
/// - Parameter x: Real input.
/// - Returns: Ai′(x).
/// - Throws: `SpecialFunctionError.parameterNotFinite` if `x` is not finite.
@inlinable public func airyAiPrime<T: BinaryFloatingPoint>(_ x: T) throws -> T {
    let dx = D(x)
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    return T(bs_airy_ai_prime(dx))
}

/// Derivative Bi′(x) of the Airy function Bi.
///
/// Uses Boost.Math’s `airy_bi_prime`.
///
/// - Parameter x: Real input.
/// - Returns: Bi′(x).
/// - Throws: `SpecialFunctionError.parameterNotFinite` if `x` is not finite.
@inlinable public func airyBiPrime<T: BinaryFloatingPoint>(_ x: T) throws -> T {
    let dx = D(x)
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    return T(bs_airy_bi_prime(dx))
}

/// Float overload of `airyAi(_:)`.
@inlinable public func airyAi(_ x: Float) throws -> Float {
    guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    return bs_airy_ai_f(x)
}

/// Float overload of `airyBi(_:)`.
@inlinable public func airyBi(_ x: Float) throws -> Float {
    guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    return bs_airy_bi_f(x)
}

/// Float overload of `airyAiPrime(_:)`.
@inlinable public func airyAiPrime(_ x: Float) throws -> Float {
    guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    return bs_airy_ai_prime_f(x)
}

/// Float overload of `airyBiPrime(_:)`.
@inlinable public func airyBiPrime(_ x: Float) throws -> Float {
    guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    return bs_airy_bi_prime_f(x)
}

/// Float80 overload of `airyAi(_:)`.
#if arch(x86_64)
@inlinable public func airyAi(_ x: Float80) throws -> Float80 {
    guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    return bs_airy_ai_l(x)
}
#endif

/// Float80 overload of `airyBi(_:)`.
#if arch(x86_64)
@inlinable public func airyBi(_ x: Float80) throws -> Float80 {
    guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    return bs_airy_bi_l(x)
}
#endif

/// Float80 overload of `airyAiPrime(_:)`.
#if arch(x86_64)
@inlinable public func airyAiPrime(_ x: Float80) throws -> Float80 {
    guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    return bs_airy_ai_prime_l(x)
}
#endif

/// Float80 overload of `airyBiPrime(_:)`.
#if arch(x86_64)
@inlinable public func airyBiPrime(_ x: Float80) throws -> Float80 {
    guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    return bs_airy_bi_prime_l(x)
}
#endif

// MARK: - Bessel (cylindrical)

/// Cylindrical Bessel function of the first kind Jᵥ(x).
///
/// Uses Boost.Math’s `cyl_bessel_j`.
///
/// - Parameters:
///   - v: Order (real).
///   - x: Argument (real).
/// - Returns: Jᵥ(x).
/// - Throws: `SpecialFunctionError.parameterNotFinite` if inputs are not finite.
@inlinable public func besselJ<T: BinaryFloatingPoint>(v: T, x: T) throws -> T {
    let dv = D(v), dx = D(x)
    guard dv.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "v") }
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    return T(bs_cyl_bessel_j(dv, dx))
}

/// Cylindrical Neumann function (Bessel Y) Yᵥ(x).
///
/// Uses Boost.Math’s `cyl_neumann`.
///
/// - Parameters:
///   - v: Order (real).
///   - x: Argument (real), must be > 0 (singularity at 0).
/// - Returns: Yᵥ(x).
/// - Throws: `SpecialFunctionError.parameterNotFinite` if inputs are not finite; `SpecialFunctionError.parameterNotPositive` if `x ≤ 0`.
@inlinable public func besselY<T: BinaryFloatingPoint>(v: T, x: T) throws -> T {
    let dv = D(v), dx = D(x)
    guard dv.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "v") }
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    guard dx > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "x") }
    return T(bs_cyl_neumann(dv, dx))
}

/// Modified Bessel function of the first kind Iᵥ(x).
///
/// Uses Boost.Math’s `cyl_bessel_i`.
///
/// - Parameters:
///   - v: Order (real).
///   - x: Argument (real).
/// - Returns: Iᵥ(x).
/// - Throws: `SpecialFunctionError.parameterNotFinite` if inputs are not finite.
@inlinable public func modifiedBesselI<T: BinaryFloatingPoint>(v: T, x: T) throws -> T {
    let dv = D(v), dx = D(x)
    guard dv.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "v") }
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    return T(bs_cyl_bessel_i(dv, dx))
}

/// Modified Bessel function of the second kind Kᵥ(x).
///
/// Uses Boost.Math’s `cyl_bessel_k`.
///
/// - Parameters:
///   - v: Order (real).
///   - x: Argument (real), must be > 0 (singularity at 0).
/// - Returns: Kᵥ(x).
/// - Throws: `SpecialFunctionError.parameterNotFinite` if inputs are not finite; `SpecialFunctionError.parameterNotPositive` if `x ≤ 0`.
@inlinable public func modifiedBesselK<T: BinaryFloatingPoint>(v: T, x: T) throws -> T {
    let dv = D(v), dx = D(x)
    guard dv.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "v") }
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    guard dx > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "x") }
    return T(bs_cyl_bessel_k(dv, dx))
}

/// Float overload of `besselJ(v:x:)`.
@inlinable public func besselJ(v: Float, x: Float) throws -> Float {
    guard v.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "v") }
    guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    return bs_cyl_bessel_j_f(v, x)
}

/// Float overload of `besselY(v:x:)`.
@inlinable public func besselY(v: Float, x: Float) throws -> Float {
    guard v.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "v") }
    guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    guard x > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "x") }
    return bs_cyl_neumann_f(v, x)
}

/// Float overload of `modifiedBesselI(v:x:)`.
@inlinable public func modifiedBesselI(v: Float, x: Float) throws -> Float {
    guard v.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "v") }
    guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    return bs_cyl_bessel_i_f(v, x)
}

/// Float overload of `modifiedBesselK(v:x:)`.
@inlinable public func modifiedBesselK(v: Float, x: Float) throws -> Float {
    guard v.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "v") }
    guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    guard x > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "x") }
    return bs_cyl_bessel_k_f(v, x)
}

/// Float80 overload of `besselJ(v:x:)`.
#if arch(x86_64)
@inlinable public func besselJ(v: Float80, x: Float80) throws -> Float80 {
    guard v.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "v") }
    guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    return bs_cyl_bessel_j_l(v, x)
}
#endif

/// Float80 overload of `besselY(v:x:)`.
#if arch(x86_64)
@inlinable public func besselY(v: Float80, x: Float80) throws -> Float80 {
    guard v.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "v") }
    guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    guard x > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "x") }
    return bs_cyl_neumann_l(v, x)
}
#endif

/// Float80 overload of `modifiedBesselI(v:x:)`.
#if arch(x86_64)
@inlinable public func modifiedBesselI(v: Float80, x: Float80) throws -> Float80 {
    guard v.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "v") }
    guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    return bs_cyl_bessel_i_l(v, x)
}
#endif

/// Float80 overload of `modifiedBesselK(v:x:)`.
#if arch(x86_64)
@inlinable public func modifiedBesselK(v: Float80, x: Float80) throws -> Float80 {
    guard v.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "v") }
    guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    guard x > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "x") }
    return bs_cyl_bessel_k_l(v, x)
}
#endif

// MARK: - Legendre

/// Legendre polynomial Pₙ(x) for integer n ≥ 0.
///
/// Uses Boost.Math’s `legendre_p`.
///
/// - Parameters:
///   - n: Degree n ≥ 0.
///   - x: Real input.
/// - Returns: Pₙ(x).
/// - Throws: `SpecialFunctionError.parameterNotPositive` if `n < 0`; `SpecialFunctionError.parameterNotFinite` if `x` is not finite.
@inlinable public func legendreP<T: BinaryFloatingPoint>(_ n: Int, _ x: T) throws -> T {
    guard n >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "n") }
    let dx = D(x)
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    return T(bs_legendre_p(Int32(n), dx))
}

/// Associated Legendre function Pₙᵐ(x) for integers n ≥ 0 and |m| ≤ n.
///
/// Uses Boost.Math’s `legendre_p(n, m, x)`.
///
/// - Parameters:
///   - n: Degree n ≥ 0.
///   - m: Order with |m| ≤ n.
///   - x: Real input.
/// - Returns: Pₙᵐ(x).
/// - Throws: `SpecialFunctionError.parameterNotPositive` if `n < 0`; `SpecialFunctionError.parameterOutOfRange` if `|m| > n`; `SpecialFunctionError.parameterNotFinite` if `x` is not finite.
@inlinable public func associatedLegendreP<T: BinaryFloatingPoint>(_ n: Int, _ m: Int, _ x: T) throws -> T {
    guard n >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "n") }
    guard abs(m) <= n else { throw SpecialFunctionError.parameterOutOfRange(name: "m", min: Double(-n), max: Double(n)) }
    let dx = D(x)
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    return T(bs_assoc_legendre_p(Int32(n), Int32(m), dx))
}

/// Float overload of `legendreP(_:_: )`.
@inlinable public func legendreP(_ n: Int, _ x: Float) throws -> Float {
    guard n >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "n") }
    guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    return bs_legendre_p_f(Int32(n), x)
}

/// Float overload of `associatedLegendreP(_:_:_:)`.
@inlinable public func associatedLegendreP(_ n: Int, _ m: Int, _ x: Float) throws -> Float {
    guard n >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "n") }
    guard abs(m) <= n else { throw SpecialFunctionError.parameterOutOfRange(name: "m", min: Double(-n), max: Double(n)) }
    guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    return bs_assoc_legendre_p_f(Int32(n), Int32(m), x)
}

/// Float80 overload of `legendreP(_:_: )`.
#if arch(x86_64)
@inlinable public func legendreP(_ n: Int, _ x: Float80) throws -> Float80 {
    guard n >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "n") }
    guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    return bs_legendre_p_l(Int32(n), x)
}
#endif

/// Float80 overload of `associatedLegendreP(_:_:_:)`.
#if arch(x86_64)
@inlinable public func associatedLegendreP(_ n: Int, _ m: Int, _ x: Float80) throws -> Float80 {
    guard n >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "n") }
    guard abs(m) <= n else { throw SpecialFunctionError.parameterOutOfRange(name: "m", min: Double(-n), max: Double(n)) }
    guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    return bs_assoc_legendre_p_l(Int32(n), Int32(m), x)
}
#endif

// MARK: - Elliptic integrals (Legendre forms)

/// Complete elliptic integral of the first kind K(k).
///
/// Uses Boost.Math’s `ellint_1(k)`.
///
/// - Parameter k: Modulus with |k| ≤ 1.
/// - Returns: K(k).
/// - Throws: `SpecialFunctionError.parameterOutOfRange` if |k| > 1; `SpecialFunctionError.parameterNotFinite` if `k` is not finite.
@inlinable public func completeEllipticIntegralK<T: BinaryFloatingPoint>(_ k: T) throws -> T {
    let dk = D(k)
    guard dk.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "k") }
    guard abs(dk) <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "k", min: -1, max: 1) }
    return T(bs_ellint_1_complete(dk))
}

/// Incomplete elliptic integral of the first kind F(k, φ).
///
/// Uses Boost.Math’s `ellint_1(k, φ)`.
///
/// - Parameters:
///   - k: Modulus with |k| ≤ 1.
///   - phi: Amplitude φ.
/// - Returns: F(k, φ).
/// - Throws: `SpecialFunctionError.parameterOutOfRange` if |k| > 1; `SpecialFunctionError.parameterNotFinite` if inputs are not finite.
@inlinable public func incompleteEllipticIntegralF<T: BinaryFloatingPoint>(_ k: T, phi: T) throws -> T {
    let dk = D(k), dphi = D(phi)
    guard dk.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "k") }
    guard dphi.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "phi") }
    guard abs(dk) <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "k", min: -1, max: 1) }
    return T(bs_ellint_1(dk, dphi))
}

/// Complete elliptic integral of the second kind E(k).
///
/// Uses Boost.Math’s `ellint_2(k)`.
///
/// - Parameter k: Modulus with |k| ≤ 1.
/// - Returns: E(k).
/// - Throws: `SpecialFunctionError.parameterOutOfRange` if |k| > 1; `SpecialFunctionError.parameterNotFinite` if `k` is not finite.
@inlinable public func completeEllipticIntegralE<T: BinaryFloatingPoint>(_ k: T) throws -> T {
    let dk = D(k)
    guard dk.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "k") }
    guard abs(dk) <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "k", min: -1, max: 1) }
    return T(bs_ellint_2_complete(dk))
}

/// Incomplete elliptic integral of the second kind E(k, φ).
///
/// Uses Boost.Math’s `ellint_2(k, φ)`.
///
/// - Parameters:
///   - k: Modulus with |k| ≤ 1.
///   - phi: Amplitude φ.
/// - Returns: E(k, φ).
/// - Throws: `SpecialFunctionError.parameterOutOfRange` if |k| > 1; `SpecialFunctionError.parameterNotFinite` if inputs are not finite.
@inlinable public func incompleteEllipticIntegralE<T: BinaryFloatingPoint>(_ k: T, phi: T) throws -> T {
    let dk = D(k), dphi = D(phi)
    guard dk.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "k") }
    guard dphi.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "phi") }
    guard abs(dk) <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "k", min: -1, max: 1) }
    return T(bs_ellint_2(dk, dphi))
}

/// Complete elliptic integral of the third kind Π(n; k).
///
/// Uses Boost.Math’s `ellint_3(k, n)`.
///
/// - Parameters:
///   - k: Modulus with |k| ≤ 1.
///   - n: Characteristic n.
/// - Returns: Π(n; k).
/// - Throws: `SpecialFunctionError.parameterOutOfRange` if |k| > 1; `SpecialFunctionError.parameterNotFinite` if inputs are not finite.
@inlinable public func completeEllipticIntegralPi<T: BinaryFloatingPoint>(_ k: T, characteristic n: T) throws -> T {
    let dk = D(k), dn = D(n)
    guard dk.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "k") }
    guard dn.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "n") }
    guard abs(dk) <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "k", min: -1, max: 1) }
    return T(bs_ellint_3_complete(dk, dn))
}

/// Incomplete elliptic integral of the third kind Π(n; k, φ).
///
/// Uses Boost.Math’s `ellint_3(k, n, φ)`.
///
/// - Parameters:
///   - k: Modulus with |k| ≤ 1.
///   - n: Characteristic n.
///   - phi: Amplitude φ.
/// - Returns: Π(n; k, φ).
/// - Throws: `SpecialFunctionError.parameterOutOfRange` if |k| > 1; `SpecialFunctionError.parameterNotFinite` if inputs are not finite.
@inlinable public func incompleteEllipticIntegralPi<T: BinaryFloatingPoint>(_ k: T, characteristic n: T, phi: T) throws -> T {
    let dk = D(k), dn = D(n), dphi = D(phi)
    guard dk.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "k") }
    guard dn.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "n") }
    guard dphi.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "phi") }
    guard abs(dk) <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "k", min: -1, max: 1) }
    return T(bs_ellint_3(dk, dn, dphi))
}

// Float overloads (Legendre forms)

/// Float overload of `completeEllipticIntegralK(_:)`.
@inlinable public func completeEllipticIntegralK(_ k: Float) throws -> Float {
    guard k.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "k") }
    guard abs(k) <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "k", min: -1, max: 1) }
    return bs_ellint_1_complete_f(k)
}

/// Float overload of `incompleteEllipticIntegralF(_:phi:)`.
@inlinable public func incompleteEllipticIntegralF(_ k: Float, phi: Float) throws -> Float {
    guard k.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "k") }
    guard phi.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "phi") }
    guard abs(k) <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "k", min: -1, max: 1) }
    return bs_ellint_1_f(k, phi)
}

/// Float overload of `completeEllipticIntegralE(_:)`.
@inlinable public func completeEllipticIntegralE(_ k: Float) throws -> Float {
    guard k.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "k") }
    guard abs(k) <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "k", min: -1, max: 1) }
    return bs_ellint_2_complete_f(k)
}

/// Float overload of `incompleteEllipticIntegralE(_:phi:)`.
@inlinable public func incompleteEllipticIntegralE(_ k: Float, phi: Float) throws -> Float {
    guard k.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "k") }
    guard phi.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "phi") }
    guard abs(k) <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "k", min: -1, max: 1) }
    return bs_ellint_2_f(k, phi)
}

/// Float overload of `completeEllipticIntegralPi(_:characteristic:)`.
@inlinable public func completeEllipticIntegralPi(_ k: Float, characteristic n: Float) throws -> Float {
    guard k.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "k") }
    guard n.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "n") }
    guard abs(k) <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "k", min: -1, max: 1) }
    return bs_ellint_3_complete_f(k, n)
}

/// Float overload of `incompleteEllipticIntegralPi(_:characteristic:phi:)`.
@inlinable public func incompleteEllipticIntegralPi(_ k: Float, characteristic n: Float, phi: Float) throws -> Float {
    guard k.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "k") }
    guard n.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "n") }
    guard phi.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "phi") }
    guard abs(k) <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "k", min: -1, max: 1) }
    return bs_ellint_3_f(k, n, phi)
}

// Float80 overloads (Legendre forms)

/// Float80 overload of `completeEllipticIntegralK(_:)`.
#if arch(x86_64)
@inlinable public func completeEllipticIntegralK(_ k: Float80) throws -> Float80 {
    guard k.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "k") }
    guard abs(k) <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "k", min: -1, max: 1) }
    return bs_ellint_1_complete_l(k)
}
#endif

/// Float80 overload of `incompleteEllipticIntegralF(_:phi:)`.
#if arch(x86_64)
@inlinable public func incompleteEllipticIntegralF(_ k: Float80, phi: Float80) throws -> Float80 {
    guard k.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "k") }
    guard phi.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "phi") }
    guard abs(k) <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "k", min: -1, max: 1) }
    return bs_ellint_1_l(k, phi)
}
#endif

/// Float80 overload of `completeEllipticIntegralE(_:)`.
#if arch(x86_64)
@inlinable public func completeEllipticIntegralE(_ k: Float80) throws -> Float80 {
    guard k.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "k") }
    guard abs(k) <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "k", min: -1, max: 1) }
    return bs_ellint_2_complete_l(k)
}
#endif

/// Float80 overload of `incompleteEllipticIntegralE(_:phi:)`.
#if arch(x86_64)
@inlinable public func incompleteEllipticIntegralE(_ k: Float80, phi: Float80) throws -> Float80 {
    guard k.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "k") }
    guard phi.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "phi") }
    guard abs(k) <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "k", min: -1, max: 1) }
    return bs_ellint_2_l(k, phi)
}
#endif

/// Float80 overload of `completeEllipticIntegralPi(_:characteristic:)`.
#if arch(x86_64)
@inlinable public func completeEllipticIntegralPi(_ k: Float80, characteristic n: Float80) throws -> Float80 {
    guard k.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "k") }
    guard n.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "n") }
    guard abs(k) <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "k", min: -1, max: 1) }
    return bs_ellint_3_complete_l(k, n)
}
#endif

/// Float80 overload of `incompleteEllipticIntegralPi(_:characteristic:phi:)`.
#if arch(x86_64)
@inlinable public func incompleteEllipticIntegralPi(_ k: Float80, characteristic n: Float80, phi: Float80) throws -> Float80 {
    guard k.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "k") }
    guard n.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "n") }
    guard phi.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "phi") }
    guard abs(k) <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "k", min: -1, max: 1) }
    return bs_ellint_3_l(k, n, phi)
}
#endif

// MARK: - Elliptic integrals (Carlson symmetric forms)

/// Carlson symmetric integral R_C(x, y).
///
/// Uses Boost.Math’s `ellint_rc`.
///
/// - Parameters:
///   - x: Real input.
///   - y: Real input, must be > 0 to avoid singularities.
/// - Returns: R_C(x, y).
/// - Throws: `SpecialFunctionError.parameterNotFinite` if inputs are not finite; `SpecialFunctionError.parameterNotPositive` if `y ≤ 0`.
@inlinable public func carlsonRC<T: BinaryFloatingPoint>(_ x: T, _ y: T) throws -> T {
    let dx = D(x), dy = D(y)
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    guard dy.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "y") }
    guard dy > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "y") }
    return T(bs_ellint_rc(dx, dy))
}

/// Carlson symmetric integral R_F(x, y, z).
///
/// Uses Boost.Math’s `ellint_rf`.
///
/// - Parameters:
///   - x: Real input (≥ 0).
///   - y: Real input (≥ 0).
///   - z: Real input (≥ 0).
/// - Returns: R_F(x, y, z).
/// - Throws: `SpecialFunctionError.parameterNotFinite` if inputs are not finite; `SpecialFunctionError.parameterOutOfRange` if any input < 0.
@inlinable public func carlsonRF<T: BinaryFloatingPoint>(_ x: T, _ y: T, _ z: T) throws -> T {
    let dx = D(x), dy = D(y), dz = D(z)
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    guard dy.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "y") }
    guard dz.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "z") }
    guard dx >= 0, dy >= 0, dz >= 0 else {
        throw SpecialFunctionError.parameterOutOfRange(name: "x/y/z", min: 0, max: Double.infinity)
    }
    return T(bs_ellint_rf(dx, dy, dz))
}

/// Carlson symmetric integral R_D(x, y, z).
///
/// Uses Boost.Math’s `ellint_rd`.
///
/// - Parameters:
///   - x: Real input (≥ 0).
///   - y: Real input (≥ 0).
///   - z: Real input (≥ 0).
/// - Returns: R_D(x, y, z).
/// - Throws: `SpecialFunctionError.parameterNotFinite` if inputs are not finite; `SpecialFunctionError.parameterOutOfRange` if any input < 0.
@inlinable public func carlsonRD<T: BinaryFloatingPoint>(_ x: T, _ y: T, _ z: T) throws -> T {
    let dx = D(x), dy = D(y), dz = D(z)
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    guard dy.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "y") }
    guard dz.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "z") }
    guard dx >= 0, dy >= 0, dz >= 0 else {
        throw SpecialFunctionError.parameterOutOfRange(name: "x/y/z", min: 0, max: Double.infinity)
    }
    return T(bs_ellint_rd(dx, dy, dz))
}

/// Carlson symmetric integral R_J(x, y, z, p).
///
/// Uses Boost.Math’s `ellint_rj`.
///
/// - Parameters:
///   - x: Real input (≥ 0).
///   - y: Real input (≥ 0).
///   - z: Real input (≥ 0).
///   - p: Real input (≥ 0).
/// - Returns: R_J(x, y, z, p).
/// - Throws: `SpecialFunctionError.parameterNotFinite` if inputs are not finite; `SpecialFunctionError.parameterOutOfRange` if any input < 0.
@inlinable public func carlsonRJ<T: BinaryFloatingPoint>(_ x: T, _ y: T, _ z: T, _ p: T) throws -> T {
    let dx = D(x), dy = D(y), dz = D(z), dp = D(p)
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    guard dy.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "y") }
    guard dz.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "z") }
    guard dp.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "p") }
    guard dx >= 0, dy >= 0, dz >= 0, dp >= 0 else {
        throw SpecialFunctionError.parameterOutOfRange(name: "x/y/z/p", min: 0, max: Double.infinity)
    }
    return T(bs_ellint_rj(dx, dy, dz, dp))
}

/// Carlson symmetric integral R_G(x, y, z).
///
/// Uses Boost.Math’s `ellint_rg`.
///
/// - Parameters:
///   - x: Real input (≥ 0).
///   - y: Real input (≥ 0).
///   - z: Real input (≥ 0).
/// - Returns: R_G(x, y, z).
/// - Throws: `SpecialFunctionError.parameterNotFinite` if inputs are not finite; `SpecialFunctionError.parameterOutOfRange` if any input < 0.
@inlinable public func carlsonRG<T: BinaryFloatingPoint>(_ x: T, _ y: T, _ z: T) throws -> T {
    let dx = D(x), dy = D(y), dz = D(z)
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    guard dy.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "y") }
    guard dz.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "z") }
    guard dx >= 0, dy >= 0, dz >= 0 else {
        throw SpecialFunctionError.parameterOutOfRange(name: "x/y/z", min: 0, max: Double.infinity)
    }
    return T(bs_ellint_rg(dx, dy, dz))
}

// Float overloads (Carlson)

/// Float overload of `carlsonRC(_:_: )`.
@inlinable public func carlsonRC(_ x: Float, _ y: Float) throws -> Float {
    guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    guard y.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "y") }
    guard y > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "y") }
    return bs_ellint_rc_f(x, y)
}

/// Float overload of `carlsonRF(_:_:_:)`.
@inlinable public func carlsonRF(_ x: Float, _ y: Float, _ z: Float) throws -> Float {
    guard x.isFinite, y.isFinite, z.isFinite else { throw SpecialFunctionError.invalidCombination(message: "non-finite inputs") }
    guard x >= 0, y >= 0, z >= 0 else { throw SpecialFunctionError.parameterOutOfRange(name: "x/y/z", min: 0, max: Double.infinity) }
    return bs_ellint_rf_f(x, y, z)
}

/// Float overload of `carlsonRD(_:_:_:)`.
@inlinable public func carlsonRD(_ x: Float, _ y: Float, _ z: Float) throws -> Float {
    guard x.isFinite, y.isFinite, z.isFinite else { throw SpecialFunctionError.invalidCombination(message: "non-finite inputs") }
    guard x >= 0, y >= 0, z >= 0 else { throw SpecialFunctionError.parameterOutOfRange(name: "x/y/z", min: 0, max: Double.infinity) }
    return bs_ellint_rd_f(x, y, z)
}

/// Float overload of `carlsonRJ(_:_:_:_:)`.
@inlinable public func carlsonRJ(_ x: Float, _ y: Float, _ z: Float, _ p: Float) throws -> Float {
    guard x.isFinite, y.isFinite, z.isFinite, p.isFinite else { throw SpecialFunctionError.invalidCombination(message: "non-finite inputs") }
    guard x >= 0, y >= 0, z >= 0, p >= 0 else { throw SpecialFunctionError.parameterOutOfRange(name: "x/y/z/p", min: 0, max: Double.infinity) }
    return bs_ellint_rj_f(x, y, z, p)
}

/// Float overload of `carlsonRG(_:_:_:)`.
@inlinable public func carlsonRG(_ x: Float, _ y: Float, _ z: Float) throws -> Float {
    guard x.isFinite, y.isFinite, z.isFinite else { throw SpecialFunctionError.invalidCombination(message: "non-finite inputs") }
    guard x >= 0, y >= 0, z >= 0 else { throw SpecialFunctionError.parameterOutOfRange(name: "x/y/z", min: 0, max: Double.infinity) }
    return bs_ellint_rg_f(x, y, z)
}

// Float80 overloads (Carlson)

/// Float80 overload of `carlsonRC(_:_: )`.
#if arch(x86_64)
@inlinable public func carlsonRC(_ x: Float80, _ y: Float80) throws -> Float80 {
    guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    guard y.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "y") }
    guard y > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "y") }
    return bs_ellint_rc_l(x, y)
}
#endif

/// Float80 overload of `carlsonRF(_:_:_:)`.
#if arch(x86_64)
@inlinable public func carlsonRF(_ x: Float80, _ y: Float80, _ z: Float80) throws -> Float80 {
    guard x.isFinite, y.isFinite, z.isFinite else { throw SpecialFunctionError.invalidCombination(message: "non-finite inputs") }
    guard x >= 0, y >= 0, z >= 0 else { throw SpecialFunctionError.parameterOutOfRange(name: "x/y/z", min: 0, max: Double.infinity) }
    return bs_ellint_rf_l(x, y, z)
}
#endif

/// Float80 overload of `carlsonRD(_:_:_:)`.
#if arch(x86_64)
@inlinable public func carlsonRD(_ x: Float80, _ y: Float80, _ z: Float80) throws -> Float80 {
    guard x.isFinite, y.isFinite, z.isFinite else { throw SpecialFunctionError.invalidCombination(message: "non-finite inputs") }
    guard x >= 0, y >= 0, z >= 0 else { throw SpecialFunctionError.parameterOutOfRange(name: "x/y/z", min: 0, max: Double.infinity) }
    return bs_ellint_rd_l(x, y, z)
}
#endif

/// Float80 overload of `carlsonRJ(_:_:_:_:)`.
#if arch(x86_64)
@inlinable public func carlsonRJ(_ x: Float80, _ y: Float80, _ z: Float80, _ p: Float80) throws -> Float80 {
    guard x.isFinite, y.isFinite, z.isFinite, p.isFinite else { throw SpecialFunctionError.invalidCombination(message: "non-finite inputs") }
    guard x >= 0, y >= 0, z >= 0, p >= 0 else { throw SpecialFunctionError.parameterOutOfRange(name: "x/y/z/p", min: 0, max: Double.infinity) }
    return bs_ellint_rj_l(x, y, z, p)
}
#endif

/// Float80 overload of `carlsonRG(_:_:_:)`.
#if arch(x86_64)
@inlinable public func carlsonRG(_ x: Float80, _ y: Float80, _ z: Float80) throws -> Float80 {
    guard x.isFinite, y.isFinite, z.isFinite else { throw SpecialFunctionError.invalidCombination(message: "non-finite inputs") }
    guard x >= 0, y >= 0, z >= 0 else { throw SpecialFunctionError.parameterOutOfRange(name: "x/y/z", min: 0, max: Double.infinity) }
    return bs_ellint_rg_l(x, y, z)
}
#endif

// MARK: - Lambert W (real branches)

/// Principal branch of the Lambert W function, W₀(x).
///
/// Uses Boost.Math’s `lambert_w0`. Domain is x ≥ −1/e.
///
/// - Parameter x: Real input (x ≥ −1/e).
/// - Returns: W₀(x).
/// - Throws: `SpecialFunctionError.parameterOutOfRange` if `x < −1/e`; `SpecialFunctionError.parameterNotFinite` if `x` is not finite.
/// - SeeAlso: `lambertWm1(_:)`
@inlinable public func lambertW0<T: BinaryFloatingPoint>(_ x: T) throws -> T {
    let dx = D(x)
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    let minX = -1.0 / bs_const_e()
    guard dx >= minX else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: minX, max: Double.infinity) }
    return T(bs_lambert_w0(dx))
}

/// Lower real branch of the Lambert W function, W₋₁(x).
///
/// Uses Boost.Math’s `lambert_wm1`. Domain is −1/e ≤ x < 0.
///
/// - Parameter x: Real input (−1/e ≤ x < 0).
/// - Returns: W₋₁(x).
/// - Throws: `SpecialFunctionError.parameterOutOfRange` if `x ∉ [−1/e, 0)`; `SpecialFunctionError.parameterNotFinite` if `x` is not finite.
/// - SeeAlso: `lambertW0(_:)`
@inlinable public func lambertWm1<T: BinaryFloatingPoint>(_ x: T) throws -> T {
    let dx = D(x)
    guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    let minX = -1.0 / bs_const_e()
    guard dx >= minX && dx < 0 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: minX, max: 0.0) }
    return T(bs_lambert_wm1(dx))
}

/// Float overload of `lambertW0(_:)`.
@inlinable public func lambertW0(_ x: Float) throws -> Float {
    guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    let minX = -1.0 as Float / bs_const_e_f()
    guard x >= minX else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: Double(minX), max: Double.infinity) }
    return bs_lambert_w0_f(x)
}

/// Float overload of `lambertWm1(_:)`.
@inlinable public func lambertWm1(_ x: Float) throws -> Float {
    guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    let minX = -1.0 as Float / bs_const_e_f()
    guard x >= minX && x < 0 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: Double(minX), max: 0.0) }
    return bs_lambert_wm1_f(x)
}

/// Float80 overload of `lambertW0(_:)`.
#if arch(x86_64)
@inlinable public func lambertW0(_ x: Float80) throws -> Float80 {
    guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    let minX = -1.0 as Float80 / bs_const_e_l()
    guard x >= minX else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: Double(minX), max: Double.infinity) }
    return bs_lambert_w0_l(x)
}
#endif

/// Float80 overload of `lambertWm1(_:)`.
#if arch(x86_64)
@inlinable public func lambertWm1(_ x: Float80) throws -> Float80 {
    guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
    let minX = -1.0 as Float80 / bs_const_e_l()
    guard x >= minX && x < 0 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: Double(minX), max: 0.0) }
    return bs_lambert_wm1_l(x)
}
#endif
