//
//  EllipticCarlson.swift
//  Math/SpecialFunctions
//
//  Swift wrappers for Carlson’s symmetric elliptic integrals RC, RF, RD, RJ and the
//  related mean RG. These functions provide type-generic and type-specific entry
//  points, validate common real-domain constraints, and delegate numerical work to
//  CBoostBridge (Boost.Math).
//
//  Carlson symmetric integrals (real-valued):
//  - RC(x, y) = ∫₀^∞ dt / ((t + x) √(t + y)), with y > 0 (wrapper enforces y > 0).
//  - RF(x, y, z) = ∫₀^∞ dt / √((t + x)(t + y)(t + z)), wrapper enforces x, y, z ≥ 0.
//  - RD(x, y, z) = ∫₀^∞ dt / ((t + z) √((t + x)(t + y)(t + z))), wrapper enforces x, y, z ≥ 0.
//  - RJ(x, y, z, p) = ∫₀^∞ dt / ((t + p) √((t + x)(t + y)(t + z))), wrapper enforces x, y, z, p ≥ 0.
//  - RG(x, y, z) is Carlson’s symmetric “mean” (closely related to RF/RD); wrapper enforces x, y, z ≥ 0.
//
//  Notes on domains:
//  - The wrappers here impose simple, sufficient real-domain constraints that match
//    typical usage and the Boost.Math C entry points used by CBoostBridge.
//  - In particular, we enforce non-negativity for RF/RD/RJ/RG arguments and y > 0
//    for RC. These are sufficient for real-valued results and avoid singularities
//    on common branches. Stricter mathematical domains exist (see DLMF §19), but
//    we keep the Swift API conservative and predictable.
//  - Float overloads use a single “non-finite inputs” invalidCombination error for
//    brevity/performance; generic overloads and Float80 variants check each input.
//
//  Error model:
//  - Non-finite inputs throw SpecialFunctionError.parameterNotFinite (generic/Float80)
//    or SpecialFunctionError.invalidCombination(message:value:) in some Float overloads, with `value` capturing the offending argument.
//  - Real-domain violations throw SpecialFunctionError.parameterOutOfRange or
//    SpecialFunctionError.parameterNotPositive as documented per function.
//
//  References:
//  - NIST DLMF §19 (Elliptic Integrals): https://dlmf.nist.gov/19
//  - B. C. Carlson, “Numerical computation of real or complex elliptic integrals”,
//    Numer. Math. 33, 1–16 (1979).
//  - Boost.Math Carlson symmetric integrals:
//    https://www.boost.org/doc/libs/release/libs/math/doc/html/math_toolkit/ellint/carlson.html
//

import SwiftyBoostPrelude
public extension SpecialFunctions {
    
    
    // MARK: - Generic BinaryFloatingPoint overloads (Double-backed)
    
    /// Carlson’s symmetric integral RC(x, y).
    ///
    /// Real-domain (this wrapper):
    /// - Requires y > 0 for a real-valued result. x may be any finite real; common
    ///   practice is x ≥ 0, but Boost.Math supports broader real x when y > 0.
    /// - This function validates finiteness and y > 0.
    ///
    /// Parameters:
    /// - x: First argument (finite real).
    /// - y: Second argument (finite real, strictly positive).
    ///
    /// Returns:
    /// - RC(x, y) as `T`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotFinite(name: "x")` if `x` is NaN or ±∞.
    /// - `SpecialFunctionError.parameterNotFinite(name: "y")` if `y` is NaN or ±∞.
    /// - `SpecialFunctionError.parameterNotPositive(name: "y", value: y)` if `y ≤ 0`.
    @inlinable static func carlsonRC<T: Real & BinaryFloatingPoint & Sendable>(_ x: T, _ y: T) throws -> T {
        let dx = D(x), dy = D(y)
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        guard dy.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "y", value: y) }
        guard dy > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "y", value: y) }
        return T(bs_ellint_rc_d(dx, dy))
    }

    /// Carlson’s symmetric integral RF(x, y, z).
    ///
    /// Real-domain (this wrapper):
    /// - Requires x ≥ 0, y ≥ 0, z ≥ 0 for a real-valued result. Singularities occur
    ///   when any argument is negative enough to cross branch cuts; we conservatively
    ///   enforce non-negativity.
    ///
    /// Parameters:
    /// - x: First argument (finite, ≥ 0).
    /// - y: Second argument (finite, ≥ 0).
    /// - z: Third argument (finite, ≥ 0).
    ///
    /// Returns:
    /// - RF(x, y, z) as `T`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotFinite(name: "x"|"y"|"z")` if any is NaN or ±∞.
    /// - `SpecialFunctionError.parameterOutOfRange(name: "x/y/z", min: 0, max: +∞)` if any < 0.
    @inlinable static func carlsonRF<T: Real & BinaryFloatingPoint & Sendable>(_ x: T, _ y: T, _ z: T) throws -> T {
        let dx = D(x), dy = D(y), dz = D(z)
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        guard dy.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "y", value: y) }
        guard dz.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "z", value: z) }
        guard dx >= 0, dy >= 0, dz >= 0 else {
            throw SpecialFunctionError.parameterOutOfRange(name: "x/y/z", min: 0, max: Double.infinity)
        }
        return T(bs_ellint_rf_d(dx, dy, dz))
    }

    /// Carlson’s symmetric integral RD(x, y, z).
    ///
    /// Real-domain (this wrapper):
    /// - Requires x ≥ 0, y ≥ 0, z ≥ 0 for a real-valued result. Additional singular
    ///   cases (e.g., z = 0) are handled by Boost.Math; the wrapper enforces basic
    ///   non-negativity and finiteness.
    ///
    /// Parameters:
    /// - x: First argument (finite, ≥ 0).
    /// - y: Second argument (finite, ≥ 0).
    /// - z: Third argument (finite, ≥ 0).
    ///
    /// Returns:
    /// - RD(x, y, z) as `T`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotFinite(name: "x"|"y"|"z")` if any is NaN or ±∞.
    /// - `SpecialFunctionError.parameterOutOfRange(name: "x/y/z", min: 0, max: +∞)` if any < 0.
    @inlinable static func carlsonRD<T: Real & BinaryFloatingPoint & Sendable>(_ x: T, _ y: T, _ z: T) throws -> T {
        let dx = D(x), dy = D(y), dz = D(z)
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        guard dy.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "y", value: y) }
        guard dz.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "z", value: z) }
        guard dx >= 0, dy >= 0, dz >= 0 else {
            throw SpecialFunctionError.parameterOutOfRange(name: "x/y/z", min: 0, max: Double.infinity)
        }
        return T(bs_ellint_rd_d(dx, dy, dz))
    }

    /// Carlson’s symmetric integral RJ(x, y, z, p).
    ///
    /// Real-domain (this wrapper):
    /// - Requires x ≥ 0, y ≥ 0, z ≥ 0, p ≥ 0 for a real-valued result. Singularities
    ///   at p = 0 or negative values are avoided by the non-negativity check.
    ///
    /// Parameters:
    /// - x: First argument (finite, ≥ 0).
    /// - y: Second argument (finite, ≥ 0).
    /// - z: Third argument (finite, ≥ 0).
    /// - p: Fourth argument (finite, ≥ 0).
    ///
    /// Returns:
    /// - RJ(x, y, z, p) as `T`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotFinite(name: "x"|"y"|"z"|"p")` if any is NaN or ±∞.
    /// - `SpecialFunctionError.parameterOutOfRange(name: "x/y/z/p", min: 0, max: +∞)` if any < 0.
    @inlinable static func carlsonRJ<T: Real & BinaryFloatingPoint & Sendable>(_ x: T, _ y: T, _ z: T, _ p: T) throws -> T {
        let dx = D(x), dy = D(y), dz = D(z), dp = D(p)
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        guard dy.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "y", value: y) }
        guard dz.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "z", value: z) }
        guard dp.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "p", value: p) }
        guard dx >= 0, dy >= 0, dz >= 0, dp >= 0 else {
            throw SpecialFunctionError.parameterOutOfRange(name: "x/y/z/p", min: 0, max: Double.infinity)
        }
        return T(bs_ellint_rj_d(dx, dy, dz, dp))
    }

    /// Carlson’s symmetric mean RG(x, y, z).
    ///
    /// Real-domain (this wrapper):
    /// - Requires x ≥ 0, y ≥ 0, z ≥ 0 for a real-valued result.
    ///
    /// Parameters:
    /// - x: First argument (finite, ≥ 0).
    /// - y: Second argument (finite, ≥ 0).
    /// - z: Third argument (finite, ≥ 0).
    ///
    /// Returns:
    /// - RG(x, y, z) as `T`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotFinite(name: "x"|"y"|"z")` if any is NaN or ±∞.
    /// - `SpecialFunctionError.parameterOutOfRange(name: "x/y/z", min: 0, max: +∞)` if any < 0.
    @inlinable static func carlsonRG<T: Real & BinaryFloatingPoint & Sendable>(_ x: T, _ y: T, _ z: T) throws -> T {
        let dx = D(x), dy = D(y), dz = D(z)
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        guard dy.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "y", value: y) }
        guard dz.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "z", value: z) }
        guard dx >= 0, dy >= 0, dz >= 0 else {
            throw SpecialFunctionError.parameterOutOfRange(name: "x/y/z", min: 0, max: Double.infinity)
        }
        return T(bs_ellint_rg_d(dx, dy, dz))
    }

    // MARK: - Float overloads
    // These overloads call directly into the Float-precision C implementations for
    // performance and avoid intermediate conversions. Domain checks are streamlined.
    
    /// RC(x, y) for `Float`. Requires y > 0.
    @inlinable static func carlsonRC(_ x: Float, _ y: Float) throws -> Float {
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        guard y.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "y", value: y) }
        guard y > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "y", value: y) }
        return bs_ellint_rc_f(x, y)
    }
    
    /// RF(x, y, z) for `Float`. Requires x, y, z ≥ 0 and finite.
    @inlinable static func carlsonRF(_ x: Float, _ y: Float, _ z: Float) throws -> Float {
        guard x.isFinite else { throw SpecialFunctionError.invalidCombination(message: "non-finite inputs", value: x) }
        guard y.isFinite else { throw SpecialFunctionError.invalidCombination(message: "non-finite inputs", value: y) }
        guard z.isFinite else { throw SpecialFunctionError.invalidCombination(message: "non-finite inputs", value: z) }
        guard x >= 0, y >= 0, z >= 0 else { throw SpecialFunctionError.parameterOutOfRange(name: "x/y/z", min: 0, max: Float.infinity) }
        return bs_ellint_rf_f(x, y, z)
    }
    
    /// RD(x, y, z) for `Float`. Requires x, y, z ≥ 0 and finite.
    @inlinable static func carlsonRD(_ x: Float, _ y: Float, _ z: Float) throws -> Float {
        guard x.isFinite else { throw SpecialFunctionError.invalidCombination(message: "non-finite inputs", value: x) }
        guard y.isFinite else { throw SpecialFunctionError.invalidCombination(message: "non-finite inputs", value: y) }
        guard z.isFinite else { throw SpecialFunctionError.invalidCombination(message: "non-finite inputs", value: z) }
        guard x >= 0, y >= 0, z >= 0 else { throw SpecialFunctionError.parameterOutOfRange(name: "x/y/z", min: 0, max: Float.infinity) }
        return bs_ellint_rd_f(x, y, z)
    }
    
    /// RJ(x, y, z, p) for `Float`. Requires x, y, z, p ≥ 0 and finite.
    @inlinable static func carlsonRJ(_ x: Float, _ y: Float, _ z: Float, _ p: Float) throws -> Float {
        guard x.isFinite else { throw SpecialFunctionError.invalidCombination(message: "non-finite inputs", value: x) }
        guard y.isFinite else { throw SpecialFunctionError.invalidCombination(message: "non-finite inputs", value: y) }
        guard z.isFinite else { throw SpecialFunctionError.invalidCombination(message: "non-finite inputs", value: z) }
        guard p.isFinite else { throw SpecialFunctionError.invalidCombination(message: "non-finite inputs", value: p) }
        guard x >= 0, y >= 0, z >= 0, p >= 0 else { throw SpecialFunctionError.parameterOutOfRange(name: "x/y/z/p", min: 0, max: Float.infinity) }
        return bs_ellint_rj_f(x, y, z, p)
    }
    
    /// RG(x, y, z) for `Float`. Requires x, y, z ≥ 0 and finite.
    @inlinable static func carlsonRG(_ x: Float, _ y: Float, _ z: Float) throws -> Float {
        guard x.isFinite else { throw SpecialFunctionError.invalidCombination(message: "non-finite inputs", value: x) }
        guard y.isFinite else { throw SpecialFunctionError.invalidCombination(message: "non-finite inputs", value: y) }
        guard z.isFinite else { throw SpecialFunctionError.invalidCombination(message: "non-finite inputs", value: z) }
        guard x >= 0, y >= 0, z >= 0 else { throw SpecialFunctionError.parameterOutOfRange(name: "x/y/z", min: 0, max: Float.infinity) }
        return bs_ellint_rg_f(x, y, z)
    }
    
    // MARK: - Float80 overloads (x86_64)
    // Extended-precision entry points for platforms that support Float80.
    
#if arch(x86_64)
    /// RC(x, y) for `Float80` (x86_64 only). Requires y > 0.
    @inlinable static func carlsonRC(_ x: Float80, _ y: Float80) throws -> Float80 {
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        guard y.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "y", value: y) }
        guard y > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "y", value: y) }
        return bs_ellint_rc_l(x, y)
    }

    /// RF(x, y, z) for `Float80` (x86_64 only). Requires x, y, z ≥ 0 and finite.
    @inlinable static func carlsonRF(_ x: Float80, _ y: Float80, _ z: Float80) throws -> Float80 {
        guard x.isFinite else { throw SpecialFunctionError.invalidCombination(message: "non-finite inputs", value: x) }
        guard y.isFinite else { throw SpecialFunctionError.invalidCombination(message: "non-finite inputs", value: y) }
        guard z.isFinite else { throw SpecialFunctionError.invalidCombination(message: "non-finite inputs", value: z) }
        guard x >= 0, y >= 0, z >= 0 else { throw SpecialFunctionError.parameterOutOfRange(name: "x/y/z", min: 0, max: Float80.infinity) }
        return bs_ellint_rf_l(x, y, z)
    }

    /// RD(x, y, z) for `Float80` (x86_64 only). Requires x, y, z ≥ 0 and finite.
    @inlinable static func carlsonRD(_ x: Float80, _ y: Float80, _ z: Float80) throws -> Float80 {
        guard x.isFinite else { throw SpecialFunctionError.invalidCombination(message: "non-finite inputs", value: x) }
        guard y.isFinite else { throw SpecialFunctionError.invalidCombination(message: "non-finite inputs", value: y) }
        guard z.isFinite else { throw SpecialFunctionError.invalidCombination(message: "non-finite inputs", value: z) }
        guard x >= 0, y >= 0, z >= 0 else { throw SpecialFunctionError.parameterOutOfRange(name: "x/y/z", min: 0, max: Float80.infinity) }
        return bs_ellint_rd_l(x, y, z)
    }

    /// RJ(x, y, z, p) for `Float80` (x86_64 only). Requires x, y, z, p ≥ 0 and finite.
    @inlinable static func carlsonRJ(_ x: Float80, _ y: Float80, _ z: Float80, _ p: Float80) throws -> Float80 {
        guard x.isFinite else { throw SpecialFunctionError.invalidCombination(message: "non-finite inputs", value: x) }
        guard y.isFinite else { throw SpecialFunctionError.invalidCombination(message: "non-finite inputs", value: y) }
        guard z.isFinite else { throw SpecialFunctionError.invalidCombination(message: "non-finite inputs", value: z) }
        guard p.isFinite else { throw SpecialFunctionError.invalidCombination(message: "non-finite inputs", value: p) }
        guard x >= 0, y >= 0, z >= 0, p >= 0 else { throw SpecialFunctionError.parameterOutOfRange(name: "x/y/z/p", min: 0, max: Float80.infinity) }
        return bs_ellint_rj_l(x, y, z, p)
    }

    /// RG(x, y, z) for `Float80` (x86_64 only). Requires x, y, z ≥ 0 and finite.
    @inlinable static func carlsonRG(_ x: Float80, _ y: Float80, _ z: Float80) throws -> Float80 {
        guard x.isFinite else { throw SpecialFunctionError.invalidCombination(message: "non-finite inputs", value: x) }
        guard y.isFinite else { throw SpecialFunctionError.invalidCombination(message: "non-finite inputs", value: y) }
        guard z.isFinite else { throw SpecialFunctionError.invalidCombination(message: "non-finite inputs", value: z) }
        guard x >= 0, y >= 0, z >= 0 else { throw SpecialFunctionError.parameterOutOfRange(name: "x/y/z", min: 0, max: Float80.infinity) }
        return bs_ellint_rg_l(x, y, z)
    }
#endif
}
