//
//  LambertW.swift
//  Math/SpecialFunctions
//
//  Swift wrappers for the real branches of the Lambert W function, defined by
//  W(x) e^{W(x)} = x. The two real branches are:
//  - Principal branch W0(x), defined for x ∈ [−1/e, +∞).
//  - Lower branch W−1(x), defined for x ∈ [−1/e, 0).
//
//  These APIs validate real-domain constraints and delegate numerical work to
//  CBoostBridge (wrapping Boost.Math).
//
//  Functions (real-valued):
//  - lambertW0(x): principal branch W0 on [−1/e, +∞).
//  - lambertWm1(x): lower branch W−1 on [−1/e, 0).
//
//  Notes:
//  - At the branch point x = −1/e, both branches meet with W0(−1/e) = W−1(−1/e) = −1.
//  - At x = 0, W0(0) = 0, while W−1 is not defined (the real branch is open at 0).
//  - For large x, W0(x) grows like ln x − ln ln x + …
//  - These wrappers enforce finiteness and proper domain intervals and surface
//    violations via SpecialFunctionError.
//
//  Overloads:
//  - Generic BinaryFloatingPoint overloads funnel through a Double-backed
//    implementation using the helper D(_:), returning the result as T.
//  - Type-specific overloads for Float and (on x86_64) Float80 call directly
//    into the corresponding C functions to avoid conversions.
//
//  Error model:
//  - Non-finite arguments throw SpecialFunctionError.parameterNotFinite.
//  - Domain violations throw SpecialFunctionError.parameterOutOfRange with the
//    appropriate interval.
//
//  References:
//  - NIST DLMF §4.13 (Lambert W): https://dlmf.nist.gov/4.13
//  - Boost.Math Lambert W: https://www.boost.org/doc/libs/release/libs/math/doc/html/math_toolkit/lambert_w.html
//

import SwiftyBoostPrelude
public extension SpecialFunctions {
    
    
    // MARK: - Generic overloads (Double-backed)
    
    /// Principal branch Lambert W, W0(x), satisfying W e^W = x, with domain x ∈ [−1/e, +∞).
    ///
    /// Domain:
    /// - Requires x ≥ −1/e for a real-valued result on the principal branch.
    ///
    /// Parameters:
    /// - x: Input value (finite real, x ≥ −1/e).
    ///
    /// Returns:
    /// - W0(x) as `T`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotFinite(name: "x")` if `x` is NaN or ±∞.
    /// - `SpecialFunctionError.parameterOutOfRange(name: "x", min: −1/e, max: +∞)` if `x < −1/e`.
    @inlinable static func lambertW0<T: Real & BinaryFloatingPoint>(_ x: T) throws -> T {
        let dx = D(x)
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        let minX = -1.0 / bs_const_e_d()
        guard dx >= minX else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: minX, max: Double.infinity) }
        return T(bs_lambert_w0(dx))
    }
    
    /// Lower branch Lambert W, W−1(x), satisfying W e^W = x, with domain x ∈ [−1/e, 0).
    ///
    /// Domain:
    /// - Requires −1/e ≤ x < 0 for a real-valued result on the lower branch.
    ///
    /// Parameters:
    /// - x: Input value (finite real, −1/e ≤ x < 0).
    ///
    /// Returns:
    /// - W−1(x) as `T`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotFinite(name: "x")` if `x` is NaN or ±∞.
    /// - `SpecialFunctionError.parameterOutOfRange(name: "x", min: −1/e, max: 0)` if `x ∉ [−1/e, 0)`.
    @inlinable static func lambertWm1<T: Real & BinaryFloatingPoint>(_ x: T) throws -> T {
        let dx = D(x)
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        let minX = -1.0 / bs_const_e_d()
        guard dx >= minX && dx < 0 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: minX, max: 0.0) }
        return T(bs_lambert_wm1(dx))
    }
    
    // MARK: - Float overloads
    
    /// W0(x) for `Float`. Requires x ≥ −1/e.
    @inlinable static func lambertW0(_ x: Float) throws -> Float {
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        let minX = -1.0 as Float / bs_const_e_f()
        guard x >= minX else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: Double(minX), max: Double.infinity) }
        return bs_lambert_w0_f(x)
    }
    
    /// W−1(x) for `Float`. Requires −1/e ≤ x < 0.
    @inlinable static func lambertWm1(_ x: Float) throws -> Float {
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        let minX = -1.0 as Float / bs_const_e_f()
        guard x >= minX && x < 0 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: Double(minX), max: 0.0) }
        return bs_lambert_wm1_f(x)
    }
    
    // MARK: - Float80 overloads (x86_64)
    
#if arch(x86_64)
    /// W0(x) for `Float80` (x86_64 only). Requires x ≥ −1/e.
    @inlinable static func lambertW0(_ x: Float80) throws -> Float80 {
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        let minX = -1.0 as Float80 / bs_const_e_l()
        guard x >= minX else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: Double(minX), max: Double.infinity) }
        return bs_lambert_w0_l(x)
    }
    
    /// W−1(x) for `Float80` (x86_64 only). Requires −1/e ≤ x < 0.
    @inlinable static func lambertWm1(_ x: Float80) throws -> Float80 {
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        let minX = -1.0 as Float80 / bs_const_e_l()
        guard x >= minX && x < 0 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: Double(minX), max: 0.0) }
        return bs_lambert_wm1_l(x)
    }
#endif
}
