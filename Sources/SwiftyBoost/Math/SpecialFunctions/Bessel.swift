//
//  Bessel.swift
//  Math/SpecialFunctions
//
//  Thin Swift wrappers around Boost.Math cylindrical and modified Bessel functions,
//  with argument validation and a Swifty error model.
//
//  Functions included:
//  - J_v(x): Cylindrical Bessel function of the first kind.
//  - Y_v(x): Cylindrical Bessel function of the second kind (Neumann function).
//  - I_v(x): Modified Bessel function of the first kind.
//  - K_v(x): Modified Bessel function of the second kind (Macdonald function).
//
//  Overloads are provided for generic BinaryFloatingPoint plus fast paths
//  for Float and (on x86_64) Float80. All implementations delegate to
//  CBoostBridge (Boost.Math) for numerical evaluation.
//
//  Domain notes (real-valued):
//  - J_v(x): defined for all finite v and x; if x < 0 and v is non-integer, the result is complex.
//            Boost.Math enforces this and will error for x < 0 with non-integer v.
//            We do not pre-validate the sign here to allow integer v (including large values).
//  - Y_v(x): requires x > 0; the real-valued function has a branch cut and is complex for x ≤ 0.
//  - I_v(x): defined for all finite v and x; if x < 0 and v is non-integer, the result is complex,
//            but Boost.Math reflects using parity for integer v. As with J_v, we do not pre-validate sign.
//  - K_v(x): requires x > 0; has a singularity at x = 0.
//
//  Error model:
//  - Inputs are checked for finiteness. Non-finite arguments throw SpecialFunctionError.parameterNotFinite.
//  - Functions with real-domain restrictions on x (Y_v, K_v) throw SpecialFunctionError.parameterNotPositive
//    when x ≤ 0.
//  - Further domain checks (e.g., J_v/I_v with x < 0 and non-integer v) are enforced by Boost.Math;
//    if violated, Boost’s policy raises a domain error which is surfaced as a NaN/inf by the C layer.
//    The wrappers here focus on common real-domain constraints and finiteness.
//
//  References:
//  - NIST DLMF: https://dlmf.nist.gov/10
//  - Boost.Math Bessel functions documentation:
//    https://www.boost.org/doc/libs/release/libs/math/doc/html/math_toolkit/special/bessel.html
//

import CBoostBridge

public extension SpecialFunctions {
    
    
    // MARK: - Generic BinaryFloatingPoint overloads
    
    /// Cylindrical Bessel function of the first kind J_v(x).
    ///
    /// Definition (real-valued):
    /// - J_v(x) is real for real order v and real x. For x < 0, a real result is only obtained when v is an integer
    ///   (odd/even parity applies). Otherwise the result is complex and outside the scope of this real-valued API.
    ///
    /// Behavior:
    /// - Delegates to Boost.Math cyl_bessel_j(v, x).
    /// - This wrapper validates that both v and x are finite.
    /// - No pre-check on the sign of x: Boost.Math will reflect for integer v, and raise a domain error
    ///   for non-integer v with x < 0.
    ///
    /// Parameters:
    /// - v: Order of the function (any finite real).
    /// - x: Argument (any finite real; see note above for x < 0).
    ///
    /// Returns:
    /// - J_v(x) as T.
    ///
    /// Throws:
    /// - SpecialFunctionError.parameterNotFinite(name: "v") if v is NaN or ±∞.
    /// - SpecialFunctionError.parameterNotFinite(name: "x") if x is NaN or ±∞.
    ///
    /// See also:
    /// - NIST DLMF §10.2–10.4 for properties and series/recurrence relations.
    @inlinable static func besselJ<T: BinaryFloatingPoint>(v: T, x: T) throws -> T {
        let dv = D(v), dx = D(x)
        guard dv.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "v") }
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return T(bs_cyl_bessel_j(dv, dx))
    }
    
    /// Cylindrical Bessel function of the second kind Y_v(x) (Neumann function).
    ///
    /// Domain:
    /// - Requires x > 0 for a real-valued result. For x ≤ 0 the value is complex or singular.
    ///   This wrapper enforces x > 0 and throws otherwise.
    ///
    /// Behavior:
    /// - Delegates to Boost.Math cyl_neumann(v, x).
    /// - Validates that v and x are finite and that x > 0.
    ///
    /// Parameters:
    /// - v: Order of the function (any finite real).
    /// - x: Argument (must satisfy x > 0).
    ///
    /// Returns:
    /// - Y_v(x) as T.
    ///
    /// Throws:
    /// - SpecialFunctionError.parameterNotFinite(name: "v") if v is NaN or ±∞.
    /// - SpecialFunctionError.parameterNotFinite(name: "x") if x is NaN or ±∞.
    /// - SpecialFunctionError.parameterNotPositive(name: "x") if x ≤ 0.
    ///
    /// Notes:
    /// - Also known as the Neumann function or Bessel function of the second kind.
    /// - Has a logarithmic singularity at x → 0⁺ for many orders.
    ///
    /// See also:
    /// - NIST DLMF §10.2–10.4 and §10.7 for behavior and relations.
    @inlinable static func besselY<T: BinaryFloatingPoint>(v: T, x: T) throws -> T {
        let dv = D(v), dx = D(x)
        guard dv.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "v") }
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        guard dx > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "x") }
        return T(bs_cyl_neumann(dv, dx))
    }
    
    /// Modified Bessel function of the first kind I_v(x).
    ///
    /// Definition (real-valued):
    /// - I_v(x) is real for real order v and real x. For x < 0 and non-integer v, the analytic continuation is complex;
    ///   Boost.Math applies parity for integer v. This wrapper does not pre-restrict the sign of x.
    ///
    /// Behavior:
    /// - Delegates to Boost.Math cyl_bessel_i(v, x).
    /// - Validates that v and x are finite.
    ///
    /// Parameters:
    /// - v: Order of the function (any finite real).
    /// - x: Argument (any finite real; see note above for x < 0).
    ///
    /// Returns:
    /// - I_v(x) as T.
    ///
    /// Throws:
    /// - SpecialFunctionError.parameterNotFinite(name: "v") if v is NaN or ±∞.
    /// - SpecialFunctionError.parameterNotFinite(name: "x") if x is NaN or ±∞.
    ///
    /// See also:
    /// - NIST DLMF §10.25–10.27.
    @inlinable static func modifiedBesselI<T: BinaryFloatingPoint>(v: T, x: T) throws -> T {
        let dv = D(v), dx = D(x)
        guard dv.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "v") }
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return T(bs_cyl_bessel_i(dv, dx))
    }
    
    /// Modified Bessel function of the second kind K_v(x) (Macdonald function).
    ///
    /// Domain:
    /// - Requires x > 0 for a real-valued result. K_v(x) has a singularity at x = 0 and is complex for x ≤ 0.
    ///
    /// Behavior:
    /// - Delegates to Boost.Math cyl_bessel_k(v, x).
    /// - Validates that v and x are finite and that x > 0.
    ///
    /// Parameters:
    /// - v: Order of the function (any finite real).
    /// - x: Argument (must satisfy x > 0).
    ///
    /// Returns:
    /// - K_v(x) as T.
    ///
    /// Throws:
    /// - SpecialFunctionError.parameterNotFinite(name: "v") if v is NaN or ±∞.
    /// - SpecialFunctionError.parameterNotFinite(name: "x") if x is NaN or ±∞.
    /// - SpecialFunctionError.parameterNotPositive(name: "x") if x ≤ 0.
    ///
    /// Notes:
    /// - Rapidly decays for large x and is often used in boundary-value problems.
    /// - Sometimes denoted as the Macdonald function.
    ///
    /// See also:
    /// - NIST DLMF §10.25–10.31.
    @inlinable static func modifiedBesselK<T: BinaryFloatingPoint>(v: T, x: T) throws -> T {
        let dv = D(v), dx = D(x)
        guard dv.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "v") }
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        guard dx > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "x") }
        return T(bs_cyl_bessel_k(dv, dx))
    }
    
    // MARK: - Float overloads
    
    /// Cylindrical Bessel function of the first kind J_v(x) for Float.
    @inlinable static func besselJ_f(v: Float, x: Float) throws -> Float {
        guard v.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "v") }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return bs_cyl_bessel_j_f(v, x)
    }
    
    /// Cylindrical Bessel function of the second kind Y_v(x) for Float.
    ///
    /// Domain: requires x > 0.
    @inlinable static func besselY_f(v: Float, x: Float) throws -> Float {
        guard v.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "v") }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        guard x > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "x") }
        return bs_cyl_neumann_f(v, x)
    }
    
    /// Modified Bessel function of the first kind I_v(x) for Float.
    @inlinable static func modifiedBessel0I_f(v: Float, x: Float) throws -> Float {
        guard v.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "v") }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return bs_cyl_bessel_i_f(v, x)
    }
    
    /// Modified Bessel function of the second kind K_v(x) for Float.
    ///
    /// Domain: requires x > 0.
    @inlinable static func modifiedBesselK_f(v: Float, x: Float) throws -> Float {
        guard v.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "v") }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        guard x > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "x") }
        return bs_cyl_bessel_k_f(v, x)
    }
    
    // MARK: - Float80 overloads (x86_64)
    
#if arch(x86_64)
    /// Cylindrical Bessel function of the first kind J_v(x) for Float80 (x86_64 only).
    @inlinable static func besselJ_l(v: Float80, x: Float80) throws -> Float80 {
        guard v.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "v") }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return bs_cyl_bessel_j_l(v, x)
    }
    
    /// Cylindrical Bessel function of the second kind Y_v(x) for Float80 (x86_64 only).
    ///
    /// Domain: requires x > 0.
    @inlinable static func besselY_l(v: Float80, x: Float80) throws -> Float80 {
        guard v.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "v") }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        guard x > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "x") }
        return bs_cyl_neumann_l(v, x)
    }
    
    /// Modified Bessel function of the first kind I_v(x) for Float80 (x86_64 only).
    @inlinable static func modifiedBesselI_l(v: Float80, x: Float80) throws -> Float80 {
        guard v.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "v") }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return bs_cyl_bessel_i_l(v, x)
    }
    
    /// Modified Bessel function of the second kind K_v(x) for Float80 (x86_64 only).
    ///
    /// Domain: requires x > 0.
    @inlinable static func modifiedBesselK_l(v: Float80, x: Float80) throws -> Float80 {
        guard v.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "v") }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        guard x > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "x") }
        return bs_cyl_bessel_k_l(v, x)
    }
#endif
    
}

