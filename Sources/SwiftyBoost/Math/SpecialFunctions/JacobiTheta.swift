//  Created by Volker Thieme 2025.
//  License: MIT (see project root)
//
//  Swift wrappers for the Jacobi theta functions θ₁–θ₄. Each helper exposes
//  both nome-based (q) and half-period ratio (τ) parameterizations and delegates
//  the evaluation to Boost.Math while mirroring SwiftyBoost’s validation rules.
//
//  Overview
//  --------
//  The Jacobi theta functions θ₁, θ₂, θ₃, and θ₄ are special functions that arise
//  in elliptic function theory, modular forms, lattice sums, and heat kernel
//  solutions. This file provides thin, type-generic Swift wrappers around the
//  corresponding Boost.Math implementations, with consistent parameter validation
//  and Swift-friendly errors.
//
//  Parameterizations
//  -----------------
//  - Nome form: θₖ(x, q) with real q constrained to 0 < q < 1.
//    This is the “convergent q-series” regime used for real-valued evaluation.
//  - Tau form: θₖ(x, τ), where Boost.Math accepts a real τ and internally relates
//    it to a nome in (0, 1). In the real-valued API, τ must be finite.
//    Notes:
//      * In complex analysis, τ is typically complex with Im(τ) > 0. These
//        wrappers expose the real-valued Boost parameterization only.
//      * Boost.Math defines a real-valued τ that corresponds to a real nome in
//        (0, 1) internally (see Boost’s jacobi theta documentation for the exact
//        mapping used in the implementation).
//
//  Domain and Validation
//  ---------------------
//  - For θₖ(x, q): x must be finite; q must be finite and satisfy 0 < q < 1.
//  - For θₖ(x, τ): x and τ must be finite (real-valued API).
//  - Violations raise SpecialFunctionError with parameterNotInDomain or
//    parameterNotFinite, mirroring SwiftyBoost’s conventions.
//
//  Precision and Architecture Notes
//  --------------------------------
//  - Generic T overloads dispatch to the double-precision Boost backend and cast
//    results back to T.
//  - Dedicated Float overloads call the float backend directly.
//  - On x86 architectures, Float80 overloads are provided and call the long
//    double backend.
//  - Numerical sensitivity increases as q approaches 1 from below (or when τ is
//    small in the τ-parameterization). Expect slower convergence and potential
//    loss of precision near these edges.
//
//  Performance
//  -----------
//  - Calls are direct shims into Boost.Math and are O(N) in the number of series
//    terms used by Boost internally. There is no additional allocation.
//  - Prefer using the overload that matches your working precision to avoid
//    unnecessary casts.
//
//  See Also
//  --------
//  - Boost.Math documentation for jacobi_theta1..4 and jacobi_theta?tau.
//  - SpecialFunctions.jacobiTheta1..4 and SpecialFunctions.jacobiTheta?Tau in
//    this module for consistent usage patterns.
//  - Elliptic integrals and modular functions where theta functions commonly
//    appear.
//

import SwiftyBoostPrelude

public extension SpecialFunctions {
    
    /// Jacobi theta function θ₁(x, q) parameterized by the nome `q`.
    ///
    /// θ₁ is an odd function of `x` and is defined here for real `x` and real
    /// nome `q` with 0 < q < 1. This overload is generic over real
    /// BinaryFloatingPoint types and delegates to the double-precision Boost
    /// backend, casting back to `T`.
    ///
    /// - Parameters:
    ///   - x: Real argument. Must be finite.
    ///   - q: Real nome in the open interval (0, 1). Must be finite.
    /// - Returns: θ₁(x, q) as `T`.
    /// - Throws: ``SpecialFunctionError.parameterNotInDomain`` when `q` falls outside `(0, 1)`,
    ///           or ``SpecialFunctionError.parameterNotFinite`` if `x` or `q` is not finite.
    ///
    /// Notes:
    /// - Numerical difficulty increases as `q → 1⁻`.
    /// - For periodicity/parity properties and series definitions, consult standard
    ///   references or Boost.Math documentation.
    ///
    /// Delegates to Boost.Math `jacobi_theta1` and enforces SwiftyBoost’s parameter checks.
    @inlinable static func jacobiTheta1<T: Real & BinaryFloatingPoint & Sendable>(_ x: T, q: T) throws -> T {
        guard q > 0 && q < 1 else {
            throw SpecialFunctionError.parameterNotInDomain(name: "q", value: q)
        }
        guard q.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "q", value: q) }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        return T(bs_jacobi_theta1_d(Double(x), Double(q)))
    }
    
    /// Jacobi theta θ₁(x, q) for `Float`.
    ///
    /// - Parameters:
    ///   - x: Real argument. Must be finite.
    ///   - q: Real nome with `0 < q < 1`. Must be finite.
    /// - Returns: θ₁(x, q) as `Float`.
    /// - Throws: See generic overload.
    ///
    /// Calls Boost.Math `jacobi_theta1` (float backend).
    @inlinable static func jacobiTheta1(_ x: Float, q: Float) throws -> Float {
        guard q > 0 && q < 1 else {
            throw SpecialFunctionError.parameterNotInDomain(name: "q", value: q)
        }
        guard q.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "q", value: q) }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        return bs_jacobi_theta1_f(Float(x), Float(q))
    }
    
#if arch(x86_64) || arch(i386)
    /// Jacobi theta θ₁(x, q) for `Float80` (x86 only).
    ///
    /// - Parameters:
    ///   - x: Real argument. Must be finite.
    ///   - q: Real nome with `0 < q < 1`. Must be finite.
    /// - Returns: θ₁(x, q) as `Float80`.
    /// - Throws: See generic overload.
    ///
    /// Calls Boost.Math `jacobi_theta1` (long double backend).
    @inlinable static func jacobiTheta1(_ x: Float, q: Float) throws -> Float80 {
        guard q > 0 && q < 1 else {
            throw SpecialFunctionError.parameterNotInDomain(name: "q", value: q)
        }
        guard q.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "q", value: q) }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        return bs_jacobi_theta1_l(Float80(x), Float80(q))
    }
#endif
    
    /// Jacobi theta function θ₁(x, τ) parameterized by the half-period ratio `τ`.
    ///
    /// This overload exposes the real-τ variant provided by Boost.Math. Both `x`
    /// and `tau` must be finite. Internally, Boost relates a real `tau` to a
    /// nome in (0, 1) in order to evaluate the function.
    ///
    /// - Parameters:
    ///   - x: Real argument. Must be finite.
    ///   - tau: Real half-period ratio τ. Must be finite.
    /// - Returns: θ₁(x, τ) as `T`.
    /// - Throws: ``SpecialFunctionError.parameterNotFinite`` when `x` or `tau` is not finite.
    ///
    /// Wraps Boost.Math `jacobi_theta1tau` with SwiftyBoost validation.
    @inlinable static func jacobiTheta1Tau<T: Real & BinaryFloatingPoint & Sendable>(_ x: T, tau: T) throws -> T {
        guard tau.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "tau", value: tau) }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        return T(bs_jacobi_theta1tau_d(Double(x), Double(tau)))
    }

    /// Jacobi theta θ₁(x, τ) for `Float`.
    ///
    /// - Parameters:
    ///   - x: Real argument. Must be finite.
    ///   - tau: Real half-period ratio τ. Must be finite.
    /// - Returns: θ₁(x, τ) as `Float`.
    /// - Throws: See generic overload.
    @inlinable static func jacobiTheta1Tau(_ x: Float, tau: Float) throws -> Float {
        guard tau.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "tau", value: tau) }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        return bs_jacobi_theta1tau_f(Float(x), Float(tau))
    }

    #if arch(x86_64) || arch(i386)
    /// Jacobi theta θ₁(x, τ) for `Float80` (x86 only).
    ///
    /// - Parameters:
    ///   - x: Real argument. Must be finite.
    ///   - tau: Real half-period ratio τ. Must be finite.
    /// - Returns: θ₁(x, τ) as `Float80`.
    /// - Throws: See generic overload.
    @inlinable static func jacobiTheta1Tau(_ x: Float, tau: Float) throws -> Float80 {
        guard tau.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "tau", value: tau) }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        return bs_jacobi_theta1tau_l(Float80(x), Float80(tau))
    }
    #endif
    /// Jacobi theta function θ₂(x, q) parameterized by the nome `q`.
    ///
    /// θ₂ is an even/odd-shifted companion to θ₁ and is defined here for real `x`
    /// and real nome `q` with 0 < q < 1. This generic overload uses the
    /// double-precision backend and casts to `T`.
    ///
    /// - Parameters:
    ///   - x: Real argument. Must be finite.
    ///   - q: Nome `q` with `0 < q < 1`. Must be finite.
    /// - Returns: θ₂(x, q) as `T`.
    /// - Throws: ``SpecialFunctionError.parameterNotInDomain`` when `q` falls outside `(0, 1)`,
    ///           or ``SpecialFunctionError.parameterNotFinite`` if `x` or `q` is not finite.
    ///
    /// Delegates to Boost.Math `jacobi_theta2`.
    @inlinable static func jacobiTheta2<T: Real & BinaryFloatingPoint & Sendable>(_ x: T, q: T) throws -> T {
        guard q > 0 && q < 1 else {
            throw SpecialFunctionError.parameterNotInDomain(name: "q", value: q)
        }
        guard q.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "q", value: q) }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        return T(bs_jacobi_theta2_d(Double(x), Double(q)))
    }
    
    /// Jacobi theta θ₂(x, q) for `Float`.
    ///
    /// - Parameters:
    ///   - x: Real argument. Must be finite.
    ///   - q: Real nome with `0 < q < 1`. Must be finite.
    /// - Returns: θ₂(x, q) as `Float`.
    /// - Throws: See generic overload.
    @inlinable static func jacobiTheta2(_ x: Float, q: Float) throws -> Float {
        guard q > 0 && q < 1 else {
            throw SpecialFunctionError.parameterNotInDomain(name: "q", value: q)
        }
        guard q.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "q", value: q) }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        return bs_jacobi_theta2_f(Float(x), Float(q))
    }
    
#if arch(x86_64) || arch(i386)
    /// Jacobi theta θ₂(x, q) for `Float80` (x86 only).
    ///
    /// - Parameters:
    ///   - x: Real argument. Must be finite.
    ///   - q: Real nome with `0 < q < 1`. Must be finite.
    /// - Returns: θ₂(x, q) as `Float80`.
    /// - Throws: See generic overload.
    @inlinable static func jacobiTheta2(_ x: Float, q: Float) throws -> Float80 {
        guard q > 0 && q < 1 else {
            throw SpecialFunctionError.parameterNotInDomain(name: "q", value: q)
        }
        guard q.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "q", value: q) }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        return bs_jacobi_theta2_l(Float80(x), Float80(q))
    }
#endif
    
    /// Jacobi theta function θ₂(x, τ) parameterized by the half-period ratio `τ`.
    ///
    /// - Parameters:
    ///   - x: Real argument. Must be finite.
    ///   - tau: Real half-period ratio τ. Must be finite.
    /// - Returns: θ₂(x, τ) as `T`.
    /// - Throws: ``SpecialFunctionError.parameterNotFinite`` when `x` or `tau` is not finite.
    ///
    /// Wraps Boost.Math `jacobi_theta2tau`.
    @inlinable static func jacobiTheta2Tau<T: Real & BinaryFloatingPoint & Sendable>(_ x: T, tau: T) throws -> T {
        guard tau.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "tau", value: tau) }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        return T(bs_jacobi_theta2tau_d(Double(x), Double(tau)))
    }

    /// Jacobi theta θ₂(x, τ) for `Float`.
    ///
    /// - Parameters:
    ///   - x: Real argument. Must be finite.
    ///   - tau: Real half-period ratio τ. Must be finite.
    /// - Returns: θ₂(x, τ) as `Float`.
    /// - Throws: See generic overload.
    @inlinable static func jacobiTheta2Tau(_ x: Float, tau: Float) throws -> Float {
        guard tau.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "tau", value: tau) }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        return bs_jacobi_theta2tau_f(Float(x), Float(tau))
    }

    #if arch(x86_64) || arch(i386)
    /// Jacobi theta θ₂(x, τ) for `Float80` (x86 only).
    ///
    /// - Parameters:
    ///   - x: Real argument. Must be finite.
    ///   - tau: Real half-period ratio τ. Must be finite.
    /// - Returns: θ₂(x, τ) as `Float80`.
    /// - Throws: See generic overload.
    @inlinable static func jacobiTheta2Tau(_ x: Float, tau: Float) throws -> Float80 {
        guard tau.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "tau", value: tau) }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        return bs_jacobi_theta2tau_l(Float80(x), Float80(tau))
    }
    #endif

    /// Jacobi theta function θ₃(x, q) parameterized by the nome `q`.
    ///
    /// θ₃ is an even function of `x`. This overload requires finite real `x`
    /// and real nome `q` with 0 < q < 1. It uses the double-precision backend
    /// and casts back to `T`.
    ///
    /// - Parameters:
    ///   - x: Real argument. Must be finite.
    ///   - q: Nome `q` with `0 < q < 1`. Must be finite.
    /// - Returns: θ₃(x, q) as `T`.
    /// - Throws: ``SpecialFunctionError.parameterNotInDomain`` when `q` falls outside `(0, 1)`,
    ///           or ``SpecialFunctionError.parameterNotFinite`` if `x` or `q` is not finite.
    ///
    /// Delegates to Boost.Math `jacobi_theta3`.
    @inlinable static func jacobiTheta3<T: Real & BinaryFloatingPoint & Sendable>(_ x: T, q: T) throws -> T {
        guard q > 0 && q < 1 else {
            throw SpecialFunctionError.parameterNotInDomain(name: "q", value: q)
        }
        guard q.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "q", value: q) }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        return T(bs_jacobi_theta3_d(Double(x), Double(q)))
    }
    
    /// Jacobi theta θ₃(x, q) for `Float`.
    ///
    /// - Parameters:
    ///   - x: Real argument. Must be finite.
    ///   - q: Real nome with `0 < q < 1`. Must be finite.
    /// - Returns: θ₃(x, q) as `Float`.
    /// - Throws: See generic overload.
    @inlinable static func jacobiTheta3(_ x: Float, q: Float) throws -> Float {
        guard q > 0 && q < 1 else {
            throw SpecialFunctionError.parameterNotInDomain(name: "q", value: q)
        }
        guard q.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "q", value: q) }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        return bs_jacobi_theta3_f(Float(x), Float(q))
    }
    
#if arch(x86_64) || arch(i386)
    /// Jacobi theta θ₃(x, q) for `Float80` (x86 only).
    ///
    /// - Parameters:
    ///   - x: Real argument. Must be finite.
    ///   - q: Real nome with `0 < q < 1`. Must be finite.
    /// - Returns: θ₃(x, q) as `Float80`.
    /// - Throws: See generic overload.
    @inlinable static func jacobiTheta3(_ x: Float, q: Float) throws -> Float80 {
        guard q > 0 && q < 1 else {
            throw SpecialFunctionError.parameterNotInDomain(name: "q", value: q)
        }
        guard q.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "q", value: q) }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        return bs_jacobi_theta3_l(Float80(x), Float80(q))
    }
#endif
    
    /// Jacobi theta function θ₃(x, τ) parameterized by the half-period ratio `τ`.
    ///
    /// - Parameters:
    ///   - x: Real argument. Must be finite.
    ///   - tau: Real half-period ratio τ. Must be finite.
    /// - Returns: θ₃(x, τ) as `T`.
    /// - Throws: ``SpecialFunctionError.parameterNotFinite`` when `x` or `tau` is not finite.
    ///
    /// Wraps Boost.Math `jacobi_theta3tau`.
    @inlinable static func jacobiTheta3Tau<T: Real & BinaryFloatingPoint & Sendable>(_ x: T, tau: T) throws -> T {
        guard tau.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "tau", value: tau) }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        return T(bs_jacobi_theta3tau_d(Double(x), Double(tau)))
    }

    /// Jacobi theta θ₃(x, τ) for `Float`.
    ///
    /// - Parameters:
    ///   - x: Real argument. Must be finite.
    ///   - tau: Real half-period ratio τ. Must be finite.
    /// - Returns: θ₃(x, τ) as `Float`.
    /// - Throws: See generic overload.
    @inlinable static func jacobiTheta3Tau(_ x: Float, tau: Float) throws -> Float {
        guard tau.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "tau", value: tau) }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        return bs_jacobi_theta3tau_f(Float(x), Float(tau))
    }

    #if arch(x86_64) || arch(i386)
    /// Jacobi theta θ₃(x, τ) for `Float80` (x86 only).
    ///
    /// - Parameters:
    ///   - x: Real argument. Must be finite.
    ///   - tau: Real half-period ratio τ. Must be finite.
    /// - Returns: θ₃(x, τ) as `Float80`.
    /// - Throws: See generic overload.
    @inlinable static func jacobiTheta3Tau(_ x: Float, tau: Float) throws -> Float80 {
        guard tau.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "tau", value: tau) }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        return bs_jacobi_theta3tau_l(Float80(x), Float80(tau))
    }
    #endif

    
    /// Jacobi theta function θ₄(x, q) parameterized by the nome `q`.
    ///
    /// θ₄ is an even function of `x`. This overload requires finite real `x` and
    /// real nome `q` with 0 < q < 1. It uses the double-precision backend and
    /// casts back to `T`.
    ///
    /// - Parameters:
    ///   - x: Real argument. Must be finite.
    ///   - q: Nome `q` with `0 < q < 1`. Must be finite.
    /// - Returns: θ₄(x, q) as `T`.
    /// - Throws: ``SpecialFunctionError.parameterNotInDomain`` when `q` falls outside `(0, 1)`,
    ///           or ``SpecialFunctionError.parameterNotFinite`` if `x` or `q` is not finite.
    ///
    /// Delegates to Boost.Math `jacobi_theta4`.
    @inlinable static func jacobiTheta4<T: Real & BinaryFloatingPoint & Sendable>(_ x: T, q: T) throws -> T {
        guard q > 0 && q < 1 else {
            throw SpecialFunctionError.parameterNotInDomain(name: "q", value: q)
        }
        guard q.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "q", value: q) }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        return T(bs_jacobi_theta4_d(Double(x), Double(q)))
    }
    
    /// Jacobi theta θ₄(x, q) for `Float`.
    ///
    /// - Parameters:
    ///   - x: Real argument. Must be finite.
    ///   - q: Real nome with `0 < q < 1`. Must be finite.
    /// - Returns: θ₄(x, q) as `Float`.
    /// - Throws: See generic overload.
    @inlinable static func jacobiTheta4(_ x: Float, q: Float) throws -> Float {
        guard q > 0 && q < 1 else {
            throw SpecialFunctionError.parameterNotInDomain(name: "q", value: q)
        }
        guard q.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "q", value: q) }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        return bs_jacobi_theta4_f(Float(x), Float(q))
    }
    
#if arch(x86_64) || arch(i386)
    /// Jacobi theta θ₄(x, q) for `Float80` (x86 only).
    ///
    /// - Parameters:
    ///   - x: Real argument. Must be finite.
    ///   - q: Real nome with `0 < q < 1`. Must be finite.
    /// - Returns: θ₄(x, q) as `Float80`.
    /// - Throws: See generic overload.
    @inlinable static func jacobiTheta4(_ x: Float, q: Float) throws -> Float80 {
        guard q > 0 && q < 1 else {
            throw SpecialFunctionError.parameterNotInDomain(name: "q", value: q)
        }
        guard q.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "q", value: q) }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        return bs_jacobi_theta4_l(Float80(x), Float80(q))
    }
#endif
    
    /// Jacobi theta function θ₄(x, τ) parameterized by the half-period ratio `τ`.
    ///
    /// - Parameters:
    ///   - x: Real argument. Must be finite.
    ///   - tau: Real half-period ratio τ. Must be finite.
    /// - Returns: θ₄(x, τ) as `T`.
    /// - Throws: ``SpecialFunctionError.parameterNotFinite`` when `x` or `tau` is not finite.
    ///
    /// Wraps Boost.Math `jacobi_theta4tau`.
    @inlinable static func jacobiTheta4Tau<T: Real & BinaryFloatingPoint & Sendable>(_ x: T, tau: T) throws -> T {
        guard tau.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "tau", value: tau) }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        return T(bs_jacobi_theta4tau_d(Double(x), Double(tau)))
    }

    /// Jacobi theta θ₄(x, τ) for `Float`.
    ///
    /// - Parameters:
    ///   - x: Real argument. Must be finite.
    ///   - tau: Real half-period ratio τ. Must be finite.
    /// - Returns: θ₄(x, τ) as `Float`.
    /// - Throws: See generic overload.
    @inlinable static func jacobiTheta4Tau(_ x: Float, tau: Float) throws -> Float {
        guard tau.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "tau", value: tau) }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        return bs_jacobi_theta4tau_f(Float(x), Float(tau))
    }

    #if arch(x86_64) || arch(i386)
    /// Jacobi theta θ₄(x, τ) for `Float80` (x86 only).
    ///
    /// - Parameters:
    ///   - x: Real argument. Must be finite.
    ///   - tau: Real half-period ratio τ. Must be finite.
    /// - Returns: θ₄(x, τ) as `Float80`.
    /// - Throws: See generic overload.
    @inlinable static func jacobiTheta4Tau(_ x: Float, tau: Float) throws -> Float80 {
        guard tau.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "tau", value: tau) }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        return bs_jacobi_theta4tau_l(Float80(x), Float80(tau))
    }
    #endif
}
