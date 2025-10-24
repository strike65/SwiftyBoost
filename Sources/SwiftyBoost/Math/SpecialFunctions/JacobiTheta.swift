//  Created by Volker Thieme 2025.
//  License: MIT (see project root)
//
//  Swift wrappers for the Jacobi theta functions θ₁–θ₄. Each helper exposes
//  both nome-based (q) and half-period ratio (τ) parameterizations and delegates
//  the evaluation to Boost.Math while mirroring SwiftyBoost’s validation rules.
//
//  Overview
//  --------
//  The Jacobi theta functions θ₁, θ₂, θ₃, and θ₄ are classical special functions
//  central to elliptic functions, modular forms, lattice sums, and heat kernels.
//  This file provides thin Swift wrappers around Boost.Math’s real-valued
//  implementations with consistent parameter validation and Swift-friendly errors.
//
//  Parameterizations
//  -----------------
//  - Nome form: θₖ(x, q) with real q constrained to 0 < q < 1.
//    This is the standard convergent q-series regime for real evaluation.
//  - Tau form: θₖ(x, τ), using Boost.Math’s real-τ API. In complex analysis,
//    τ is typically complex with Im(τ) > 0. Boost’s real-τ API treats τ as the
//    imaginary component (see Boost’s header notes). For real τ, Boost relates
//    the parameters by:
//      q = exp(-π · τ)
//    These wrappers expose the real-valued τ API and require finiteness checks
//    only; domain restrictions beyond finiteness are delegated to Boost.
//
//  Domains and validation
//  ----------------------
//  - θₖ(x, q): requires x finite; q finite and strictly 0 < q < 1.
//    • On violation, throws SpecialFunctionError.parameterNotInDomain("q", q).
//    • Non-finite x or q throws SpecialFunctionError.parameterNotFinite.
//  - θₖ(x, τ): requires x and τ finite (real-valued API).
//    • Non-finite x or τ throws SpecialFunctionError.parameterNotFinite.
//    • Additional domain handling (e.g., τ ≤ 0) is performed by Boost; these
//      wrappers intentionally do not preempt Boost’s behavior to preserve parity
//      with the C bridge shims and tests.
//
//  Precision and architecture
//  --------------------------
//  - Generic T overloads call the double-precision backend and cast to T.
//  - Float overloads call the float backend directly.
//  - Float80 overloads (when available on x86 targets) call the long double backend.
//  - As q → 1⁻ (or for small τ in the τ-parameterization), series convergence
//    slows and sensitivity increases; expect potential precision loss near edges.
//
//  Performance
//  -----------
//  - Each call delegates directly to Boost.Math with no extra allocation.
//  - Prefer the overload matching your working precision (Float, Double, Float80)
//    to avoid unnecessary casts.
//
//  Usage examples (Double)
//  -----------------------
//  - Using the nome q (recommended when you already have q):
//      let x: Double = 0.75
//      let q: Double = 0.25
//      let t3 = try SpecialFunctions.jacobiTheta3(x, q: q)
//  - Using τ (recommended when q is naturally expressed as exp(-π·τ)):
//      let x: Double = 0.75
//      let tau: Double = 1.2
//      let t4 = try SpecialFunctions.jacobiTheta4Tau(x, tau: tau)
//
//  See also
//  --------
//  - Boost.Math jacobi_theta1..4 and jacobi_theta?tau reference
//  - DLMF §20 (https://dlmf.nist.gov/20) and MathWorld “JacobiThetaFunctions”
//  - Elliptic integrals and Jacobi elliptic functions (relations via θₖ)
//
//  Note: These wrappers are real-valued. For complex arguments or full complex τ,
//  consult Boost’s complex-number APIs directly.
//

import SwiftyBoostPrelude

public extension SpecialFunctions {
    
    /// Jacobi theta function θ₁(x, q) parameterized by the nome `q`.
    ///
    /// θ₁ is an odd function of `x`. This overload is generic over real
    /// BinaryFloatingPoint types. It dispatches to the double-precision Boost
    /// backend and casts back to `T`.
    ///
    /// - Parameters:
    ///   - x: Real argument. Must be finite.
    ///   - q: Real nome in the open interval (0, 1). Must be finite.
    /// - Returns: θ₁(x, q) as `T`.
    /// - Throws:
    ///   - `SpecialFunctionError.parameterNotInDomain("q", q)` if `q ∉ (0, 1)`.
    ///   - `SpecialFunctionError.parameterNotFinite` if `x` or `q` is not finite.
    ///
    /// Notes:
    /// - Numerical difficulty increases as `q → 1⁻`.
    /// - For series definitions and identities, see DLMF §20.2 and §20.7.
    ///
    /// Delegates to Boost.Math `jacobi_theta1`.
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
    /// - Throws: See the generic overload for error behavior.
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
    /// - Throws: See the generic overload for error behavior.
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
    /// Boost’s real-τ API relates `q` and `τ` by `q = exp(-π · τ)` for real τ.
    /// This wrapper enforces finiteness only; any additional domain handling is
    /// delegated to Boost to preserve behavior parity with the bridge.
    ///
    /// - Parameters:
    ///   - x: Real argument. Must be finite.
    ///   - tau: Real half-period ratio τ. Must be finite.
    /// - Returns: θ₁(x, τ) as `T`.
    /// - Throws: `SpecialFunctionError.parameterNotFinite` when `x` or `tau` is not finite.
    ///
    /// Wraps Boost.Math `jacobi_theta1tau`.
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
    /// - Throws: See the generic overload for error behavior.
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
    /// - Throws: See the generic overload for error behavior.
    @inlinable static func jacobiTheta1Tau(_ x: Float, tau: Float) throws -> Float80 {
        guard tau.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "tau", value: tau) }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        return bs_jacobi_theta1tau_l(Float80(x), Float80(tau))
    }
    #endif

    /// Jacobi theta function θ₂(x, q) parameterized by the nome `q`.
    ///
    /// θ₂ is closely related to θ₁ via half-period shifts. This overload is
    /// generic over real BinaryFloatingPoint types and dispatches to the
    /// double-precision backend.
    ///
    /// - Parameters:
    ///   - x: Real argument. Must be finite.
    ///   - q: Nome `q` with `0 < q < 1`. Must be finite.
    /// - Returns: θ₂(x, q) as `T`.
    /// - Throws:
    ///   - `SpecialFunctionError.parameterNotInDomain("q", q)` if `q ∉ (0, 1)`.
    ///   - `SpecialFunctionError.parameterNotFinite` if `x` or `q` is not finite.
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
    /// - Throws: See the generic overload for error behavior.
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
    /// - Throws: See the generic overload for error behavior.
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
    /// Uses Boost’s real-τ API (with `q = exp(-π · τ)` internally). Finiteness
    /// is enforced here; further domain handling is delegated to Boost.
    ///
    /// - Parameters:
    ///   - x: Real argument. Must be finite.
    ///   - tau: Half-period ratio τ. Must be finite.
    /// - Returns: θ₂(x, τ) as `T`.
    /// - Throws: `SpecialFunctionError.parameterNotFinite` when `x` or `tau` is not finite.
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
    /// - Throws: See the generic overload for error behavior.
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
    /// - Throws: See the generic overload for error behavior.
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
    /// - Throws:
    ///   - `SpecialFunctionError.parameterNotInDomain("q", q)` if `q ∉ (0, 1)`.
    ///   - `SpecialFunctionError.parameterNotFinite` if `x` or `q` is not finite.
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
    /// - Throws: See the generic overload for error behavior.
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
    /// - Throws: See the generic overload for error behavior.
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
    /// Uses Boost’s real-τ API (internally `q = exp(-π · τ)`). Finiteness is
    /// enforced here; further domain handling is delegated to Boost.
    ///
    /// - Parameters:
    ///   - x: Real argument. Must be finite.
    ///   - tau: Real half-period ratio τ. Must be finite.
    /// - Returns: θ₃(x, τ) as `T`.
    /// - Throws: `SpecialFunctionError.parameterNotFinite` when `x` or `tau` is not finite.
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
    /// - Throws: See the generic overload for error behavior.
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
    /// - Throws: See the generic overload for error behavior.
    @inlinable static func jacobiTheta3Tau(_ x: Float, tau: Float) throws -> Float80 {
        guard tau.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "tau", value: tau) }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        return bs_jacobi_theta3tau_l(Float80(x), Float80(tau))
    }
    #endif

    
    /// Jacobi theta function θ₄(x, q) parameterized by the nome `q`.
    ///
    /// θ₄ is an even function of `x`. This overload requires finite real `x`
    /// and real nome `q` with 0 < q < 1. It uses the double-precision backend
    /// and casts back to `T`.
    ///
    /// - Parameters:
    ///   - x: Real argument. Must be finite.
    ///   - q: Nome `q` with `0 < q < 1`. Must be finite.
    /// - Returns: θ₄(x, q) as `T`.
    /// - Throws:
    ///   - `SpecialFunctionError.parameterNotInDomain("q", q)` if `q ∉ (0, 1)`.
    ///   - `SpecialFunctionError.parameterNotFinite` if `x` or `q` is not finite.
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
    /// - Throws: See the generic overload for error behavior.
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
    /// - Throws: See the generic overload for error behavior.
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
    /// Uses Boost’s real-τ API (internally `q = exp(-π · τ)`). Finiteness is
    /// enforced here; further domain handling is delegated to Boost.
    ///
    /// - Parameters:
    ///   - x: Real argument. Must be finite.
    ///   - tau: Real half-period ratio τ. Must be finite.
    /// - Returns: θ₄(x, τ) as `T`.
    /// - Throws: `SpecialFunctionError.parameterNotFinite` when `x` or `tau` is not finite.
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
    /// - Throws: See the generic overload for error behavior.
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
    /// - Throws: See the generic overload for error behavior.
    @inlinable static func jacobiTheta4Tau(_ x: Float, tau: Float) throws -> Float80 {
        guard tau.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "tau", value: tau) }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        return bs_jacobi_theta4tau_l(Float80(x), Float80(tau))
    }
    #endif
}
