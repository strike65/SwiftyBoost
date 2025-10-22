//
//  Error.swift
//  Math/SpecialFunctions
//
//  Swift wrappers for the Gaussian error function erf(x) and the complementary
//  error function erfc(x). These APIs validate common real-domain constraints
//  (finiteness) and delegate numerical work to CBoostBridge (wrapping Boost.Math).
//
//  Functions (real-valued):
//  - errorFunction(x) = erf(x)  = 2 / √π ∫₀ˣ e^(−t²) dt
//  - complementaryErrorFunction(x) = erfc(x) = 1 − erf(x) = 2 / √π ∫ₓ^∞ e^(−t²) dt
//
//  Notes:
//  - For all real x, erf(x) ∈ (−1, 1) and erfc(x) ∈ (0, 2).
//  - These functions are commonly used in probability, statistics, and diffusion
//    problems (e.g., normal CDF and tail probabilities).
//
//  Overloads:
//  - Generic BinaryFloatingPoint overloads funnel through a Double-backed
//    implementation using the helper D(_:), returning the result as T.
//  - Type-specific overloads for Float and (on x86_64) Float80 call directly
//    into the corresponding C functions to avoid conversions.
//
//  Error model:
//  - Non-finite arguments throw SpecialFunctionError.parameterNotFinite.
//
//  References:
//  - NIST DLMF §7 (Error Functions): https://dlmf.nist.gov/7
//  - Boost.Math erf/erfc: https://www.boost.org/doc/libs/release/libs/math/doc/html/math_toolkit/sf_erf/error_function.html
//

import SwiftyBoostPrelude
public extension SpecialFunctions {
    
    
    // MARK: - Generic overloads (Double-backed)
    
    /// Compute the Gaussian error function erf(x).
    ///
    /// Definition:
    /// - erf(x) = 2 / √π ∫₀ˣ e^(−t²) dt
    ///
    /// Range:
    /// - erf(x) ∈ (−1, 1) for all real x.
    ///
    /// Parameters:
    /// - x: The input value (finite real).
    ///
    /// Returns:
    /// - erf(x) as `T`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotFinite(name: "x")` if `x` is NaN or ±∞.
    @inlinable static func errorFunction<T: Real & BinaryFloatingPoint>(_ x: T) throws -> T {
        let dx = D(x)
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return T(bs_erf_d(dx))
    }
    
    /// Compute the complementary error function erfc(x) = 1 − erf(x).
    ///
    /// Range:
    /// - erfc(x) ∈ (0, 2) for all real x.
    ///
    /// Par/Users/volker/Documents/Develop_current/personal/SwiftyBoost/SwiftyBoost/Sources/CBoostBridge/include/CBoostBridge.hameters:
    /// - x: The input value (finite real).
    ///
    /// Returns:
    /// - erfc(x) as `T`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotFinite(name: "x")` if `x` is NaN or ±∞.
    @inlinable static func complementaryErrorFunction<T: Real & BinaryFloatingPoint>(_ x: T) throws -> T {
        let dx = D(x)
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return T(bs_erfc_d(dx))
    }
    
    /// Compute the inverse Gaussian error function erf⁻¹(z).
    ///
    /// Definition:
    /// - Given z = erf(x), this returns x = erf⁻¹(z).
    ///
    /// Domain (real-valued result):
    /// - Requires |z| < 1. At z = ±1 the inverse diverges (±∞), and for |z| > 1 there is no real solution.
    ///
    /// Parameters:
    /// - z: The target value of erf(x) (finite real).
    ///
    /// Returns:
    /// - x such that erf(x) = z, as `T`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotFinite(name: "z")` if `z` is NaN or ±∞.
    /// - `SpecialFunctionError.parameterOutOfRange(name: "z", min: -1.nextUp, max: 1.nextDown)` if `|z| ≥ 1`.
    ///
    /// Notes:
    /// - This generic overload funnels through a Double-backed backend and converts the result back to `T`.
    @inlinable static func inverseErrorFunction<T: Real & BinaryFloatingPoint>(_ z: T) throws -> T {
        let dp = D(z)
        guard dp.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "z") }
        guard abs(dp) < 1 else {
            throw SpecialFunctionError.parameterOutOfRange(
                name: "z",
                min: (-1.0).nextUp,
                max: (1.0).nextDown
            )
        }
        return T(bs_erf_inv_d(dp))
    }
    
    
    // MARK: - Float overloads
    
    /// erf(x) for `Float`.
    ///
    /// Throws if `x` is not finite.
    @inlinable static func errorFunction(_ x: Float) throws -> Float {
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return bs_erf_f(x)
    }
    
    /// erfc(x) for `Float`.
    ///
    /// Throws if `x` is not finite.
    @inlinable static func complementaryErrorFunction(_ x: Float) throws -> Float {
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return bs_erfc_f(x)
    }
    
    /// Inverse error function erf⁻¹(z) for `Float`.
    ///
    /// Domain:
    /// - Requires |z| < 1 for a finite real result.
    ///
    /// Parameters:
    /// - z: Target value of erf(x) (finite real).
    ///
    /// Returns:
    /// - x such that erf(x) = z, as `Float`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotFinite(name: "z")` if `z` is NaN or ±∞.
    /// - `SpecialFunctionError.parameterOutOfRange(name: "z", min: -1.nextUp, max: 1.nextDown)` if `|z| ≥ 1`.
    @inlinable static func inverseErrorFunction(_ z: Float) throws -> Float {
        guard z.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "z") }
        guard abs(z) < 1 else {
            throw SpecialFunctionError.parameterOutOfRange(
                name: "z",
                min: Double(Float(-1).nextUp),
                max: Double(Float(1).nextDown)
            )
        }
        return bs_erf_inv_f(z)
    }
    // MARK: - Float80 overloads (x86_64 only)
    
#if arch(x86_64)
    /// erf(x) for `Float80` (x86_64 only).
    ///
    /// Throws if `x` is not finite.
    @inlinable static func errorFunction(_ x: Float80) throws -> Float80 {
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return bs_erf_l(x)
    }
    
    /// erfc(x) for `Float80` (x86_64 only).
    ///
    /// Throws if `x` is not finite.
    @inlinable static func complementaryErrorFunction(_ x: Float80) throws -> Float80 {
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return bs_erfc_l(x)
    }
    
    /// Inverse error function erf⁻¹(z) for `Float80` (x86_64 only).
    ///
    /// Domain:
    /// - Requires |z| < 1 for a finite real result.
    ///
    /// Parameters:
    /// - z: Target value of erf(x) (finite real).
    ///
    /// Returns:
    /// - x such that erf(x) = z, as `Float80`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotFinite(name: "z")` if `z` is NaN or ±∞.
    /// - `SpecialFunctionError.parameterOutOfRange(name: "z", min: -1.nextUp, max: 1.nextDown)` if `|z| ≥ 1`.
    @inlinable static func inverseErrorFunction(_ z: Float80) throws -> Float80 {
        guard z.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "z") }
        guard abs(z) < 1 else {
            throw SpecialFunctionError.parameterOutOfRange(
                name: "z",
                min: Double((Float80(-1)).nextUp),
                max: Double((Float80(1)).nextDown)
            )
        }
        return bs_erf_inv_l(z)
    }

#endif
}
