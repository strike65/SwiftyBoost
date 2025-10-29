//
//  Created by Volker Thieme 2025.
//  Copyright © 2025 Volker Thieme.
//  License: MIT (see project root)
//
//  Chebyshev.swift
//  Math/SpecialFunctions
//
//  Swift wrappers for Chebyshev polynomials of the first and second kind,
//  Tₙ(x) and Uₙ(x). These validate common real-domain constraints and
//  delegate numerical work to CBoostBridge (wrapping Boost.Math).
//
//  Overview
//  - The Chebyshev polynomials of the first kind Tₙ and second kind Uₙ are
//    classic families of orthogonal polynomials on [−1, 1] with numerous uses
//    in approximation theory, spectral methods, quadrature, and signal processing.
//  - Definitions (trigonometric):
//      • Tₙ(cos θ) = cos(nθ)
//      • Uₙ(cos θ) = sin((n + 1)θ) / sin θ
//    valid for real θ. These identities imply stable evaluation for |x| ≤ 1.
//  - Three-term recurrences:
//      • T₀(x) = 1, T₁(x) = x,      Tₙ₊₁(x) = 2x Tₙ(x) − Tₙ₋₁(x)
//      • U₀(x) = 1, U₁(x) = 2x,     Uₙ₊₁(x) = 2x Uₙ(x) − Uₙ₋₁(x)
//  - Generating functions:
//      • Σₙ≥0 Tₙ(x) tⁿ = (1 − x t) / (1 − 2 x t + t²), |t| < 1
//      • Σₙ≥0 Uₙ(x) tⁿ = 1 / (1 − 2 x t + t²),         |t| < 1
//  - Values at endpoints and parity:
//      • Tₙ(1) = 1,  Tₙ(−1) = (−1)ⁿ,  Tₙ is even/odd with n
//      • Uₙ(1) = n + 1, Uₙ(−1) = (−1)ⁿ (n + 1)
//  - Orthogonality on [−1, 1]:
//      • Tₙ ⟂ Tₘ with weight w(x) = (1 − x²)^{−1/2}
//      • Uₙ ⟂ Uₘ with weight w(x) = (1 − x²)^{+1/2}
//  - Behavior for |x| > 1:
//      • Tₙ(x) relates to cosh(n arcosh x) and grows rapidly with n and |x|.
//      • Uₙ(x) relates to sinh((n+1) arcosh x) / √(x² − 1).
//    The backend (Boost.Math) selects numerically stable forms across regions.
//
//  API Summary (real-valued):
//  - chebyshevT(n, x): Tₙ(x) — first kind
//  - chebyshevU(n, x): Uₙ(x) — second kind
//  - chebyshev_next(Tn, Tn1, x): next term via the three-term recurrence
//  - chebyshevClenshawRecurrence(coefficients, x, halfWeightC0:): Chebyshev T-series evaluation via Clenshaw
//
//  Domain (real-valued use):
//  - Degree n must be an integer with n ≥ 0.
//  - Argument x must be a finite real (NaN/±∞ rejected).
//  - No additional range restriction on x; results are defined for all finite x,
//    but note potential rapid growth when |x| > 1.
//
//  Overloads:
//  - Generic BinaryFloatingPoint overloads convert to Double using D(_:) and
//    return results as the same generic type T.
//  - Type-specific overloads for Float and (on x86_64) Float80 call directly
//    into the matching C functions for performance and to avoid conversions.
//  - chebyshev_next is implemented directly in Swift using the standard recurrence.
//  - chebyshevClenshawRecurrence is implemented directly in Swift using the
//    numerically stable Clenshaw algorithm.
//
//  Error model:
//  - n < 0 throws SpecialFunctionError.parameterNotPositive(name: "n", value: n).
//  - Non-finite x throws SpecialFunctionError.parameterNotFinite(name: "x").
//  - chebyshev_next throws SpecialFunctionError.parameterNotFinite if any input is non-finite.
//  - chebyshevClenshawRecurrence throws SpecialFunctionError.parameterNotFinite if x or any coefficient is non-finite.
//
//  Performance and numerical notes:
//  - The underlying Boost.Math implementation uses region-appropriate formulas
//    (e.g., trigonometric for |x| ≤ 1, hyperbolic for |x| > 1) and stable
//    recurrences to reduce cancellation and overflow risk.
//  - The Clenshaw algorithm evaluates Chebyshev T-series in O(N) time with
//    superior numerical stability compared to naive summation.
//  - Extremely large n and/or |x| may still overflow the target floating-point
//    type; these wrappers do not pre-check magnitude, they rely on backend
//    algorithms and IEEE behavior of the selected precision.
//  - For tight loops at fixed precision, prefer the Float/Float80 overloads
//    to avoid generic conversions.
//
//  Examples:
//    // Evaluate Tₙ and Uₙ in Double
//    let t3 = try SpecialFunctions.chebyshevT(3, 0.5 as Double) // T₃(0.5) = 4·0.5³ − 3·0.5 = −0.5
//    let u2 = try SpecialFunctions.chebyshevU(2, 0.5 as Double) // U₂(0.5) = 4·0.5² − 1 = 0
//
//    // Float overload
//    let tf: Float = try SpecialFunctions.chebyshevT(4, 1.0 as Float) // 1.0
//
//    // Relation to cos/arccos for |x| ≤ 1
//    // T₇(x) = cos(7 arccos x)
//    let x = 0.3
//    let t7 = try SpecialFunctions.chebyshevT(7, x)
//    let t7Alt = cos(7 * acos(x)) // numerically close for |x| ≤ 1
//
//    // Recurrence stepping using chebyshev_next (works for Tₙ or Uₙ)
//    // Build T₂(x) from T₀, T₁ at x = 0.5:
//    // T₀=1, T₁=x, then T₂ = 2x·T₁ − T₀ = 2·0.5·0.5 − 1 = −0.5
//    let xrec = 0.5 as Double
//    let T0 = 1.0, T1 = xrec
//    let T2 = try SpecialFunctions.chebyshev_next(T1, T0, xrec) // −0.5
//
//    // Build U₂(x) from U₀, U₁ at x = 0.5:
//    // U₀=1, U₁=2x, then U₂ = 2x·U₁ − U₀ = 2·0.5·1.0 − 1 = 0
//    let U0 = 1.0, U1 = 2.0 * xrec
//    let U2 = try SpecialFunctions.chebyshev_next(U1, U0, xrec) // 0
//
//  References:
//  - NIST DLMF §18.3 (Chebyshev Polynomials): https://dlmf.nist.gov/18.3
//  - Boost.Math chebyshev: https://www.boost.org/doc/libs/release/libs/math/doc/html/math_toolkit/sf_poly/chebyshev.html
//

import SwiftyBoostPrelude

public extension SpecialFunctions {
    // MARK: - Generic overloads (Double-backed)

    /// Compute the Chebyshev polynomial of the first kind Tₙ(x).
    ///
    /// Mathematical definition:
    /// - Tₙ(cos θ) = cos(nθ)
    ///
    /// Domain:
    /// - Requires integer degree n ≥ 0 and finite real x.
    ///
    /// Parameters:
    /// - n: Polynomial degree (n ≥ 0).
    /// - x: Evaluation point (finite real).
    ///
    /// Returns:
    /// - Tₙ(x) as `T`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotPositive(name: "n", value: n)` if `n < 0`.
    /// - `SpecialFunctionError.parameterNotFinite(name: "x")` if `x` is NaN or ±∞.
    ///
    /// Notes:
    /// - For |x| ≤ 1, Tₙ(x) ∈ [−1, 1]. For |x| > 1, magnitude can grow quickly with n.
    @inlinable static func chebyshevT<T: Real & BinaryFloatingPoint & Sendable>(_ n: Int, _ x: T) throws -> T {
        guard n >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "n", value: T(n)) }
        let dx = D(x)
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        return T(bs_chebyshev_T_d(UInt32(n), dx))
    }

    /// Compute the Chebyshev polynomial of the second kind Uₙ(x).
    ///
    /// Mathematical definition:
    /// - Uₙ(cos θ) = sin((n + 1)θ) / sin θ
    ///
    /// Domain:
    /// - Requires integer degree n ≥ 0 and finite real x.
    ///
    /// Parameters:
    /// - n: Polynomial degree (n ≥ 0).
    /// - x: Evaluation point (finite real).
    ///
    /// Returns:
    /// - Uₙ(x) as `T`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotPositive(name: "n", value: n)` if `n < 0`.
    /// - `SpecialFunctionError.parameterNotFinite(name: "x")` if `x` is NaN or ±∞.
    ///
    /// Notes:
    /// - Uₙ has (n) simple zeros in (−1, 1). Uₙ(±1) = (±1)ⁿ (n + 1).
    @inlinable static func chebyshevU<T: Real & BinaryFloatingPoint & Sendable>(_ n: Int, _ x: T) throws -> T {
        guard n >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "n", value: T(n)) }
        let dx = D(x)
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        return T(bs_chebyshev_U_d(UInt32(n), dx))
    }

    // MARK: - Chebyshev T-series via Clenshaw recurrence
    
    /// Evaluate a Chebyshev T-series at x using the Clenshaw recurrence.
    ///
    /// Series convention:
    /// - By default, coefficients are interpreted with the standard half-weight on c₀:
    ///   S(x) = c₀/2 + Σ_{k=1}^{N} c_k · T_k(x).
    /// - To evaluate the series without the half-weight, set `halfWeightC0` to `false`, i.e.:
    ///   S(x) = Σ_{k=0}^{N} c_k · T_k(x).
    ///
    /// Implementation:
    /// - Calls bridged Boost functions:
    ///     • Double:     bs_chebyshev_clenshaw
    ///     • Float:      bs_chebyshev_clenshaw_f
    ///     • Float80:    bs_chebyshev_clenshaw_l (x86_64)
    ///
    /// - If `halfWeightC0 == false`, we pass a modified coefficient vector with c₀' = 2·c₀
    ///   so that Boost’s half-weight interpretation yields the desired unweighted c₀.
    @inlinable static func chebyshevClenshawRecurrence<T: Real & BinaryFloatingPoint & Sendable>(_ coefficients: [T], x: T, halfWeightC0: Bool = true) throws -> T {
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        for c in coefficients {
            guard c.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "coefficients", value: c) }
        }
        if coefficients.isEmpty { return .zero }

        // Generic path: convert to Double and call the Double backend.
        // Adjust c0 if halfWeightC0 is false.
        var cD = coefficients.map { Double($0) }
        if !halfWeightC0, !cD.isEmpty {
            cD[0] *= 2.0
        }
        let result: Double = cD.withUnsafeBufferPointer { bp in
            bs_chebyshev_clenshaw_d(bp.baseAddress, bp.count, Double(x))
        }
        return T(result)
    }

    /// Float overload for Chebyshev T-series via Clenshaw recurrence (bridged).
    @inlinable static func chebyshevClenshawRecurrence(_ coefficients: [Float], x: Float, halfWeightC0: Bool = true) throws -> Float {
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        for c in coefficients {
            guard c.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "coefficients", value: c) }
        }
        if coefficients.isEmpty { return 0 }

        if halfWeightC0 {
            return coefficients.withUnsafeBufferPointer { bp in
                bs_chebyshev_clenshaw_f(bp.baseAddress, bp.count, x)
            }
        } else {
            var coeffs = coefficients
            coeffs[0] *= 2
            return coeffs.withUnsafeBufferPointer { bp in
                bs_chebyshev_clenshaw_f(bp.baseAddress, bp.count, x)
            }
        }
    }

    #if arch(x86_64)
    /// Float80 overload for Chebyshev T-series via Clenshaw recurrence (bridged, x86_64 only).
    @inlinable static func chebyshevClenshawRecurrence(_ coefficients: [Float80], x: Float80, halfWeightC0: Bool = true) throws -> Float80 {
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        for c in coefficients {
            guard c.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "coefficients", value: c) }
        }
        if coefficients.isEmpty { return 0 }

        if halfWeightC0 {
            return coefficients.withUnsafeBufferPointer { bp in
                bs_chebyshev_clenshaw_l(bp.baseAddress, bp.count, x)
            }
        } else {
            var coeffs = coefficients
            coeffs[0] *= 2
            return coeffs.withUnsafeBufferPointer { bp in
                bs_chebyshev_clenshaw_l(bp.baseAddress, bp.count, x)
            }
        }
    }

#endif

    // MARK: - Recurrence helper

    /// Compute the next term in a Chebyshev-style three-term recurrence:
    /// Pₙ₊₁(x) = 2x · Pₙ(x) − Pₙ₋₁(x).
    ///
    /// Use with either Chebyshev Tₙ or Uₙ sequences by supplying the appropriate
    /// initial values:
    /// - For Tₙ: T₀ = 1, T₁ = x.
    /// - For Uₙ: U₀ = 1, U₁ = 2x.
    ///
    /// Parameters:
    /// - Pn:   The current term Pₙ(x) (finite).
    /// - Pn1:  The previous term Pₙ₋₁(x) (finite).
    /// - x:    The evaluation point (finite).
    ///
    /// Returns:
    /// - Pₙ₊₁(x) computed as 2x·Pₙ − Pₙ₋₁.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotFinite(name: "...")` if any input is NaN or ±∞.
    ///
    /// Example:
    /// ```swift
    /// // Build T₂(x) at x = 0.5 from T₀ and T₁:
    /// let x = 0.5 as Double
    /// var T0 = 1.0, T1 = x
    /// let T2 = try SpecialFunctions.chebyshev_next(T1, T0, x) // −0.5
    ///
    /// // Build U₂(x) at x = 0.5 from U₀ and U₁:
    /// var U0 = 1.0, U1 = 2.0 * x
    /// let U2 = try SpecialFunctions.chebyshev_next(U1, U0, x) // 0.0
    /// ```
    @inlinable static func chebyshev_next<T: Real & BinaryFloatingPoint & Sendable>(_ Pn: T, _ Pn1: T, _ x: T) throws -> T {
        guard Pn.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "Pn", value: Pn) }
        guard Pn1.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "Pn1", value: Pn1) }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        return 2 * x * Pn - Pn1
    }

    /// Float overload of chebyshev_next.
    @inlinable static func chebyshev_next(_ Pn: Float, _ Pn1: Float, _ x: Float) throws -> Float {
        guard Pn.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "Pn", value: Pn) }
        guard Pn1.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "Pn1", value: Pn1) }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        return 2 * x * Pn - Pn1
    }

    #if arch(x86_64)
    /// Float80 overload of chebyshev_next (x86_64 only).
    @inlinable static func chebyshev_next(_ Pn: Float80, _ Pn1: Float80, _ x: Float80) throws -> Float80 {
        guard Pn.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "Pn", value: Pn) }
        guard Pn1.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "Pn1", value: Pn1) }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x", value: x) }
        return 2 * x * Pn - Pn1
    }
    #endif
}
