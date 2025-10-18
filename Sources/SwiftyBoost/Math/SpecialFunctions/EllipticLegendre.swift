//
//  EllipticLegendre.swift
//  Math/SpecialFunctions
//
//  Swift wrappers for the classical Legendre elliptic integrals of the first,
//  second, and third kind, in complete and incomplete forms. The APIs validate
//  common real-domain constraints and delegate numerical work to CBoostBridge
//  (wrapping Boost.Math).
//
//  Functions (real-valued):
//  - completeEllipticIntegralK(k):    K(k)      = F(π/2 | k)    — Legendre’s complete integral of the 1st kind.
//  - incompleteEllipticIntegralF(k, phi): F(φ | k)             — Legendre’s incomplete integral of the 1st kind.
//  - completeEllipticIntegralE(k):    E(k)      = E(π/2 | k)    — Legendre’s complete integral of the 2nd kind.
//  - incompleteEllipticIntegralE(k, phi): E(φ | k)             — Legendre’s incomplete integral of the 2nd kind.
//  - completeEllipticIntegralPi(k, characteristic: n): Π(n | k) — Legendre’s complete integral of the 3rd kind.
//  - incompleteEllipticIntegralPi(k, characteristic: n, phi): Π(n; φ | k) — Legendre’s incomplete integral of the 3rd kind.
//
//  Parameterization:
//  - These wrappers use the “characteristic” form for the 3rd-kind integral Π(n; φ | k)
//    and Π(n | k), consistent with Boost.Math’s API. Here k is the modulus (not m = k^2).
//
//  Domain notes (real-valued):
//  - For the real-valued functions exposed here, we require |k| ≤ 1 (k ∈ [−1, 1]).
//  - φ (phi) is a real amplitude with no additional restriction in these wrappers.
//  - The complete/third-kind functions have additional singular curves for some (n, k);
//    this Swift layer enforces common constraints (|k| ≤ 1) and finiteness, while Boost.Math
//    performs the full evaluation and handles edge behavior.
//  - All generic overloads validate finiteness of inputs where appropriate.
//
//  Overloads:
//  - Generic BinaryFloatingPoint overloads convert to Double using the helper D(_:)
//    for evaluation and return results as the same generic type T.
//  - Type-specific overloads for Float and (on x86_64) Float80 call directly into
//    the matching C functions for performance.
//
//  Error model:
//  - Non-finite arguments throw SpecialFunctionError.parameterNotFinite.
//  - |k| > 1 throws SpecialFunctionError.parameterOutOfRange(name: "k", min: -1, max: 1).
//
//  References:
//  - NIST DLMF §19 (Elliptic Integrals): https://dlmf.nist.gov/19
//  - Boost.Math elliptic integrals:
//    https://www.boost.org/doc/libs/release/libs/math/doc/html/math_toolkit/ellint.html
//

import CBoostBridge
public extension SpecialFunctions {
    
    
    /// Compute Legendre’s complete elliptic integral of the first kind K(k).
    ///
    /// Definition (modulus k):
    /// - K(k) = F(π/2 | k) = ∫₀^{π/2} dθ / √(1 − k² sin²θ)
    ///
    /// Domain:
    /// - Requires |k| ≤ 1 for real-valued results in this wrapper.
    ///
    /// Parameters:
    /// - k: The modulus (not the parameter m = k²).
    ///
    /// Returns:
    /// - K(k) as `T`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotFinite(name: "k")` if `k` is NaN or ±∞.
    /// - `SpecialFunctionError.parameterOutOfRange(name: "k", min: -1, max: 1)` if `|k| > 1`.
    @inlinable static func completeEllipticIntegralK<T: BinaryFloatingPoint>(_ k: T) throws -> T {
        let dk = D(k)
        guard dk.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "k") }
        guard abs(dk) <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "k", min: -1, max: 1) }
        return T(bs_ellint_1_complete(dk))
    }
    
    /// Compute Legendre’s incomplete elliptic integral of the first kind F(φ | k).
    ///
    /// Definition (modulus k, amplitude φ):
    /// - F(φ | k) = ∫₀^{φ} dθ / √(1 − k² sin²θ)
    ///
    /// Domain:
    /// - Requires |k| ≤ 1. φ is any finite real (no restriction here).
    ///
    /// Parameters:
    /// - k: The modulus (|k| ≤ 1).
    /// - phi: The amplitude φ (finite real).
    ///
    /// Returns:
    /// - F(φ | k) as `T`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotFinite(name: "k")` if `k` is NaN or ±∞.
    /// - `SpecialFunctionError.parameterNotFinite(name: "phi")` if `phi` is NaN or ±∞.
    /// - `SpecialFunctionError.parameterOutOfRange(name: "k", min: -1, max: 1)` if `|k| > 1`.
    @inlinable static func incompleteEllipticIntegralF<T: BinaryFloatingPoint>(_ k: T, phi: T) throws -> T {
        let dk = D(k), dphi = D(phi)
        guard dk.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "k") }
        guard dphi.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "phi") }
        guard abs(dk) <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "k", min: -1, max: 1) }
        return T(bs_ellint_1(dk, dphi))
    }

    // Mixed-precision promotions (Float ↔ Double) → Double
    /// F(φ | k) with mixed `Float`/`Double`; returns `Double`.
    @inlinable static func incompleteEllipticIntegralF(_ k: Float, phi: Double) throws -> Double { try incompleteEllipticIntegralF(Double(k), phi: phi) }
    /// F(φ | k) with mixed `Double`/`Float`; returns `Double`.
    @inlinable static func incompleteEllipticIntegralF(_ k: Double, phi: Float) throws -> Double { try incompleteEllipticIntegralF(k, phi: Double(phi)) }
    
    /// Compute Legendre’s complete elliptic integral of the second kind E(k).
    ///
    /// Definition (modulus k):
    /// - E(k) = E(π/2 | k) = ∫₀^{π/2} √(1 − k² sin²θ) dθ
    ///
    /// Domain:
    /// - Requires |k| ≤ 1.
    ///
    /// Parameters:
    /// - k: The modulus (|k| ≤ 1).
    ///
    /// Returns:
    /// - E(k) as `T`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotFinite(name: "k")` if `k` is NaN or ±∞.
    /// - `SpecialFunctionError.parameterOutOfRange(name: "k", min: -1, max: 1)` if `|k| > 1`.
    @inlinable static func completeEllipticIntegralE<T: BinaryFloatingPoint>(_ k: T) throws -> T {
        let dk = D(k)
        guard dk.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "k") }
        guard abs(dk) <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "k", min: -1, max: 1) }
        return T(bs_ellint_2_complete(dk))
    }
    
    /// Compute Legendre’s incomplete elliptic integral of the second kind E(φ | k).
    ///
    /// Definition (modulus k, amplitude φ):
    /// - E(φ | k) = ∫₀^{φ} √(1 − k² sin²θ) dθ
    ///
    /// Domain:
    /// - Requires |k| ≤ 1. φ is any finite real (no restriction here).
    ///
    /// Parameters:
    /// - k: The modulus (|k| ≤ 1).
    /// - phi: The amplitude φ (finite real).
    ///
    /// Returns:
    /// - E(φ | k) as `T`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotFinite(name: "k")` if `k` is NaN or ±∞.
    /// - `SpecialFunctionError.parameterNotFinite(name: "phi")` if `phi` is NaN or ±∞.
    /// - `SpecialFunctionError.parameterOutOfRange(name: "k", min: -1, max: 1)` if `|k| > 1`.
    @inlinable static func incompleteEllipticIntegralE<T: BinaryFloatingPoint>(_ k: T, phi: T) throws -> T {
        let dk = D(k), dphi = D(phi)
        guard dk.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "k") }
        guard dphi.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "phi") }
        guard abs(dk) <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "k", min: -1, max: 1) }
        return T(bs_ellint_2(dk, dphi))
    }

    // Mixed-precision promotions (Float ↔ Double) → Double
    /// E(φ | k) with mixed `Float`/`Double`; returns `Double`.
    @inlinable static func incompleteEllipticIntegralE(_ k: Float, phi: Double) throws -> Double { try incompleteEllipticIntegralE(Double(k), phi: phi) }
    /// E(φ | k) with mixed `Double`/`Float`; returns `Double`.
    @inlinable static func incompleteEllipticIntegralE(_ k: Double, phi: Float) throws -> Double { try incompleteEllipticIntegralE(k, phi: Double(phi)) }
    
    /// Compute Legendre’s incomplete elliptic integral of the third kind Π(n; φ | k) in its characteristic form.
    ///
    /// Definitions and notation:
    /// - This API follows the “characteristic” parameterization used by Boost.Math:
    ///   Π(n; φ | k) = ∫₀^{φ} dθ / [(1 − n sin²θ) √(1 − k² sin²θ)]
    ///   where k is the modulus (not the parameter m = k²), φ is the amplitude, and n is the characteristic (often denoted ν).
    ///
    /// Key identities:
    /// - Π(0; φ | k) = F(φ | k) (reduces to the elliptic integral of the first kind).
    /// - When φ = 0, Π(n; 0 | k) = 0 for any k, n in the regular real domain.
    ///
    /// Domain and singular behavior (real-valued use):
    /// - For many applications one restricts |k| ≤ 1 and uses real φ, n.
    /// - The integrand has potential singularities where 1 − n sin²θ = 0 or 1 − k² sin²θ = 0 within [0, φ].
    ///   Depending on (n, k, φ), Π can diverge or require interpretation via principal value/analytic continuation.
    ///
    /// Validation and error handling:
    /// - This particular wrapper does not perform argument validation and does not throw.
    ///   It forwards directly to the Boost backend via CBoostBridge (bs_ellint_3), which implements the numerical evaluation
    ///   and determines how to handle edge cases, singularities, and analytic continuation.
    /// - If you need stricter real-domain checking (e.g. enforcing |k| ≤ 1 or finiteness), perform it at the call site.
    ///
    /// Parameters:
    /// - k: Modulus (not m = k²). In many real-valued contexts, |k| ≤ 1 is typical, but not enforced here.
    /// - nu: Characteristic n (also written ν), appearing as 1 − n sin²θ in the denominator.
    /// - phi: Amplitude φ (in radians).
    ///
    /// Returns:
    /// - The value Π(n; φ | k) as `T`, computed by the backend and converted back to `T`.
    ///
    /// References:
    /// - NIST DLMF §19.2, §19.7 (Legendre normal forms and properties)
    /// - Boost.Math ellint_3 (characteristic form)
    @inlinable static func incompleteEllipticIntegralPi<T: BinaryFloatingPoint>(_ k: T, _ nu: T, _ phi: T) -> T {
        let dk = D(k), dphi = D(phi), dnu = D(nu)
        return T(bs_ellint_3(dk, dnu, dphi))
    }

    // Mixed-precision promotions for Π(n; φ | k) (Float ↔ Double) → Double
    /// Π with mixed argument types; returns `Double`.
    @inlinable static func incompleteEllipticIntegralPi(_ k: Float, _ nu: Double, _ phi: Double) -> Double { incompleteEllipticIntegralPi(Double(k), nu, phi) }
    @inlinable static func incompleteEllipticIntegralPi(_ k: Double, _ nu: Float, _ phi: Double) -> Double { incompleteEllipticIntegralPi(k, Double(nu), phi) }
    @inlinable static func incompleteEllipticIntegralPi(_ k: Double, _ nu: Double, _ phi: Float) -> Double { incompleteEllipticIntegralPi(k, nu, Double(phi)) }
    
    /// Compute Legendre’s complete elliptic integral of the third kind Π(n | k) in its characteristic form.
    ///
    /// Definition (characteristic n, modulus k):
    /// - Π(n | k) = Π(n; π/2 | k) = ∫₀^{π/2} dθ / [(1 − n sin²θ) √(1 − k² sin²θ)]
    ///
    /// Identities and special cases:
    /// - Relation to the incomplete form: Π(n | k) = Π(n; π/2 | k).
    /// - Reduction to K(k): Π(0 | k) = K(k), i.e. setting n = 0 removes the (1 − n sin²θ) factor.
    ///
    /// Domain and singular behavior (real-valued use):
    /// - In many real-valued applications, |k| ≤ 1 is assumed (k ∈ [−1, 1]); this wrapper does not enforce it here.
    /// - Singularities can occur for parameter curves where 1 − n sin²θ = 0 within θ ∈ [0, π/2], e.g. n ≥ 1 intersects the integration domain.
    ///   Along such curves Π(n | k) may diverge or require principal-value interpretation; the backend governs the behavior.
    ///
    /// Numerical considerations:
    /// - Different evaluation paths (complete vs. incomplete at π/2) can yield tiny roundoff differences at Double precision.
    ///   When comparing Π(n | k) to Π(n; π/2 | k), prefer an absolute tolerance rather than exact equality.
    ///
    /// Validation and error handling:
    /// - This wrapper performs no Swift-side validation and does not throw.
    /// - The CBoostBridge/Boost.Math backend (bs_ellint_3_complete) performs the numerical evaluation and handles edge cases.
    /// - If your use requires explicit domain checks (e.g., |k| ≤ 1, finite inputs), validate at the call site before invoking this API.
    ///
    /// Parameters:
    /// - k: Modulus (not m = k²).
    /// - nu: Characteristic n (ν), appearing as 1 − n sin²θ in the denominator.
    ///
    /// Returns:
    /// - Π(n | k) as `T`, computed by the backend and converted back to `T`.
    ///
    /// References:
    /// - NIST DLMF §19.2, §19.7 (Legendre normal forms and properties)
    /// - Boost.Math ellint_3_complete (characteristic form)
    @inlinable static func completeEllipticIntegralPi<T: BinaryFloatingPoint>(_ k: T, _ nu: T) -> T {
        let dk = D(k), dnu = D(nu)
        return T(bs_ellint_3_complete(dk, dnu))

    }

    // Mixed-precision promotions for Π(n | k) (Float ↔ Double) → Double
    @inlinable static func completeEllipticIntegralPi(_ k: Float, _ nu: Double) -> Double { completeEllipticIntegralPi(Double(k), nu) }
    @inlinable static func completeEllipticIntegralPi(_ k: Double, _ nu: Float) -> Double { completeEllipticIntegralPi(k, Double(nu)) }
    
    // MARK: - Float overloads
    // Direct Float-precision entry points that avoid generic conversions and call
    // the corresponding C implementations.
    
    /// K(k) for `Float`. Requires |k| ≤ 1.
    @inlinable static func completeEllipticIntegralK(_ k: Float) throws -> Float {
        guard k.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "k") }
        guard abs(k) <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "k", min: -1, max: 1) }
        return bs_ellint_1_complete_f(k)
    }
    
    /// F(φ | k) for `Float`. Requires |k| ≤ 1 and finite φ.
    @inlinable static func incompleteEllipticIntegralF(_ k: Float, phi: Float) throws -> Float {
        guard k.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "k") }
        guard phi.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "phi") }
        guard abs(k) <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "k", min: -1, max: 1) }
        return bs_ellint_1_f(k, phi)
    }
    
    /// E(k) for `Float`. Requires |k| ≤ 1.
    @inlinable static func completeEllipticIntegralE(_ k: Float) throws -> Float {
        guard k.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "k") }
        guard abs(k) <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "k", min: -1, max: 1) }
        return bs_ellint_2_complete_f(k)
    }
    
    /// E(φ | k) for `Float`. Requires |k| ≤ 1 and finite φ.
    @inlinable static func incompleteEllipticIntegralE(_ k: Float, phi: Float) throws -> Float {
        guard k.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "k") }
        guard phi.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "phi") }
        guard abs(k) <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "k", min: -1, max: 1) }
        return bs_ellint_2_f(k, phi)
    }

    /// Π(n; φ | k) for `Float` in characteristic form.
    ///
    /// See the generic `incompleteEllipticIntegralPi` for full discussion. This overload:
    /// - Performs no argument validation and does not throw.
    /// - Delegates directly to the CBoostBridge backend (bs_ellint_3_f).
    ///
    /// Parameters:
    /// - k: Modulus (not m = k²). Typically |k| ≤ 1 for real-valued contexts.
    /// - nu: Characteristic n (ν).
    /// - phi: Amplitude φ (radians).
    ///
    /// Returns:
    /// - Π(n; φ | k) as `Float`.
    @inlinable static func incompleteEllipticIntegralPi(_ k: Float, _ nu: Float, _ phi: Float) -> Float {
        return bs_ellint_3_f(k, nu, phi)
    }

    /// Π(n | k) for `Float` in characteristic form (complete third-kind integral).
    ///
    /// Definition and notes:
    /// - Π(n | k) = Π(n; π/2 | k).
    /// - Reduction: Π(0 | k) = K(k).
    /// - No argument validation is performed here; see the generic overload docs for domain notes.
    ///
    /// Returns:
    /// - Π(n | k) as `Float`.
    @inlinable static func completeEllipticIntegralPi(_ k: Float, _ nu: Float) -> Float {
        return bs_ellint_3_complete_f(k, nu)
    }

   
    // MARK: - Float80 overloads (x86_64)
    // Extended-precision versions for platforms that support Float80.
    
#if arch(x86_64)
    /// K(k) for `Float80` (x86_64 only). Requires |k| ≤ 1.
    @inlinable static func completeEllipticIntegralK(_ k: Float80) throws -> Float80 {
        guard k.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "k") }
        guard abs(k) <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "k", min: -1, max: 1) }
        return bs_ellint_1_complete_l(k)
    }
    
    /// F(φ | k) for `Float80` (x86_64 only). Requires |k| ≤ 1 and finite φ.
    @inlinable static func incompleteEllipticIntegralF(_ k: Float80, phi: Float80) throws -> Float80 {
        guard k.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "k") }
        guard phi.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "phi") }
        guard abs(k) <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "k", min: -1, max: 1) }
        return bs_ellint_1_l(k, phi)
    }

    // Mixed promotions with Float80 → Float80
    @inlinable static func incompleteEllipticIntegralF(_ k: Float80, phi: Double) throws -> Float80 { try incompleteEllipticIntegralF(k, phi: Float80(phi)) }
    @inlinable static func incompleteEllipticIntegralF(_ k: Double, phi: Float80) throws -> Float80 { try incompleteEllipticIntegralF(Float80(k), phi: phi) }
    @inlinable static func incompleteEllipticIntegralF(_ k: Float80, phi: Float) throws -> Float80 { try incompleteEllipticIntegralF(k, phi: Float80(phi)) }
    @inlinable static func incompleteEllipticIntegralF(_ k: Float, phi: Float80) throws -> Float80 { try incompleteEllipticIntegralF(Float80(k), phi: phi) }
    
    /// E(k) for `Float80` (x86_64 only). Requires |k| ≤ 1.
    @inlinable static func completeEllipticIntegralE(_ k: Float80) throws -> Float80 {
        guard k.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "k") }
        guard abs(k) <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "k", min: -1, max: 1) }
        return bs_ellint_2_complete_l(k)
    }
    
    /// E(φ | k) for `Float80` (x86_64 only). Requires |k| ≤ 1 and finite φ.
    @inlinable static func incompleteEllipticIntegralE(_ k: Float80, phi: Float80) throws -> Float80 {
        guard k.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "k") }
        guard phi.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "phi") }
        guard abs(k) <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "k", min: -1, max: 1) }
        return bs_ellint_2_l(k, phi)
    }

    // Mixed promotions with Float80 → Float80
    @inlinable static func incompleteEllipticIntegralE(_ k: Float80, phi: Double) throws -> Float80 { try incompleteEllipticIntegralE(k, phi: Float80(phi)) }
    @inlinable static func incompleteEllipticIntegralE(_ k: Double, phi: Float80) throws -> Float80 { try incompleteEllipticIntegralE(Float80(k), phi: phi) }
    @inlinable static func incompleteEllipticIntegralE(_ k: Float80, phi: Float) throws -> Float80 { try incompleteEllipticIntegralE(k, phi: Float80(phi)) }
    @inlinable static func incompleteEllipticIntegralE(_ k: Float, phi: Float80) throws -> Float80 { try incompleteEllipticIntegralE(Float80(k), phi: phi) }

    /// Π(n; φ | k) for `Float80` (x86_64 only) in characteristic form.
    ///
    /// Behavior mirrors the `Float` overload:
    /// - No argument validation or throwing at this layer.
    /// - Delegates to the backend (bs_ellint_3_f_l) for evaluation and handling of edge cases.
    ///
    /// Parameters:
    /// - k: Modulus (not m = k²).
    /// - nu: Characteristic n (ν).
    /// - phi: Amplitude φ (radians).
    ///
    /// Returns:
    /// - Π(n; φ | k) as `Float80`.
    @inlinable static func incompleteEllipticIntegralPi(_ k: Float80, _ nu: Float80, _ phi: Float80) -> Float80 {
        return bs_ellint_3_f_l(k, nu, phi)
    }

    // Mixed promotions with Float80 for Π(n; φ | k) → Float80
    @inlinable static func incompleteEllipticIntegralPi(_ k: Float80, _ nu: Double, _ phi: Double) -> Float80 { incompleteEllipticIntegralPi(k, Float80(nu), Float80(phi)) }
    @inlinable static func incompleteEllipticIntegralPi(_ k: Double, _ nu: Float80, _ phi: Double) -> Float80 { incompleteEllipticIntegralPi(Float80(k), nu, Float80(phi)) }
    @inlinable static func incompleteEllipticIntegralPi(_ k: Double, _ nu: Double, _ phi: Float80) -> Float80 { incompleteEllipticIntegralPi(Float80(k), Float80(nu), phi) }
    @inlinable static func incompleteEllipticIntegralPi(_ k: Float80, _ nu: Float, _ phi: Float) -> Float80 { incompleteEllipticIntegralPi(k, Float80(nu), Float80(phi)) }
    @inlinable static func incompleteEllipticIntegralPi(_ k: Float, _ nu: Float80, _ phi: Float) -> Float80 { incompleteEllipticIntegralPi(Float80(k), nu, Float80(phi)) }
    @inlinable static func incompleteEllipticIntegralPi(_ k: Float, _ nu: Float, _ phi: Float80) -> Float80 { incompleteEllipticIntegralPi(Float80(k), Float80(nu), phi) }
    
    /// Π(n | k) for `Float80` (x86_64 only) in characteristic form (complete third-kind integral).
    ///
    /// Definition and notes:
    /// - Π(n | k) = Π(n; π/2 | k).
    /// - Reduction: Π(0 | k) = K(k).
    /// - This overload performs no argument validation and does not throw.
    ///
    /// Returns:
    /// - Π(n | k) as `Float80`.
    @inlinable static func completeEllipticIntegralPi(_ k: Float80, _ nu: Float80) -> Float80 {
        return bs_ellint_3_complete_l(k, nu)
    }

    // Mixed promotions with Float80 for Π(n | k) → Float80
    @inlinable static func completeEllipticIntegralPi(_ k: Float80, _ nu: Double) -> Float80 { completeEllipticIntegralPi(k, Float80(nu)) }
    @inlinable static func completeEllipticIntegralPi(_ k: Double, _ nu: Float80) -> Float80 { completeEllipticIntegralPi(Float80(k), nu) }
    @inlinable static func completeEllipticIntegralPi(_ k: Float80, _ nu: Float) -> Float80 { completeEllipticIntegralPi(k, Float80(nu)) }
    @inlinable static func completeEllipticIntegralPi(_ k: Float, _ nu: Float80) -> Float80 { completeEllipticIntegralPi(Float80(k), nu) }

#endif
    
}
