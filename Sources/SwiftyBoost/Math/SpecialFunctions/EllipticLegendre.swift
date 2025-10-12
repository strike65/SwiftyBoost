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
    
    /// Compute Legendre’s complete elliptic integral of the third kind Π(n | k).
    ///
    /// Definition (modulus k, characteristic n):
    /// - Π(n | k) = Π(n; π/2 | k) = ∫₀^{π/2} dθ / ((1 − n sin²θ) √(1 − k² sin²θ))
    ///
    /// Domain:
    /// - Requires |k| ≤ 1 for this wrapper. Additional singularities may occur for
    ///   some (n, k) when 1 − n sin²θ crosses zero; these are handled by Boost.Math.
    ///   We validate finiteness of inputs and the |k| bound.
    ///
    /// Parameters:
    /// - k: The modulus (|k| ≤ 1).
    /// - n: The characteristic (finite real).
    ///
    /// Returns:
    /// - Π(n | k) as `T`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotFinite(name: "k")` if `k` is NaN or ±∞.
    /// - `SpecialFunctionError.parameterNotFinite(name: "n")` if `n` is NaN or ±∞.
    /// - `SpecialFunctionError.parameterOutOfRange(name: "k", min: -1, max: 1)` if `|k| > 1`.
    @inlinable static func completeEllipticIntegralPi<T: BinaryFloatingPoint>(_ k: T, characteristic n: T) throws -> T {
        let dk = D(k), dn = D(n)
        guard dk.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "k") }
        guard dn.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "n") }
        guard abs(dk) <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "k", min: -1, max: 1) }
        return T(bs_ellint_3_complete(dk, dn))
    }
    
    /// Compute Legendre’s incomplete elliptic integral of the third kind Π(n; φ | k).
    ///
    /// Definition (modulus k, characteristic n, amplitude φ):
    /// - Π(n; φ | k) = ∫₀^{φ} dθ / ((1 − n sin²θ) √(1 − k² sin²θ))
    ///
    /// Domain:
    /// - Requires |k| ≤ 1 for this wrapper. Additional singularities may occur for
    ///   some (n, k, φ) when 1 − n sin²θ crosses zero; these are handled by Boost.Math.
    ///   We validate finiteness of inputs and the |k| bound.
    ///
    /// Parameters:
    /// - k: The modulus (|k| ≤ 1).
    /// - n: The characteristic (finite real).
    /// - phi: The amplitude φ (finite real).
    ///
    /// Returns:
    /// - Π(n; φ | k) as `T`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotFinite(name: "k")` if `k` is NaN or ±∞.
    /// - `SpecialFunctionError.parameterNotFinite(name: "n")` if `n` is NaN or ±∞.
    /// - `SpecialFunctionError.parameterNotFinite(name: "phi")` if `phi` is NaN or ±∞.
    /// - `SpecialFunctionError.parameterOutOfRange(name: "k", min: -1, max: 1)` if `|k| > 1`.
    @inlinable static func incompleteEllipticIntegralPi<T: BinaryFloatingPoint>(_ k: T, characteristic n: T, phi: T) throws -> T {
        let dk = D(k), dn = D(n), dphi = D(phi)
        guard dk.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "k") }
        guard dn.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "n") }
        guard dphi.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "phi") }
        guard abs(dk) <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "k", min: -1, max: 1) }
        // Boost/C bridge signature is ellint_3(k, n, phi): pass (k, n, phi) in that order.
        return T(bs_ellint_3(dk, dn, dphi))
    }
    
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
    
    /// Π(n | k) for `Float`. Requires |k| ≤ 1 and finite n.
    @inlinable static func completeEllipticIntegralPi(_ k: Float, characteristic n: Float) throws -> Float {
        guard k.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "k") }
        guard n.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "n") }
        guard abs(k) <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "k", min: -1, max: 1) }
        return bs_ellint_3_complete_f(k, n)
    }
    
    /// Π(n; φ | k) for `Float`. Requires |k| ≤ 1 and finite n, φ.
    @inlinable static func incompleteEllipticIntegralPi(_ k: Float, characteristic n: Float, phi: Float) throws -> Float {
        guard k.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "k") }
        guard n.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "n") }
        guard phi.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "phi") }
        guard abs(k) <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "k", min: -1, max: 1) }
        // Pass (k, n, phi)
        return bs_ellint_3_f(k, n, phi)
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
    
    /// Π(n | k) for `Float80` (x86_64 only). Requires |k| ≤ 1 and finite n.
    @inlinable static func completeEllipticIntegralPi(_ k: Float80, characteristic n: Float80) throws -> Float80 {
        guard k.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "k") }
        guard n.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "n") }
        guard abs(k) <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "k", min: -1, max: 1) }
        return bs_ellint_3_complete_l(k, n)
    }
    
    /// Π(n; φ | k) for `Float80` (x86_64 only). Requires |k| ≤ 1 and finite n, φ.
    @inlinable static func incompleteEllipticIntegralPi(_ k: Float80, characteristic n: Float80, phi: Float80) throws -> Float80 {
        guard k.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "k") }
        guard n.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "n") }
        guard phi.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "phi") }
        guard abs(k) <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "k", min: -1, max: 1) }
        // Pass (k, n, phi)
        return bs_ellint_3_l(k, n, phi)
    }
#endif
    
}
