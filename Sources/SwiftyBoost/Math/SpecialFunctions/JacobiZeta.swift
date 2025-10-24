//
//  JacobiZeta.swift
//  Math/SpecialFunctions
//
//  Created by Volker Thieme 2025.
//  License: MIT (see project root)
//
//  Swift wrappers for the Jacobi Zeta function Z(φ, k), which relates the
//  incomplete elliptic integral of the third kind to Legendre’s normal form.
//  The implementation delegates to Boost.Math and mirrors SwiftyBoost’s
//  validation conventions (finite parameters, architecture-specific overloads).
//

import SwiftyBoostPrelude

public extension SpecialFunctions {

    /// Jacobi Zeta function Z(φ, k) = E(φ, k) − (E(k) / K(k)) F(φ, k),
    /// where F/E denote the incomplete elliptic integrals of the first and second kind.
    ///
    /// - Parameters:
    ///   - phi: Amplitude φ (in radians). Must be finite.
    ///   - k: Modulus k (real). Must be finite.
    /// - Returns: `Z(φ, k)` as `T`.
    ///
    /// Boost’s implementation handles the typical branch structure and sign symmetry.
    /// This wrapper simply enforces finiteness before delegating to the bridge.
    @inlinable static func jacobiZeta<T: Real & BinaryFloatingPoint & Sendable>(_ phi: T, modulus k: T) throws -> T {
        let dphi = D(phi), dk = D(k)
        guard dphi.isFinite else { throw SpecialFunctionError<T>.parameterNotFinite(name: "phi", value: phi) }
        guard dk.isFinite else { throw SpecialFunctionError<T>.parameterNotFinite(name: "k", value: k) }
        return T(bs_jacobi_zeta_d(dk, dphi))
    }

    /// Jacobi Zeta Z(φ, k) for `Float`.
    @inlinable static func jacobiZeta(_ phi: Float, modulus k: Float) throws -> Float {
        guard phi.isFinite else { throw SpecialFunctionError<Float>.parameterNotFinite(name: "phi", value: phi) }
        guard k.isFinite else { throw SpecialFunctionError<Float>.parameterNotFinite(name: "k", value: k) }
        return bs_jacobi_zeta_f(k, phi)
    }

    #if arch(x86_64) || arch(i386)
    /// Jacobi Zeta Z(φ, k) for `Float80` (x86 only).
    @inlinable static func jacobiZeta(_ phi: Float80, modulus k: Float80) throws -> Float80 {
        guard phi.isFinite else { throw SpecialFunctionError<Float80>.parameterNotFinite(name: "phi", value: phi) }
        guard k.isFinite else { throw SpecialFunctionError<Float80>.parameterNotFinite(name: "k", value: k) }
        return bs_jacobi_zeta_l(k, phi)
    }
    #endif
}
