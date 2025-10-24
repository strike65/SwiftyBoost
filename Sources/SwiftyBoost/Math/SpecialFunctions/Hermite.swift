//
//  Hermite.swift
//  Math/SpecialFunctions
//
//  Created by Volker Thieme 2025.
//  License: MIT (see project root)
//
//  Swift wrappers for the physicists' Hermite polynomials Hₙ(x) and the
//  three-term recurrence helper exposed by Boost.Math. These wrappers validate
//  domain constraints (n ≥ 0, finite inputs) before delegating to the C bridge.
//

import SwiftyBoostPrelude

public extension SpecialFunctions {

    // MARK: - Generic overloads (Double-backed)

    /// Evaluate the physicists' Hermite polynomial Hₙ(x).
    ///
    /// - Parameters:
    ///   - n: Non-negative order.
    ///   - x: Evaluation point (finite real).
    /// - Returns: `Hₙ(x)` as `T`.
    @inlinable static func hermite<T: Real & BinaryFloatingPoint & Sendable>(_ n: Int, _ x: T) throws -> T {
        let nu = try validateHermiteOrder(n, as: T.self)
        let dx = D(x)
        guard dx.isFinite else { throw SpecialFunctionError<T>.parameterNotFinite(name: "x", value: x) }
        return T(bs_hermite_d(nu, dx))
    }

    /// Compute the next Hermite value using the three-term recurrence.
    ///
    /// Returns Hₙ₊₁(x) given Hₙ(x) and Hₙ₋₁(x).
    ///
    /// - Parameters:
    ///   - n: Order corresponding to `Hn` (requires `n ≥ 1`).
    ///   - x: Evaluation point (finite real).
    ///   - hn: `Hₙ(x)` (finite real).
    ///   - hnm1: `Hₙ₋₁(x)` (finite real).
    /// - Returns: `Hₙ₊₁(x)` as `T`.
    @inlinable static func hermiteNext<T: Real & BinaryFloatingPoint & Sendable>(_ n: Int, _ x: T, hn: T, hnm1: T) throws -> T {
        let nu = try validateHermiteNextOrder(n, as: T.self)
        let dx = D(x)
        guard dx.isFinite else { throw SpecialFunctionError<T>.parameterNotFinite(name: "x", value: x) }
        let dhn = D(hn)
        guard dhn.isFinite else { throw SpecialFunctionError<T>.parameterNotFinite(name: "hn", value: hn) }
        let dhnm1 = D(hnm1)
        guard dhnm1.isFinite else { throw SpecialFunctionError<T>.parameterNotFinite(name: "hnm1", value: hnm1) }
        return T(bs_hermite_next_d(nu, dx, dhn, dhnm1))
    }

    // MARK: - Float overloads

    @inlinable static func hermite(_ n: Int, _ x: Float) throws -> Float {
        let nu = try validateHermiteOrder(n, as: Float.self)
        guard x.isFinite else { throw SpecialFunctionError<Float>.parameterNotFinite(name: "x", value: x) }
        return bs_hermite_f(nu, x)
    }

    @inlinable static func hermiteNext(_ n: Int, _ x: Float, hn: Float, hnm1: Float) throws -> Float {
        let nu = try validateHermiteNextOrder(n, as: Float.self)
        guard x.isFinite else { throw SpecialFunctionError<Float>.parameterNotFinite(name: "x", value: x) }
        guard hn.isFinite else { throw SpecialFunctionError<Float>.parameterNotFinite(name: "hn", value: hn) }
        guard hnm1.isFinite else { throw SpecialFunctionError<Float>.parameterNotFinite(name: "hnm1", value: hnm1) }
        return bs_hermite_next_f(nu, x, hn, hnm1)
    }

    // MARK: - Float80 overloads (x86 only)

    #if arch(x86_64) || arch(i386)
    @inlinable static func hermite(_ n: Int, _ x: Float80) throws -> Float80 {
        let nu = try validateHermiteOrder(n, as: Float80.self)
        guard x.isFinite else { throw SpecialFunctionError<Float80>.parameterNotFinite(name: "x", value: x) }
        return bs_hermite_l(nu, x)
    }

    @inlinable static func hermiteNext(_ n: Int, _ x: Float80, hn: Float80, hnm1: Float80) throws -> Float80 {
        let nu = try validateHermiteNextOrder(n, as: Float80.self)
        guard x.isFinite else { throw SpecialFunctionError<Float80>.parameterNotFinite(name: "x", value: x) }
        guard hn.isFinite else { throw SpecialFunctionError<Float80>.parameterNotFinite(name: "hn", value: hn) }
        guard hnm1.isFinite else { throw SpecialFunctionError<Float80>.parameterNotFinite(name: "hnm1", value: hnm1) }
        return bs_hermite_next_l(nu, x, hn, hnm1)
    }
    #endif

    // MARK: - Validation helpers

    @usableFromInline internal static func validateHermiteOrder<T: Real & BinaryFloatingPoint & Sendable>(_ n: Int, as _: T.Type) throws -> UInt32 {
        guard n >= 0 else { throw SpecialFunctionError<T>.parameterNotPositive(name: "n", value: T(n)) }
        guard let nu = UInt32(exactly: n) else {
            throw SpecialFunctionError<T>.parameterExceedsMaximumIntegerValue(name: "n", max: Int(UInt32.max))
        }
        return nu
    }

    @usableFromInline internal static func validateHermiteNextOrder<T: Real & BinaryFloatingPoint & Sendable>(_ n: Int, as _: T.Type) throws -> UInt32 {
        guard n >= 1 else { throw SpecialFunctionError<T>.parameterNotPositive(name: "n", value: T(n)) }
        guard let nu = UInt32(exactly: n) else {
            throw SpecialFunctionError<T>.parameterExceedsMaximumIntegerValue(name: "n", max: Int(UInt32.max))
        }
        return nu
    }
}
