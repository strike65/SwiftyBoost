//
//  JacobiElliptic.swift
//  Math/SpecialFunctions
//
//  Created by Volker Thieme 2025.
//  License: MIT (see project root)
//
//  Swift wrappers for the Jacobi elliptic functions sn, cn, dn and the ten
//  quotients cd/cs/dc/ds/nc/nd/ns/sc/sd. These delegate to Boost.Math via the
//  C bridge and mirror the standard SwiftyBoost validation conventions
//  (finiteness of modulus k and argument θ).
//

import SwiftyBoostPrelude

public extension SpecialFunctions {

    // MARK: - Generic (Double-backed) overloads

    /// Computes the three primary Jacobi elliptic functions sn(θ|k), cn(θ|k), and dn(θ|k) simultaneously.
    ///
    /// - Parameters:
    ///   - k: The elliptic modulus (not the parameter m = k²). Must be finite.
    ///   - theta: The argument θ (amplitude or phase) in radians. Must be finite.
    /// - Returns: A tuple `(sn, cn, dn)` where:
    ///   - `sn` = sn(θ|k)
    ///   - `cn` = cn(θ|k)
    ///   - `dn` = dn(θ|k)
    /// - Throws: `SpecialFunctionError.parameterNotFinite` if either `k` or `theta` is not finite.
    /// - Discussion:
    ///   - This generic overload accepts any `BinaryFloatingPoint` conforming to `Real & Sendable`.
    ///   - Internally computes using `Double` precision and converts back to `T`.
    ///   - Periodicity:
    ///     - sn and cn are 4K(k)-periodic in θ; dn is 2K(k)-periodic (K is the complete elliptic integral of the first kind).
    /// - SeeAlso: `jacobiEllipticSn(_:theta:)`, `jacobiEllipticCn(_:theta:)`, `jacobiEllipticDn(_:theta:)`
    @inlinable static func jacobiElliptic<T: Real & BinaryFloatingPoint & Sendable>(_ k: T, theta: T) throws -> (sn: T, cn: T, dn: T) {
        var cn: Double = .zero
        var dn: Double = .zero
        let (dk, dtheta) = try checkedInputs(k, theta, as: T.self)
        let sn = T(bs_jacobi_elliptic_d(dk, dtheta, &cn, &dn))
        return (sn, T(cn), T(dn))
    }

    /// Computes the Jacobi elliptic function sn(θ|k).
    ///
    /// - Parameters:
    ///   - k: Elliptic modulus (finite).
    ///   - theta: Argument θ in radians (finite).
    /// - Returns: sn(θ|k).
    /// - Throws: `SpecialFunctionError.parameterNotFinite` if inputs are not finite.
    @inlinable static func jacobiEllipticSn<T: Real & BinaryFloatingPoint & Sendable>(_ k: T, theta: T) throws -> T {
        let (dk, dtheta) = try checkedInputs(k, theta, as: T.self)
        return T(bs_jacobi_elliptic_sn_d(dk, dtheta))
    }

    /// Computes the Jacobi elliptic function cn(θ|k).
    ///
    /// - Parameters:
    ///   - k: Elliptic modulus (finite).
    ///   - theta: Argument θ in radians (finite).
    /// - Returns: cn(θ|k).
    /// - Throws: `SpecialFunctionError.parameterNotFinite` if inputs are not finite.
    @inlinable static func jacobiEllipticCn<T: Real & BinaryFloatingPoint & Sendable>(_ k: T, theta: T) throws -> T {
        let (dk, dtheta) = try checkedInputs(k, theta, as: T.self)
        return T(bs_jacobi_elliptic_cn_d(dk, dtheta))
    }

    /// Computes the Jacobi elliptic function dn(θ|k).
    ///
    /// - Parameters:
    ///   - k: Elliptic modulus (finite).
    ///   - theta: Argument θ in radians (finite).
    /// - Returns: dn(θ|k).
    /// - Throws: `SpecialFunctionError.parameterNotFinite` if inputs are not finite.
    @inlinable static func jacobiEllipticDn<T: Real & BinaryFloatingPoint & Sendable>(_ k: T, theta: T) throws -> T {
        let (dk, dtheta) = try checkedInputs(k, theta, as: T.self)
        return T(bs_jacobi_elliptic_dn_d(dk, dtheta))
    }

    /// Computes the Jacobi quotient sd(θ|k) = sn(θ|k)/dn(θ|k).
    ///
    /// - Parameters:
    ///   - k: Elliptic modulus (finite).
    ///   - theta: Argument θ in radians (finite).
    /// - Returns: sd(θ|k).
    /// - Throws: `SpecialFunctionError.parameterNotFinite` if inputs are not finite.
    @inlinable static func jacobiEllipticSd<T: Real & BinaryFloatingPoint & Sendable>(_ k: T, theta: T) throws -> T {
        let (dk, dtheta) = try checkedInputs(k, theta, as: T.self)
        return T(bs_jacobi_elliptic_sd_d(dk, dtheta))
    }

    /// Computes the Jacobi quotient sc(θ|k) = sn(θ|k)/cn(θ|k).
    ///
    /// - Parameters:
    ///   - k: Elliptic modulus (finite).
    ///   - theta: Argument θ in radians (finite).
    /// - Returns: sc(θ|k).
    /// - Throws: `SpecialFunctionError.parameterNotFinite` if inputs are not finite.
    @inlinable static func jacobiEllipticSc<T: Real & BinaryFloatingPoint & Sendable>(_ k: T, theta: T) throws -> T {
        let (dk, dtheta) = try checkedInputs(k, theta, as: T.self)
        return T(bs_jacobi_elliptic_sc_d(dk, dtheta))
    }

    /// Computes the Jacobi quotient ns(θ|k) = 1/sn(θ|k).
    ///
    /// - Parameters:
    ///   - k: Elliptic modulus (finite).
    ///   - theta: Argument θ in radians (finite).
    /// - Returns: ns(θ|k).
    /// - Throws: `SpecialFunctionError.parameterNotFinite` if inputs are not finite.
    /// - Note: May overflow or underflow near the zeros of sn.
    @inlinable static func jacobiEllipticNs<T: Real & BinaryFloatingPoint & Sendable>(_ k: T, theta: T) throws -> T {
        let (dk, dtheta) = try checkedInputs(k, theta, as: T.self)
        return T(bs_jacobi_elliptic_ns_d(dk, dtheta))
    }

    /// Computes the Jacobi quotient nd(θ|k) = 1/dn(θ|k).
    ///
    /// - Parameters:
    ///   - k: Elliptic modulus (finite).
    ///   - theta: Argument θ in radians (finite).
    /// - Returns: nd(θ|k).
    /// - Throws: `SpecialFunctionError.parameterNotFinite` if inputs are not finite.
    /// - Note: May overflow or underflow near the zeros of dn.
    @inlinable static func jacobiEllipticNd<T: Real & BinaryFloatingPoint & Sendable>(_ k: T, theta: T) throws -> T {
        let (dk, dtheta) = try checkedInputs(k, theta, as: T.self)
        return T(bs_jacobi_elliptic_nd_d(dk, dtheta))
    }

    /// Computes the Jacobi quotient nc(θ|k) = 1/cn(θ|k).
    ///
    /// - Parameters:
    ///   - k: Elliptic modulus (finite).
    ///   - theta: Argument θ in radians (finite).
    /// - Returns: nc(θ|k).
    /// - Throws: `SpecialFunctionError.parameterNotFinite` if inputs are not finite.
    /// - Note: May overflow or underflow near the zeros of cn.
    @inlinable static func jacobiEllipticNc<T: Real & BinaryFloatingPoint & Sendable>(_ k: T, theta: T) throws -> T {
        let (dk, dtheta) = try checkedInputs(k, theta, as: T.self)
        return T(bs_jacobi_elliptic_nc_d(dk, dtheta))
    }

    /// Computes the Jacobi quotient ds(θ|k) = dn(θ|k)/sn(θ|k).
    ///
    /// - Parameters:
    ///   - k: Elliptic modulus (finite).
    ///   - theta: Argument θ in radians (finite).
    /// - Returns: ds(θ|k).
    /// - Throws: `SpecialFunctionError.parameterNotFinite` if inputs are not finite.
    @inlinable static func jacobiEllipticDs<T: Real & BinaryFloatingPoint & Sendable>(_ k: T, theta: T) throws -> T {
        let (dk, dtheta) = try checkedInputs(k, theta, as: T.self)
        return T(bs_jacobi_elliptic_ds_d(dk, dtheta))
    }

    /// Computes the Jacobi quotient dc(θ|k) = dn(θ|k)/cn(θ|k).
    ///
    /// - Parameters:
    ///   - k: Elliptic modulus (finite).
    ///   - theta: Argument θ in radians (finite).
    /// - Returns: dc(θ|k).
    /// - Throws: `SpecialFunctionError.parameterNotFinite` if inputs are not finite.
    @inlinable static func jacobiEllipticDc<T: Real & BinaryFloatingPoint & Sendable>(_ k: T, theta: T) throws -> T {
        let (dk, dtheta) = try checkedInputs(k, theta, as: T.self)
        return T(bs_jacobi_elliptic_dc_d(dk, dtheta))
    }

    /// Computes the Jacobi quotient cs(θ|k) = cn(θ|k)/sn(θ|k).
    ///
    /// - Parameters:
    ///   - k: Elliptic modulus (finite).
    ///   - theta: Argument θ in radians (finite).
    /// - Returns: cs(θ|k).
    /// - Throws: `SpecialFunctionError.parameterNotFinite` if inputs are not finite.
    @inlinable static func jacobiEllipticCs<T: Real & BinaryFloatingPoint & Sendable>(_ k: T, theta: T) throws -> T {
        let (dk, dtheta) = try checkedInputs(k, theta, as: T.self)
        return T(bs_jacobi_elliptic_cs_d(dk, dtheta))
    }

    /// Computes the Jacobi quotient cd(θ|k) = cn(θ|k)/dn(θ|k).
    ///
    /// - Parameters:
    ///   - k: Elliptic modulus (finite).
    ///   - theta: Argument θ in radians (finite).
    /// - Returns: cd(θ|k).
    /// - Throws: `SpecialFunctionError.parameterNotFinite` if inputs are not finite.
    @inlinable static func jacobiEllipticCd<T: Real & BinaryFloatingPoint & Sendable>(_ k: T, theta: T) throws -> T {
        let (dk, dtheta) = try checkedInputs(k, theta, as: T.self)
        return T(bs_jacobi_elliptic_cd_d(dk, dtheta))
    }

    // MARK: - Float overloads

    /// Float-specialized overload computing sn(θ|k), cn(θ|k), dn(θ|k) simultaneously.
    ///
    /// - Parameters:
    ///   - k: Elliptic modulus (finite).
    ///   - theta: Argument θ in radians (finite).
    /// - Returns: `(sn, cn, dn)` for `Float`.
    /// - Throws: `SpecialFunctionError.parameterNotFinite` if inputs are not finite.
    /// - Performance: Calls the native `float` C bridge for Boost.Math.
    @inlinable static func jacobiElliptic(_ k: Float, theta: Float) throws -> (sn: Float, cn: Float, dn: Float) {
        try checkedFloat(k, theta)
        var cn: Float = .zero
        var dn: Float = .zero
        let sn = bs_jacobi_elliptic_f(k, theta, &cn, &dn)
        return (sn, cn, dn)
    }

    /// Float-specialized sn(θ|k).
    ///
    /// - Throws: `SpecialFunctionError.parameterNotFinite` if inputs are not finite.
    @inlinable static func jacobiEllipticSn(_ k: Float, theta: Float) throws -> Float {
        try checkedFloat(k, theta)
        return bs_jacobi_elliptic_sn_f(k, theta)
    }

    /// Float-specialized cn(θ|k).
    ///
    /// - Throws: `SpecialFunctionError.parameterNotFinite` if inputs are not finite.
    @inlinable static func jacobiEllipticCn(_ k: Float, theta: Float) throws -> Float {
        try checkedFloat(k, theta)
        return bs_jacobi_elliptic_cn_f(k, theta)
    }

    /// Float-specialized dn(θ|k).
    ///
    /// - Throws: `SpecialFunctionError.parameterNotFinite` if inputs are not finite.
    @inlinable static func jacobiEllipticDn(_ k: Float, theta: Float) throws -> Float {
        try checkedFloat(k, theta)
        return bs_jacobi_elliptic_dn_f(k, theta)
    }

    /// Float-specialized sd(θ|k) = sn/dn.
    ///
    /// - Throws: `SpecialFunctionError.parameterNotFinite` if inputs are not finite.
    @inlinable static func jacobiEllipticSd(_ k: Float, theta: Float) throws -> Float {
        try checkedFloat(k, theta)
        return bs_jacobi_elliptic_sd_f(k, theta)
    }

    /// Float-specialized sc(θ|k) = sn/cn.
    ///
    /// - Throws: `SpecialFunctionError.parameterNotFinite` if inputs are not finite.
    @inlinable static func jacobiEllipticSc(_ k: Float, theta: Float) throws -> Float {
        try checkedFloat(k, theta)
        return bs_jacobi_elliptic_sc_f(k, theta)
    }

    /// Float-specialized ns(θ|k) = 1/sn.
    ///
    /// - Throws: `SpecialFunctionError.parameterNotFinite` if inputs are not finite.
    @inlinable static func jacobiEllipticNs(_ k: Float, theta: Float) throws -> Float {
        try checkedFloat(k, theta)
        return bs_jacobi_elliptic_ns_f(k, theta)
    }

    /// Float-specialized nd(θ|k) = 1/dn.
    ///
    /// - Throws: `SpecialFunctionError.parameterNotFinite` if inputs are not finite.
    @inlinable static func jacobiEllipticNd(_ k: Float, theta: Float) throws -> Float {
        try checkedFloat(k, theta)
        return bs_jacobi_elliptic_nd_f(k, theta)
    }

    /// Float-specialized nc(θ|k) = 1/cn.
    ///
    /// - Throws: `SpecialFunctionError.parameterNotFinite` if inputs are not finite.
    @inlinable static func jacobiEllipticNc(_ k: Float, theta: Float) throws -> Float {
        try checkedFloat(k, theta)
        return bs_jacobi_elliptic_nc_f(k, theta)
    }

    /// Float-specialized ds(θ|k) = dn/sn.
    ///
    /// - Throws: `SpecialFunctionError.parameterNotFinite` if inputs are not finite.
    @inlinable static func jacobiEllipticDs(_ k: Float, theta: Float) throws -> Float {
        try checkedFloat(k, theta)
        return bs_jacobi_elliptic_ds_f(k, theta)
    }

    /// Float-specialized dc(θ|k) = dn/cn.
    ///
    /// - Throws: `SpecialFunctionError.parameterNotFinite` if inputs are not finite.
    @inlinable static func jacobiEllipticDc(_ k: Float, theta: Float) throws -> Float {
        try checkedFloat(k, theta)
        return bs_jacobi_elliptic_dc_f(k, theta)
    }

    /// Float-specialized cs(θ|k) = cn/sn.
    ///
    /// - Throws: `SpecialFunctionError.parameterNotFinite` if inputs are not finite.
    @inlinable static func jacobiEllipticCs(_ k: Float, theta: Float) throws -> Float {
        try checkedFloat(k, theta)
        return bs_jacobi_elliptic_cs_f(k, theta)
    }

    /// Float-specialized cd(θ|k) = cn/dn.
    ///
    /// - Throws: `SpecialFunctionError.parameterNotFinite` if inputs are not finite.
    @inlinable static func jacobiEllipticCd(_ k: Float, theta: Float) throws -> Float {
        try checkedFloat(k, theta)
        return bs_jacobi_elliptic_cd_f(k, theta)
    }

    // MARK: - Float80 overloads (x86 only)

#if arch(x86_64) || arch(i386)
    /// Float80-specialized overload computing sn(θ|k), cn(θ|k), dn(θ|k) simultaneously (x86 only).
    ///
    /// - Parameters:
    ///   - k: Elliptic modulus (finite).
    ///   - theta: Argument θ in radians (finite).
    /// - Returns: `(sn, cn, dn)` for `Float80`.
    /// - Throws: `SpecialFunctionError.parameterNotFinite` if inputs are not finite.
    /// - Availability: Only available on x86 architectures that support `Float80`.
    @inlinable static func jacobiElliptic(_ k: Float80, theta: Float80) throws -> (sn: Float80, cn: Float80, dn: Float80) {
        try checkedFloat80(k, theta)
        var cn: Float80 = .zero
        var dn: Float80 = .zero
        let sn = bs_jacobi_elliptic_l(k, theta, &cn, &dn)
        return (sn, cn, dn)
    }

    /// Float80-specialized sn(θ|k) (x86 only).
    ///
    /// - Throws: `SpecialFunctionError.parameterNotFinite` if inputs are not finite.
    /// - Availability: Only available on x86 architectures that support `Float80`.
    @inlinable static func jacobiEllipticSn(_ k: Float80, theta: Float80) throws -> Float80 {
        try checkedFloat80(k, theta)
        return bs_jacobi_elliptic_sn_l(k, theta)
    }

    /// Float80-specialized cn(θ|k) (x86 only).
    ///
    /// - Throws: `SpecialFunctionError.parameterNotFinite` if inputs are not finite.
    /// - Availability: Only available on x86 architectures that support `Float80`.
    @inlinable static func jacobiEllipticCn(_ k: Float80, theta: Float80) throws -> Float80 {
        try checkedFloat80(k, theta)
        return bs_jacobi_elliptic_cn_l(k, theta)
    }

    /// Float80-specialized dn(θ|k) (x86 only).
    ///
    /// - Throws: `SpecialFunctionError.parameterNotFinite` if inputs are not finite.
    /// - Availability: Only available on x86 architectures that support `Float80`.
    @inlinable static func jacobiEllipticDn(_ k: Float80, theta: Float80) throws -> Float80 {
        try checkedFloat80(k, theta)
        return bs_jacobi_elliptic_dn_l(k, theta)
    }

    /// Float80-specialized sd(θ|k) = sn/dn (x86 only).
    ///
    /// - Throws: `SpecialFunctionError.parameterNotFinite` if inputs are not finite.
    /// - Availability: Only available on x86 architectures that support `Float80`.
    @inlinable static func jacobiEllipticSd(_ k: Float80, theta: Float80) throws -> Float80 {
        try checkedFloat80(k, theta)
        return bs_jacobi_elliptic_sd_l(k, theta)
    }

    /// Float80-specialized sc(θ|k) = sn/cn (x86 only).
    ///
    /// - Throws: `SpecialFunctionError.parameterNotFinite` if inputs are not finite.
    /// - Availability: Only available on x86 architectures that support `Float80`.
    @inlinable static func jacobiEllipticSc(_ k: Float80, theta: Float80) throws -> Float80 {
        try checkedFloat80(k, theta)
        return bs_jacobi_elliptic_sc_l(k, theta)
    }

    /// Float80-specialized ns(θ|k) = 1/sn (x86 only).
    ///
    /// - Throws: `SpecialFunctionError.parameterNotFinite` if inputs are not finite.
    /// - Availability: Only available on x86 architectures that support `Float80`.
    @inlinable static func jacobiEllipticNs(_ k: Float80, theta: Float80) throws -> Float80 {
        try checkedFloat80(k, theta)
        return bs_jacobi_elliptic_ns_l(k, theta)
    }

    /// Float80-specialized nd(θ|k) = 1/dn (x86 only).
    ///
    /// - Throws: `SpecialFunctionError.parameterNotFinite` if inputs are not finite.
    /// - Availability: Only available on x86 architectures that support `Float80`.
    @inlinable static func jacobiEllipticNd(_ k: Float80, theta: Float80) throws -> Float80 {
        try checkedFloat80(k, theta)
        return bs_jacobi_elliptic_nd_l(k, theta)
    }

    /// Float80-specialized nc(θ|k) = 1/cn (x86 only).
    ///
    /// - Throws: `SpecialFunctionError.parameterNotFinite` if inputs are not finite.
    /// - Availability: Only available on x86 architectures that support `Float80`.
    @inlinable static func jacobiEllipticNc(_ k: Float80, theta: Float80) throws -> Float80 {
        try checkedFloat80(k, theta)
        return bs_jacobi_elliptic_nc_l(k, theta)
    }

    /// Float80-specialized ds(θ|k) = dn/sn (x86 only).
    ///
    /// - Throws: `SpecialFunctionError.parameterNotFinite` if inputs are not finite.
    /// - Availability: Only available on x86 architectures that support `Float80`.
    @inlinable static func jacobiEllipticDs(_ k: Float80, theta: Float80) throws -> Float80 {
        try checkedFloat80(k, theta)
        return bs_jacobi_elliptic_ds_l(k, theta)
    }

    /// Float80-specialized dc(θ|k) = dn/cn (x86 only).
    ///
    /// - Throws: `SpecialFunctionError.parameterNotFinite` if inputs are not finite.
    /// - Availability: Only available on x86 architectures that support `Float80`.
    @inlinable static func jacobiEllipticDc(_ k: Float80, theta: Float80) throws -> Float80 {
        try checkedFloat80(k, theta)
        return bs_jacobi_elliptic_dc_l(k, theta)
    }

    /// Float80-specialized cs(θ|k) = cn/sn (x86 only).
    ///
    /// - Throws: `SpecialFunctionError.parameterNotFinite` if inputs are not finite.
    /// - Availability: Only available on x86 architectures that support `Float80`.
    @inlinable static func jacobiEllipticCs(_ k: Float80, theta: Float80) throws -> Float80 {
        try checkedFloat80(k, theta)
        return bs_jacobi_elliptic_cs_l(k, theta)
    }

    /// Float80-specialized cd(θ|k) = cn/dn (x86 only).
    ///
    /// - Throws: `SpecialFunctionError.parameterNotFinite` if inputs are not finite.
    /// - Availability: Only available on x86 architectures that support `Float80`.
    @inlinable static func jacobiEllipticCd(_ k: Float80, theta: Float80) throws -> Float80 {
        try checkedFloat80(k, theta)
        return bs_jacobi_elliptic_cd_l(k, theta)
    }
    #endif

    // MARK: - Validation helpers

    /// Validates generic inputs and converts them to `Double` for the C bridge.
    ///
    /// - Parameters:
    ///   - k: Elliptic modulus as generic `T`.
    ///   - theta: Argument θ as generic `T`.
    ///   - as: Type witness (unused).
    /// - Returns: `(Double(k), Double(theta))` if both are finite.
    /// - Throws: `SpecialFunctionError.parameterNotFinite` if either is not finite.
    @usableFromInline internal static func checkedInputs<T: Real & BinaryFloatingPoint & Sendable>(_ k: T, _ theta: T, as _: T.Type) throws -> (Double, Double) {
        let dk = D(k), dtheta = D(theta)
        guard dk.isFinite else { throw SpecialFunctionError<T>.parameterNotFinite(name: "k", value: k) }
        guard dtheta.isFinite else { throw SpecialFunctionError<T>.parameterNotFinite(name: "theta", value: theta) }
        return (dk, dtheta)
    }

    /// Validates `Float` inputs for finiteness.
    ///
    /// - Throws: `SpecialFunctionError.parameterNotFinite` if either argument is not finite.
    @usableFromInline internal static func checkedFloat(_ k: Float, _ theta: Float) throws {
        guard k.isFinite else { throw SpecialFunctionError<Float>.parameterNotFinite(name: "k", value: k) }
        guard theta.isFinite else { throw SpecialFunctionError<Float>.parameterNotFinite(name: "theta", value: theta) }
    }

    #if arch(x86_64) || arch(i386)
    /// Validates `Float80` inputs for finiteness (x86 only).
    ///
    /// - Throws: `SpecialFunctionError.parameterNotFinite` if either argument is not finite.
    @usableFromInline internal static func checkedFloat80(_ k: Float80, _ theta: Float80) throws {
        guard k.isFinite else { throw SpecialFunctionError<Float80>.parameterNotFinite(name: "k", value: k) }
        guard theta.isFinite else { throw SpecialFunctionError<Float80>.parameterNotFinite(name: "theta", value: theta) }
    }
    #endif
}
