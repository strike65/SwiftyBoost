//
//  Airy.swift
//  Math/SpecialFunctions
//
//  Created by Volker Thieme 2025.
//  License: MIT (see project root)
//
//  Swift wrappers for the Airy functions Ai(x), Bi(x) and their first derivatives.
//  These wrappers validate finite inputs then delegate to Boost.Math via the C bridge.
//

import SwiftyBoostPrelude

public extension SpecialFunctions {

    // MARK: - Generic (Double-backed) overloads

    /// Airy Ai(x).
    @inlinable static func airyAi<T: Real & BinaryFloatingPoint & Sendable>(_ x: T) throws -> T {
        let dx = D(x)
        guard dx.isFinite else { throw SpecialFunctionError<T>.parameterNotFinite(name: "x", value: x) }
        return T(bs_airy_ai_d(dx))
    }

    /// Airy Bi(x).
    @inlinable static func airyBi<T: Real & BinaryFloatingPoint & Sendable>(_ x: T) throws -> T {
        let dx = D(x)
        guard dx.isFinite else { throw SpecialFunctionError<T>.parameterNotFinite(name: "x", value: x) }
        return T(bs_airy_bi_d(dx))
    }

    /// Airy Ai′(x).
    @inlinable static func airyAiPrime<T: Real & BinaryFloatingPoint & Sendable>(_ x: T) throws -> T {
        let dx = D(x)
        guard dx.isFinite else { throw SpecialFunctionError<T>.parameterNotFinite(name: "x", value: x) }
        return T(bs_airy_ai_prime_d(dx))
    }

    /// Airy Bi′(x).
    @inlinable static func airyBiPrime<T: Real & BinaryFloatingPoint & Sendable>(_ x: T) throws -> T {
        let dx = D(x)
        guard dx.isFinite else { throw SpecialFunctionError<T>.parameterNotFinite(name: "x", value: x) }
        return T(bs_airy_bi_prime_d(dx))
    }

    /// n-th zero of Airy Ai(x) (Zero-based index: `n = 0` returns the first zero).
    @inlinable static func airyAiZero<T: Real & BinaryFloatingPoint & Sendable>(_ n: Int) throws -> T {
        let idx = try validateAiryIndex(n, as: T.self)
        return T(bs_airy_ai_zero_d(idx))
    }

    /// n-th zero of Airy Bi(x) (Zero-based index).
    @inlinable static func airyBiZero<T: Real & BinaryFloatingPoint & Sendable>(_ n: Int) throws -> T {
        let idx = try validateAiryIndex(n, as: T.self)
        return T(bs_airy_bi_zero_d(idx))
    }

    /// Consecutive zeros of Airy Ai starting at `startIndex`.
    @inlinable static func airyAiZeros<T: Real & BinaryFloatingPoint & Sendable>(startIndex: Int, count: Int) throws -> [T] {
        let si = try validateAiryIndex(startIndex, as: T.self)
        guard count >= 0 else { throw SpecialFunctionError<T>.parameterNotPositive(name: "count", value: T(count)) }
        guard count > 0 else { return [] }
        guard let uCount = UInt32(exactly: count) else {
            throw SpecialFunctionError<T>.parameterExceedsMaximumIntegerValue(name: "count", max: Int(UInt32.max))
        }
        var tmp = Array<Double>(repeating: .zero, count: count)
        tmp.withUnsafeMutableBufferPointer { buf in
            bs_airy_ai_zeros_d(si, uCount, buf.baseAddress)
        }
        return tmp.map(T.init)
    }

    /// Consecutive zeros of Airy Bi starting at `startIndex`.
    @inlinable static func airyBiZeros<T: Real & BinaryFloatingPoint & Sendable>(startIndex: Int, count: Int) throws -> [T] {
        let si = try validateAiryIndex(startIndex, as: T.self)
        guard count >= 0 else { throw SpecialFunctionError<T>.parameterNotPositive(name: "count", value: T(count)) }
        guard count > 0 else { return [] }
        guard let uCount = UInt32(exactly: count) else {
            throw SpecialFunctionError<T>.parameterExceedsMaximumIntegerValue(name: "count", max: Int(UInt32.max))
        }
        var tmp = Array<Double>(repeating: .zero, count: count)
        tmp.withUnsafeMutableBufferPointer { buf in
            bs_airy_bi_zeros_d(si, uCount, buf.baseAddress)
        }
        return tmp.map(T.init)
    }

    // MARK: - Float overloads

    @inlinable static func airyAi(_ x: Float) throws -> Float {
        guard x.isFinite else { throw SpecialFunctionError<Float>.parameterNotFinite(name: "x", value: x) }
        return bs_airy_ai_f(x)
    }

    @inlinable static func airyBi(_ x: Float) throws -> Float {
        guard x.isFinite else { throw SpecialFunctionError<Float>.parameterNotFinite(name: "x", value: x) }
        return bs_airy_bi_f(x)
    }

    @inlinable static func airyAiPrime(_ x: Float) throws -> Float {
        guard x.isFinite else { throw SpecialFunctionError<Float>.parameterNotFinite(name: "x", value: x) }
        return bs_airy_ai_prime_f(x)
    }

    @inlinable static func airyBiPrime(_ x: Float) throws -> Float {
        guard x.isFinite else { throw SpecialFunctionError<Float>.parameterNotFinite(name: "x", value: x) }
        return bs_airy_bi_prime_f(x)
    }

    @inlinable static func airyAiZero(_ n: Int) throws -> Float {
        let idx = try validateAiryIndex(n, as: Float.self)
        return bs_airy_ai_zero_f(idx)
    }

    @inlinable static func airyBiZero(_ n: Int) throws -> Float {
        let idx = try validateAiryIndex(n, as: Float.self)
        return bs_airy_bi_zero_f(idx)
    }

    @inlinable static func airyAiZeros(startIndex: Int, count: Int) throws -> [Float] {
        let si = try validateAiryIndex(startIndex, as: Float.self)
        guard count >= 0 else { throw SpecialFunctionError<Float>.parameterNotPositive(name: "count", value: Float(count)) }
        guard count > 0 else { return [] }
        guard let uCount = UInt32(exactly: count) else {
            throw SpecialFunctionError<Float>.parameterExceedsMaximumIntegerValue(name: "count", max: Int(UInt32.max))
        }
        var tmp = Array<Float>(repeating: .zero, count: count)
        tmp.withUnsafeMutableBufferPointer { buf in
            bs_airy_ai_zeros_f(si, uCount, buf.baseAddress)
        }
        return tmp
    }

    @inlinable static func airyBiZeros(startIndex: Int, count: Int) throws -> [Float] {
        let si = try validateAiryIndex(startIndex, as: Float.self)
        guard count >= 0 else { throw SpecialFunctionError<Float>.parameterNotPositive(name: "count", value: Float(count)) }
        guard count > 0 else { return [] }
        guard let uCount = UInt32(exactly: count) else {
            throw SpecialFunctionError<Float>.parameterExceedsMaximumIntegerValue(name: "count", max: Int(UInt32.max))
        }
        var tmp = Array<Float>(repeating: .zero, count: count)
        tmp.withUnsafeMutableBufferPointer { buf in
            bs_airy_bi_zeros_f(si, uCount, buf.baseAddress)
        }
        return tmp
    }

    // MARK: - Float80 overloads (x86 only)

    #if arch(x86_64) || arch(i386)
    @inlinable static func airyAi(_ x: Float80) throws -> Float80 {
        guard x.isFinite else { throw SpecialFunctionError<Float80>.parameterNotFinite(name: "x", value: x) }
        return bs_airy_ai_l(x)
    }

    @inlinable static func airyBi(_ x: Float80) throws -> Float80 {
        guard x.isFinite else { throw SpecialFunctionError<Float80>.parameterNotFinite(name: "x", value: x) }
        return bs_airy_bi_l(x)
    }

    @inlinable static func airyAiPrime(_ x: Float80) throws -> Float80 {
        guard x.isFinite else { throw SpecialFunctionError<Float80>.parameterNotFinite(name: "x", value: x) }
        return bs_airy_ai_prime_l(x)
    }

    @inlinable static func airyBiPrime(_ x: Float80) throws -> Float80 {
        guard x.isFinite else { throw SpecialFunctionError<Float80>.parameterNotFinite(name: "x", value: x) }
        return bs_airy_bi_prime_l(x)
    }

    @inlinable static func airyAiZero(_ n: Int) throws -> Float80 {
        let idx = try validateAiryIndex(n, as: Float80.self)
        return bs_airy_ai_zero_l(idx)
    }

    @inlinable static func airyBiZero(_ n: Int) throws -> Float80 {
        let idx = try validateAiryIndex(n, as: Float80.self)
        return bs_airy_bi_zero_l(idx)
    }

    @inlinable static func airyAiZeros(startIndex: Int, count: Int) throws -> [Float80] {
        let si = try validateAiryIndex(startIndex, as: Float80.self)
        guard count >= 0 else { throw SpecialFunctionError<Float80>.parameterNotPositive(name: "count", value: Float80(count)) }
        guard count > 0 else { return [] }
        guard let uCount = UInt32(exactly: count) else {
            throw SpecialFunctionError<Float80>.parameterExceedsMaximumIntegerValue(name: "count", max: Int(UInt32.max))
        }
        var tmp = Array<Float80>(repeating: .zero, count: count)
        tmp.withUnsafeMutableBufferPointer { buf in
            bs_airy_ai_zeros_l(si, uCount, buf.baseAddress)
        }
        return tmp
    }

    @inlinable static func airyBiZeros(startIndex: Int, count: Int) throws -> [Float80] {
        let si = try validateAiryIndex(startIndex, as: Float80.self)
        guard count >= 0 else { throw SpecialFunctionError<Float80>.parameterNotPositive(name: "count", value: Float80(count)) }
        guard count > 0 else { return [] }
        guard let uCount = UInt32(exactly: count) else {
            throw SpecialFunctionError<Float80>.parameterExceedsMaximumIntegerValue(name: "count", max: Int(UInt32.max))
        }
        var tmp = Array<Float80>(repeating: .zero, count: count)
        tmp.withUnsafeMutableBufferPointer { buf in
            bs_airy_bi_zeros_l(si, uCount, buf.baseAddress)
        }
        return tmp
    }
    #endif

    // MARK: - Validation

    @usableFromInline internal static func validateAiryIndex<T: Real & BinaryFloatingPoint & Sendable>(_ n: Int, as _: T.Type) throws -> Int32 {
        guard n >= 0 else { throw SpecialFunctionError<T>.parameterNotPositive(name: "n", value: T(n)) }
        guard let idx = Int32(exactly: n) else {
            throw SpecialFunctionError<T>.parameterExceedsMaximumIntegerValue(name: "n", max: Int(Int32.max))
        }
        return idx
    }
}
