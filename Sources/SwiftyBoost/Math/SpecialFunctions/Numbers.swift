//
//  Numbers.swift
//  Math/SpecialFunctions
//
//  Swift wrappers for Bernoulli numbers B_{2n} and Tangent numbers T_{2n}.
//  These APIs provide scalar and vector entry points, validate simple argument
//  expectations, and delegate the numerical work to CBoostBridge (Boost.Math).
//
//  Overview
//  - Bernoulli numbers (even index):
//    • bernoulli_b2n(_:)   → B_{2n} as Double
//    • bernoulli_b2n_f(_:) → B_{2n} as Float
//    • bernoulli_b2n_l(_:) → B_{2n} as Float80 (x86_64/i386)
//  - Tangent numbers (even index; aka “zag numbers”, coefficients in tan series):
//    • tangent_t2n(_:)     → T_{2n} as Double
//    • tangent_t2n_f(_:)   → T_{2n} as Float
//    • tangent_t2n_l(_:)   → T_{2n} as Float80 (x86_64/i386)
//    • tangent_t2n(count:) → [T2, T4, …, T_{2·count}] as [Double]
//    • tangent_t2n_f(count:) → same as [Float]
//    • tangent_t2n_l(count:) → same as [Float80] (x86_64/i386)
//  - Prime numbers:
//    • prime(_:)           → nth prime (bridged from CBoostBridge)
//  - Fibonacci (optional convenience):
//    • fibonacci(_:)       → nth Fibonacci number as UInt64
//
//  Indexing convention
//  - Both families use the half-index n to denote the even-indexed values:
//    B_{2n} and T_{2n}. In particular, tangent numbers appear in the Maclaurin
//    series tan(x) = Σ_{n≥1} T_{2n} x^{2n−1} / (2n−1)!
//
//  Notes
//  - These wrappers assume finite integer indices and leave numeric details to Boost.Math.
//  - Vector functions are efficient: Swift preallocates the array and the C bridge fills it.
//  - For Float80 variants, availability is limited to x86_64/i386 architectures.
//
//  References
//  - NIST DLMF: §24.2 (Bernoulli numbers), §4.23 (Tangent numbers)
//  - Boost.Math: bernoulli_b2n, tangent_t2n
//

import CBoostBridge

public extension SpecialFunctions {
    // MARK: - Bernoulli numbers B_{2n}
    
    /// Returns the even-indexed Bernoulli number B_{2n} converted to the generic
    /// floating-point type `T`.
    ///
    /// This is a convenience wrapper that computes `B_{2n}` in `Double` precision via
    /// the C bridge (`bs_bernoulli_b2n`) and then converts the result to `T`.
    /// It is most appropriate when you are already working in `Double` or when the
    /// small cost/precision trade-off of converting from `Double` to `Float` (or
    /// to `Float80` on supported architectures) is acceptable.
    ///
    /// Parameters:
    /// - n: The half-index n (n ≥ 0 in typical usage). This function returns `B_{2n}`.
    ///
    /// Returns:
    /// - The Bernoulli number `B_{2n}` represented as the generic floating-point type `T`.
    ///
    /// Discussion:
    /// - Precision: The computation is performed in `Double` and then converted to `T`.
    ///   For best-possible precision in `Float` or `Float80`, prefer the specialized
    ///   overloads `bernoulli_b2n_f(_:)` and `bernoulli_b2n_l(_:)` (the latter only on
    ///   x86_64/i386).
    /// - Range/overflow: Very large indices can lead to overflows or infinities
    ///   depending on `T`. Behavior mirrors the underlying Boost.Math routine.
    ///
    /// Example:
    /// - let b0: Double = SpecialFunctions.bernoulli_b2n(0)   // 1.0
    /// - let b1: Float  = SpecialFunctions.bernoulli_b2n(1)   // ≈ 0.16666667
    /// - let b2 = SpecialFunctions.bernoulli_b2n(2) as Double // -1/30
    @inlinable static func bernoulli_b2n<T: BinaryFloatingPoint>(_ n: Int32) -> T {
        return T(bs_bernoulli_b2n(Int32(n)))
    }

    /// Returns the even-indexed Bernoulli number B_{2n} as Double.
    ///
    /// Parameters:
    /// - n: The half-index n (n ≥ 0 in typical usage). This function returns B_{2n}.
    ///
    /// Returns:
    /// - The Bernoulli number B_{2n} represented as `Double`.
    ///
    /// Discussion:
    /// - Odd-index Bernoulli numbers B_{2n+1} are zero for n ≥ 1; this API targets even indices directly.
    /// - The computation is delegated to Boost.Math via the C bridge.
    ///
    /// Example:
    /// - B0 = 1 → bernoulli_b2n(0) = 1
    /// - B2 = 1/6 → bernoulli_b2n(1) ≈ 0.166666...
    @inlinable static func bernoulli_b2n(_ n: Int32) -> Double {
        return bs_bernoulli_b2n(Int32(n))
    }
    
    /// Returns the even-indexed Bernoulli number B_{2n} as Float.
    ///
    /// See `bernoulli_b2n(_:)` for semantics and discussion.
    ///
    /// Parameter:
    /// - n: The half-index n (n ≥ 0 in typical usage).
    ///
    /// Returns:
    /// - The Bernoulli number B_{2n} represented as `Float`.
    @inlinable static func bernoulli_b2n_f(_ n: Int32) -> Float {
        return bs_bernoulli_b2n_f(n)
    }
    
    #if arch(i386) || arch(x86_64)
    /// Returns the even-indexed Bernoulli number B_{2n} as Float80 (x86_64/i386 only).
    ///
    /// Parameters:
    /// - n: The half-index n (n ≥ 0 in typical usage).
    ///
    /// Returns:
    /// - The Bernoulli number B_{2n} represented as `Float80`.
    ///
    /// Notes:
    /// - Availability limited to x86_64/i386 architectures.
    @inlinable static func bernoulli_b2n_l(_ n: Int) -> Float80 {
        return bs_bernoulli_b2n_l(Int32(n))
    }
    #endif
    
    // MARK: - Tangent numbers T_{2n} (scalar)
    
    /// Returns the (even) tangent number T_{2n} converted to the generic
    /// floating-point type `T`.
    ///
    /// This is a convenience wrapper that computes `B_{2n}` in `Double` precision via
    /// the C bridge (`bs_bernoulli_b2n`) and then converts the result to `T`.
    /// It is most appropriate when you are already working in `Double` or when the
    /// small cost/precision trade-off of converting from `Double` to `Float` (or
    /// to `Float80` on supported architectures) is acceptable.
    ///
    /// Definition:
    /// - Tangent numbers T_{2n} arise in the Maclaurin expansion
    ///   tan(x) = Σ_{n≥1} T_{2n} x^{2n−1} / (2n−1)!.
    ///
    /// Parameters:
    /// - n: The half-index n (n ≥ 1 in typical usage). This function returns T_{2n}.
    ///
    /// Returns:
    /// - T_{2n} as `Double`.
    ///
    /// Complexity:
    /// - Delegated to Boost.Math; consult its documentation for complexity and range behavior.
    @inlinable static func tangent_t2n<T: BinaryFloatingPoint>(_ n: Int32) -> T {
        return T(bs_tangent_t2n(Int32(n)))
    }

    /// Returns the (even) tangent number T_{2n} as Double.
    ///
    /// Definition:
    /// - Tangent numbers T_{2n} arise in the Maclaurin expansion
    ///   tan(x) = Σ_{n≥1} T_{2n} x^{2n−1} / (2n−1)!.
    ///
    /// Parameters:
    /// - n: The half-index n (n ≥ 1 in typical usage). This function returns T_{2n}.
    ///
    /// Returns:
    /// - T_{2n} as `Double`.
    ///
    /// Complexity:
    /// - Delegated to Boost.Math; consult its documentation for complexity and range behavior.
    @inlinable static func tangent_t2n(_ n: Int32) -> Double {
        return bs_tangent_t2n(Int32(n))
    }
    
    /// Returns the (even) tangent number T_{2n} as Float.
    ///
    /// See `tangent_t2n(_:)` for semantics and discussion.
    ///
    /// Parameters:
    /// - n: The half-index n (n ≥ 1 in typical usage).
    ///
    /// Returns:
    /// - T_{2n} as `Float`.
    @inlinable static func tangent_t2n_f(_ n: Int32) -> Float {
        return bs_tangent_t2n_f(n)
    }
    
    #if arch(i386) || arch(x86_64)
    /// Returns the (even) tangent number T_{2n} as Float80 (x86_64/i386 only).
    ///
    /// Parameters:
    /// - n: The half-index n (n ≥ 1 in typical usage).
    ///
    /// Returns:
    /// - T_{2n} as `Float80`.
    ///
    /// Notes:
    /// - Availability limited to x86_64/i386 architectures.
    @inlinable static func tangent_t2n_l(_ n: Int) -> Float80 {
        return bs_tangent_t2n_l(Int32(n))
    }
    #endif
    
    // MARK: - Tangent numbers T_{2n} (bulk sequence)
    //
    // These functions compute a contiguous sequence:
    // [T_{2·startIndex}, T_{2·(startIndex+1)}, ..., count elements]
    // Defaults to startIndex = 1, i.e. [T2, T4, ..., T_{2·count}].
    
    /// Computes a sequence of tangent numbers as Double.
    ///
    /// Parameters:
    /// - count: Number of values to compute. Must be ≥ 0. If 0, returns an empty array.
    /// - startIndex: The half-index to start with (default 1), producing
    ///   T_{2·startIndex}, T_{2·(startIndex+1)}, ... for `count` elements.
    ///
    /// Returns:
    /// - Array of length `count` with consecutive tangent numbers.
    ///
    /// Example:
    /// - tangent_t2n(count: 3, startIndex: 1) → [T2, T4, T6]
    @inlinable static func tangent_t2n(count: Int, startIndex: Int32 = 1) -> [Double] {
        guard count > 0 else { return [] }
        var result = Array<Double>(repeating: .zero, count: count)
        result.withUnsafeMutableBufferPointer { buf in
            bs_tangent_t2n_seq(Int32(startIndex), UInt32(count), buf.baseAddress!)
        }
        return result
    }
    
    /// Computes a sequence of tangent numbers as Float.
    ///
    /// See `tangent_t2n(count:startIndex:)` for semantics and discussion.
    ///
    /// Parameters:
    /// - count: Number of values to compute. Must be ≥ 0.
    /// - startIndex: The half-index to start with (default 1).
    ///
    /// Returns:
    /// - Array of length `count` with consecutive tangent numbers (`Float`).
    @inlinable static func tangent_t2n_f(count: Int, startIndex: Int32 = 1) -> [Float] {
        guard count > 0 else { return [] }
        var result = Array<Float>(repeating: .zero, count: count)
        result.withUnsafeMutableBufferPointer { buf in
            bs_tangent_t2n_seq_f(Int32(startIndex), UInt32(count), buf.baseAddress!)
        }
        return result
    }
    
    #if arch(i386) || arch(x86_64)
    /// Computes a sequence of tangent numbers as Float80 (x86_64/i386 only).
    ///
    /// See `tangent_t2n(count:startIndex:)` for semantics and discussion.
    ///
    /// Parameters:
    /// - count: Number of values to compute. Must be ≥ 0.
    /// - startIndex: The half-index to start with (default 1).
    ///
    /// Returns:
    /// - Array of length `count` with consecutive tangent numbers (`Float80`).
    ///
    /// Notes:
    /// - Availability limited to x86_64/i386 architectures.
    @inlinable static func tangent_t2n_l(count: Int, startIndex: Int32 = 1) -> [Float80] {
        guard count > 0 else { return [] }
        var result = Array<Float80>(repeating: .zero, count: count)
        result.withUnsafeMutableBufferPointer { buf in
            bs_tangent_t2n_seq_l(Int32(startIndex), UInt32(count), buf.baseAddress!)
        }
        return result
    }
    #endif
    
    // MARK: - Prime numbers
    
    /// Returns the nth prime (1-based) as `UInt32`.
    ///
    /// Parameters:
    /// - n: 1-based index of the prime to obtain. For example, n=1 → 2, n=2 → 3, n=3 → 5, etc.
    ///
    /// Returns:
    /// - The nth prime as `UInt32`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterExceedsMaximumIntegerValue` if `n` exceeds the supported limit.
    ///
    /// Notes:
    /// - This wrapper forwards to `bs_prime_d` in the C bridge. A simple upper-bound check is applied
    ///   to guard the underlying implementation (currently 10000).
    ///
    /// Example:
    /// - prime(1) → 2
    /// - prime(5) → 11
    @inlinable static func prime(_ n: UInt32) throws -> UInt32 {
        guard n < 10000 else {throw SpecialFunctionError.parameterExceedsMaximumIntegerValue(name: "n", max: 10000)}
        return bs_prime(n)
    }
    
    // MARK: - Fibonacci numbers
    
    /// Returns the nth Fibonacci number (F(0)=0, F(1)=1) as `UInt64`.
    ///
    /// Parameters:
    /// - n: The index in the Fibonacci sequence.
    ///
    /// Returns:
    /// - F(n) as `UInt64`.
    ///
    /// Notes:
    /// - Delegates to Boost.Math via the C bridge (`bs_fibonacci_ull`).
    /// - Beware of overflow for large `n` as the result is bounded by `UInt64`.
    ///
    /// Example:
    /// - fibonacci(0) → 0
    /// - fibonacci(1) → 1
    /// - fibonacci(10) → 55
    @inlinable static func fibonacci(_ n: UInt64) -> UInt64 {
        return bs_fibonacci_ull(n)
    }
}
