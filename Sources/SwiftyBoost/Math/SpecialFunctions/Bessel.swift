//
//  Bessel.swift
//  Math/SpecialFunctions
//
//  Thin Swift wrappers around Boost.Math cylindrical and modified Bessel functions,
//  with argument validation and a Swifty error model.
//
//  Functions included:
//  - J_v(x): Cylindrical Bessel function of the first kind.
//  - Y_v(x): Cylindrical Bessel function of the second kind (Neumann function).
//  - I_v(x): Modified Bessel function of the first kind.
//  - K_v(x): Modified Bessel function of the second kind (Macdonald function).
//
//  Overloads are provided for generic BinaryFloatingPoint plus fast paths
//  for Float and (on x86_64) Float80. All implementations delegate to
//  CBoostBridge (Boost.Math) for numerical evaluation.
//
//  Domain notes (real-valued):
//  - J_v(x): defined for all finite v and x; if x < 0 and v is non-integer, the result is complex.
//            Boost.Math enforces this and will error for x < 0 with non-integer v.
//            We do not pre-validate the sign here to allow integer v (including large values).
//  - Y_v(x): requires x > 0; the real-valued function has a branch cut and is complex for x ≤ 0.
//  - I_v(x): defined for all finite v and x; if x < 0 and v is non-integer, the result is complex,
//            but Boost.Math reflects using parity for integer v. As with J_v, we do not pre-validate sign.
//  - K_v(x): requires x > 0; has a singularity at x = 0.
//
//  Error model:
//  - Inputs are checked for finiteness. Non-finite arguments throw SpecialFunctionError.parameterNotFinite.
//  - Functions with real-domain restrictions on x (Y_v, K_v) throw SpecialFunctionError.parameterNotPositive
//    when x ≤ 0.
//  - Further domain checks (e.g., J_v/I_v with x < 0 and non-integer v) are enforced by Boost.Math;
//    if violated, Boost’s policy raises a domain error which is surfaced as a NaN/inf by the C layer.
//    The wrappers here focus on common real-domain constraints and finiteness.
//
//  References:
//  - NIST DLMF: https://dlmf.nist.gov/10
//  - Boost.Math Bessel functions documentation:
//    https://www.boost.org/doc/libs/release/libs/math/doc/html/math_toolkit/special/bessel.html
//

import CBoostBridge

public extension SpecialFunctions {
    
    
    // MARK: - Generic BinaryFloatingPoint overloads
    
    /// Cylindrical Bessel function of the first kind J_v(x).
    ///
    /// Definition (real-valued):
    /// - J_v(x) is real for real order v and real x. For x < 0, a real result is only obtained when v is an integer
    ///   (odd/even parity applies). Otherwise the result is complex and outside the scope of this real-valued API.
    ///
    /// Behavior:
    /// - Delegates to Boost.Math cyl_bessel_j(v, x).
    /// - This wrapper validates that both v and x are finite.
    /// - No pre-check on the sign of x: Boost.Math will reflect for integer v, and raise a domain error
    ///   for non-integer v with x < 0.
    ///
    /// Parameters:
    /// - v: Order of the function (any finite real).
    /// - x: Argument (any finite real; see note above for x < 0).
    ///
    /// Returns:
    /// - J_v(x) as T.
    ///
    /// Throws:
    /// - SpecialFunctionError.parameterNotFinite(name: "v") if v is NaN or ±∞.
    /// - SpecialFunctionError.parameterNotFinite(name: "x") if x is NaN or ±∞.
    ///
    /// See also:
    /// - NIST DLMF §10.2–10.4 for properties and series/recurrence relations.
    @inlinable static func besselJ<T: BinaryFloatingPoint>(v: T, x: T) throws -> T {
        let dv = D(v), dx = D(x)
        guard dv.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "v") }
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return T(bs_cyl_bessel_j_d(dv, dx))
    }
    
    /// Cylindrical Bessel function of the second kind Y_v(x) (Neumann function).
    ///
    /// Domain:
    /// - Requires x > 0 for a real-valued result. For x ≤ 0 the value is complex or singular.
    ///   This wrapper enforces x > 0 and throws otherwise.
    ///
    /// Behavior:
    /// - Delegates to Boost.Math cyl_neumann(v, x).
    /// - Validates that v and x are finite and that x > 0.
    ///
    /// Parameters:
    /// - v: Order of the function (any finite real).
    /// - x: Argument (must satisfy x > 0).
    ///
    /// Returns:
    /// - Y_v(x) as T.
    ///
    /// Throws:
    /// - SpecialFunctionError.parameterNotFinite(name: "v") if v is NaN or ±∞.
    /// - SpecialFunctionError.parameterNotFinite(name: "x") if x is NaN or ±∞.
    /// - SpecialFunctionError.parameterNotPositive(name: "x") if x ≤ 0.
    ///
    /// Notes:
    /// - Also known as the Neumann function or Bessel function of the second kind.
    /// - Has a logarithmic singularity at x → 0⁺ for many orders.
    ///
    /// See also:
    /// - NIST DLMF §10.2–10.4 and §10.7 for behavior and relations.
    @inlinable static func besselY<T: BinaryFloatingPoint>(v: T, x: T) throws -> T {
        let dv = D(v), dx = D(x)
        guard dv.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "v") }
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        guard dx > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "x") }
        return T(bs_cyl_neumann_d(dv, dx))
    }
    
    /// Modified Bessel function of the first kind I_v(x).
    ///
    /// Definition (real-valued):
    /// - I_v(x) is real for real order v and real x. For x < 0 and non-integer v, the analytic continuation is complex;
    ///   Boost.Math applies parity for integer v. This wrapper does not pre-restrict the sign of x.
    ///
    /// Behavior:
    /// - Delegates to Boost.Math cyl_bessel_i(v, x).
    /// - Validates that v and x are finite.
    ///
    /// Parameters:
    /// - v: Order of the function (any finite real).
    /// - x: Argument (any finite real; see note above for x < 0).
    ///
    /// Returns:
    /// - I_v(x) as T.
    ///
    /// Throws:
    /// - SpecialFunctionError.parameterNotFinite(name: "v") if v is NaN or ±∞.
    /// - SpecialFunctionError.parameterNotFinite(name: "x") if x is NaN or ±∞.
    ///
    /// See also:
    /// - NIST DLMF §10.25–10.27.
    @inlinable static func modifiedBesselI<T: BinaryFloatingPoint>(v: T, x: T) throws -> T {
        let dv = D(v), dx = D(x)
        guard dv.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "v") }
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return T(bs_cyl_bessel_i_d(dv, dx))
    }
    
    /// Modified Bessel function of the second kind K_v(x) (Macdonald function).
    ///
    /// Domain:
    /// - Requires x > 0 for a real-valued result. K_v(x) has a singularity at x = 0 and is complex for x ≤ 0.
    ///
    /// Behavior:
    /// - Delegates to Boost.Math cyl_bessel_k(v, x).
    /// - Validates that v and x are finite and that x > 0.
    ///
    /// Parameters:
    /// - v: Order of the function (any finite real).
    /// - x: Argument (must satisfy x > 0).
    ///
    /// Returns:
    /// - K_v(x) as T.
    ///
    /// Throws:
    /// - SpecialFunctionError.parameterNotFinite(name: "v") if v is NaN or ±∞.
    /// - SpecialFunctionError.parameterNotFinite(name: "x") if x is NaN or ±∞.
    /// - SpecialFunctionError.parameterNotPositive(name: "x") if x ≤ 0.
    ///
    /// Notes:
    /// - Rapidly decays for large x and is often used in boundary-value problems.
    /// - Sometimes denoted as the Macdonald function.
    ///
    /// See also:
    /// - NIST DLMF §10.25–10.31.
    @inlinable static func modifiedBesselK<T: BinaryFloatingPoint>(v: T, x: T) throws -> T {
        let dv = D(v), dx = D(x)
        guard dv.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "v") }
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        guard dx > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "x") }
        return T(bs_cyl_bessel_k_d(dv, dx))
    }

    // MARK: - Mixed-precision promotions for cylindrical Bessel (Float ↔ Double)

    /// J_v(x) with mixed `Float`/`Double` arguments; returns `Double`.
    @inlinable static func besselJ(v: Float, x: Double) throws -> Double { try besselJ(v: Double(v), x: x) }
    /// J_v(x) with mixed `Double`/`Float` arguments; returns `Double`.
    @inlinable static func besselJ(v: Double, x: Float) throws -> Double { try besselJ(v: v, x: Double(x)) }

    /// Y_v(x) with mixed `Float`/`Double` arguments; returns `Double`.
    @inlinable static func besselY(v: Float, x: Double) throws -> Double { try besselY(v: Double(v), x: x) }
    /// Y_v(x) with mixed `Double`/`Float` arguments; returns `Double`.
    @inlinable static func besselY(v: Double, x: Float) throws -> Double { try besselY(v: v, x: Double(x)) }

    /// I_v(x) with mixed `Float`/`Double` arguments; returns `Double`.
    @inlinable static func modifiedBesselI(v: Float, x: Double) throws -> Double { try modifiedBesselI(v: Double(v), x: x) }
    /// I_v(x) with mixed `Double`/`Float` arguments; returns `Double`.
    @inlinable static func modifiedBesselI(v: Double, x: Float) throws -> Double { try modifiedBesselI(v: v, x: Double(x)) }

    /// K_v(x) with mixed `Float`/`Double` arguments; returns `Double`.
    @inlinable static func modifiedBesselK(v: Float, x: Double) throws -> Double { try modifiedBesselK(v: Double(v), x: x) }
    /// K_v(x) with mixed `Double`/`Float` arguments; returns `Double`.
    @inlinable static func modifiedBesselK(v: Double, x: Float) throws -> Double { try modifiedBesselK(v: v, x: Double(x)) }

    // MARK: - Mixed-precision promotions for cylindrical Bessel derivatives (Float ↔ Double)

    /// J′_v(x) with mixed `Float`/`Double` arguments; returns `Double`.
    @inlinable static func besselJPrime(v: Float, x: Double) throws -> Double { try besselJPrime(v: Double(v), x: x) }
    /// J′_v(x) with mixed `Double`/`Float` arguments; returns `Double`.
    @inlinable static func besselJPrime(v: Double, x: Float) throws -> Double { try besselJPrime(v: v, x: Double(x)) }

    /// I′_v(x) with mixed `Float`/`Double` arguments; returns `Double`.
    @inlinable static func modifiedBesselIPrime(v: Float, x: Double) throws -> Double { try modifiedBesselIPrime(v: Double(v), x: x) }
    /// I′_v(x) with mixed `Double`/`Float` arguments; returns `Double`.
    @inlinable static func modifiedBesselIPrime(v: Double, x: Float) throws -> Double { try modifiedBesselIPrime(v: v, x: Double(x)) }

    /// K′_v(x) with mixed `Float`/`Double` arguments; returns `Double`.
    @inlinable static func modifiedBesselKPrime(v: Float, x: Double) throws -> Double { try modifiedBesselKPrime(v: Double(v), x: x) }
    /// K′_v(x) with mixed `Double`/`Float` arguments; returns `Double`.
    @inlinable static func modifiedBesselKPrime(v: Double, x: Float) throws -> Double { try modifiedBesselKPrime(v: v, x: Double(x)) }

    // MARK: - Integer order convenience overloads (result type follows x)

    /// J_v(x) where `v` is integer; result type matches `x`.
    @inlinable static func besselJ<V: BinaryInteger>(v: V, x: Double) throws -> Double { try besselJ(v: Double(v), x: x) }
    @inlinable static func besselJ<V: BinaryInteger>(v: V, x: Float) throws -> Float { try besselJ_f(v: Float(v), x: x) }
    #if arch(x86_64)
    @inlinable static func besselJ<V: BinaryInteger>(v: V, x: Float80) throws -> Float80 { try besselJ_l(v: Float80(v), x: x) }
    #endif

    /// Y_v(x) where `v` is integer; result type matches `x`.
    @inlinable static func besselY<V: BinaryInteger>(v: V, x: Double) throws -> Double { try besselY(v: Double(v), x: x) }
    @inlinable static func besselY<V: BinaryInteger>(v: V, x: Float) throws -> Float { try besselY_f(v: Float(v), x: x) }
    #if arch(x86_64)
    @inlinable static func besselY<V: BinaryInteger>(v: V, x: Float80) throws -> Float80 { try besselY_l(v: Float80(v), x: x) }
    #endif

    /// I_v(x) where `v` is integer; result type matches `x`.
    @inlinable static func modifiedBesselI<V: BinaryInteger>(v: V, x: Double) throws -> Double { try modifiedBesselI(v: Double(v), x: x) }
    @inlinable static func modifiedBesselI<V: BinaryInteger>(v: V, x: Float) throws -> Float { try modifiedBesselI_f(v: Float(v), x: x) }
    #if arch(x86_64)
    @inlinable static func modifiedBesselI<V: BinaryInteger>(v: V, x: Float80) throws -> Float80 { try modifiedBesselI_l(v: Float80(v), x: x) }
    #endif

    /// K_v(x) where `v` is integer; result type matches `x`.
    @inlinable static func modifiedBesselK<V: BinaryInteger>(v: V, x: Double) throws -> Double { try modifiedBesselK(v: Double(v), x: x) }
    @inlinable static func modifiedBesselK<V: BinaryInteger>(v: V, x: Float) throws -> Float { try modifiedBesselK_f(v: Float(v), x: x) }
    #if arch(x86_64)
    @inlinable static func modifiedBesselK<V: BinaryInteger>(v: V, x: Float80) throws -> Float80 { try modifiedBesselK_l(v: Float80(v), x: x) }
    #endif

    /// J′_v(x) where `v` is integer; result type matches `x`.
    @inlinable static func besselJPrime<V: BinaryInteger>(v: V, x: Double) throws -> Double { try besselJPrime(v: Double(v), x: x) }
    @inlinable static func besselJPrime<V: BinaryInteger>(v: V, x: Float) throws -> Float { try besselJPrime(v: Float(v), x: x) }
    #if arch(x86_64)
    @inlinable static func besselJPrime<V: BinaryInteger>(v: V, x: Float80) throws -> Float80 { try besselJPrime(v: Float80(v), x: x) }
    #endif

    /// I′_v(x) where `v` is integer; result type matches `x`.
    @inlinable static func modifiedBesselIPrime<V: BinaryInteger>(v: V, x: Double) throws -> Double { try modifiedBesselIPrime(v: Double(v), x: x) }
    @inlinable static func modifiedBesselIPrime<V: BinaryInteger>(v: V, x: Float) throws -> Float { try modifiedBesselIPrime(v: Float(v), x: x) }
    #if arch(x86_64)
    @inlinable static func modifiedBesselIPrime<V: BinaryInteger>(v: V, x: Float80) throws -> Float80 { try modifiedBesselIPrime(v: Float80(v), x: x) }
    #endif

    /// K′_v(x) where `v` is integer; result type matches `x`.
    @inlinable static func modifiedBesselKPrime<V: BinaryInteger>(v: V, x: Double) throws -> Double { try modifiedBesselKPrime(v: Double(v), x: x) }
    @inlinable static func modifiedBesselKPrime<V: BinaryInteger>(v: V, x: Float) throws -> Float { try modifiedBesselKPrime(v: Float(v), x: x) }
    #if arch(x86_64)
    @inlinable static func modifiedBesselKPrime<V: BinaryInteger>(v: V, x: Float80) throws -> Float80 { try modifiedBesselKPrime(v: Float80(v), x: x) }
    #endif
    
    @inlinable static func besselJZero<T: BinaryFloatingPoint>(v: T, m: Int32) -> T {
        return T(bs_cyl_bessel_j_zero_d(Double(v), m))
    }
    
    /// Compute multiple zeros of J_v(x), starting at a given index, using the C bridge.
    ///
    /// Parameters:
    /// - v: Order of the Bessel function (precision is driven by `v`’s type).
    /// - start: 1-based start index m for the first zero (as in Boost). Some orders allow m = 0 per Boost rules.
    /// - count: Number of consecutive zeros to compute (must be ≥ 0 and fit into UInt32).
    ///
    /// Returns:
    /// - An array of `count` zeros: [j_{v, start}, j_{v, start+1}, ...].
    ///
    /// Notes:
    /// - This uses the plural bridge function to fill results in one call.
    /// - Domain specifics (e.g., invalid start for certain orders) are handled by the backend.
    @inlinable static func besselJZeros<T: BinaryFloatingPoint>(v: T, start: Int32, count: Int) -> [T] {
        precondition(count >= 0, "count must be non-negative")
        if count == 0 { return [] }
        precondition(count <= Int(UInt32.max), "count exceeds UInt32 capacity")
        
        let dv = Double(v)
        var buf = Array<Double>(repeating: 0, count: count)
        buf.withUnsafeMutableBufferPointer { bp in
            // Bridge call: fills bp.baseAddress with `count` zeros starting at `start`.
            bs_cyl_bessel_j_zeros_d(dv, start, UInt32(count), bp.baseAddress)
        }
        // Map the Double results back to T.
        return buf.map { T($0) }
    }
    
    // MARK: - Spherical Bessel functions (generic)
    
    /// Spherical Bessel function of the first kind j_n(x).
    ///
    /// Behavior:
    /// - Delegates to Boost.Math sph_bessel(n, x) via the C bridge.
    /// - Validates that n ≥ 0, n ≤ UInt32.max, and x is finite.
    ///
    /// Domain:
    /// - Real-valued for all real x.
    ///
    /// - Parameters:
    ///   - n: Non-negative integer order.
    ///   - x: Real argument (finite).
    ///
    /// - Returns: j_n(x) as T.
    @inlinable static func sphericalBesselJ<T: BinaryFloatingPoint>(n: Int, x: T) throws -> T {
        let dx = D(x)
        guard n >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "n") }
        guard n <= Int(UInt32.max) else { throw SpecialFunctionError.parameterExceedsMaximumIntegerValue(name: "n", max: Int(UInt32.max)) }
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return T(bs_sph_bessel_d(UInt32(n), dx))
    }
    
    /// Spherical Bessel function of the second kind y_n(x) (spherical Neumann).
    ///
    /// Behavior:
    /// - Delegates to Boost.Math sph_neumann(n, x) via the C bridge.
    /// - Validates that n ≥ 0, n ≤ UInt32.max, x is finite, and x > 0 for the real-valued branch.
    ///
    /// Domain:
    /// - Requires x > 0 for real-valued result; singular at x = 0 and complex for x ≤ 0.
    ///
    /// - Parameters:
    ///   - n: Non-negative integer order.
    ///   - x: Real argument (must satisfy x > 0).
    ///
    /// - Returns: y_n(x) as T.
    @inlinable static func sphericalBesselY<T: BinaryFloatingPoint>(n: Int, x: T) throws -> T {
        let dx = D(x)
        guard n >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "n") }
        guard n <= Int(UInt32.max) else { throw SpecialFunctionError.parameterExceedsMaximumIntegerValue(name: "n", max: Int(UInt32.max)) }
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        guard dx > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "x") }
        return T(bs_sph_neumann_d(UInt32(n), dx))
    }
    
    // MARK: - Bessel derivatives (generic)
    
    /// Derivative with respect to x of the cylindrical Bessel function J_v(x).
    ///
    /// Definition:
    /// - J′_v(x) = d/dx J_v(x).
    ///
    /// Behavior:
    /// - Delegates to Boost.Math cyl_bessel_j_prime(v, x).
    /// - Validates that v and x are finite.
    /// - No pre-check on the sign of x: for x < 0 and non-integer v, the analytic continuation is complex;
    ///   Boost handles integer-order reflection and raises a domain error otherwise.
    ///
    /// Parameters:
    /// - v: Order of the function (any finite real).
    /// - x: Argument (any finite real; see note above for x < 0).
    ///
    /// Returns:
    /// - J′_v(x) as T.
    ///
    /// Throws:
    /// - SpecialFunctionError.parameterNotFinite(name: "v") if v is NaN or ±∞.
    /// - SpecialFunctionError.parameterNotFinite(name: "x") if x is NaN or ±∞.
    @inlinable static func besselJPrime<T: BinaryFloatingPoint>(v: T, x: T) throws -> T {
        let dv = D(v), dx = D(x)
        guard dv.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "v") }
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return T(bs_cyl_bessel_j_prime_d(dv, dx))
    }
    
    /// Derivative with respect to x of the modified Bessel function I_v(x).
    ///
    /// Definition:
    /// - I′_v(x) = d/dx I_v(x).
    ///
    /// Behavior:
    /// - Delegates to Boost.Math cyl_bessel_i_prime(v, x).
    /// - Validates that v and x are finite.
    ///
    /// Parameters:
    /// - v: Order of the function (any finite real).
    /// - x: Argument (any finite real; for x < 0 and non-integer v the analytic continuation is complex).
    ///
    /// Returns:
    /// - I′_v(x) as T.
    ///
    /// Throws:
    /// - SpecialFunctionError.parameterNotFinite(name: "v") if v is NaN or ±∞.
    /// - SpecialFunctionError.parameterNotFinite(name: "x") if x is NaN or ±∞.
    @inlinable static func modifiedBesselIPrime<T: BinaryFloatingPoint>(v: T, x: T) throws -> T {
        let dv = D(v), dx = D(x)
        guard dv.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "v") }
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return T(bs_cyl_bessel_i_prime_d(dv, dx))
    }
    
    /// Derivative with respect to x of the modified Bessel function K_v(x).
    ///
    /// Definition:
    /// - K′_v(x) = d/dx K_v(x).
    ///
    /// Domain:
    /// - Requires x > 0 for a real-valued result (singular at x = 0).
    ///
    /// Behavior:
    /// - Delegates to Boost.Math cyl_bessel_k_prime(v, x).
    /// - Validates that v and x are finite and that x > 0.
    ///
    /// Parameters:
    /// - v: Order of the function (any finite real).
    /// - x: Argument (must satisfy x > 0).
    ///
    /// Returns:
    /// - K′_v(x) as T.
    ///
    /// Throws:
    /// - SpecialFunctionError.parameterNotFinite(name: "v") if v is NaN or ±∞.
    /// - SpecialFunctionError.parameterNotFinite(name: "x") if x is NaN or ±∞.
    /// - SpecialFunctionError.parameterNotPositive(name: "x") if x ≤ 0.
    @inlinable static func modifiedBesselKPrime<T: BinaryFloatingPoint>(v: T, x: T) throws -> T {
        let dv = D(v), dx = D(x)
        guard dv.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "v") }
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        guard dx > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "x") }
        return T(bs_cyl_bessel_k_prime_d(dv, dx))
    }
    
    /// Derivative with respect to x of the spherical Bessel function j_n(x).
    ///
    /// Definition:
    /// - j′_n(x) = d/dx j_n(x).
    ///
    /// Domain:
    /// - Real-valued branch requires x ≥ 0. This wrapper enforces x ≥ 0.
    ///
    /// Behavior:
    /// - Delegates to Boost.Math sph_bessel_prime(n, x).
    /// - Validates that n ≥ 0, n ≤ UInt32.max, x is finite, and x ≥ 0.
    ///
    /// Parameters:
    /// - n: Non-negative integer order.
    /// - x: Argument (must satisfy x ≥ 0).
    ///
    /// Returns:
    /// - j′_n(x) as T.
    ///
    /// Throws:
    /// - SpecialFunctionError.parameterNotPositive(name: "n") if n < 0.
    /// - SpecialFunctionError.parameterExceedsMaximumIntegerValue(name: "n", max: ...) if n > UInt32.max.
    /// - SpecialFunctionError.parameterNotFinite(name: "x") if x is NaN or ±∞.
    /// - SpecialFunctionError.parameterOutOfRange(name: "x", min: 0, max: +∞) if x < 0.
    @inlinable static func sphericalBesselJPrime<T: BinaryFloatingPoint>(n: Int, x: T) throws -> T {
        let dx = D(x)
        guard n >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "n") }
        guard n <= Int(UInt32.max) else { throw SpecialFunctionError.parameterExceedsMaximumIntegerValue(name: "n", max: Int(UInt32.max)) }
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        guard dx >= 0 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 0.0, max: Double.infinity) }
        return T(bs_sph_bessel_prime_d(UInt32(n), dx))
    }
    
    /// Derivative with respect to x of the spherical Neumann function y_n(x).
    ///
    /// Definition:
    /// - y′_n(x) = d/dx y_n(x).
    ///
    /// Domain:
    /// - Requires x > 0 for a real-valued result (singular at x = 0).
    ///
    /// Behavior:
    /// - Delegates to Boost.Math sph_neumann_prime(n, x).
    /// - Validates that n ≥ 0, n ≤ UInt32.max, x is finite, and x > 0.
    ///
    /// Parameters:
    /// - n: Non-negative integer order.
    /// - x: Argument (must satisfy x > 0).
    ///
    /// Returns:
    /// - y′_n(x) as T.
    ///
    /// Throws:
    /// - SpecialFunctionError.parameterNotPositive(name: "n") if n < 0.
    /// - SpecialFunctionError.parameterExceedsMaximumIntegerValue(name: "n", max: ...) if n > UInt32.max.
    /// - SpecialFunctionError.parameterNotFinite(name: "x") if x is NaN or ±∞.
    /// - SpecialFunctionError.parameterNotPositive(name: "x") if x ≤ 0.
    @inlinable static func sphericalBesselYPrime<T: BinaryFloatingPoint>(n: Int, x: T) throws -> T {
        let dx = D(x)
        guard n >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "n") }
        guard n <= Int(UInt32.max) else { throw SpecialFunctionError.parameterExceedsMaximumIntegerValue(name: "n", max: Int(UInt32.max)) }
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        guard dx > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "x") }
        return T(bs_sph_neumann_prime_d(UInt32(n), dx))
    }
    
    
    // MARK: - Float overloads
    
    /// Cylindrical Bessel function of the first kind J_v(x) for Float.
    @inlinable static func besselJ_f(v: Float, x: Float) throws -> Float {
        guard v.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "v") }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return bs_cyl_bessel_j_f(v, x)
    }
    
    /// Cylindrical Bessel function of the second kind Y_v(x) for Float.
    ///
    /// Domain: requires x > 0.
    @inlinable static func besselY_f(v: Float, x: Float) throws -> Float {
        guard v.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "v") }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        guard x > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "x") }
        return bs_cyl_neumann_f(v, x)
    }
    
    /// Modified Bessel function of the first kind I_v(x) for Float.
    @inlinable static func modifiedBesselI_f(v: Float, x: Float) throws -> Float {
        guard v.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "v") }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return bs_cyl_bessel_i_f(v, x)
    }
    
    /// Modified Bessel function of the second kind K_v(x) for Float.
    ///
    /// Domain: requires x > 0.
    @inlinable static func modifiedBesselK_f(v: Float, x: Float) throws -> Float {
        guard v.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "v") }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        guard x > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "x") }
        return bs_cyl_bessel_k_f(v, x)
    }
    
    /// Compute multiple zeros of J_v(x), starting at a given index, using the C bridge - Float
    @inlinable static func besselJZeros(v: Float, start: Int32, count: Int) -> [Float] {
        precondition(count >= 0, "count must be non-negative")
        if count == 0 { return [] }
        precondition(count <= Int(UInt32.max), "count exceeds UInt32 capacity")
        
        var buf = Array<Float>(repeating: 0, count: count)
        buf.withUnsafeMutableBufferPointer { bp in
            // Bridge call: fills bp.baseAddress with `count` zeros starting at `start`.
            bs_cyl_bessel_j_zeros_f(v, start, UInt32(count), bp.baseAddress)
        }
        return buf
    }
    
    // Spherical Bessel (Float)
    
    /// Spherical Bessel j_n(x) for Float.
    @inlinable static func sphericalBesselJ(n: Int, x: Float) throws -> Float {
        guard n >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "n") }
        guard n <= Int(UInt32.max) else { throw SpecialFunctionError.parameterExceedsMaximumIntegerValue(name: "n", max: Int(UInt32.max)) }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return bs_sph_bessel_f(UInt32(n), x)
    }
    
    /// Spherical Neumann y_n(x) for Float. Domain: x > 0.
    @inlinable static func sphericalBesselY(n: Int, x: Float) throws -> Float {
        guard n >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "n") }
        guard n <= Int(UInt32.max) else { throw SpecialFunctionError.parameterExceedsMaximumIntegerValue(name: "n", max: Int(UInt32.max)) }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        guard x > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "x") }
        return bs_sph_neumann_f(UInt32(n), x)
    }
    
    // Bessel derivatives (Float)
    
    /// J′_v(x) for Float.
    ///
    /// See ``SpecialFunctions/besselJPrime(v:x:)->T`` for definition, domain, and behavior.
    @inlinable static func besselJPrime(v: Float, x: Float) throws -> Float {
        guard v.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "v") }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return bs_cyl_bessel_j_prime_f(v, x)
    }
    
    /// I′_v(x) for Float.
    ///
    /// See ``SpecialFunctions/modifiedBesselIPrime(v:x:)->T`` for definition, domain, and behavior.
    @inlinable static func modifiedBesselIPrime(v: Float, x: Float) throws -> Float {
        guard v.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "v") }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return bs_cyl_bessel_i_prime_f(v, x)
    }
    
    /// K′_v(x) for Float.
    ///
    /// Domain:
    /// - Requires x > 0 for a real-valued result (singular at x = 0).
    ///
    /// See ``SpecialFunctions/modifiedBesselKPrime(v:x:)->T`` for additional details.
    @inlinable static func modifiedBesselKPrime(v: Float, x: Float) throws -> Float {
        guard v.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "v") }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        guard x > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "x") }
        return bs_cyl_bessel_k_prime_f(v, x)
    }
    
    /// j′_n(x) for Float.
    ///
    /// Domain:
    /// - Real-valued branch requires x ≥ 0.
    ///
    /// See ``SpecialFunctions/sphericalBesselJPrime(n:x:)->T`` for additional details.
    @inlinable static func sphericalBesselJPrime(n: Int, x: Float) throws -> Float {
        guard n >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "n") }
        guard n <= Int(UInt32.max) else { throw SpecialFunctionError.parameterExceedsMaximumIntegerValue(name: "n", max: Int(UInt32.max)) }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        guard x >= 0 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 0.0, max: Double.infinity) }
        return bs_sph_bessel_prime_f(UInt32(n), x)
    }
    
    /// y′_n(x) for Float.
    ///
    /// Domain:
    /// - Requires x > 0 for a real-valued result (singular at x = 0).
    ///
    /// See ``SpecialFunctions/sphericalBesselYPrime(n:x:)->T`` for additional details.
    @inlinable static func sphericalBesselYPrime(n: Int, x: Float) throws -> Float {
        guard n >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "n") }
        guard n <= Int(UInt32.max) else { throw SpecialFunctionError.parameterExceedsMaximumIntegerValue(name: "n", max: Int(UInt32.max)) }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        guard x > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "x") }
        return bs_sph_neumann_prime_f(UInt32(n), x)
    }


    // MARK: - Float80 overloads (x86_64)
    
#if arch(x86_64)
    // MARK: - Mixed-precision promotions with Float80 (result → Float80)

    /// J_v(x) with one `Float80` argument; returns `Float80`.
    @inlinable static func besselJ(v: Float80, x: Double) throws -> Float80 { try besselJ_l(v: v, x: Float80(x)) }
    @inlinable static func besselJ(v: Double, x: Float80) throws -> Float80 { try besselJ_l(v: Float80(v), x: x) }
    @inlinable static func besselJ(v: Float80, x: Float) throws -> Float80 { try besselJ_l(v: v, x: Float80(x)) }
    @inlinable static func besselJ(v: Float, x: Float80) throws -> Float80 { try besselJ_l(v: Float80(v), x: x) }

    /// Y_v(x) with one `Float80` argument; returns `Float80`.
    @inlinable static func besselY(v: Float80, x: Double) throws -> Float80 { try besselY_l(v: v, x: Float80(x)) }
    @inlinable static func besselY(v: Double, x: Float80) throws -> Float80 { try besselY_l(v: Float80(v), x: x) }
    @inlinable static func besselY(v: Float80, x: Float) throws -> Float80 { try besselY_l(v: v, x: Float80(x)) }
    @inlinable static func besselY(v: Float, x: Float80) throws -> Float80 { try besselY_l(v: Float80(v), x: x) }

    /// I_v(x) with one `Float80` argument; returns `Float80`.
    @inlinable static func modifiedBesselI(v: Float80, x: Double) throws -> Float80 { try modifiedBesselI_l(v: v, x: Float80(x)) }
    @inlinable static func modifiedBesselI(v: Double, x: Float80) throws -> Float80 { try modifiedBesselI_l(v: Float80(v), x: x) }
    @inlinable static func modifiedBesselI(v: Float80, x: Float) throws -> Float80 { try modifiedBesselI_l(v: v, x: Float80(x)) }
    @inlinable static func modifiedBesselI(v: Float, x: Float80) throws -> Float80 { try modifiedBesselI_l(v: Float80(v), x: x) }

    /// K_v(x) with one `Float80` argument; returns `Float80`.
    @inlinable static func modifiedBesselK(v: Float80, x: Double) throws -> Float80 { try modifiedBesselK_l(v: v, x: Float80(x)) }
    @inlinable static func modifiedBesselK(v: Double, x: Float80) throws -> Float80 { try modifiedBesselK_l(v: Float80(v), x: x) }
    @inlinable static func modifiedBesselK(v: Float80, x: Float) throws -> Float80 { try modifiedBesselK_l(v: v, x: Float80(x)) }
    @inlinable static func modifiedBesselK(v: Float, x: Float80) throws -> Float80 { try modifiedBesselK_l(v: Float80(v), x: x) }

    /// J′_v(x) with one `Float80` argument; returns `Float80`.
    @inlinable static func besselJPrime(v: Float80, x: Double) throws -> Float80 { try besselJPrime(v: v, x: Float80(x)) }
    @inlinable static func besselJPrime(v: Double, x: Float80) throws -> Float80 { try besselJPrime(v: Float80(v), x: x) }
    @inlinable static func besselJPrime(v: Float80, x: Float) throws -> Float80 { try besselJPrime(v: v, x: Float80(x)) }
    @inlinable static func besselJPrime(v: Float, x: Float80) throws -> Float80 { try besselJPrime(v: Float80(v), x: x) }

    /// I′_v(x) with one `Float80` argument; returns `Float80`.
    @inlinable static func modifiedBesselIPrime(v: Float80, x: Double) throws -> Float80 { try modifiedBesselIPrime(v: v, x: Float80(x)) }
    @inlinable static func modifiedBesselIPrime(v: Double, x: Float80) throws -> Float80 { try modifiedBesselIPrime(v: Float80(v), x: x) }
    @inlinable static func modifiedBesselIPrime(v: Float80, x: Float) throws -> Float80 { try modifiedBesselIPrime(v: v, x: Float80(x)) }
    @inlinable static func modifiedBesselIPrime(v: Float, x: Float80) throws -> Float80 { try modifiedBesselIPrime(v: Float80(v), x: x) }

    /// K′_v(x) with one `Float80` argument; returns `Float80`.
    @inlinable static func modifiedBesselKPrime(v: Float80, x: Double) throws -> Float80 { try modifiedBesselKPrime(v: v, x: Float80(x)) }
    @inlinable static func modifiedBesselKPrime(v: Double, x: Float80) throws -> Float80 { try modifiedBesselKPrime(v: Float80(v), x: x) }
    @inlinable static func modifiedBesselKPrime(v: Float80, x: Float) throws -> Float80 { try modifiedBesselKPrime(v: v, x: Float80(x)) }
    @inlinable static func modifiedBesselKPrime(v: Float, x: Float80) throws -> Float80 { try modifiedBesselKPrime(v: Float80(v), x: x) }

    /// Cylindrical Bessel function of the first kind J_v(x) for Float80 (x86_64 only).
    @inlinable static func besselJ_l(v: Float80, x: Float80) throws -> Float80 {
        guard v.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "v") }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return bs_cyl_bessel_j_l(v, x)
    }
    
    /// Cylindrical Bessel function of the second kind Y_v(x) for Float80 (x86_64 only).
    ///
    /// Domain: requires x > 0.
    @inlinable static func besselY_l(v: Float80, x: Float80) throws -> Float80 {
        guard v.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "v") }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        guard x > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "x") }
        return bs_cyl_neumann_l(v, x)
    }
    
    /// Modified Bessel function of the first kind I_v(x) for Float80 (x86_64 only).
    @inlinable static func modifiedBesselI_l(v: Float80, x: Float80) throws -> Float80 {
        guard v.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "v") }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return bs_cyl_bessel_i_l(v, x)
    }
    
    /// Modified Bessel function of the second kind K_v(x) for Float80 (x86_64 only).
    ///
    /// Domain: requires x > 0.
    @inlinable static func modifiedBesselK_l(v: Float80, x: Float80) throws -> Float80 {
        guard v.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "v") }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        guard x > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "x") }
        return bs_cyl_bessel_k_l(v, x)
    }

    /// Compute multiple zeros of J_v(x), starting at a given index, using the C bridge - Float
    @inlinable static func besselJZeros(v: Float80, start: Int32, count: Int) -> [Float80] {
        precondition(count >= 0, "count must be non-negative")
        if count == 0 { return [] }
        precondition(count <= Int(UInt32.max), "count exceeds UInt32 capacity")
        
        var buf = Array<Float80>(repeating: 0, count: count)
        buf.withUnsafeMutableBufferPointer { bp in
            // Bridge call: fills bp.baseAddress with `count` zeros starting at `start`.
            bs_cyl_bessel_j_zeros_l(v, start, UInt32(count), bp.baseAddress)
        }
        return buf
    }

    // Spherical Bessel (Float80)
    
    /// Spherical Bessel j_n(x) for Float80 (x86_64 only).
    @inlinable static func sphericalBesselJ(n: Int, x: Float80) throws -> Float80 {
        guard n >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "n") }
        guard n <= Int(UInt32.max) else { throw SpecialFunctionError.parameterExceedsMaximumIntegerValue(name: "n", max: Int(UInt32.max)) }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return bs_sph_bessel_l(UInt32(n), x)
    }
    
    /// Spherical Neumann y_n(x) for Float80 (x86_64 only). Domain: x > 0.
    @inlinable static func sphericalBesselY(n: Int, x: Float80) throws -> Float80 {
        guard n >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "n") }
        guard n <= Int(UInt32.max) else { throw SpecialFunctionError.parameterExceedsMaximumIntegerValue(name: "n", max: Int(UInt32.max)) }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        guard x > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "x") }
        return bs_sph_neumann_l(UInt32(n), x)
    }

    // Bessel derivatives (Float80)
    
    /// J′_v(x) for Float80 (x86_64 only).
    ///
    /// See ``SpecialFunctions/besselJPrime(v:x:)->T`` for definition, domain, and behavior.
    @inlinable static func besselJPrime(v: Float80, x: Float80) throws -> Float80 {
        guard v.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "v") }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return bs_cyl_bessel_j_prime_l(v, x)
    }
    
    /// I′_v(x) for Float80 (x86_64 only).
    ///
    /// See ``SpecialFunctions/modifiedBesselIPrime(v:x:)->T`` for definition, domain, and behavior.
    @inlinable static func modifiedBesselIPrime(v: Float80, x: Float80) throws -> Float80 {
        guard v.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "v") }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return bs_cyl_bessel_i_prime_l(v, x)
    }
    
    /// K′_v(x) for Float80 (x86_64 only).
    ///
    /// Domain:
    /// - Requires x > 0 for a real-valued result (singular at x = 0).
    ///
    /// See ``SpecialFunctions/modifiedBesselKPrime(v:x:)->T`` for additional details.
    @inlinable static func modifiedBesselKPrime(v: Float80, x: Float80) throws -> Float80 {
        guard v.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "v") }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        guard x > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "x") }
        return bs_cyl_bessel_k_prime_l(v, x)
    }
    
    /// j′_n(x) for Float80 (x86_64 only).
    ///
    /// Domain:
    /// - Real-valued branch requires x ≥ 0.
    ///
    /// See ``SpecialFunctions/sphericalBesselJPrime(n:x:)->T`` for additional details.
    @inlinable static func sphericalBesselJPrime(n: Int, x: Float80) throws -> Float80 {
        guard n >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "n") }
        guard n <= Int(UInt32.max) else { throw SpecialFunctionError.parameterExceedsMaximumIntegerValue(name: "n", max: Int(UInt32.max)) }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        guard x >= 0 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 0.0, max: Double.infinity) }
        return bs_sph_bessel_prime_l(UInt32(n), x)
    }
    
    /// y′_n(x) for Float80 (x86_64 only).
    ///
    /// Domain:
    /// - Requires x > 0 for a real-valued result (singular at x = 0).
    ///
    /// See ``SpecialFunctions/sphericalBesselYPrime(n:x:)->T`` for additional details.
    @inlinable static func sphericalBesselYPrime(n: Int, x: Float80) throws -> Float80 {
        guard n >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "n") }
        guard n <= Int(UInt32.max) else { throw SpecialFunctionError.parameterExceedsMaximumIntegerValue(name: "n", max: Int(UInt32.max)) }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        guard x > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "x") }
        return bs_sph_neumann_prime_l(UInt32(n), x)
    }

#endif
    
}
