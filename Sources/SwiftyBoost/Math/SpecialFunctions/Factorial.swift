//
//  Factorial.swift
//  Math/SpecialFunctions
//
//  Swift wrappers for factorials and combinatorics-related special functions.
//  These APIs validate inputs (range, finiteness) and surface numeric issues
//  via SpecialFunctionError, while delegating the numerical work to the
//  CBoostBridge (Boost.Math) backend.
//
//  Functions included (real-valued):
//  - Factorial n!:
//    • factorial(_:) as Double, factorial_f(_:) as Float, factorial_l(_:) as Float80 (x86 only)
//    • factorial<T>(_:) generic overload that honors target type limits
//  - Rising factorial / Pochhammer (x)_n = x (x+1) … (x+n−1):
//    • rising_factorial(_:_:), rising_factorial_f(_:_:), rising_factorial_l(_:_:)
//    • pochhammer aliases for the same
//  - Binomial coefficient C(n, k):
//    • binomial_coeff(_:_:), binomial_coeff_f(_:_:), binomial_coeff_l(_:_:)
//  - Double factorial n!!:
//    • double_factorial(_:) as Double, double_factorial_f(_:) as Float, double_factorial_l(_:) as Float80 (x86 only)
//    • double_factorial<T>(_:) generic overload
//
//  Safety and numeric behavior:
//  - Factorial: These wrappers enforce the largest admissible n for each IEEE type
//    and throw parameterExceedsMaximumIntegerValue when exceeded. The backend returns
//    +∞ beyond those maxima; throwing avoids silent overflow.
//  - Rising factorial: Performs finiteness checks on x, detects exact zeros when
//    x is a non-positive integer within the multiplicative range, and uses ln Γ
//    differences to pre-check overflow/underflow for the target type. On underflow
//    (magnitude below least nonzero), returns 0. Throws invalidCombination on overflow.
//  - Binomial coefficient: Uses ln Γ to pre-check overflow for the target type. Returns
//    0 when k > n, 1 when k == 0 or k == n. Throws invalidCombination on overflow.
//  - Double factorial: Uses a closed-form ln(n!!) for magnitude checks before calling
//    the backend. Throws invalidCombination on overflow in the target type.
//
//  Architecture notes:
//  - Float80 variants are provided only on i386/x86_64 architectures.
//  - Generic overloads funnel computation through Double-backed implementations and
//    convert the result to the requested T, but their pre-check limits and bounds
//    are chosen according to T to respect the caller’s target type.
//
//  References:
//  - NIST DLMF §5 (Gamma function and related), §26 (Combinatorial functions)
//  - Boost.Math documentation
//

import CBoostBridge
#if canImport(Darwin)
import Darwin
#else
#if canImport(Glibc)
import Glibc
#endif
#endif


public extension SpecialFunctions {
    
    // MARK: - Maximum admissible n for n! per floating-point type
    
    /// Maximum n such that n! is finite in IEEE-754 single precision (Float).
    ///
    /// Discussion:
    /// - For n > 34, Float factorial overflows; the backend would return +∞.
    ///   The Float and generic wrappers enforce this limit and throw instead.
    @inlinable static var maxFactorialInputFloat: UInt32 { 34 }
    
    /// Maximum n such that n! is finite in IEEE-754 double precision (Double).
    ///
    /// Discussion:
    /// - For n > 170, Double factorial overflows; the backend would return +∞.
    ///   The Double and generic wrappers enforce this limit and throw instead.
    @inlinable static var maxFactorialInputDouble: UInt32 { 170 }
#if arch(i386) || arch(x86_64)
    /// Maximum n such that n! is finite in 80-bit extended precision (Float80, x86 only).
    ///
    /// Discussion:
    /// - For n > 1754, Float80 factorial overflows; the backend would return +∞.
    ///   The Float80 wrapper enforces this limit and throws instead.
    @inlinable static var maxFactorialInputFloat80: UInt32 { 1754 }
#endif
    
    /// Returns the appropriate factorial limit for the requested floating-point type T.
    ///
    /// Behavior:
    /// - For generic calls that evaluate via Double internally, this limit is still
    ///   selected based on `T` to honor the caller’s target type. The actual computation
    ///   is performed by the Double-backed backend and converted to `T`.
    ///
    /// Parameters:
    /// - for: The concrete floating-point type.
    ///
    /// Returns:
    /// - The largest n such that n! is finite in type `T`.
    @inlinable static func maxFactorialInput<T: BinaryFloatingPoint>(for _: T.Type) -> UInt32 {
        if T.self == Float.self { return maxFactorialInputFloat }
#if arch(i386) || arch(x86_64)
        if T.self == Float80.self { return maxFactorialInputFloat80 }
#endif
        return maxFactorialInputDouble
    }
    
    // MARK: - Factorial n! (throwing on overflow)
    
    /// Compute the factorial n! returned as `T`. Throws if `n` exceeds the representable limit for `T`.
    ///
    /// Limits:
    /// - Float: 34
    /// - Double: 170
    /// - Float80 (x86 only): 1754
    ///
    /// Parameters:
    /// - n: Non-negative integer for which to compute n!.
    ///
    /// Returns:
    /// - n! as `T`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterExceedsMaximumIntegerValue(name: "n", max: …)` if `n` exceeds the limit for `T`.
    ///
    /// Notes:
    /// - The computation is delegated to a Double-backed backend and then converted to `T`.
    @inlinable static func factorial<T: BinaryFloatingPoint> (_ n: UInt32) throws -> T {
        let limit = maxFactorialInput(for: T.self)
        guard n <= limit else {
            throw SpecialFunctionError.parameterExceedsMaximumIntegerValue(name: "n", max: Int(limit))
        }
        // Evaluate via Double-backed backend and cast to T.
        return T(bs_factorial(n))
    }
    
    /// Compute the factorial n! as `Double`. Throws if `n` > 170.
    ///
    /// Parameters:
    /// - n: Non-negative integer for which to compute n!.
    ///
    /// Returns:
    /// - n! as `Double`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterExceedsMaximumIntegerValue(name: "n", max: 170)` if `n` > 170.
    @inlinable static func factorial (_ n: UInt32) throws -> Double {
        let limit = maxFactorialInputDouble
        guard n <= limit else {
            throw SpecialFunctionError.parameterExceedsMaximumIntegerValue(name: "n", max: Int(limit))
        }
        return bs_factorial(n)
    }
    
    /// Compute the factorial n! as `Float`. Throws if `n` > 34.
    ///
    /// Parameters:
    /// - n: Non-negative integer for which to compute n!.
    ///
    /// Returns:
    /// - n! as `Float`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterExceedsMaximumIntegerValue(name: "n", max: 34)` if `n` > 34.
    @inlinable static func factorial_f (_ n: UInt32) throws -> Float {
        let limit = maxFactorialInputFloat
        guard n <= limit else {
            throw SpecialFunctionError.parameterExceedsMaximumIntegerValue(name: "n", max: Int(limit))
        }
        return bs_factorial_f(n)
    }
    
#if arch(i386) || arch(x86_64)
    /// Compute the factorial n! as `Float80` (x86 only). Throws if `n` > 1754.
    ///
    /// Parameters:
    /// - n: Non-negative integer for which to compute n!.
    ///
    /// Returns:
    /// - n! as `Float80`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterExceedsMaximumIntegerValue(name: "n", max: 1754)` if `n` > 1754.
    @inlinable static func factorial_l (_ n: UInt32) throws -> Float80 {
        let limit = maxFactorialInputFloat80
        guard n <= limit else {
            throw SpecialFunctionError.parameterExceedsMaximumIntegerValue(name: "n", max: Int(limit))
        }
        return bs_factorial_l(n)
    }
#endif
    
    // MARK: - Internal helpers for rising factorial bounds
    
    /// Natural log of the greatest finite magnitude for each floating-point type, expressed as `Double`.
    ///
    /// Purpose:
    /// - Used to pre-check overflow for rising factorial and related computations by
    ///   comparing ln-magnitudes against ln(max finite).
    @inline(__always)  static func logMaxFinite<T: BinaryFloatingPoint>(for _: T.Type) -> Double {
        if T.self == Float.self {
            return 88.72283905206835 // log(Float.greatestFiniteMagnitude)
        }
#if arch(i386) || arch(x86_64)
        if T.self == Float80.self {
            return 11356.523406294143 // log(1.189731495357231765e4932) for 80-bit extended
        }
#endif
        return 709.782712893384 // log(Double.greatestFiniteMagnitude)
    }
    
    /// Natural log of the least positive nonzero magnitude for each floating-point type, expressed as `Double`.
    ///
    /// Purpose:
    /// - Used as an optional underflow-to-zero policy for rising factorial. If the
    ///   predicted ln-magnitude is below this value, the function returns 0.
    @inline(__always)  static func logLeastPosNonzero<T: BinaryFloatingPoint>(for _: T.Type) -> Double {
        if T.self == Float.self {
            return -103.27892990343185 // log(Float.leastNonzeroMagnitude)
        }
#if arch(i386) || arch(x86_64)
        if T.self == Float80.self {
            // Use a conservative value; underflow pre-check for Float80 can be skipped if undesired.
            return -11399.0
        }
#endif
        return -744.4400719213812 // log(Double.leastNonzeroMagnitude)
    }
    
    /// Exact-zero check for rising factorial: if there exists an integer k in [0, n−1] with x + k == 0, then (x)_n == 0.
    ///
    /// Notes:
    /// - Implemented in Double for robust integer detection, even when T is Float/Float80.
    /// - Short-circuits the computation and returns zero when the multiplicative sequence
    ///   contains a factor equal to zero.
    @inline(__always) static func risingFactorialIsExactlyZero(_ x: Double, n: UInt32) -> Bool {
        guard n > 0 else { return false }
        if x <= 0 {
            let negx = -x
            let k = negx.rounded(.towardZero)
            if negx == k {
                // k is a nonnegative integer; zero occurs if k ∈ [0, n-1]
                return k >= 0 && k < Double(n)
            }
        }
        return false
    }
    
    // MARK: - Rising factorial / Pochhammer (x)_{n} with overflow/underflow checks
    
    /// Compute the rising factorial (x)_n = x (x+1) … (x+n−1) returned as `T`, with magnitude checks.
    ///
    /// Behavior:
    /// - Validates finiteness of x and throws if not finite.
    /// - Detects exact zeros when x is a non-positive integer within the multiplicative range.
    /// - Uses S = ln Γ(x+n) − ln Γ(x) to pre-check overflow/underflow for the target type.
    ///   Throws invalidCombination on overflow; returns 0 on underflow.
    ///
    /// Parameters:
    /// - x: The base value (finite).
    /// - n: The non-negative integer order (n ≥ 0).
    ///
    /// Returns:
    /// - (x)_n as `T`. If n == 0, returns 1 for all finite x.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotFinite(name: "x")` if x is NaN or ±∞.
    /// - `SpecialFunctionError.invalidCombination(message: ...)` if the result would overflow `T`.
    @inlinable static func rising_factorial<T: BinaryFloatingPoint> (_ x: T, _ n: UInt32) throws -> T {
        let dx = D(x)
        // Finiteness check
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        
        // Exact zero short-circuit
        if risingFactorialIsExactlyZero(dx, n: n) { return T(0) }
        
        // Compute S = ln Γ(x+n) − ln Γ(x) in Double
        let lnGammaXN = bs_lgamma(dx + Double(n))
        let lnGammaX  = bs_lgamma(dx)
        let S = lnGammaXN - lnGammaX
        
        // Overflow check
        let lnMax = logMaxFinite(for: T.self)
        if S > lnMax {
            throw SpecialFunctionError.invalidCombination(message: "rising_factorial overflows for given x and n in the target type")
        }
        
        // Optional underflow-to-zero
        let lnMin = logLeastPosNonzero(for: T.self)
        if S < lnMin {
            return T(0)
        }
        
        return T(bs_rising_factorial(dx, n))
    }
    
    /// Rising factorial (x)_n for `Double` with magnitude checks.
    ///
    /// Behavior:
    /// - Throws on overflow based on ln-magnitude check.
    /// - Returns 0 on exact zero (factor equals zero in the product) or underflow.
    ///
    /// Parameters:
    /// - x: Finite base value.
    /// - n: Non-negative integer order.
    ///
    /// Returns:
    /// - (x)_n as `Double`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotFinite(name: "x")` if x is NaN or ±∞.
    /// - `SpecialFunctionError.invalidCombination(message: ...)` if the result would overflow Double.
    @inlinable static func rising_factorial (_ x: Double,_ n: UInt32) throws -> Double {
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        if risingFactorialIsExactlyZero(x, n: n) { return 0 }
        
        let S = bs_lgamma(x + Double(n)) - bs_lgamma(x)
        if S > 709.782712893384 { // log(Double.greatestFiniteMagnitude)
            throw SpecialFunctionError.invalidCombination(message: "rising_factorial overflows for given x and n in Double")
        }
        if S < -744.4400719213812 { // log(Double.leastNonzeroMagnitude)
            return 0
        }
        
        return bs_rising_factorial(x, n)
    }
    
    /// Rising factorial (x)_n for `Float` with magnitude checks.
    ///
    /// Behavior:
    /// - Throws on overflow based on ln-magnitude check.
    /// - Returns 0 on exact zero or underflow.
    ///
    /// Parameters:
    /// - x: Finite base value.
    /// - n: Non-negative integer order.
    ///
    /// Returns:
    /// - (x)_n as `Float`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotFinite(name: "x")` if x is NaN or ±∞.
    /// - `SpecialFunctionError.invalidCombination(message: ...)` if the result would overflow Float.
    @inlinable static func rising_factorial_f (_ x: Float, _ n: UInt32) throws -> Float {
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        // Do exact-zero check in Float via Double for integer detection robustness
        let dx = Double(x)
        if risingFactorialIsExactlyZero(dx, n: n) { return 0 }
        
        // Compute in Float precision for log-gamma
        let S = Double(bs_lgamma_f(x + Float(n)) - bs_lgamma_f(x))
        if S > 88.72283905206835 { // log(Float.greatestFiniteMagnitude)
            throw SpecialFunctionError.invalidCombination(message: "rising_factorial overflows for given x and n in Float")
        }
        if S < -103.27892990343185 { // log(Float.leastNonzeroMagnitude)
            return 0
        }
        
        return bs_rising_factorial_f(x, n)
    }
    
#if arch(i386) || arch(x86_64)
    /// Rising factorial (x)_n for `Float80` (x86 only) with magnitude checks.
    ///
    /// Behavior:
    /// - Throws on overflow based on ln-magnitude check.
    /// - Exact zero is detected via a Double-based check.
    /// - Underflow pre-check is conservative and may be omitted.
    ///
    /// Parameters:
    /// - x: Finite base value.
    /// - n: Non-negative integer order.
    ///
    /// Returns:
    /// - (x)_n as `Float80`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotFinite(name: "x")` if x is NaN or ±∞.
    /// - `SpecialFunctionError.invalidCombination(message: ...)` if the result would overflow Float80.
    @inlinable static func rising_factorial_l (_ x: Float80, _ n: UInt32) throws -> Float80 {
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        // Exact-zero check via Double
        if risingFactorialIsExactlyZero(Double(x), n: n) { return 0 }
        
        // Compute in extended precision for log-gamma
        let S = Double(bs_lgamma_l(x + Float80(n)) - bs_lgamma_l(x))
        if S > 11356.523406294143 { // log(max Float80)
            throw SpecialFunctionError.invalidCombination(message: "rising_factorial overflows for given x and n in Float80")
        }
        // Underflow pre-check for Float80 is conservative; skip or set a very negative threshold.
        // if S < -11399 { return 0 }
        
        return bs_rising_factorial_l(x, n)
    }
#endif
    
    /// Alias for the Pochhammer symbol (x)_{n} = rising_factorial(x, n).
    ///
    /// Parameters:
    /// - x: Finite base value.
    /// - n: Non-negative integer order.
    ///
    /// Returns:
    /// - (x)_n as `T`, with the same checks and throwing behavior as rising_factorial(_:_:).
    @inlinable static func pochhammer<T: BinaryFloatingPoint> (_ x: T, _ n: UInt32) throws -> T {
        try rising_factorial(x, n)
    }
    
    /// Alias for the Pochhammer symbol (x)_{n} for `Double`.
    @inlinable static func pochhammer(_ x: Double,_ n: UInt32) throws -> Double {
        try rising_factorial(x, n)
    }
    
    /// Alias for the Pochhammer symbol (x)_{n} for `Float`.
    @inlinable static func pochhammer_f (_ x: Float, _ n: UInt32) throws -> Float {
        try rising_factorial_f(x, n)
    }
    
#if arch(i386) || arch(x86_64)
    /// Alias for the Pochhammer symbol (x)_{n} for `Float80` (x86 only).
    @inlinable static func pochhammer_l (_ x: Float80, _ n: UInt32) throws -> Float80 {
        try rising_factorial_l(x, n)
    }
#endif
    
    // MARK: - Binomial coefficient with overflow checks
    
    /// Compute the binomial coefficient C(n, k) returned as `T`.
    ///
    /// Behavior:
    /// - Returns 0 if `k > n`.
    /// - Returns 1 if `k == 0` or `k == n`.
    /// - Pre-checks overflow using S = ln Γ(n+1) − ln Γ(k+1) − ln Γ(n−k+1) against the
    ///   ln(max finite) for the target type `T`. Throws invalidCombination on overflow.
    ///
    /// Parameters:
    /// - n: Non-negative integer.
    /// - k: Non-negative integer.
    ///
    /// Returns:
    /// - C(n, k) as `T`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.invalidCombination(message: ...)` if the result would overflow `T`.
    @inlinable static func binomial_coeff<T: BinaryFloatingPoint> (_ n: UInt32, _ k: UInt32) throws -> T {
        if k > n { return T(0) }
        if k == 0 || k == n { return T(1) }
        
        // Compute S = ln Γ(n+1) − ln Γ(k+1) − ln Γ(n−k+1) in Double.
        let dn = Double(n), dk = Double(k), dnk = Double(n - k)
        let S = bs_lgamma(dn + 1) - bs_lgamma(dk + 1) - bs_lgamma(dnk + 1)
        
        let lnMax = logMaxFinite(for: T.self)
        if S > lnMax {
            throw SpecialFunctionError.invalidCombination(message: "binomial_coeff overflows for given n and k in the target type")
        }
        
        return T(bs_binomial_coefficient(n, k))
    }
    
    /// Compute the binomial coefficient C(n, k) as `Double`.
    ///
    /// Behavior:
    /// - Returns 0 if `k > n`. Returns 1 if `k == 0` or `k == n`.
    /// - Throws invalidCombination if the result would overflow Double.
    ///
    /// Parameters:
    /// - n: Non-negative integer.
    /// - k: Non-negative integer.
    ///
    /// Returns:
    /// - C(n, k) as `Double`.
    @inlinable static func binomial_coeff(_ n: UInt32,_ k: UInt32) throws -> Double {
        if k > n { return 0 }
        if k == 0 || k == n { return 1 }
        
        let dn = Double(n), dk = Double(k), dnk = Double(n - k)
        let S = bs_lgamma(dn + 1) - bs_lgamma(dk + 1) - bs_lgamma(dnk + 1)
        if S > 709.782712893384 { // log(Double.greatestFiniteMagnitude)
            throw SpecialFunctionError.invalidCombination(message: "binomial_coeff overflows for given n and k in Double")
        }
        
        return bs_binomial_coefficient(n, k)
    }
    
    /// Compute the binomial coefficient C(n, k) as `Float`.
    ///
    /// Behavior:
    /// - Returns 0 if `k > n`. Returns 1 if `k == 0` or `k == n`.
    /// - Throws invalidCombination if the result would overflow Float.
    ///
    /// Parameters:
    /// - n: Non-negative integer.
    /// - k: Non-negative integer.
    ///
    /// Returns:
    /// - C(n, k) as `Float`.
    @inlinable static func binomial_coeff_f(_ n: UInt32, _ k: UInt32) throws -> Float {
        if k > n { return 0 }
        if k == 0 || k == n { return 1 }
        
        // Use Float lgamma for consistency with rising_factorial_f and cast to Double for the bound check.
        let fn = Float(n), fk = Float(k), fnk = Float(n - k)
        let S = Double(bs_lgamma_f(fn + 1) - bs_lgamma_f(fk + 1) - bs_lgamma_f(fnk + 1))
        if S > 88.72283905206835 { // log(Float.greatestFiniteMagnitude)
            throw SpecialFunctionError.invalidCombination(message: "binomial_coeff overflows for given n and k in Float")
        }
        
        return bs_binomial_coefficient_f(n, k)
    }
    
#if arch(i386) || arch(x86_64)
    /// Compute the binomial coefficient C(n, k) as `Float80` (x86 only).
    ///
    /// Behavior:
    /// - Returns 0 if `k > n`. Returns 1 if `k == 0` or `k == n`.
    /// - Throws invalidCombination if the result would overflow Float80.
    ///
    /// Parameters:
    /// - n: Non-negative integer.
    /// - k: Non-negative integer.
    ///
    /// Returns:
    /// - C(n, k) as `Float80`.
    @inlinable static func binomial_coeff_l (_ n: UInt32, _ k: UInt32) throws -> Float80 {
        if k > n { return 0 }
        if k == 0 || k == n { return 1 }
        
        let ln = Float80(n), lk = Float80(k), lnk = Float80(n - k)
        let S = Double(bs_lgamma_l(ln + 1) - bs_lgamma_l(lk + 1) - bs_lgamma_l(lnk + 1))
        if S > 11356.523406294143 { // log(max Float80)
            throw SpecialFunctionError.invalidCombination(message: "binomial_coeff overflows for given n and k in Float80")
        }
        
        return bs_binomial_coefficient_l(n, k)
    }
#endif

    // MARK: - Double factorial n!! with overflow checks
    
    /// Compute the natural logarithm of n!! in `Double` for magnitude checks.
    ///
    /// Definitions:
    /// - Even n = 2k: (2k)!! = 2^k · k!
    /// - Odd  n = 2k − 1: (2k − 1)!! = 2^k · Γ(k + 1/2) / √π
    ///
    /// Purpose:
    /// - Used to pre-check overflow for double factorial before calling the backend.
    ///
    /// Parameters:
    /// - n: Non-negative integer.
    ///
    /// Returns:
    /// - ln(n!!) as `Double`.
    @inline(__always) static func logDoubleFactorial(_ n: UInt32) -> Double {
        if n == 0 || n == 1 { return 0 }
        if n % 2 == 0 {
            // n = 2k: (2k)!! = 2^k * k!
            let k = Double(n / 2)
            return k * log(2.0) + bs_lgamma(k + 1.0)
        } else {
            // n = 2k - 1: (2k-1)!! = 2^k * Γ(k + 1/2) / √π
            let k = Double((n + 1) / 2) // integer ceil(n/2)
            return k * log(2.0) + bs_lgamma(k + 0.5) - 0.5 * log(Double.pi)
        }
    }
    
    /// Compute the double factorial n!! returned as `T`, with overflow checking.
    ///
    /// Behavior:
    /// - Uses a logarithmic magnitude estimate via logDoubleFactorial(_:) to ensure
    ///   the result fits in the target type `T`. Throws invalidCombination on overflow.
    ///
    /// Parameters:
    /// - n: Non-negative integer.
    ///
    /// Returns:
    /// - n!! as `T`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.invalidCombination(message: ...)` if the result would overflow `T`.
    @inlinable static func double_factorial<T: BinaryFloatingPoint>(_ n: UInt32) throws -> T {
        // Magnitude check via log
        let S = logDoubleFactorial(n)
        let lnMax = logMaxFinite(for: T.self)
        if S > lnMax {
            throw SpecialFunctionError.invalidCombination(message: "double_factorial overflows for given n in the target type")
        }
        return T(bs_double_factorial(n))
    }
    
    /// Compute the double factorial n!! as `Double`, with overflow checking.
    ///
    /// Parameters:
    /// - n: Non-negative integer.
    ///
    /// Returns:
    /// - n!! as `Double`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.invalidCombination(message: ...)` if the result would overflow Double.
    @inlinable static func double_factorial(_ n: UInt32) throws -> Double {
        let S = logDoubleFactorial(n)
        if S > 709.782712893384 { // log(Double.greatestFiniteMagnitude)
            throw SpecialFunctionError.invalidCombination(message: "double_factorial overflows for given n in Double")
        }
        return bs_double_factorial(n)
    }
    
    /// Compute the double factorial n!! as `Float`, with overflow checking.
    ///
    /// Parameters:
    /// - n: Non-negative integer.
    ///
    /// Returns:
    /// - n!! as `Float`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.invalidCombination(message: ...)` if the result would overflow Float.
    @inlinable static func double_factorial_f(_ n: UInt32) throws -> Float {
        // For Float we can compute the log using Double (sufficient and safe).
        let S = logDoubleFactorial(n)
        if S > 88.72283905206835 { // log(Float.greatestFiniteMagnitude)
            throw SpecialFunctionError.invalidCombination(message: "double_factorial overflows for given n in Float")
        }
        return bs_double_factorial_f(n)
    }
    
#if arch(i386) || arch(x86_64)
    /// Compute the double factorial n!! as `Float80` (x86 only), with overflow checking.
    ///
    /// Parameters:
    /// - n: Non-negative integer.
    ///
    /// Returns:
    /// - n!! as `Float80`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.invalidCombination(message: ...)` if the result would overflow Float80.
    @inlinable static func double_factorial_l(_ n: UInt32) throws -> Float80 {
        let S = logDoubleFactorial(n)
        if S > 11356.523406294143 { // log(max Float80)
            throw SpecialFunctionError.invalidCombination(message: "double_factorial overflows for given n in Float80")
        }
        return bs_double_factorial_l(n)
    }
#endif
}
