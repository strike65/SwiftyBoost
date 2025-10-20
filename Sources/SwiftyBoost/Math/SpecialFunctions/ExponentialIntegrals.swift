// Sources/Math/SpecialFunctions/ExponentialIntegrals.swift

import CBoostBridge

public extension SpecialFunctions {
    
    // MARK: - Exponential integral Ei(x)
    
    /// Exponential integral Ei(x) (principal value).
    ///
    /// Definition:
    /// - Ei(x) = PV ∫_{−∞}^x (e^t / t) dt, where PV denotes the Cauchy principal value.
    ///
    /// Properties:
    /// - Ei(x) is real-valued for real x ≠ 0, with a logarithmic singularity at x = 0.
    /// - As x → +∞, Ei(x) ~ e^x / x · (1 + 1/x + 2!/x^2 + …).
    /// - As x → −∞, Ei(x) → 0⁻.
    ///
    /// Domain and validation:
    /// - Input must be finite (NaN and ±∞ are rejected).
    /// - No attempt is made to mask the singularity at x = 0; very near zero, large magnitude values are expected.
    ///
    /// Numeric behavior:
    /// - Delegates to Boost.Math’s expint(x) overload for robust evaluation across the real line.
    ///
    /// Parameters:
    /// - x: Real argument (finite).
    ///
    /// Returns:
    /// - Ei(x) as `T`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotFinite(name: "x")` if `x` is NaN or ±∞.
    ///
    /// References:
    /// - NIST DLMF §6.2 (Exponential Integral Ei)
    /// - Boost.Math expint(x)
    @inlinable static func exponentialIntegralEi<T: BinaryFloatingPoint>(_ x: T) throws -> T {
        let dx = D(x)
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return T(bs_expint_Ei_d(dx))
    }
    
    /// Exponential integral Ei(x) for `Float`.
    @inlinable static func exponentialIntegralEi(_ x: Float) throws -> Float {
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return bs_expint_Ei_f(x)
    }
    
#if arch(x86_64)
    /// Exponential integral Ei(x) for `Float80` (x86_64 only).
    @inlinable static func exponentialIntegralEi(_ x: Float80) throws -> Float80 {
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return bs_expint_Ei_l(x)
    }
#endif
    
    // MARK: - Exponential integral E_n(x)
    
    /// Generalized exponential integral E_n(x).
    ///
    /// Definition:
    /// - For integer n ≥ 0 and real x, the generalized exponential integral is
    ///   E_n(x) = ∫₁^∞ e^{−xt} t^{−n} dt.
    ///
    /// Special cases and relations:
    /// - E_1(x) = E₁(x) is sometimes denoted as the exponential integral of order one.
    /// - For x > 0 and n ≥ 1, E_n(x) is finite and positive.
    /// - Recurrence relations connect E_{n+1}, E_n, and derivatives with respect to x.
    ///
    /// Domain and validation:
    /// - Requires `n ≥ 0` (integer).
    /// - Requires `x` to be finite (NaN and ±∞ are rejected).
    ///
    /// Numeric behavior:
    /// - Delegates to Boost.Math’s expint(n, x), which selects stable algorithms across parameter ranges.
    ///
    /// Parameters:
    /// - n: Nonnegative integer order.
    /// - x: Real argument (finite).
    ///
    /// Returns:
    /// - E_n(x) as `T`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotPositive(name: "n")` if `n < 0`.
    /// - `SpecialFunctionError.parameterNotFinite(name: "x")` if `x` is NaN or ±∞.
    ///
    /// References:
    /// - NIST DLMF §8.20 (Exponential Integrals E_n)
    /// - Boost.Math expint(n, x)
    @inlinable static func exponentialIntegralEn<T: BinaryFloatingPoint>(_ n: Int, _ x: T) throws -> T {
        guard n >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "n") }
        let dx = D(x)
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return T(bs_expint_En_d(Int32(n), dx))
    }
    
    /// Generalized exponential integral E_n(x) for `Float`.
    @inlinable static func exponentialIntegralEn(_ n: Int, _ x: Float) throws -> Float {
        guard n >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "n") }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return bs_expint_En_f(Int32(n), x)
    }
    
#if arch(x86_64)
    /// Generalized exponential integral E_n(x) for `Float80` (x86_64 only).
    @inlinable static func exponentialIntegralEn(_ n: Int, _ x: Float80) throws -> Float80 {
        guard n >= 0 else { throw SpecialFunctionError.parameterNotPositive(name: "n") }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        return bs_expint_En_l(Int32(n), x)
    }
#endif
}
