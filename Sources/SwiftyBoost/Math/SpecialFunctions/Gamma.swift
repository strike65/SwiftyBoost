//
//  GammaAndError.swift
//  Math/SpecialFunctions
//
//  This file exposes Gamma-related special functions with type-generic
//  and type-specific overloads. The implementations are backed by the
//  CBoostBridge C shim (wrapping Boost.Math), while performing Swift-side
//  argument validation and domain checks.
//
//  All functions throw Swift errors for invalid inputs (e.g., non-finite
//  values) and at poles where the function is not defined (e.g., non-positive
//  integers for Gamma). See each function’s documentation for details.
//

import CBoostBridge
public extension SpecialFunctions {
    
    
    // MARK: - Gamma (generic Double-backed)
    
    /// Compute the Gamma function Γ(x).
    ///
    /// Overview:
    /// - This generic overload accepts any `BinaryFloatingPoint` argument and returns
    ///   the value as the same generic type `T`.
    /// - Internally, the argument is converted to `Double` and evaluated using the
    ///   C shim provided by `CBoostBridge` (wrapping Boost.Math), then converted
    ///   back to `T`.
    ///
    /// Mathematical domain:
    /// - Γ(x) has simple poles at non-positive integers (x ∈ {0, -1, -2, ...}).
    ///   Calling this function at those points will throw a `poleAtNonPositiveInteger`.
    ///
    /// Numeric behavior:
    /// - For finite, non-pole inputs, the result is finite where defined.
    /// - Very large positive inputs can overflow depending on `T`.
    /// - Subnormal and very small magnitudes are supported but can lose precision.
    /// - If you require extended precision on x86_64, see the `Float80` overloads.
    ///
    /// Thread-safety:
    /// - This function is pure and thread-safe (no shared mutable state).
    ///
    /// Parameters:
    /// - x: The input value `x`.
    ///
    /// Returns:
    /// - The value of Γ(x) as `T`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotFinite(name: "x")` if `x` is NaN or ±∞.
    /// - `SpecialFunctionError.poleAtNonPositiveInteger(name: "x")` if `x` is a non-positive integer.
    ///
    /// Example:
    /// ```swift
    /// do {
    ///     let valueD: Double = try gamma(5.0)     // 24.0
    ///     let valueF: Float  = try gamma(5 as Float) // 24.0
    /// } catch {
    ///     // Handle invalid inputs or poles
    /// }
    /// ```
    @inlinable static func gamma<T: BinaryFloatingPoint>(_ x: T) throws -> T {
        let dx = D(x)
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        if dx <= 0, dx == dx.rounded(.towardZero) {
            throw SpecialFunctionError.poleAtNonPositiveInteger(name: "x")
        }
        return T(bs_tgamma(dx))
    }
    
    /// Compute the natural logarithm of the Gamma function ln Γ(x).
    ///
    /// Overview:
    /// - This generic overload accepts any `BinaryFloatingPoint` argument and returns
    ///   the value as the same generic type `T`.
    /// - Internally, the argument is converted to `Double` and evaluated using the
    ///   C shim provided by `CBoostBridge` (wrapping Boost.Math), then converted
    ///   back to `T`.
    ///
    /// Why log-gamma:
    /// - `logGamma` is often numerically more stable than `gamma` for large `x`
    ///   because it avoids overflow.
    ///
    /// Mathematical domain:
    /// - ln Γ(x) inherits the poles of Γ(x) at non-positive integers.
    ///   Calling this function at those points will throw.
    ///
    /// Thread-safety:
    /// - This function is pure and thread-safe.
    ///
    /// Parameters:
    /// - x: The input value `x`.
    ///
    /// Returns:
    /// - The value of ln Γ(x) as `T`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotFinite(name: "x")` if `x` is NaN or ±∞.
    /// - `SpecialFunctionError.poleAtNonPositiveInteger(name: "x")` if `x` is a non-positive integer.
    ///
    /// Example:
    /// ```swift
    /// do {
    ///     let lg10: Double = try logGamma(10.0) // ≈ 12.801827...
    /// } catch {
    ///     // Handle invalid inputs or poles
    /// }
    /// ```
    @inlinable static func logGamma<T: BinaryFloatingPoint>(_ x: T) throws -> T {
        let dx = D(x)
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        if dx <= 0, dx == dx.rounded(.towardZero) {
            throw SpecialFunctionError.poleAtNonPositiveInteger(name: "x")
        }
        return T(bs_lgamma(dx))
    }
    
    // MARK: - Gamma ratios
    
    /// Compute the ratio Γ(a) / Γ(b) in a numerically stable way.
    ///
    /// Discussion:
    /// - Directly computing `gamma(a) / gamma(b)` is often numerically unstable and can
    ///   overflow/underflow even when the ratio itself is well-scaled.
    /// - This function delegates to Boost.Math’s specialized routine to evaluate the ratio
    ///   more stably across a broad range of inputs.
    ///
    /// Domain:
    /// - Requires `a` and `b` to be finite reals.
    /// - Γ has simple poles at non-positive integers. If either `a` or `b` is a non-positive
    ///   integer (i.e., ≤ 0 and integer-valued), this function throws `poleAtNonPositiveInteger`.
    ///
    /// Numeric behavior:
    /// - Stable across a wide range of inputs.
    ///
    /// Parameters:
    /// - a: Numerator argument to Γ(·).
    /// - b: Denominator argument to Γ(·).
    ///
    /// Returns:
    /// - The value of Γ(a) / Γ(b) as `T`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotFinite(name: "a")` if `a` is NaN or ±∞.
    /// - `SpecialFunctionError.parameterNotFinite(name: "b")` if `b` is NaN or ±∞.
    /// - `SpecialFunctionError.poleAtNonPositiveInteger(name: "a")` if `a` is a non‑positive integer.
    /// - `SpecialFunctionError.poleAtNonPositiveInteger(name: "b")` if `b` is a non‑positive integer.
    ///
    /// SeeAlso:
    /// - ``SpecialFunctions/gammaDeltaRatio(_:delta:)->T`` for Γ(a)/Γ(a+δ).
    /// - ``SpecialFunctions/gamma(_:)->T`` and ``SpecialFunctions/logGamma(_:)->T`` for base functions.
    ///
    /// References:
    /// - Boost.Math `tgamma_ratio`
    @inlinable static func gammaRatio<T: BinaryFloatingPoint>(_ a: T, _ b: T) throws -> T {
        let da = D(a), db = D(b)
        guard da.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a") }
        guard db.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "b") }
        if da <= 0, da == da.rounded(.towardZero) { throw SpecialFunctionError.poleAtNonPositiveInteger(name: "a") }
        if db <= 0, db == db.rounded(.towardZero) { throw SpecialFunctionError.poleAtNonPositiveInteger(name: "b") }
        return T(bs_tgamma_ratio(da, db))
    }
    
    /// Compute the ratio Γ(a) / Γ(a + δ) in a numerically stable way.
    ///
    /// Discussion:
    /// - This “delta ratio” occurs frequently in probability, statistics, and combinatorics.
    ///   Specialized evaluation avoids forming Γ values explicitly, improving stability and
    ///   performance compared to naive evaluation.
    ///
    /// Domain:
    /// - Requires `a` and `a + δ` to be finite reals.
    /// - Γ has simple poles at non-positive integers. If `a` or `a+δ` is a non-positive
    ///   integer, this function throws `poleAtNonPositiveInteger`.
    ///
    /// Numeric behavior:
    /// - Stable across a wide range of inputs; may overflow for extreme magnitudes depending on `T`.
    ///
    /// Parameters:
    /// - a: Base argument to Γ(·) in the numerator.
    /// - delta: Offset applied to `a` in the denominator (Γ(a + δ)).
    ///
    /// Returns:
    /// - The value of Γ(a) / Γ(a + δ) as `T`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotFinite(name: "a")` if `a` is NaN or ±∞.
    /// - `SpecialFunctionError.parameterNotFinite(name: "delta")` if `delta` is NaN or ±∞.
    /// - `SpecialFunctionError.poleAtNonPositiveInteger(name: "a")` if `a` is a non‑positive integer.
    /// - `SpecialFunctionError.poleAtNonPositiveInteger(name: "a + delta")` if `a + δ` is a non‑positive integer.
    ///
    /// SeeAlso:
    /// - ``SpecialFunctions/gammaRatio(_:_:)`` for Γ(a)/Γ(b) with independent parameters.
    ///
    /// References:
    /// - Boost.Math `tgamma_delta_ratio`
    @inlinable static func gammaDeltaRatio<T: BinaryFloatingPoint>(_ a: T, delta: T) throws -> T {
        let da = D(a), dd = D(delta)
        guard da.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a") }
        guard dd.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "delta") }
        let adb = da + dd
        if da <= 0, da == da.rounded(.towardZero) { throw SpecialFunctionError.poleAtNonPositiveInteger(name: "a") }
        if adb <= 0, adb == adb.rounded(.towardZero) { throw SpecialFunctionError.poleAtNonPositiveInteger(name: "a + delta") }
        return T(bs_tgamma_delta_ratio(da, dd))
    }
    
    // MARK: - Incomplete gamma (lower/upper, regularized, and inverses)
    
    /// Compute the lower incomplete gamma function γ(a, x).
    ///
    /// Definition:
    /// - γ(a, x) = ∫₀ˣ t^{a−1} e^{−t} dt, for a > 0 and x ≥ 0.
    ///
    /// Domain:
    /// - Requires `a > 0` and `x ≥ 0`.
    ///
    /// Discussion:
    /// - This function is monotone increasing in `x` for fixed `a > 0`.
    /// - Related identity: Γ(a) = γ(a, x) + Γ(a, x).
    ///
    /// Numeric behavior:
    /// - Delegates to a robust Boost.Math routine with good coverage for small and large `x`.
    ///
    /// Parameters:
    /// - a: Shape parameter `a` (strictly positive).
    /// - x: Upper integration limit `x` (non-negative).
    ///
    /// Returns:
    /// - The value of γ(a, x) as `T`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotFinite(name: "a"|"x")` if inputs are NaN or ±∞.
    /// - `SpecialFunctionError.parameterNotPositive(name: "a")` if `a ≤ 0`.
    /// - `SpecialFunctionError.parameterOutOfRange(name: "x", min: 0, max: +∞)` if `x < 0`.
    ///
    /// SeeAlso:
    /// - ``SpecialFunctions/incompleteGammaUpper(_:x:)->T``
    /// - ``SpecialFunctions/regularizedGammaP(_:x:)->T`` and ``SpecialFunctions/regularizedGammaQ(_:x:)->T``
    @inlinable static func incompleteGammaLower<T: BinaryFloatingPoint>(_ a: T, x: T) throws -> T {
        let da = D(a), dx = D(x)
        guard da.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a") }
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        guard da > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
        guard dx >= 0 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 0.0, max: Double.infinity) }
        return T(bs_tgamma_lower(da, dx))
    }
    
    /// Compute the upper incomplete gamma function Γ(a, x).
    ///
    /// Definition:
    /// - Γ(a, x) = ∫ₓ^∞ t^{a−1} e^{−t} dt, for a > 0 and x ≥ 0.
    ///
    /// Relation:
    /// - Γ(a) = γ(a, x) + Γ(a, x).
    ///
    /// Domain:
    /// - Requires `a > 0` and `x ≥ 0`.
    ///
    /// Discussion:
    /// - This function is monotone decreasing in `x` for fixed `a > 0`.
    ///
    /// Parameters:
    /// - a: Shape parameter `a` (strictly positive).
    /// - x: Lower integration limit `x` (non-negative).
    ///
    /// Returns:
    /// - The value of Γ(a, x) as `T`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotFinite(name: "a"|"x")` if inputs are NaN or ±∞.
    /// - `SpecialFunctionError.parameterNotPositive(name: "a")` if `a ≤ 0`.
    /// - `SpecialFunctionError.parameterOutOfRange(name: "x", min: 0, max: +∞)` if `x < 0`.
    ///
    /// SeeAlso:
    /// - ``SpecialFunctions/incompleteGammaLower(_:x:)->T``
    /// - ``SpecialFunctions/regularizedGammaP(_:x:)->T`` and ``SpecialFunctions/regularizedGammaQ(_:x:)->T``
    @inlinable static func incompleteGammaUpper<T: BinaryFloatingPoint>(_ a: T, x: T) throws -> T {
        let da = D(a), dx = D(x)
        guard da.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a") }
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        guard da > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
        guard dx >= 0 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 0.0, max: Double.infinity) }
        return T(bs_tgamma_upper(da, dx))
    }
    
    /// Compute the regularized lower incomplete gamma function P(a, x).
    ///
    /// Definition:
    /// - P(a, x) = γ(a, x) / Γ(a), for a > 0 and x ≥ 0.
    ///
    /// Range:
    /// - P(a, x) ∈ [0, 1] for all valid inputs.
    ///
    /// Domain:
    /// - Requires `a > 0` and `x ≥ 0`.
    ///
    /// Discussion:
    /// - Monotone increasing in `x` for fixed `a > 0`.
    /// - Related identity: Q(a, x) = 1 − P(a, x).
    ///
    /// Parameters:
    /// - a: Shape parameter `a` (strictly positive).
    /// - x: Evaluation point (non-negative).
    ///
    /// Returns:
    /// - The value of P(a, x) as `T`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotFinite(name: "a"|"x")` if inputs are NaN or ±∞.
    /// - `SpecialFunctionError.parameterNotPositive(name: "a")` if `a ≤ 0`.
    /// - `SpecialFunctionError.parameterOutOfRange(name: "x", min: 0, max: +∞)` if `x < 0`.
    ///
    /// SeeAlso:
    /// - ``SpecialFunctions/regularizedGammaQ(_:x:)->T``
    /// - ``SpecialFunctions/regularizedGammaPInv(_:p:)->T``
    @inlinable static func regularizedGammaP<T: BinaryFloatingPoint>(_ a: T, x: T) throws -> T {
        let da = D(a), dx = D(x)
        guard da.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a") }
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        guard da > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
        guard dx >= 0 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 0.0, max: Double.infinity) }
        return T(bs_gamma_p(da, dx))
    }
    
    /// Compute the regularized upper incomplete gamma function Q(a, x).
    ///
    /// Definition:
    /// - Q(a, x) = Γ(a, x) / Γ(a) = 1 − P(a, x), for a > 0 and x ≥ 0.
    ///
    /// Range:
    /// - Q(a, x) ∈ [0, 1] for all valid inputs.
    ///
    /// Domain:
    /// - Requires `a > 0` and `x ≥ 0`.
    ///
    /// Discussion:
    /// - Monotone decreasing in `x` for fixed `a > 0`.
    ///
    /// Parameters:
    /// - a: Shape parameter `a` (strictly positive).
    /// - x: Evaluation point (non-negative).
    ///
    /// Returns:
    /// - The value of Q(a, x) as `T`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotFinite(name: "a"|"x")` if inputs are NaN or ±∞.
    /// - `SpecialFunctionError.parameterNotPositive(name: "a")` if `a ≤ 0`.
    /// - `SpecialFunctionError.parameterOutOfRange(name: "x", min: 0, max: +∞)` if `x < 0`.
    ///
    /// SeeAlso:
    /// - ``SpecialFunctions/regularizedGammaP(_:x:)->T``
    /// - ``SpecialFunctions/regularizedGammaQInv(_:q:)->T``
    @inlinable static func regularizedGammaQ<T: BinaryFloatingPoint>(_ a: T, x: T) throws -> T {
        let da = D(a), dx = D(x)
        guard da.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a") }
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        guard da > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
        guard dx >= 0 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 0.0, max: Double.infinity) }
        return T(bs_gamma_q(da, dx))
    }
    
    /// Invert the regularized lower incomplete gamma: find x such that P(a, x) = p.
    ///
    /// Domain:
    /// - Requires `a > 0` and `p ∈ [0, 1]`.
    ///
    /// Discussion:
    /// - Returns the unique `x ≥ 0` satisfying the equation for valid `(a, p)`.
    ///
    /// Parameters:
    /// - a: Shape parameter `a` (strictly positive).
    /// - p: Target probability in `[0, 1]`.
    ///
    /// Returns:
    /// - The unique `x ≥ 0` such that P(a, x) = p.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotFinite(name: "a"|"p")` if inputs are NaN or ±∞.
    /// - `SpecialFunctionError.parameterNotPositive(name: "a")` if `a ≤ 0`.
    /// - `SpecialFunctionError.parameterOutOfRange(name: "p", min: 0, max: 1)` if `p ∉ [0, 1]`.
    ///
    /// SeeAlso:
    /// - ``SpecialFunctions/regularizedGammaP(_:x:)->T``
    @inlinable static func regularizedGammaPInv<T: BinaryFloatingPoint>(_ a: T, p: T) throws -> T {
        let da = D(a), dp = D(p)
        guard da.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a") }
        guard dp.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "p") }
        guard da > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
        guard dp >= 0 && dp <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "p", min: 0.0, max: 1.0) }
        return T(bs_gamma_p_inv(da, dp))
    }
    
    /// Invert the regularized upper incomplete gamma: find x such that Q(a, x) = q.
    ///
    /// Domain:
    /// - Requires `a > 0` and `q ∈ [0, 1]`.
    ///
    /// Discussion:
    /// - Returns the unique `x ≥ 0` satisfying the equation for valid `(a, q)`.
    ///
    /// Parameters:
    /// - a: Shape parameter `a` (strictly positive).
    /// - q: Target tail probability in `[0, 1]`.
    ///
    /// Returns:
    /// - The unique `x ≥ 0` such that Q(a, x) = q.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotFinite(name: "a"|"q")` if inputs are NaN or ±∞.
    /// - `SpecialFunctionError.parameterNotPositive(name: "a")` if `a ≤ 0`.
    /// - `SpecialFunctionError.parameterOutOfRange(name: "q", min: 0, max: 1)` if `q ∉ [0, 1]`.
    ///
    /// SeeAlso:
    /// - ``SpecialFunctions/regularizedGammaQ(_:x:)->T``
    @inlinable static func regularizedGammaQInv<T: BinaryFloatingPoint>(_ a: T, q: T) throws -> T {
        let da = D(a), dq = D(q)
        guard da.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a") }
        guard dq.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "q") }
        guard da > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
        guard dq >= 0 && dq <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "q", min: 0.0, max: 1.0) }
        return T(bs_gamma_q_inv(da, dq))
    }
    
    // MARK: - Derivatives of regularized incomplete gamma (w.r.t. x)
    
    /// Derivative with respect to x of the regularized lower incomplete gamma P(a, x).
    ///
    /// Definition:
    /// - For `a > 0` and `x ≥ 0`, the derivative is:
    ///   d/dx P(a, x) = exp(−x) x^{a−1} / Γ(a).
    ///
    /// Properties:
    /// - Nonnegative for all valid inputs.
    /// - Integrates to P(a, x) up to a constant of integration.
    /// - Since Q(a, x) = 1 − P(a, x), it follows that d/dx Q(a, x) = − d/dx P(a, x).
    ///
    /// Domain:
    /// - Requires `a > 0` and `x ≥ 0`.
    ///
    /// Numeric behavior:
    /// - Stable across a wide range via Boost.Math; extremely small or large x may
    ///   lead to underflow/overflow depending on the floating-point type.
    ///
    /// Parameters:
    /// - a: Shape parameter (strictly positive).
    /// - x: Evaluation point (nonnegative).
    ///
    /// Returns:
    /// - The value of d/dx P(a, x) as `T`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotFinite(name: "a"|"x")` if inputs are NaN or ±∞.
    /// - `SpecialFunctionError.parameterNotPositive(name: "a")` if `a ≤ 0`.
    /// - `SpecialFunctionError.parameterOutOfRange(name: "x", min: 0, max: +∞)` if `x < 0`.
    ///
    /// SeeAlso:
    /// - ``SpecialFunctions/regularizedGammaP(_:x:)->T``
    /// - ``SpecialFunctions/regularizedGammaQ(_:x:)->T``
    @inlinable static func regularizedGammaPDerivative<T: BinaryFloatingPoint>(_ a: T, x: T) throws -> T {
        let da = D(a), dx = D(x)
        guard da.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a") }
        guard dx.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        guard da > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
        guard dx >= 0 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 0.0, max: Double.infinity) }
        return T(bs_gamma_p_derivative(da, dx))
    }
    

    // MARK: - Float overloads
    
    /// Compute the Gamma function Γ(x) for `Float`.
    ///
    /// This overload evaluates directly in `Float` precision via `CBoostBridge`.
    ///
    /// Parameters:
    /// - x: The input value `x` as `Float`.
    ///
    /// Returns:
    /// - The value of Γ(x) as `Float`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotFinite(name: "x")` if `x` is NaN or ±∞.
    /// - `SpecialFunctionError.poleAtNonPositiveInteger(name: "x")` if `x` is a non-positive integer.
    ///
    /// Example:
    /// ```swift
    /// let value = try gamma(3.5 as Float)
    /// ```
    @inlinable static func gamma(_ x: Float) throws -> Float {
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        if x <= 0, x == x.rounded(.towardZero) { throw SpecialFunctionError.poleAtNonPositiveInteger(name: "x") }
        return bs_tgamma_f(x)
    }
    
    /// Compute the natural logarithm of the Gamma function ln Γ(x) for `Float`.
    ///
    /// This overload evaluates directly in `Float` precision via `CBoostBridge`.
    ///
    /// Parameters:
    /// - x: The input value `x` as `Float`.
    ///
    /// Returns:
    /// - The value of ln Γ(x) as `Float`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotFinite(name: "x")` if `x` is NaN or ±∞.
    /// - `SpecialFunctionError.poleAtNonPositiveInteger(name: "x")` if `x` is a non-positive integer.
    ///
    /// Example:
    /// ```swift
    /// let value = try logGamma(10 as Float)
    /// ```
    @inlinable static func logGamma(_ x: Float) throws -> Float {
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        if x <= 0, x == x.rounded(.towardZero) { throw SpecialFunctionError.poleAtNonPositiveInteger(name: "x") }
        return bs_lgamma_f(x)
    }
    
    /// Γ(a) / Γ(b) for `Float`.
    ///
    /// Discussion:
    /// - Direct `gamma(a)/gamma(b)` is unstable; this wrapper calls the Boost routine
    ///   designed for stable evaluation.
    ///
    /// See ``SpecialFunctions/gammaRatio(_:_:)`` for domain, behavior, and references.
    @inlinable static func gammaRatio(_ a: Float, _ b: Float) throws -> Float {
        guard a.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a") }
        guard b.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "b") }
        if a <= 0, a == a.rounded(.towardZero) { throw SpecialFunctionError.poleAtNonPositiveInteger(name: "a") }
        if b <= 0, b == b.rounded(.towardZero) { throw SpecialFunctionError.poleAtNonPositiveInteger(name: "b") }
        return bs_tgamma_ratio_f(a, b)
    }
    
    /// Γ(a) / Γ(a + δ) for `Float`.
    ///
    /// Discussion:
    /// - Specialized Boost routine for the delta ratio.
    ///
    /// See ``SpecialFunctions/gammaDeltaRatio(_:delta:)->T`` for domain, behavior, and references.
    @inlinable static func gammaDeltaRatio(_ a: Float, delta: Float) throws -> Float {
        guard a.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a") }
        guard delta.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "delta") }
        let adb = a + delta
        if a <= 0, a == a.rounded(.towardZero) { throw SpecialFunctionError.poleAtNonPositiveInteger(name: "a") }
        if adb <= 0, adb == adb.rounded(.towardZero) { throw SpecialFunctionError.poleAtNonPositiveInteger(name: "a + delta") }
        return bs_tgamma_delta_ratio_f(a, delta)
    }
    
    /// Lower incomplete gamma γ(a, x) for `Float`.
    ///
    /// See ``SpecialFunctions/incompleteGammaLower(_:x:)->T`` for definition, domain, and discussion.
    @inlinable static func incompleteGammaLower(_ a: Float, x: Float) throws -> Float {
        guard a.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a") }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
        guard x >= 0 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 0.0, max: Double.infinity) }
        return bs_tgamma_lower_f(a, x)
    }
    
    /// Upper incomplete gamma Γ(a, x) for `Float`.
    ///
    /// See ``SpecialFunctions/incompleteGammaUpper(_:x:)->T`` for definition, domain, and discussion.
    @inlinable static func incompleteGammaUpper(_ a: Float, x: Float) throws -> Float {
        guard a.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a") }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
        guard x >= 0 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 0.0, max: Double.infinity) }
        return bs_tgamma_upper_f(a, x)
    }
    
    /// Regularized lower incomplete gamma P(a, x) for `Float`.
    ///
    /// See ``SpecialFunctions/regularizedGammaP(_:x:)->T`` for definition, domain, and discussion.
    @inlinable static func regularizedGammaP(_ a: Float, x: Float) throws -> Float {
        guard a.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a") }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
        guard x >= 0 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 0.0, max: Double.infinity) }
        return bs_gamma_p_f(a, x)
    }
    
    /// Regularized upper incomplete gamma Q(a, x) for `Float`.
    ///
    /// See ``SpecialFunctions/regularizedGammaQ(_:x:)->T`` for definition, domain, and discussion.
    @inlinable static func regularizedGammaQ(_ a: Float, x: Float) throws -> Float {
        guard a.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a") }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
        guard x >= 0 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 0.0, max: Double.infinity) }
        return bs_gamma_q_f(a, x)
    }
    
    /// Inverse regularized lower incomplete gamma x = P⁻¹(a, p) for `Float`.
    ///
    /// See ``SpecialFunctions/regularizedGammaPInv(_:p:)->T`` for domain and discussion.
    @inlinable static func regularizedGammaPInv(_ a: Float, p: Float) throws -> Float {
        guard a.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a") }
        guard p.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "p") }
        guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
        guard p >= 0 && p <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "p", min: 0.0, max: 1.0) }
        return bs_gamma_p_inv_f(a, p)
    }
    
    /// Inverse regularized upper incomplete gamma x = Q⁻¹(a, q) for `Float`.
    ///
    /// See ``SpecialFunctions/regularizedGammaQInv(_:q:)->T`` for domain and discussion.
    @inlinable static func regularizedGammaQInv(_ a: Float, q: Float) throws -> Float {
        guard a.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a") }
        guard q.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "q") }
        guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
        guard q >= 0 && q <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "q", min: 0.0, max: 1.0) }
        return bs_gamma_q_inv_f(a, q)
    }
    
    /// Derivative d/dx P(a, x) for `Float`.
    ///
    /// See ``SpecialFunctions/regularizedGammaPDerivative(_:x:)->T`` for definition, domain, and discussion.
    @inlinable static func regularizedGammaPDerivative(_ a: Float, x: Float) throws -> Float {
        guard a.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a") }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
        guard x >= 0 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 0.0, max: Double.infinity) }
        return bs_gamma_p_derivative_f(a, x)
    }
    
    
    // MARK: - Float80 overloads (x86_64 only)
    
#if arch(x86_64)
    /// Compute the Gamma function Γ(x) for `Float80` (x86_64 only).
    ///
    /// This overload evaluates directly in extended precision (`Float80`) via `CBoostBridge`.
    /// Prefer this overload when you need more precision than `Double` on x86_64.
    ///
    /// Parameters:
    /// - x: The input value `x` as `Float80`.
    ///
    /// Returns:
    /// - The value of Γ(x) as `Float80`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotFinite(name: "x")` if `x` is NaN or ±∞.
    /// - `SpecialFunctionError.poleAtNonPositiveInteger(name: "x")` if `x` is a non-positive integer.
    ///
    /// Availability:
    /// - Only on x86_64 architectures where `Float80` is available.
    @inlinable static func gamma(_ x: Float80) throws -> Float80 {
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        if x <= 0, x == x.rounded(.towardZero) { throw SpecialFunctionError.poleAtNonPositiveInteger(name: "x") }
        return bs_tgamma_l(x)
    }
    
    /// Compute the natural logarithm of the Gamma function ln Γ(x) for `Float80` (x86_64 only).
    ///
    /// This overload evaluates directly in extended precision (`Float80`) via `CBoostBridge`.
    ///
    /// Parameters:
    /// - x: The input value `x` as `Float80`.
    ///
    /// Returns:
    /// - The value of ln Γ(x) as `Float80`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotFinite(name: "x")` if `x` is NaN or ±∞.
    /// - `SpecialFunctionError.poleAtNonPositiveInteger(name: "x")` if `x` is a non-positive integer.
    ///
    /// Availability:
    /// - Only on x86_64 architectures where `Float80` is available.
    @inlinable static func logGamma(_ x: Float80) throws -> Float80 {
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        if x <= 0, x == x.rounded(.towardZero) { throw SpecialFunctionError.poleAtNonPositiveInteger(name: "x") }
        return bs_lgamma_l(x)
    }
    
    /// Γ(a) / Γ(b) for `Float80` (x86_64 only).
    ///
    /// See ``SpecialFunctions/gammaRatio(_:_: )`` for discussion of domain and motivation.
    @inlinable static func gammaRatio(_ a: Float80, _ b: Float80) throws -> Float80 {
        guard a.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a") }
        guard b.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "b") }
        if a <= 0, a == a.rounded(.towardZero) { throw SpecialFunctionError.poleAtNonPositiveInteger(name: "a") }
        if b <= 0, b == b.rounded(.towardZero) { throw SpecialFunctionError.poleAtNonPositiveInteger(name: "b") }
        return bs_tgamma_ratio_l(a, b)
    }
    
    /// Γ(a) / Γ(a + δ) for `Float80` (x86_64 only).
    ///
    /// See ``SpecialFunctions/gammaDeltaRatio(_:delta:)`` for discussion of domain and motivation.
    @inlinable static func gammaDeltaRatio(_ a: Float80, delta: Float80) throws -> Float80 {
        guard a.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a") }
        guard delta.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "delta") }
        let adb = a + delta
        if a <= 0, a == a.rounded(.towardZero) { throw SpecialFunctionError.poleAtNonPositiveInteger(name: "a") }
        if adb <= 0, adb == adb.rounded(.towardZero) { throw SpecialFunctionError.poleAtNonPositiveInteger(name: "a + delta") }
        return bs_tgamma_delta_ratio_l(a, delta)
    }
    
    /// Lower incomplete gamma γ(a, x) for `Float80` (x86_64 only).
    ///
    /// See ``SpecialFunctions/incompleteGammaLower(_:x:)`` for definition, domain, and discussion.
    @inlinable static func incompleteGammaLower(_ a: Float80, x: Float80) throws -> Float80 {
        guard a.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a") }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
        guard x >= 0 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 0.0, max: Double.infinity) }
        return bs_tgamma_lower_l(a, x)
    }
    
    /// Upper incomplete gamma Γ(a, x) for `Float80` (x86_64 only).
    ///
    /// See ``SpecialFunctions/incompleteGammaUpper(_:x:)`` for definition, domain, and discussion.
    @inlinable static func incompleteGammaUpper(_ a: Float80, x: Float80) throws -> Float80 {
        guard a.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a") }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
        guard x >= 0 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 0.0, max: Double.infinity) }
        return bs_tgamma_upper_l(a, x)
    }
    
    /// Regularized lower incomplete gamma P(a, x) for `Float80` (x86_64 only).
    ///
    /// See ``SpecialFunctions/regularizedGammaP(_:x:)`` for definition, domain, and discussion.
    @inlinable static func regularizedGammaP(_ a: Float80, x: Float80) throws -> Float80 {
        guard a.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a") }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
        guard x >= 0 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 0.0, max: Double.infinity) }
        return bs_gamma_p_l(a, x)
    }
    
    /// Regularized upper incomplete gamma Q(a, x) for `Float80` (x86_64 only).
    ///
    /// See ``SpecialFunctions/regularizedGammaQ(_:x:)`` for definition, domain, and discussion.
    @inlinable static func regularizedGammaQ(_ a: Float80, x: Float80) throws -> Float80 {
        guard a.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a") }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
        guard x >= 0 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 0.0, max: Double.infinity) }
        return bs_gamma_q_l(a, x)
    }
    
    /// Inverse regularized lower incomplete gamma x = P⁻¹(a, p) for `Float80` (x86_64 only).
    ///
    /// See ``SpecialFunctions/regularizedGammaPInv(_:p:)`` for domain and discussion.
    @inlinable static func regularizedGammaPInv(_ a: Float80, p: Float80) throws -> Float80 {
        guard a.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a") }
        guard p.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "p") }
        guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
        guard p >= 0 && p <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "p", min: 0.0, max: 1.0) }
        return bs_gamma_p_inv_l(a, p)
    }
    
    /// Inverse regularized upper incomplete gamma x = Q⁻¹(a, q) for `Float80` (x86_64 only).
    ///
    /// See ``SpecialFunctions/regularizedGammaQInv(_:q:)`` for domain and discussion.
    @inlinable static func regularizedGammaQInv(_ a: Float80, q: Float80) throws -> Float80 {
        guard a.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a") }
        guard q.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "q") }
        guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
        guard q >= 0 && q <= 1 else { throw SpecialFunctionError.parameterOutOfRange(name: "q", min: 0.0, max: 1.0) }
        return bs_gamma_q_inv_l(a, q)
    }
    
    /// Derivative d/dx P(a, x) for `Float80` (x86_64 only).
    ///
    /// See ``SpecialFunctions/regularizedGammaPDerivative(_:x:)`` for definition, domain, and discussion.
    @inlinable static func regularizedGammaPDerivative(_ a: Float80, x: Float80) throws -> Float80 {
        guard a.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a") }
        guard x.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "x") }
        guard a > 0 else { throw SpecialFunctionError.parameterNotPositive(name: "a") }
        guard x >= 0 else { throw SpecialFunctionError.parameterOutOfRange(name: "x", min: 0.0, max: Double.infinity) }
        return bs_gamma_p_derivative_l(a, x)
    }
    
#endif
    
}
