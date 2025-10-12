//
//  OwensT.swift
//  Math/SpecialFunctions
//
//  This file exposes Owen’s T function T(h, a) with type-generic and
//  type-specific overloads. The implementations are backed by the
//  CBoostBridge C shim (wrapping Boost.Math), while performing Swift-side
//  argument validation (finiteness checks).
//
//  All functions throw Swift errors for invalid inputs (e.g., non-finite
//  values). See each function’s documentation for details.
//

import CBoostBridge
public extension SpecialFunctions {
    
    
    // MARK: - Owen’s T (generic Double-backed)
    
    /// Compute Owen’s T function T(h, a).
    ///
    /// Definition:
    /// - Owen’s T is the integral
    ///   T(h, a) = (1 / 2π) ∫₀^{a} exp(-½ h² (1 + t²)) / (1 + t²) dt
    ///   It frequently appears in statistics (e.g., bivariate normal probabilities).
    ///
    /// Overview:
    /// - This generic overload accepts any `BinaryFloatingPoint` arguments and returns
    ///   a value of the same generic type `T`.
    /// - Internally, arguments are converted to `Double` and evaluated using the
    ///   C shim provided by `CBoostBridge` (wrapping Boost.Math), then converted
    ///   back to `T`.
    ///
    /// Domain and behavior:
    /// - T(h, a) is well-defined for all finite real h and a.
    /// - The function is smooth in both parameters and numerically stable for a wide range
    ///   of magnitudes when delegated to Boost.Math.
    /// - Very large |h| or |a| may push the result toward limiting values with potential
    ///   precision loss depending on `T`.
    ///
    /// Thread-safety:
    /// - This function is pure and thread-safe (no shared mutable state).
    ///
    /// Parameters:
    /// - h: The “height” parameter h.
    /// - a: The slope/ratio parameter a.
    ///
    /// Returns:
    /// - The value of Owen’s T(h, a) as `T`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotFinite(name: "h")` if `h` is NaN or ±∞.
    /// - `SpecialFunctionError.parameterNotFinite(name: "a")` if `a` is NaN or ±∞.
    ///
    /// Notes:
    /// - This overload converts inputs to `Double` for evaluation and then casts the result back to `T`.
    ///   If you require extended precision on x86_64 with `Float80`, prefer the dedicated `Float80`
    ///   overload below.
    ///
    /// Example:
    /// ```swift
    /// do {
    ///     let t: Double = try owensT(h: 0.5, a: 1.0) // ≈ 0.085...
    /// } catch {
    ///     // Handle invalid inputs
    /// }
    /// ```
    @inlinable public static func owensT<T: BinaryFloatingPoint>(h: T, a: T) throws -> T {
        // Convert to Double for the C backend. Keep a single conversion point to
        // minimize rounding steps.
        let dh = D(h), da = D(a)
        // Validate finiteness before calling into the C layer.
        guard dh.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "h") }
        guard da.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a") }
        // Delegate to Boost-backed implementation and convert back to T.
        return T(bs_owens_t(dh, da))
    }
    
    // MARK: - Float overloads
    
    /// Compute Owen’s T function T(h, a) for `Float`.
    ///
    /// This overload evaluates directly in `Float` precision via `CBoostBridge`.
    ///
    /// Parameters:
    /// - h: The “height” parameter h as `Float`.
    /// - a: The slope/ratio parameter a as `Float`.
    ///
    /// Returns:
    /// - The value of Owen’s T(h, a) as `Float`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotFinite(name: "h")` if `h` is NaN or ±∞.
    /// - `SpecialFunctionError.parameterNotFinite(name: "a")` if `a` is NaN or ±∞.
    ///
    /// Example:
    /// ```swift
    /// let tf = try owensT(h: 0.5 as Float, a: 1.0 as Float)
    /// ```
    @inlinable public static func owensT(h: Float, a: Float) throws -> Float {
        guard h.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "h") }
        guard a.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a") }
        return bs_owens_t_f(h, a)
    }
    
    // MARK: - Float80 overloads (x86_64 only)
    
#if arch(x86_64)
    /// Compute Owen’s T function T(h, a) for `Float80` (x86_64 only).
    ///
    /// This overload evaluates directly in extended precision (`Float80`) via `CBoostBridge`.
    /// Prefer this overload when you need more precision than `Double` on x86_64.
    ///
    /// Parameters:
    /// - h: The “height” parameter h as `Float80`.
    /// - a: The slope/ratio parameter a as `Float80`.
    ///
    /// Returns:
    /// - The value of Owen’s T(h, a) as `Float80`.
    ///
    /// Throws:
    /// - `SpecialFunctionError.parameterNotFinite(name: "h")` if `h` is NaN or ±∞.
    /// - `SpecialFunctionError.parameterNotFinite(name: "a")` if `a` is NaN or ±∞.
    ///
    /// Availability:
    /// - Only on x86_64 architectures where `Float80` is available.
    ///
    /// Example:
    /// ```swift
    /// #if arch(x86_64)
    /// let tl = try owensT(h: 0.5 as Float80, a: 1.0 as Float80)
    /// #endif
    /// ```
    @inlinable public static func owensT(h: Float80, a: Float80) throws -> Float80 {
        guard h.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "h") }
        guard a.isFinite else { throw SpecialFunctionError.parameterNotFinite(name: "a") }
        return bs_owens_t_l(h, a)
    }
#endif
    
}
