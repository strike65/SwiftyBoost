//
//  Created by Volker Thieme 2025.
//  Copyright © 2025 Volker Thieme.
//
//  Permission is hereby granted, free of charge, to any person obtaining a copy
//  of this software and associated documentation files (the "Software"), to deal
//  in the Software without restriction, including without limitation the rights
//  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//  copies of the Software, and to permit persons to whom the Software is
//  furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included in
//  all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
//  THE SOFTWARE.
//
import SwiftyBoostPrelude

// MARK: - Public API

/// Utilities to evaluate complex and real-valued spherical harmonics Y_l^m(θ, φ).
///
/// Overview:
/// - Provides overloads for `Double`, `Float`, and (on x86_64) `Float80`.
/// - Returns complex results as project-specific complex types: ``ComplexD``, ``ComplexF``, and ``ComplexL``.
/// - Convenience functions expose the real and imaginary parts directly.
///
/// Mathematical conventions:
/// - Uses the physics convention with the Condon–Shortley phase.
/// - Angles are in radians: θ ∈ [0, π] is the polar (colatitude) angle, φ ∈ [0, 2π) is the azimuth.
/// - Degree `l` is a non-negative integer and order `m` satisfies −l ≤ m ≤ l.
///
/// Implementation notes:
/// - Values are computed via a C bridge to Boost.Math (`boost::math::spherical_harmonic`) through `CBoostBridge`.
/// - The `Float` and `Float80` variants compute via `Double` and convert the result to the target precision.
///
/// SeeAlso:
/// - ``ComplexD``
/// - ``ComplexF``
/// - ``ComplexL``
public extension SpecialFunctions {

    /// Computes the complex spherical harmonic Y_l^m(θ, φ) for `Double`.
    ///
    /// Parameters:
    /// - l: The degree (n) with `l ≥ 0`.
    /// - m: The order with −l ≤ m ≤ l.
    /// - theta: The polar (colatitude) angle θ in radians, typically in [0, π].
    /// - phi: The azimuthal angle φ in radians, typically in [0, 2π).
    ///
    /// Returns:
    /// - The complex value Y_l^m(θ, φ) as a ``ComplexD``.
    ///
    /// Precondition:
    /// - `m` must be in the closed range −Int(l)...Int(l).
    ///
    /// Notes:
    /// - Follows the physics convention with the Condon–Shortley phase and Boost.Math normalization.
    @inlinable
    static func sphericalHarmonic(n: UInt, m: Int, theta: Double, phi: Double) -> ComplexD {
        guard m <= n else { return ComplexD(.zero, .zero)}
        let z = bs_spherical_harmonic_d(UInt32(n), Int32(m), theta, phi)
        return ComplexD(z.re, z.im)
    }

    /// Y_l^m(θ, φ) with mixed `Float`/`Double` angles; returns `ComplexD`.
    @inlinable
    static func sphericalHarmonic(n: UInt, m: Int, theta: Float, phi: Double) -> ComplexD {
        sphericalHarmonic(n: n, m: m, theta: Double(theta), phi: phi)
    }
    /// Y_l^m(θ, φ) with mixed `Double`/`Float` angles; returns `ComplexD`.
    @inlinable
    static func sphericalHarmonic(n: UInt, m: Int, theta: Double, phi: Float) -> ComplexD {
        sphericalHarmonic(n: n, m: m, theta: theta, phi: Double(phi))
    }

    /// Returns the real part of the spherical harmonic Y_l^m(θ, φ) for `Double`.
    ///
    /// Parameters:
    /// - l: The degree (n) with `l ≥ 0`.
    /// - m: The order with −l ≤ m ≤ l.
    /// - theta: The polar (colatitude) angle θ in radians, typically in [0, π].
    /// - phi: The azimuthal angle φ in radians, typically in [0, 2π).
    ///
    /// Returns:
    /// - The real component Re{Y_l^m(θ, φ)} as a `Double`.
    @inlinable
    static func sperhicalHarmonic_Real(n: UInt, m: Int, theta: Double, phi: Double) -> Double {
        return sphericalHarmonic(n: n, m: m, theta: theta, phi: phi).real
    }

    /// Real part with mixed `Float`/`Double` angles; returns `Double`.
    @inlinable
    static func sphericalHarmonic_Real(n: UInt, m: Int, theta: Float, phi: Double) -> Double {
        sphericalHarmonic(n: n, m: m, theta: Double(theta), phi: phi).real
    }
    /// Real part with mixed `Double`/`Float` angles; returns `Double`.
    @inlinable
    static func sphericalHarmonic_Real(n: UInt, m: Int, theta: Double, phi: Float) -> Double {
        sphericalHarmonic(n: n, m: m, theta: theta, phi: Double(phi)).real
    }

    /// Returns the imaginary part of the spherical harmonic Y_l^m(θ, φ) for `Double`.
    ///
    /// Parameters:
    /// - l: The degree (n) with `l ≥ 0`.
    /// - m: The order with −l ≤ m ≤ l.
    /// - theta: The polar (colatitude) angle θ in radians, typically in [0, π].
    /// - phi: The azimuthal angle φ in radians, typically in [0, 2π).
    ///
    /// Returns:
    /// - The imaginary component Im{Y_l^m(θ, φ)} as a `Double`.
    @inlinable
    static func sphericalHarmonic_Imag(n: UInt, m: Int, theta: Double, phi: Double) -> Double {
        return sphericalHarmonic(n: n, m: m, theta: theta, phi: phi).imag
    }

    /// Imag part with mixed `Float`/`Double` angles; returns `Double`.
    @inlinable
    static func sphericalHarmonic_Imag(n: UInt, m: Int, theta: Float, phi: Double) -> Double {
        sphericalHarmonic(n: n, m: m, theta: Double(theta), phi: phi).imag
    }
    /// Imag part with mixed `Double`/`Float` angles; returns `Double`.
    @inlinable
    static func sphericalHarmonic_Imag(n: UInt, m: Int, theta: Double, phi: Float) -> Double {
        sphericalHarmonic(n: n, m: m, theta: theta, phi: Double(phi)).imag
    }

    /// Computes the complex spherical harmonic Y_l^m(θ, φ) for `Float`.
    ///
    /// Discussion:
    /// - This overload computes in `Double` precision internally and converts the result to `Float`.
    ///
    /// Parameters:
    /// - l: The degree (n) with `l ≥ 0`.
    /// - m: The order with −l ≤ m ≤ l.
    /// - theta: The polar (colatitude) angle θ in radians, typically in [0, π].
    /// - phi: The azimuthal angle φ in radians, typically in [0, 2π).
    ///
    /// Returns:
    /// - The complex value Y_l^m(θ, φ) as a ``ComplexF``.
    @inlinable
    static func sphericalHarmonic(n: UInt, m: Int, theta: Float, phi: Float) -> ComplexF {
        let r = bs_spherical_harmonic_f(UInt32(n), Int32(m), theta, phi)
        return ComplexF(r.re, r.im)
    }

    /// Returns the real part of the spherical harmonic Y_l^m(θ, φ) for `Float`.
    ///
    /// Parameters:
    /// - l: The degree (n) with `l ≥ 0`.
    /// - m: The order with −l ≤ m ≤ l.
    /// - theta: The polar (colatitude) angle θ in radians, typically in [0, π].
    /// - phi: The azimuthal angle φ in radians, typically in [0, 2π).
    ///
    /// Returns:
    /// - The real component Re{Y_l^m(θ, φ)} as a `Float`.
    @inlinable
    static func Y_real(n: UInt, m: Int, theta: Float, phi: Float) -> Float {
        return sphericalHarmonic(n: n, m: m, theta: theta, phi: phi).real
    }

    /// Returns the imaginary part of the spherical harmonic Y_l^m(θ, φ) for `Float`.
    ///
    /// Parameters:
    /// - l: The degree (n) with `l ≥ 0`.
    /// - m: The order with −l ≤ m ≤ l.
    /// - theta: The polar (colatitude) angle θ in radians, typically in [0, π].
    /// - phi: The azimuthal angle φ in radians, typically in [0, 2π).
    ///
    /// Returns:
    /// - The imaginary component Im{Y_l^m(θ, φ)} as a `Float`.
    @inlinable
    static func Y_imag(n: UInt, m: Int, theta: Float, phi: Float) -> Float {
        return sphericalHarmonic(n: n, m: m, theta: theta, phi: phi).imag
    }

    #if arch(x86_64)
    /// Computes the complex spherical harmonic Y_l^m(θ, φ) for `Float80` (x86_64 only).
    ///
    /// Discussion:
    /// - This overload computes in `Double` precision internally and converts the result to `Float80`.
    ///
    /// Parameters:
    /// - l: The degree (n) with `l ≥ 0`.
    /// - m: The order with −l ≤ m ≤ l.
    /// - theta: The polar (colatitude) angle θ in radians, typically in [0, π].
    /// - phi: The azimuthal angle φ in radians, typically in [0, 2π).
    ///
    /// Returns:
    /// - The complex value Y_l^m(θ, φ) as a ``ComplexL``.
    @inlinable
    static func sphericalHarmonic(n: UInt, m: Int, theta: Float80, phi: Float80) -> ComplexL {
        let z = bs_spherical_harmonic_l(UInt32(n), Int32(m), theta, phi)
        return ComplexL(z.re, z.im)
    }

    /// Returns the real part of the spherical harmonic Y_l^m(θ, φ) for `Float80` (x86_64 only).
    ///
    /// Parameters:
    /// - l: The degree (n) with `l ≥ 0`.
    /// - m: The order with −l ≤ m ≤ l.
    /// - theta: The polar (colatitude) angle θ in radians, typically in [0, π].
    /// - phi: The azimuthal angle φ in radians, typically in [0, 2π).
    ///
    /// Returns:
    /// - The real component Re{Y_l^m(θ, φ)} as a `Float80`.
    @inlinable
    static func sphericalHarmonic_Real(n l: UInt, m: Int, theta: Float80, phi: Float80) -> Float80 {
        return sphericalHarmonic(n: l, m: m, theta: theta, phi: phi).real
    }


    /// Returns the imaginary part of the spherical harmonic Y_l^m(θ, φ) for `Float80` (x86_64 only).
    ///
    /// Parameters:
    /// - l: The degree (n) with `l ≥ 0`.
    /// - m: The order with −l ≤ m ≤ l.
    /// - theta: The polar (colatitude) angle θ in radians, typically in [0, π].
    /// - phi: The azimuthal angle φ in radians, typically in [0, 2π).
    ///
    /// Returns:
    /// - The imaginary component Im{Y_l^m(θ, φ)} as a `Float80`.
    @inlinable
    static func sphericalHarmonic_Imag(n: UInt, m: Int, theta: Float80, phi: Float80) -> Float80 {
        return sphericalHarmonic(n: n, m: m, theta: theta, phi: phi).imag
    }
    #endif
}
