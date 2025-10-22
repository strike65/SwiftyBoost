//
//  Created by Volker Thieme 2025.
//  Copyright © 2025 Volker Thieme.
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

import CBoostBridge

/// A namespace-like generic that exposes high-precision mathematical constants
/// specialized for concrete floating-point types.
///
/// Overview
/// - This type provides a set of well-known mathematical constants (π, e, √2, Euler–Mascheroni, Catalan, etc.)
///   and several derived variants (e.g. 2π, π/2, π/3, 2π/3, 3π/4, π/6, π², 1/π, 2/π, √π, 1/√π, ln(2), ln(10), ln(ln(2)), φ).
/// - The values are provided via Boost.Math and bridged from C using `CBoostBridge`.
/// - You access constants through the generic `Constants<T>` with `T` constrained to `BinaryFloatingPoint`.
///   Separate specializations exist for `Double`, `Float`, and (on x86_64) `Float80`.
///
/// Precision and performance
/// - Each specialization calls a dedicated C function that returns a constant computed or represented
///   in that precision. This avoids conversion overhead and preserves as much precision as the underlying
///   Boost.Math implementation provides for that type.
/// - Accessors are `static` computed properties with trivial cost.
///
/// Availability
/// - `Float80` specializations are available only on x86_64 architectures. The extension is conditionally
///   compiled with `#if arch(x86_64)`.
///
/// Usage
/// - let tau = Constants<Double>.twoPi
/// - let halfPiF = Constants<Float>.halfPi
/// - #if arch(x86_64)
///     let sqrtPi80 = Constants<Float80>.rootPi
///   #endif
///
/// Notes
/// - All values are pure constants with no side effects.
/// - Thread-safe by construction (no mutable state).

// MARK: - Double specializations

extension Constants where T == Double {

    /// The mathematical constant π (pi) ≈ 3.14159...
    ///
    /// Backed by: `bs_const_pi_d()`
    public static var pi: Double { bs_const_pi_d() }

    /// Euler’s number e ≈ 2.71828...
    ///
    /// Backed by: `bs_const_e_d()`
    public static var e: Double { bs_const_e_d() }

    /// 2π (tau), often used as a full rotation in radians.
    ///
    /// Backed by: `bs_const_two_pi_d()`
    public static var twoPi: Double { bs_const_two_pi_d() }

    /// π/2 (half-pi), i.e. 90 degrees in radians.
    ///
    /// Backed by: `bs_const_half_pi_d()`
    public static var halfPi: Double { bs_const_half_pi_d() }

    /// π/4 (quarter-pi), i.e. 45 degrees in radians.
    ///
    /// Backed by: `bs_const_quarter_pi_d()`
    public static var quarterPi: Double { bs_const_quarter_pi_d() }

    /// π/3 (third-pi), i.e. 60 degrees in radians.
    ///
    /// Backed by: `bs_const_third_pi_d()`
    public static var thirdPi: Double { bs_const_third_pi_d() }

    /// 2π/3 (two-thirds-pi), i.e. 120 degrees in radians.
    ///
    /// Backed by: `bs_const_two_thirds_pi_d()`
    public static var twoThirdsPi: Double {
        bs_const_two_thirds_pi_d()
    }

    /// 3π/4 (three-quarters-pi), i.e. 135 degrees in radians.
    ///
    /// Backed by: `bs_const_three_quarters_pi_d()`
    public static var threeQuartersPi: Double {
        bs_const_three_quarters_pi_d()
    }

    /// π/6 (sixth-pi), i.e. 30 degrees in radians.
    ///
    /// Backed by: `bs_const_sixth_pi_d()`
    public static var sixthPi: Double { bs_const_sixth_pi_d() }

    /// π² (pi squared).
    ///
    /// Backed by: `bs_const_pi_sqr_d()`
    public static var piSqr: Double { bs_const_pi_sqr_d() }

    /// √2 (square root of 2).
    ///
    /// Backed by: `bs_const_root_two_d()`
    public static var rootTwo: Double { bs_const_root_two_d() }

    /// √3 (square root of 3).
    ///
    /// Backed by: `bs_const_root_three_d()`
    public static var rootThree: Double { bs_const_root_three_d() }

    /// √π (square root of pi).
    ///
    /// Backed by: `bs_const_root_pi_d()`
    public static var rootPi: Double { bs_const_root_pi_d() }

    /// 1/π (reciprocal of pi).
    ///
    /// Backed by: `bs_const_one_div_pi_d()`
    public static var oneDivPi: Double { bs_const_one_div_pi_d() }

    /// 1/(2π) (one over two-pi).
    ///
    /// Backed by: `bs_const_one_div_two_pi_d()`
    public static var oneDivTwoPi: Double {
        bs_const_one_div_two_pi_d()
    }

    /// 1/√π (reciprocal of root-pi).
    ///
    /// Backed by: `bs_const_one_div_root_pi_d()`
    public static var oneDivRootPi: Double {
        bs_const_one_div_root_pi_d()
    }

    /// 2/π (two over pi).
    ///
    /// Backed by: `bs_const_two_div_pi_d()`
    public static var twoDivPi: Double { bs_const_two_div_pi_d() }

    /// 2/√π (two over root-pi).
    ///
    /// Backed by: `bs_const_two_div_root_pi_d()`
    public static var twoDivRootPi: Double {
        bs_const_two_div_root_pi_d()
    }

    /// ln(2) (natural logarithm of 2).
    ///
    /// Backed by: `bs_const_ln_two_d()`
    public static var lnTwo: Double { bs_const_ln_two_d() }

    /// ln(10) (natural logarithm of 10).
    ///
    /// Backed by: `bs_const_ln_ten_d()`
    public static var lnTen: Double { bs_const_ln_ten_d() }

    /// ln(ln(2)) (natural logarithm of ln(2)).
    ///
    /// Backed by: `bs_const_ln_ln_two_d()`
    public static var lnLnTwo: Double { bs_const_ln_ln_two_d() }

    /// Euler–Mascheroni constant γ ≈ 0.57721...
    ///
    /// Backed by: `bs_const_euler_d()`
    public static var euler: Double { bs_const_euler_d() }

    /// Catalan’s constant G ≈ 0.91596...
    ///
    /// Backed by: `bs_const_catalan_d()`
    public static var catalan: Double { bs_const_catalan_d() }

    /// Apery’s constant ζ(3) ≈ 1.2020569...
    ///
    /// Backed by: `bs_const_zeta_three_d()`
    public static var zetaThree: Double { bs_const_zeta_three_d() }

    /// Golden ratio φ = (1 + √5)/2 ≈ 1.61803...
    ///
    /// Backed by: `bs_const_phi_d()`
    public static var phi: Double { bs_const_phi_d() }
}

// MARK: - Float specializations

extension Constants where T == Float {

    /// The mathematical constant π (pi) as `Float`.
    ///
    /// Backed by: `bs_const_pi_f()`
    public static var pi: Float { bs_const_pi_f() }

    /// Euler’s number e as `Float`.
    ///
    /// Backed by: `bs_const_e_f()`
    public static var e: Float { bs_const_e_f() }

    /// 2π (tau) as `Float`.
    ///
    /// Backed by: `bs_const_two_pi_f()`
    public static var twoPi: Float { bs_const_two_pi_f() }

    /// π/2 (half-pi) as `Float`.
    ///
    /// Backed by: `bs_const_half_pi_f()`
    public static var halfPi: Float { bs_const_half_pi_f() }

    /// π/4 (quarter-pi) as `Float`.
    ///
    /// Backed by: `bs_const_quarter_pi_f()`
    public static var quarterPi: Float { bs_const_quarter_pi_f() }

    /// π/3 (third-pi) as `Float`.
    ///
    /// Backed by: `bs_const_third_pi_f()`
    public static var thirdPi: Float { bs_const_third_pi_f() }

    /// 2π/3 (two-thirds-pi) as `Float`.
    ///
    /// Backed by: `bs_const_two_thirds_pi_f()`
    public static var twoThirdsPi: Float { bs_const_two_thirds_pi_f() }

    /// 3π/4 (three-quarters-pi) as `Float`.
    ///
    /// Backed by: `bs_const_three_quarters_pi_f()`
    public static var threeQuartersPi: Float {
        bs_const_three_quarters_pi_f()
    }

    /// π/6 (sixth-pi) as `Float`.
    ///
    /// Backed by: `bs_const_sixth_pi_f()`
    public static var sixthPi: Float { bs_const_sixth_pi_f() }

    /// π² (pi squared) as `Float`.
    ///
    /// Backed by: `bs_const_pi_sqr_f()`
    public static var piSqr: Float { bs_const_pi_sqr_f() }

    /// √2 as `Float`.
    ///
    /// Backed by: `bs_const_root_two_f()`
    public static var rootTwo: Float { bs_const_root_two_f() }

    /// √3 as `Float`.
    ///
    /// Backed by: `bs_const_root_three_f()`
    public static var rootThree: Float { bs_const_root_three_f() }

    /// √π as `Float`.
    ///
    /// Backed by: `bs_const_root_pi_f()`
    public static var rootPi: Float { bs_const_root_pi_f() }

    /// 1/π as `Float`.
    ///
    /// Backed by: `bs_const_one_div_pi_f()`
    public static var oneDivPi: Float { bs_const_one_div_pi_f() }

    /// 1/(2π) as `Float`.
    ///
    /// Backed by: `bs_const_one_div_two_pi_f()`
    public static var oneDivTwoPi: Float {
        bs_const_one_div_two_pi_f()
    }

    /// 1/√π as `Float`.
    ///
    /// Backed by: `bs_const_one_div_root_pi_f()`
    public static var oneDivRootPi: Float {
        bs_const_one_div_root_pi_f()
    }

    /// 2/π as `Float`.
    ///
    /// Backed by: `bs_const_two_div_pi_f()`
    public static var twoDivPi: Float { bs_const_two_div_pi_f() }

    /// 2/√π as `Float`.
    ///
    /// Backed by: `bs_const_two_div_root_pi_f()`
    public static var twoDivRootPi: Float {
        bs_const_two_div_root_pi_f()
    }

    /// ln(2) as `Float`.
    ///
    /// Backed by: `bs_const_ln_two_f()`
    public static var lnTwo: Float { bs_const_ln_two_f() }

    /// ln(10) as `Float`.
    ///
    /// Backed by: `bs_const_ln_ten_f()`
    public static var lnTen: Float { bs_const_ln_ten_f() }

    /// ln(ln(2)) as `Float`.
    ///
    /// Backed by: `bs_const_ln_ln_two_f()`
    public static var lnLnTwo: Float { bs_const_ln_ln_two_f() }

    /// Euler–Mascheroni constant γ as `Float`.
    ///
    /// Backed by: `bs_const_euler_f()`
    public static var euler: Float { bs_const_euler_f() }

    /// Catalan’s constant G as `Float`.
    ///
    /// Backed by: `bs_const_catalan_f()`
    public static var catalan: Float { bs_const_catalan_f() }

    /// Apery’s constant ζ(3) as `Float`.
    ///
    /// Backed by: `bs_const_zeta_three_f()`
    public static var zetaThree: Float { bs_const_zeta_three_f() }

    /// Golden ratio φ as `Float`.
    ///
    /// Backed by: `bs_const_phi_f()`
    public static var phi: Float { bs_const_phi_f() }
}

// MARK: - Float80 specializations (x86_64 only)

#if arch(x86_64)
    extension Constants where T == Float80 {

        /// The mathematical constant π (pi) as `Float80` (x86_64 only).
        ///
        /// Backed by: `bs_const_pi_l()`
        public static var pi: Float80 {
            bs_const_pi_l()
        }

        /// Euler’s number e as `Float80` (x86_64 only).
        ///
        /// Backed by: `bs_const_e_l()`
        public static var e: Float80 {
            bs_const_e_l()
        }

        /// 2π (tau) as `Float80` (x86_64 only).
        ///
        /// Backed by: `bs_const_two_pi_l()`
        public static var twoPi: Float80 {
            bs_const_two_pi_l()
        }

        /// π/2 (half-pi) as `Float80` (x86_64 only).
        ///
        /// Backed by: `bs_const_half_pi_l()`
        public static var halfPi: Float80 { bs_const_half_pi_l() }

        /// π/4 (quarter-pi) as `Float80` (x86_64 only).
        ///
        /// Backed by: `bs_const_quarter_pi_l()`
        public static var quarterPi: Float80 { bs_const_quarter_pi_l() }

        /// π/3 (third-pi) as `Float80` (x86_64 only).
        ///
        /// Backed by: `bs_const_third_pi_l()`
        public static var thirdPi: Float80 { bs_const_third_pi_l() }

        /// 2π/3 (two-thirds-pi) as `Float80` (x86_64 only).
        ///
        /// Backed by: `bs_const_two_thirds_pi_l()`
        public static var twoThirdsPi: Float80 { bs_const_two_thirds_pi_l() }

        /// 3π/4 (three-quarters-pi) as `Float80` (x86_64 only).
        ///
        /// Backed by: `bs_const_three_quarters_pi_l()`
        public static var threeQuartersPi: Float80 { bs_const_three_quarters_pi_l() }

        /// π/6 (sixth-pi) as `Float80` (x86_64 only).
        ///
        /// Backed by: `bs_const_sixth_pi_l()`
        public static var sixthPi: Float80 { bs_const_sixth_pi_l() }

        /// π² (pi squared) as `Float80` (x86_64 only).
        ///
        /// Backed by: `bs_const_pi_sqr_l()`
        public static var piSqr: Float80 {
            bs_const_pi_sqr_l()
        }

        /// √2 as `Float80` (x86_64 only).
        ///
        /// Backed by: `bs_const_root_two_l()`
        public static var rootTwo: Float80 { bs_const_root_two_l() }

        /// √3 as `Float80` (x86_64 only).
        ///
        /// Backed by: `bs_const_root_three_l()`
        public static var rootThree: Float80 { bs_const_root_three_l() }

        /// √π as `Float80` (x86_64 only).
        ///
        /// Backed by: `bs_const_root_pi_l()`
        public static var rootPi: Float80 { bs_const_root_pi_l() }

        /// 1/π as `Float80` (x86_64 only).
        ///
        /// Backed by: `bs_const_one_div_pi_l()`
        public static var oneDivPi: Float80 { bs_const_one_div_pi_l() }

        /// 1/(2π) as `Float80` (x86_64 only).
        ///
        /// Backed by: `bs_const_one_div_two_pi_l()`
        public static var oneDivTwoPi: Float80 { bs_const_one_div_two_pi_l() }

        /// 1/√π as `Float80` (x86_64 only).
        ///
        /// Backed by: `bs_const_one_div_root_pi_l()`
        public static var oneDivRootPi: Float80 { bs_const_one_div_root_pi_l() }

        /// 2/π as `Float80` (x86_64 only).
        ///
        /// Backed by: `bs_const_two_div_pi_l()`
        public static var twoDivPi: Float80 { bs_const_two_div_pi_l() }

        /// 2/√π as `Float80` (x86_64 only).
        ///
        /// Backed by: `bs_const_two_div_root_pi_l()`
        public static var twoDivRootPi: Float80 { bs_const_two_div_root_pi_l() }

        /// ln(2) as `Float80` (x86_64 only).
        ///
        /// Backed by: `bs_const_ln_two_l()`
        public static var lnTwo: Float80 {
            bs_const_ln_two_l()
        }

        /// ln(10) as `Float80` (x86_64 only).
        ///
        /// Backed by: `bs_const_ln_ten_l()`
        public static var lnTen: Float80 {
            bs_const_ln_ten_l()
        }

        /// ln(ln(2)) as `Float80` (x86_64 only).
        ///
        /// Backed by: `bs_const_ln_ln_two_l()`
        public static var lnLnTwo: Float80 { bs_const_ln_ln_two_l() }

        /// Euler–Mascheroni constant γ as `Float80` (x86_64 only).
        ///
        /// Backed by: `bs_const_euler_l()`
        public static var euler: Float80 {
            bs_const_euler_l()
        }

        /// Catalan’s constant G as `Float80` (x86_64 only).
        ///
        /// Backed by: `bs_const_catalan_l()`
        public static var catalan: Float80 { bs_const_catalan_l() }

        /// Apery’s constant ζ(3) as `Float80` (x86_64 only).
        ///
        /// Backed by: `bs_const_zeta_three_l()`
        public static var zetaThree: Float80 { bs_const_zeta_three_l() }

        /// Golden ratio φ as `Float80` (x86_64 only).
        ///
        /// Backed by: `bs_const_phi_l()`
        public static var phi: Float80 {
            bs_const_phi_l()
        }
    }
#endif
