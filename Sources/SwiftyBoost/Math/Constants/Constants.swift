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

import SwiftyBoostPrelude

/// High-precision mathematical constants with on-demand conversion to floating-point.
///
/// - Source: Values are taken from Boost.Math big constant tables and stored as decimal strings,
///   so they can be converted at runtime to any `Real & BinaryFloatingPoint` type (`Float`,
///   `Double`, `Float80`, or user-defined types that conform).
/// - Precision: Strings typically contain more digits than any standard binary float can represent
///   to minimize conversion error.
/// - Performance: Conversion prefers `LosslessStringConvertible` when available on the target type,
///   falling back to a lightweight custom parser for general `BinaryFloatingPoint`.
///
/// Usage examples:
/// - `let piAsDouble: Double = Constants.pi()`
/// - `let eAsFloat: Float = Constants.e()`
/// - `let tauAsFloat80: Float80 = Constants.twoPi()`
/// - `let gammaAsDouble = Constants.value(.euler, as: Double.self)`
///
/// Mathematical insight:
/// - Angles: Many constants are multiples/fractions of π and arise in trigonometry and Fourier analysis.
/// - Growth/decay: e and its reciprocal are fundamental to exponentials, logarithms, and probability.
/// - Special numbers: Euler–Mascheroni, Catalan, Apery’s ζ(3), and φ (golden ratio) appear in number
///   theory, combinatorics, geometry, and analysis.
public enum Constants {

    /// Enumerates the constant catalogue.
    ///
    /// For each case you’ll find:
    /// - Symbol: Common mathematical symbol
    /// - Definition: Canonical definition
    /// - Appears in: Typical areas where the constant naturally occurs
    public enum Identifier: CaseIterable, Sendable {
        /// π (pi)
        /// - Definition: The ratio of a circle’s circumference to its diameter.
        /// - Appears in: Geometry, trigonometry, Fourier analysis, probability (normal distribution),
        ///   complex analysis (Euler’s identity).
        case pi

        /// 2π (tau)
        /// - Definition: The full angle in radians for one revolution; 2π = circumference/radius.
        /// - Appears in: Periodic phenomena, Fourier transforms, frequency-domain analysis.
        case twoPi

        /// π/2 (half pi)
        /// - Definition: Right angle in radians; principal quarter-period for many trig functions.
        /// - Appears in: Trigonometric identities, orthogonality of sin/cos at quadrature points.
        case halfPi

        /// π/4 (quarter pi)
        /// - Definition: 45° in radians; tan(π/4) = 1; sin(π/4) = cos(π/4) = √2/2.
        /// - Appears in: Rotations, symmetry arguments, Gaussian integrals.
        case quarterPi

        /// π/3 (third pi)
        /// - Definition: 60° in radians; cos(π/3) = 1/2; sin(π/3) = √3/2.
        /// - Appears in: Equilateral triangles, hexagonal lattices.
        case thirdPi

        /// 2π/3 (two thirds pi)
        /// - Definition: 120° in radians; key angle in cubic roots of unity geometry.
        /// - Appears in: Symmetries of regular polygons, phasors.
        case twoThirdsPi

        /// 3π/4 (three quarters pi)
        /// - Definition: 135° in radians; sin and cos have equal magnitude with opposite signs.
        /// - Appears in: Rotational symmetries, signal phase shifts.
        case threeQuartersPi

        /// π/6 (sixth pi)
        /// - Definition: 30° in radians; sin(π/6) = 1/2; cos(π/6) = √3/2.
        /// - Appears in: Triangulations, classical trig problems.
        case sixthPi

        /// π² (pi squared)
        /// - Definition: Square of π.
        /// - Appears in: Basel problem (ζ(2) = π²/6), Fourier series, integrals over circles/spheres.
        case piSqr

        /// √2 (square root of two)
        /// - Definition: Hypotenuse of a unit right triangle; first known irrational.
        /// - Appears in: Rotations by 45°, normalization, lattice geometry.
        case rootTwo

        /// √3 (square root of three)
        /// - Definition: Height of an equilateral triangle of side 2.
        /// - Appears in: Hexagonal close packing, 60° rotations, power systems (√3 factors).
        case rootThree

        /// √π (square root of pi)
        /// - Definition: Principal square root of π.
        /// - Appears in: Gaussian integrals (∫ e^{-x²} dx = √π), normal distribution normalization.
        case rootPi

        /// 1/π (reciprocal of pi)
        /// - Definition: Reciprocal of π.
        /// - Appears in: Normalization constants, series for 1/π (Ramanujan-type).
        case oneDivPi

        /// 1/(2π) (reciprocal of two pi)
        /// - Definition: 1 divided by 2π.
        /// - Appears in: Fourier transform conventions, spectral densities.
        case oneDivTwoPi

        /// 1/√π (reciprocal of square root of pi)
        /// - Definition: Reciprocal of √π.
        /// - Appears in: Normalization of Gaussian kernels, error function derivatives.
        case oneDivRootPi

        /// 2/π (two over pi)
        /// - Definition: 2 divided by π.
        /// - Appears in: Fourier series of square waves, sinc integrals.
        case twoDivPi

        /// 2/√π (two over square root of pi)
        /// - Definition: 2 divided by √π.
        /// - Appears in: Derivative of the error function: d/dx erf(x) = (2/√π) e^{-x²}.
        case twoDivRootPi

        /// ln 2 (natural logarithm of two)
        /// - Definition: The natural log of 2.
        /// - Appears in: Information theory (bits ↔ nats), doubling times, binary scaling.
        case lnTwo

        /// ln 10 (natural logarithm of ten)
        /// - Definition: The natural log of 10.
        /// - Appears in: Base conversion between natural logs and common logs (log10).
        case lnTen

        /// ln ln 2 (natural logarithm of ln 2)
        /// - Definition: The natural log of ln 2.
        /// - Appears in: Asymptotics, iterated logarithms in analysis and number theory.
        case lnLnTwo

        /// e (Euler’s number)
        /// - Definition: Unique base where the derivative of a^x equals a^x at x=0; limit (1+1/n)^n.
        /// - Appears in: Exponential growth/decay, calculus, probability (Poisson, normal).
        case e

        /// 1/e (reciprocal of Euler’s number)
        /// - Definition: Reciprocal of e.
        /// - Appears in: Memoryless processes, exponential waiting times, optimization heuristics.
        case oneDivE

        /// γ (Euler–Mascheroni constant)
        /// - Definition: Limit of (H_n − ln n) as n → ∞ where H_n is the nth harmonic number.
        /// - Appears in: Analytic number theory, integrals and series involving logs and harmonic sums.
        case euler

        /// G (Catalan’s constant)
        /// - Definition: Σ_{n≥0} (-1)^n / (2n+1)^2.
        /// - Appears in: Combinatorics, geometry of lattice sums, integrals of inverse trigs.
        case catalan

        /// ζ(3) (Apéry’s constant)
        /// - Definition: Riemann zeta function at 3: Σ_{n≥1} 1/n^3.
        /// - Appears in: Quantum field theory, knot theory, multiple zeta values; irrational (Apéry).
        case zetaThree

        /// φ (golden ratio)
        /// - Definition: (1 + √5)/2; unique positive root of x^2 = x + 1.
        /// - Appears in: Recursive sequences (Fibonacci), geometry, continued fractions.
        case phi

        /// h (Holtsmark distribution differential entropy)
        /// - Definition: Differential entropy of the Holtsmark distribution (a stable law with α=3/2).
        /// - Appears in: Astrophysics (gravitational fields of random star distributions), stable laws.
        /// - Source of the numeric value: ``https://github.com/tk-yoshimura/HoltsmarkDistribution``
        case holtsmarkEntropy
    }

    /// Stores the raw decimal string for each constant.
    /// - Note: Strings are intentionally long to minimize rounding when converting to binary formats.
    private static let catalogue: [Identifier: String] = [
        .pi: "3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651e+00",
        .twoPi: "6.28318530717958647692528676655900576839433879875021164194988918461563281257241799725606965068423413596429617303e+00",
        .halfPi: "1.57079632679489661923132169163975144209858469968755291048747229615390820314310449931401741267105853399107404326e+00",
        .quarterPi: "0.785398163397448309615660845819875721049292349843776455243736148076954101571552249657008706335529266995537021628320576661773",
        .thirdPi: "1.04719755119659774615421446109316762806572313312503527365831486410260546876206966620934494178070568932738269550e+00",
        .twoThirdsPi: "2.09439510239319549230842892218633525613144626625007054731662972820521093752413933241868988356141137865476539101e+00",
        .threeQuartersPi: "2.35619449019234492884698253745962716314787704953132936573120844423086230471465674897102611900658780098661106488e+00",
        .sixthPi: "5.23598775598298873077107230546583814032861566562517636829157432051302734381034833104672470890352844663691347752e-01",
        .piSqr: "9.86960440108935861883449099987615113531369940724079062641334937622004482241920524300177340371855223182402591377e+00",
        .rootTwo: "1.41421356237309504880168872420969807856967187537694807317667973799073247846210703885038753432764157273501384623e+00",
        .rootThree: "1.73205080756887729352744634150587236694280525381038062805580697945193301690880003708114618675724857567562614142e+00",
        .rootPi: "1.77245385090551602729816748334114518279754945612238712821380778985291128459103218137495065673854466541622682362e+00",
        .oneDivPi: "0.31830988618379067153776752674502872406891929148091289749533468811779359526845307018022760553250617191214568545351",
        .oneDivTwoPi: "1.59154943091895335768883763372514362034459645740456448747667344058896797634226535090113802766253085956072842727e-01",
        .oneDivRootPi: "5.64189583547756286948079451560772585844050629328998856844085721710642468441493414486743660202107363443028347906e-01",
        .twoDivPi: "6.36619772367581343075535053490057448137838582961825794990669376235587190536906140360455211065012343824291370907e-01",
        .twoDivRootPi: "1.12837916709551257389615890312154517168810125865799771368817144342128493688298682897348732040421472688605669581272",
        .lnTwo: "6.93147180559945309417232121458176568075500134360255254120680009493393621969694715605863326996418687542001481021e-01",
        .lnTen: "2.30258509299404568401799145468436420760110148862877297603332790096757260967735248023599720508959829834196778404e+00",
        .lnLnTwo: "-3.66512920581664327012439158232669469454263447837105263053677713670561615319352738549455822856698908358302523045e-01",
        .e: "2.71828182845904523536028747135266249775724709369995957496696762772407663035354759457138217852516642742746639193e+00",
        .oneDivE:  "3.67879441171442321595523770161460867445811131031767834507836801697461495744899803357147274345919643746627325277e-01",
        .euler: "5.77215664901532860606512090082402431042159335939923598805767234884867726777664670936947063291746749514631447250e-01",
        .catalan: "9.15965594177219015054603514932384110774149374281672134266498119621763019776254769479356512926115106248574422619e-01",
        .zetaThree: "1.20205690315959428539973816151144999076498629234049888179227155534183820578631309018645587360933525814619915780e+00",
        .phi: "1.61803398874989484820458683436563811772030917980576286213544862270526046281890244970720720418939113748475408808e+00",
        .holtsmarkEntropy: "2.06944850513462440031558003845421663807675255168165702483694991535648437010785015575051452878363155526878657629"
    ]

    // MARK: - Public API

    /// Returns the requested constant converted to the supplied floating-point type.
    ///
    /// - Parameters:
    ///   - identifier: The constant to fetch from the catalogue.
    ///   - type: The target floating-point type to convert into.
    /// - Returns: The value as `T`.
    /// - Notes:
    ///   - If `T` conforms to `LosslessStringConvertible`, parsing delegates to the type’s initializer.
    ///   - Otherwise, a small numeric parser reads the decimal string (including scientific notation)
    ///     and constructs a `T` via repeated scaling by powers of ten.
    public static func value<T: Real & BinaryFloatingPoint>(
        _ identifier: Identifier,
        as type: T.Type = T.self
    ) -> T {
        guard let string = catalogue[identifier] else {
            preconditionFailure("Unsupported constant identifier: \(identifier)")
        }
        return convert(decimal: string, as: type)
    }

    // Convenience overloads mirror the existing static constants.

    /// π (pi). See `Identifier.pi` for mathematical context.
    public static func pi<T: Real & BinaryFloatingPoint>(_ type: T.Type = T.self) -> T { value(.pi, as: type) }

    /// 2π (tau). See `Identifier.twoPi` for mathematical context.
    public static func twoPi<T: Real & BinaryFloatingPoint>(_ type: T.Type = T.self) -> T { value(.twoPi, as: type) }

    /// π/2. See `Identifier.halfPi` for mathematical context.
    public static func halfPi<T: Real & BinaryFloatingPoint>(_ type: T.Type = T.self) -> T { value(.halfPi, as: type) }

    /// π/4. See `Identifier.quarterPi` for mathematical context.
    public static func quarterPi<T: Real & BinaryFloatingPoint>(_ type: T.Type = T.self) -> T { value(.quarterPi, as: type) }

    /// π/3. See `Identifier.thirdPi` for mathematical context.
    public static func thirdPi<T: Real & BinaryFloatingPoint>(_ type: T.Type = T.self) -> T { value(.thirdPi, as: type) }

    /// 2π/3. See `Identifier.twoThirdsPi` for mathematical context.
    public static func twoThirdsPi<T: Real & BinaryFloatingPoint>(_ type: T.Type = T.self) -> T { value(.twoThirdsPi, as: type) }

    /// 3π/4. See `Identifier.threeQuartersPi` for mathematical context.
    public static func threeQuartersPi<T: Real & BinaryFloatingPoint>(_ type: T.Type = T.self) -> T { value(.threeQuartersPi, as: type) }

    /// π/6. See `Identifier.sixthPi` for mathematical context.
    public static func sixthPi<T: Real & BinaryFloatingPoint>(_ type: T.Type = T.self) -> T { value(.sixthPi, as: type) }

    /// π². See `Identifier.piSqr` for mathematical context.
    public static func piSqr<T: Real & BinaryFloatingPoint>(_ type: T.Type = T.self) -> T { value(.piSqr, as: type) }

    /// √2. See `Identifier.rootTwo` for mathematical context.
    public static func rootTwo<T: Real & BinaryFloatingPoint>(_ type: T.Type = T.self) -> T { value(.rootTwo, as: type) }

    /// √3. See `Identifier.rootThree` for mathematical context.
    public static func rootThree<T: Real & BinaryFloatingPoint>(_ type: T.Type = T.self) -> T { value(.rootThree, as: type) }

    /// √π. See `Identifier.rootPi` for mathematical context.
    public static func rootPi<T: Real & BinaryFloatingPoint>(_ type: T.Type = T.self) -> T { value(.rootPi, as: type) }

    /// 1/π. See `Identifier.oneDivPi` for mathematical context.
    public static func oneDivPi<T: Real & BinaryFloatingPoint>(_ type: T.Type = T.self) -> T { value(.oneDivPi, as: type) }

    /// 1/(2π). See `Identifier.oneDivTwoPi` for mathematical context.
    public static func oneDivTwoPi<T: Real & BinaryFloatingPoint>(_ type: T.Type = T.self) -> T { value(.oneDivTwoPi, as: type) }

    /// 1/√π. See `Identifier.oneDivRootPi` for mathematical context.
    public static func oneDivRootPi<T: Real & BinaryFloatingPoint>(_ type: T.Type = T.self) -> T { value(.oneDivRootPi, as: type) }

    /// 2/π. See `Identifier.twoDivPi` for mathematical context.
    public static func twoDivPi<T: Real & BinaryFloatingPoint>(_ type: T.Type = T.self) -> T { value(.twoDivPi, as: type) }

    /// 2/√π. See `Identifier.twoDivRootPi` for mathematical context.
    public static func twoDivRootPi<T: Real & BinaryFloatingPoint>(_ type: T.Type = T.self) -> T { value(.twoDivRootPi, as: type) }

    /// ln 2. See `Identifier.lnTwo` for mathematical context.
    public static func lnTwo<T: Real & BinaryFloatingPoint>(_ type: T.Type = T.self) -> T { value(.lnTwo, as: type) }

    /// ln 10. See `Identifier.lnTen` for mathematical context.
    public static func lnTen<T: Real & BinaryFloatingPoint>(_ type: T.Type = T.self) -> T { value(.lnTen, as: type) }

    /// ln ln 2. See `Identifier.lnLnTwo` for mathematical context.
    public static func lnLnTwo<T: Real & BinaryFloatingPoint>(_ type: T.Type = T.self) -> T { value(.lnLnTwo, as: type) }

    /// e. See `Identifier.e` for mathematical context.
    public static func e<T: Real & BinaryFloatingPoint>(_ type: T.Type = T.self) -> T { value(.e, as: type) }

    /// 1/e. See `Identifier.oneDivE` for mathematical context.
    public static func oneDivE<T: Real & BinaryFloatingPoint>(_ type: T.Type = T.self) -> T { value(.oneDivE, as: type) }

    /// Euler–Mascheroni constant γ. See `Identifier.euler` for mathematical context.
    public static func euler<T: Real & BinaryFloatingPoint>(_ type: T.Type = T.self) -> T { value(.euler, as: type) }

    /// Catalan’s constant G. See `Identifier.catalan` for mathematical context.
    public static func catalan<T: Real & BinaryFloatingPoint>(_ type: T.Type = T.self) -> T { value(.catalan, as: type) }

    /// Apéry’s constant ζ(3). See `Identifier.zetaThree` for mathematical context.
    public static func zetaThree<T: Real & BinaryFloatingPoint>(_ type: T.Type = T.self) -> T { value(.zetaThree, as: type) }

    /// Golden ratio φ. See `Identifier.phi` for mathematical context.
    public static func phi<T: Real & BinaryFloatingPoint>(_ type: T.Type = T.self) -> T { value(.phi, as: type) }

    /// Holtsmark distribution differential entropy. See `Identifier.holtsmarkEntropy` for context.
    public static func holtsmarkEntropy<T: Real & BinaryFloatingPoint>(_ type: T.Type = T.self) -> T { value(.holtsmarkEntropy, as: type) }

    // MARK: - Conversion helpers

    /// Converts a decimal string (optionally in scientific notation) into a `BinaryFloatingPoint` value.
    ///
    /// - Parameters:
    ///   - decimal: A decimal literal, possibly containing a leading sign, a decimal point, and an
    ///     exponent part (`e` or `E`).
    ///   - as: The destination floating-point type.
    /// - Returns: The parsed value as `T`.
    /// - Parsing strategy:
    ///   - Try `LosslessStringConvertible.init` when available on `T` for fast, locale-independent parsing.
    ///   - Otherwise, parse digits and exponent manually and construct the value using scaling by powers of 10.
    /// - Preconditions:
    ///   - `decimal` must be a valid decimal literal; invalid characters trigger a precondition failure.
    static func convert<T: Real & BinaryFloatingPoint>(
        decimal: String,
        as _: T.Type
    ) -> T {
        let trimmed = decimal.trimmingCharacters(in: .whitespacesAndNewlines)
        precondition(!trimmed.isEmpty, "Empty decimal constant string.")

        if let convertibleType = T.self as? LosslessStringConvertible.Type,
            let parsed = convertibleType.init(trimmed) as? T
        {
            return parsed
        }

        var index = trimmed.startIndex
        let isNegative: Bool
        if trimmed[index] == "-" {
            isNegative = true
            index = trimmed.index(after: index)
        } else if trimmed[index] == "+" {
            isNegative = false
            index = trimmed.index(after: index)
        } else {
            isNegative = false
        }

        var exponent = 0
        var digits: [UInt8] = []
        var seenDecimalPoint = false

        while index < trimmed.endIndex {
            let character = trimmed[index]
            if character == "e" || character == "E" {
                let expStart = trimmed.index(after: index)
                let exponentSlice = trimmed[expStart..<trimmed.endIndex]
                guard let parsedExponent = Int(exponentSlice) else {
                    preconditionFailure("Invalid exponent in decimal literal: \(decimal)")
                }
                exponent += parsedExponent
                break
            } else if character == "." {
                seenDecimalPoint = true
            } else if let value = character.wholeNumberValue {
                digits.append(UInt8(value))
                if seenDecimalPoint {
                    exponent -= 1
                }
            } else {
                preconditionFailure("Unsupported character in decimal literal: \(character)")
            }

            index = trimmed.index(after: index)
        }

        guard !digits.isEmpty else { return .zero }

        var firstNonZero = 0
        while firstNonZero < digits.count && digits[firstNonZero] == 0 {
            firstNonZero += 1
        }

        if firstNonZero == digits.count {
            return .zero
        }

        var magnitude: T = .zero
        let ten: T = 10
        for digit in digits[firstNonZero...] {
            magnitude = magnitude * ten + T(digit)
        }

        if exponent > 0 {
            magnitude *= powTen(of: exponent)
        } else if exponent < 0 {
            magnitude /= powTen(of: -exponent)
        }

        return isNegative ? -magnitude : magnitude
    }

    /// Computes 10^exponent for `BinaryFloatingPoint` types using exponentiation by squaring.
    ///
    /// - Parameter exponent: Non-negative integer exponent.
    /// - Returns: 10 raised to `exponent` as `T`.
    /// - Complexity: O(log exponent) multiplications.
    static func powTen<T: BinaryFloatingPoint>(of exponent: Int) -> T {
        precondition(exponent >= 0, "Exponent must be non-negative.")
        if exponent == 0 { return 1 }

        var result: T = 1
        var base: T = 10
        var remaining = exponent

        while remaining > 0 {
            if remaining & 1 == 1 {
                result *= base
            }
            remaining >>= 1
            if remaining > 0 {
                base *= base
            }
        }
        return result
    }
}
