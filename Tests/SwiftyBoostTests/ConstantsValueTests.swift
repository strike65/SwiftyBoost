import Testing
@testable import SwiftyBoost

private let matlabCatalogue: [Constants.Identifier: String] = [
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
    .oneDivE: "3.67879441171442321595523770161460867445811131031767834507836801697461495744899803357147274345919643746627325277e-01",
    .euler: "5.77215664901532860606512090082402431042159335939923598805767234884867726777664670936947063291746749514631447250e-01",
    .catalan: "9.15965594177219015054603514932384110774149374281672134266498119621763019776254769479356512926115106248574422619e-01",
    .zetaThree: "1.20205690315959428539973816151144999076498629234049888179227155534183820578631309018645587360933525814619915780e+00",
    .phi: "1.61803398874989484820458683436563811772030917980576286213544862270526046281890244970720720418939113748475408808e+00",
    .holtsmarkEntropy: "2.06944850513462440031558003845421663807675255168165702483694991535648437010785015575051452878363155526878657629e+00"
]

private func approxEqual<T: BinaryFloatingPoint>(_ lhs: T, _ rhs: T, tol: T) -> Bool {
    (lhs - rhs).magnitude <= tol
}

private func matlabValue<T: BinaryFloatingPoint>(
    _ identifier: Constants.Identifier,
    as type: T.Type
) -> T {
    guard let decimal = matlabCatalogue[identifier] else {
        fatalError("Missing MATLAB constant for \(identifier)")
    }
    return parseMatlabDecimal(decimal, as: type)
}

private func parseMatlabDecimal<T: BinaryFloatingPoint>(
    _ decimal: String,
    as _: T.Type
) -> T {
    let trimmed = decimal.trimmingCharacters(in: .whitespacesAndNewlines)
    precondition(!trimmed.isEmpty, "Empty decimal literal.")

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

private func powTen<T: BinaryFloatingPoint>(of exponent: Int) -> T {
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

private func convenienceValueDouble(_ identifier: Constants.Identifier) -> Double {
    switch identifier {
    case .pi: return Constants.pi()
    case .twoPi: return Constants.twoPi()
    case .halfPi: return Constants.halfPi()
    case .quarterPi: return Constants.quarterPi()
    case .thirdPi: return Constants.thirdPi()
    case .twoThirdsPi: return Constants.twoThirdsPi()
    case .threeQuartersPi: return Constants.threeQuartersPi()
    case .sixthPi: return Constants.sixthPi()
    case .piSqr: return Constants.piSqr()
    case .rootTwo: return Constants.rootTwo()
    case .rootThree: return Constants.rootThree()
    case .rootPi: return Constants.rootPi()
    case .oneDivPi: return Constants.oneDivPi()
    case .oneDivTwoPi: return Constants.oneDivTwoPi()
    case .oneDivRootPi: return Constants.oneDivRootPi()
    case .twoDivPi: return Constants.twoDivPi()
    case .twoDivRootPi: return Constants.twoDivRootPi()
    case .lnTwo: return Constants.lnTwo()
    case .lnTen: return Constants.lnTen()
    case .lnLnTwo: return Constants.lnLnTwo()
    case .e: return Constants.e()
    case .oneDivE: return Constants.oneDivE()
    case .euler: return Constants.euler()
    case .catalan: return Constants.catalan()
    case .zetaThree: return Constants.zetaThree()
    case .phi: return Constants.phi()
    case .holtsmarkEntropy: return Constants.holtsmarkEntropy()
    }
}

private func convenienceValueFloat(_ identifier: Constants.Identifier) -> Float {
    switch identifier {
    case .pi: return Constants.pi(Float.self)
    case .twoPi: return Constants.twoPi(Float.self)
    case .halfPi: return Constants.halfPi(Float.self)
    case .quarterPi: return Constants.quarterPi(Float.self)
    case .thirdPi: return Constants.thirdPi(Float.self)
    case .twoThirdsPi: return Constants.twoThirdsPi(Float.self)
    case .threeQuartersPi: return Constants.threeQuartersPi(Float.self)
    case .sixthPi: return Constants.sixthPi(Float.self)
    case .piSqr: return Constants.piSqr(Float.self)
    case .rootTwo: return Constants.rootTwo(Float.self)
    case .rootThree: return Constants.rootThree(Float.self)
    case .rootPi: return Constants.rootPi(Float.self)
    case .oneDivPi: return Constants.oneDivPi(Float.self)
    case .oneDivTwoPi: return Constants.oneDivTwoPi(Float.self)
    case .oneDivRootPi: return Constants.oneDivRootPi(Float.self)
    case .twoDivPi: return Constants.twoDivPi(Float.self)
    case .twoDivRootPi: return Constants.twoDivRootPi(Float.self)
    case .lnTwo: return Constants.lnTwo(Float.self)
    case .lnTen: return Constants.lnTen(Float.self)
    case .lnLnTwo: return Constants.lnLnTwo(Float.self)
    case .e: return Constants.e(Float.self)
    case .oneDivE: return Constants.oneDivE(Float.self)
    case .euler: return Constants.euler(Float.self)
    case .catalan: return Constants.catalan(Float.self)
    case .zetaThree: return Constants.zetaThree(Float.self)
    case .phi: return Constants.phi(Float.self)
    case .holtsmarkEntropy: return Constants.holtsmarkEntropy(Float.self)
    }
}

#if arch(x86_64)
private func convenienceValueFloat80(_ identifier: Constants.Identifier) -> Float80 {
    switch identifier {
    case .pi: return Constants.pi(Float80.self)
    case .twoPi: return Constants.twoPi(Float80.self)
    case .halfPi: return Constants.halfPi(Float80.self)
    case .quarterPi: return Constants.quarterPi(Float80.self)
    case .thirdPi: return Constants.thirdPi(Float80.self)
    case .twoThirdsPi: return Constants.twoThirdsPi(Float80.self)
    case .threeQuartersPi: return Constants.threeQuartersPi(Float80.self)
    case .sixthPi: return Constants.sixthPi(Float80.self)
    case .piSqr: return Constants.piSqr(Float80.self)
    case .rootTwo: return Constants.rootTwo(Float80.self)
    case .rootThree: return Constants.rootThree(Float80.self)
    case .rootPi: return Constants.rootPi(Float80.self)
    case .oneDivPi: return Constants.oneDivPi(Float80.self)
    case .oneDivTwoPi: return Constants.oneDivTwoPi(Float80.self)
    case .oneDivRootPi: return Constants.oneDivRootPi(Float80.self)
    case .twoDivPi: return Constants.twoDivPi(Float80.self)
    case .twoDivRootPi: return Constants.twoDivRootPi(Float80.self)
    case .lnTwo: return Constants.lnTwo(Float80.self)
    case .lnTen: return Constants.lnTen(Float80.self)
    case .lnLnTwo: return Constants.lnLnTwo(Float80.self)
    case .e: return Constants.e(Float80.self)
    case .oneDivE: return Constants.oneDivE(Float80.self)
    case .euler: return Constants.euler(Float80.self)
    case .catalan: return Constants.catalan(Float80.self)
    case .zetaThree: return Constants.zetaThree(Float80.self)
    case .phi: return Constants.phi(Float80.self)
    case .holtsmarkEntropy: return Constants.holtsmarkEntropy(Float80.self)
    }
}
#endif

@Suite("Constants<Double>")
struct HighPrecisionDoubleConstantsTests {
    let tol: Double = 1e-15

    @Test
    func valueMatchesMatlabCatalogue() {
        for identifier in Constants.Identifier.allCases {
            let expected = matlabValue(identifier, as: Double.self)
            let got = Constants.value(identifier, as: Double.self)
            #expect(approxEqual(got, expected, tol: tol), "Mismatch for \(identifier)")
        }
    }

    @Test
    func convenienceMatchesValue() {
        for identifier in Constants.Identifier.allCases {
            let expected = Constants.value(identifier, as: Double.self)
            let got = convenienceValueDouble(identifier)
            #expect(approxEqual(got, expected, tol: tol), "Mismatch for \(identifier)")
        }
    }
}

@Suite("Constants<Float>")
struct HighPrecisionFloatConstantsTests {
    let tol: Float = 1e-6

    @Test
    func valueMatchesMatlabCatalogue() {
        for identifier in Constants.Identifier.allCases {
            let expected: Float = matlabValue(identifier, as: Float.self)
            let got = Constants.value(identifier, as: Float.self)
            #expect(approxEqual(got, expected, tol: tol), "Mismatch for \(identifier)")
        }
    }

    @Test
    func convenienceMatchesValue() {
        for identifier in Constants.Identifier.allCases {
            let expected = Constants.value(identifier, as: Float.self)
            let got = convenienceValueFloat(identifier)
            #expect(approxEqual(got, expected, tol: tol), "Mismatch for \(identifier)")
        }
    }
}

#if arch(x86_64)
@Suite("Constants<Float80>")
struct HighPrecisionFloat80ConstantsTests {
    let tol: Float80 = 1e-18

    @Test
    func valueMatchesMatlabCatalogue() {
        for identifier in Constants.Identifier.allCases {
            let expected: Float80 = matlabValue(identifier, as: Float80.self)
            let got = Constants.value(identifier, as: Float80.self)
            #expect(approxEqual(got, expected, tol: tol), "Mismatch for \(identifier)")
        }
    }

    @Test
    func convenienceMatchesValue() {
        for identifier in Constants.Identifier.allCases {
            let expected = Constants.value(identifier, as: Float80.self)
            let got = convenienceValueFloat80(identifier)
            #expect(approxEqual(got, expected, tol: tol), "Mismatch for \(identifier)")
        }
    }
}
#endif
