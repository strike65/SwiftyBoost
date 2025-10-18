import Testing
@testable import SwiftyBoost

@Suite("Mixed-Precision Promotions (Gamma, Bessel)")
struct MixedPrecisionPromotionTests {

    // MARK: - Helpers
    @inlinable
    func near(_ a: Double, _ b: Double, _ tol: Double = 1e-12, _ msg: String = "") {
        #expect(Swift.abs(a - b) <= tol)
    }

#if arch(x86_64)
    @inlinable
    func near80(_ a: Float80, _ b: Float80, _ tol: Float80 = 1e-16, _ msg: String = "") {
        #expect(Swift.abs(a - b) <= tol)
    }
#endif

    // MARK: - Gamma mixed Float ↔ Double

    @Test("gamma: Float promoted to Double in binary APIs")
    func gammaFloatToDouble() throws {
        // gammaRatio
        do {
            let aF: Float = 3.5, bD: Double = 2.25
            let got = try SpecialFunctions.gammaRatio(aF, bD)
            let exp = try SpecialFunctions.gammaRatio(Double(aF), bD)
            near(got, exp, 1e-12, "gammaRatio(Float, Double)")
        }
        // gammaDeltaRatio
        do {
            let aF: Float = 4.2, dD: Double = 0.75
            let got = try SpecialFunctions.gammaDeltaRatio(aF, delta: dD)
            let exp = try SpecialFunctions.gammaDeltaRatio(Double(aF), delta: dD)
            near(got, exp, 1e-12, "gammaDeltaRatio(Float, Double)")
        }
        // Incomplete gamma (lower/upper)
        do {
            let aF: Float = 2.5, xD: Double = 1.3
            let gotL = try SpecialFunctions.incompleteGammaLower(aF, x: xD)
            let expL = try SpecialFunctions.incompleteGammaLower(Double(aF), x: xD)
            near(gotL, expL, 1e-12, "incompleteGammaLower(Float, Double)")

            let gotU = try SpecialFunctions.incompleteGammaUpper(aF, x: xD)
            let expU = try SpecialFunctions.incompleteGammaUpper(Double(aF), x: xD)
            near(gotU, expU, 1e-12, "incompleteGammaUpper(Float, Double)")
        }
        // Regularized P/Q and inverses + derivative
        do {
            let aF: Float = 2.5, xD: Double = 1.3
            let pD: Double = 0.35, qD: Double = 0.65
            let gotP = try SpecialFunctions.regularizedGammaP(aF, x: xD)
            let expP = try SpecialFunctions.regularizedGammaP(Double(aF), x: xD)
            near(gotP, expP, 1e-12, "regularizedGammaP(Float, Double)")

            let gotQ = try SpecialFunctions.regularizedGammaQ(aF, x: xD)
            let expQ = try SpecialFunctions.regularizedGammaQ(Double(aF), x: xD)
            near(gotQ, expQ, 1e-12, "regularizedGammaQ(Float, Double)")

            let gotPinv = try SpecialFunctions.regularizedGammaPInv(aF, p: pD)
            let expPinv = try SpecialFunctions.regularizedGammaPInv(Double(aF), p: pD)
            near(gotPinv, expPinv, 1e-12, "regularizedGammaPInv(Float, Double)")

            let gotQinv = try SpecialFunctions.regularizedGammaQInv(aF, q: qD)
            let expQinv = try SpecialFunctions.regularizedGammaQInv(Double(aF), q: qD)
            near(gotQinv, expQinv, 1e-12, "regularizedGammaQInv(Float, Double)")

            let gotPp = try SpecialFunctions.regularizedGammaPDerivative(aF, x: xD)
            let expPp = try SpecialFunctions.regularizedGammaPDerivative(Double(aF), x: xD)
            near(gotPp, expPp, 1e-12, "regularizedGammaPDerivative(Float, Double)")
        }
    }

    // MARK: - Bessel mixed Float ↔ Double

    @Test("bessel: Float + Double → Double (functions and derivatives)")
    func besselFloatToDouble() throws {
        let vF: Float = 2.5, xD: Double = 1.0

        // J
        let jFD = try SpecialFunctions.besselJ(v: vF, x: xD)
        let jDD = try SpecialFunctions.besselJ(v: Double(vF), x: xD)
        near(jFD, jDD, 1e-12, "J_v mixed Float/Double")

        // Y (x > 0)
        let yFD = try SpecialFunctions.besselY(v: vF, x: xD)
        let yDD = try SpecialFunctions.besselY(v: Double(vF), x: xD)
        near(yFD, yDD, 1e-12, "Y_v mixed Float/Double")

        // I
        let iFD = try SpecialFunctions.modifiedBesselI(v: vF, x: xD)
        let iDD = try SpecialFunctions.modifiedBesselI(v: Double(vF), x: xD)
        near(iFD, iDD, 1e-12, "I_v mixed Float/Double")

        // K (x > 0)
        let kFD = try SpecialFunctions.modifiedBesselK(v: vF, x: xD)
        let kDD = try SpecialFunctions.modifiedBesselK(v: Double(vF), x: xD)
        near(kFD, kDD, 1e-12, "K_v mixed Float/Double")

        // Derivatives
        let jpFD = try SpecialFunctions.besselJPrime(v: vF, x: xD)
        let jpDD = try SpecialFunctions.besselJPrime(v: Double(vF), x: xD)
        near(jpFD, jpDD, 1e-12, "J'_v mixed Float/Double")

        let ipFD = try SpecialFunctions.modifiedBesselIPrime(v: vF, x: xD)
        let ipDD = try SpecialFunctions.modifiedBesselIPrime(v: Double(vF), x: xD)
        near(ipFD, ipDD, 1e-12, "I'_v mixed Float/Double")

        let kpFD = try SpecialFunctions.modifiedBesselKPrime(v: vF, x: xD)
        let kpDD = try SpecialFunctions.modifiedBesselKPrime(v: Double(vF), x: xD)
        near(kpFD, kpDD, 1e-12, "K'_v mixed Float/Double")
    }

    // MARK: - Integer order convenience (result type follows x)

    @Test("bessel: integer order overload result type")
    func besselIntegerOrder() throws {
        // Double
        do {
            let v = 3
            let x: Double = 1.1
            let got = try SpecialFunctions.besselJ(v: v, x: x)
            let exp = try SpecialFunctions.besselJ(v: Double(v), x: x)
            near(got, exp, 1e-12, "J_v integer order (Double)")
        }
        // Float (compare to Float-typed call via conversion)
        do {
            let v = 2
            let x: Float = 1.0
            let got = try SpecialFunctions.besselY(v: v, x: x)
            let exp = try SpecialFunctions.besselY_f(v: Float(v), x: x)
            #expect(Swift.abs(got - exp) <= 1e-5, "Y_v integer order (Float)")
        }
        // Derivative, Double
        do {
            let v = 4
            let x: Double = 1.2
            let got = try SpecialFunctions.besselJPrime(v: v, x: x)
            let exp = try SpecialFunctions.besselJPrime(v: Double(v), x: x)
            near(got, exp, 1e-12, "J'_v integer order (Double)")
        }
    }

    // MARK: - Float80 mixes (x86_64)

    #if arch(x86_64)
    @Test("gamma: Float/Double with Float80 → Float80")
    func gammaToFloat80() throws {
        let aL: Float80 = 3.5
        let bD: Double = 2.25
        let got = try SpecialFunctions.gammaRatio(aL, bD)
        let exp = try SpecialFunctions.gammaRatio(aL, Float80(bD))
        near80(got, exp, 1e-16, "gammaRatio(Float80, Double)")

        let dF: Float = 0.75
        let gotDelta = try SpecialFunctions.gammaDeltaRatio(aL, delta: dF)
        let expDelta = try SpecialFunctions.gammaDeltaRatio(aL, delta: Float80(dF))
        near80(gotDelta, expDelta, 1e-16, "gammaDeltaRatio(Float80, Float)")
    }

    @Test("bessel: Float/Double with Float80 → Float80")
    func besselToFloat80() throws {
        let vL: Float80 = 2.5
        let xD: Double = 1.0
        let jLD = try SpecialFunctions.besselJ(v: vL, x: xD)
        let jLL = try SpecialFunctions.besselJ_l(v: vL, x: Float80(xD))
        near80(jLD, jLL, 1e-16, "J_v(Float80, Double)")

        let kLF = try SpecialFunctions.modifiedBesselK(v: vL, x: Float(1.0))
        let kLL = try SpecialFunctions.modifiedBesselK_l(v: vL, x: Float80(1.0))
        near80(kLF, kLL, 1e-16, "K_v(Float80, Float)")
    }
    #endif
}
