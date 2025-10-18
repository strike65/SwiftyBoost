import Testing
@testable import SwiftyBoost

@Suite("Mixed-Precision Promotions (Beta, Elliptic, Splines)")
struct MixedPrecisionPromotionBetaEllipticSplinesTests {

    // Common helper
    @inlinable
    func near(_ a: Double, _ b: Double, _ tol: Double = 1e-12) {
        #expect(Swift.abs(a - b) <= tol)
    }

    // MARK: - Beta

    @Test("beta: Float + Double → Double")
    func betaFloatDouble() throws {
        let aF: Float = 2.25
        let bD: Double = 1.75
        let got = try SpecialFunctions.beta(aF, bD)
        let exp = try SpecialFunctions.beta(Double(aF), bD)
        near(got, exp)
    }

    @Test("regularizedIncompleteBeta: Float + Double → Double")
    func regIncBetaFloatDouble() throws {
        let aF: Float = 2.5
        let bD: Double = 3.0
        let xD: Double = 0.4
        let got = try SpecialFunctions.regularizedIncompleteBeta(aF, bD, x: xD)
        let exp = try SpecialFunctions.regularizedIncompleteBeta(Double(aF), bD, x: xD)
        near(got, exp)
    }

    @Test("inverseRegularizedIncompleteBeta: Float + Double → Double")
    func invRegIncBetaFloatDouble() throws {
        let aF: Float = 1.5
        let bD: Double = 2.0
        let pD: Double = 0.3
        let got = try SpecialFunctions.inverseRegularizedIncompleteBeta(aF, bD, p: pD)
        let exp = try SpecialFunctions.inverseRegularizedIncompleteBeta(Double(aF), bD, p: pD)
        near(got, exp)
    }

    @Test("solveAForRegularizedIncompleteBeta: mixed → Double")
    func solveARegIncBetaMixed() {
        let bF: Float = 2
        let xD: Double = 0.4
        let pD: Double = 0.3
        let got = SpecialFunctions.solveAForRegularizedIncompleteBeta(b: bF, x: xD, p: pD)
        let exp = SpecialFunctions.solveAForRegularizedIncompleteBeta(b: Double(bF), x: xD, p: pD)
        near(got, exp, 1e-10)
    }

    // MARK: - Elliptic (Legendre)

    @Test("elliptic F: Float + Double → Double")
    func ellipticF_FloatDouble() throws {
        let kF: Float = 0.6
        let phiD: Double = 0.9
        let got = try SpecialFunctions.incompleteEllipticIntegralF(kF, phi: phiD)
        let exp = try SpecialFunctions.incompleteEllipticIntegralF(Double(kF), phi: phiD)
        near(got, exp)
    }

    @Test("elliptic E: Double + Float → Double")
    func ellipticE_DoubleFloat() throws {
        let kD: Double = 0.5
        let phiF: Float = 0.7
        let got = try SpecialFunctions.incompleteEllipticIntegralE(kD, phi: phiF)
        let exp = try SpecialFunctions.incompleteEllipticIntegralE(kD, phi: Double(phiF))
        near(got, exp)
    }

    @Test("elliptic Π (complete): Float + Double → Double")
    func ellipticPiComplete_FloatDouble() {
        let kF: Float = 0.4
        let nuD: Double = 0.2
        let got = SpecialFunctions.completeEllipticIntegralPi(kF, nuD)
        let exp = SpecialFunctions.completeEllipticIntegralPi(Double(kF), nuD)
        near(got, exp, 1e-12)
    }

    @Test("elliptic Π (incomplete): mixed → Double")
    func ellipticPiIncomplete_Mixed() {
        let kD: Double = 0.3
        let nuF: Float = 0.1
        let phiD: Double = 0.75
        let got = SpecialFunctions.incompleteEllipticIntegralPi(kD, nuF, phiD)
        let exp = SpecialFunctions.incompleteEllipticIntegralPi(kD, Double(nuF), phiD)
        near(got, exp)
    }

    // MARK: - Elliptic (Carlson symmetric)

    @Test("carlsonRC: Float + Double → Double")
    func carlsonRC_FloatDouble() throws {
        let xF: Float = 0.25
        let yD: Double = 0.5 // must be > 0
        let got = try SpecialFunctions.carlsonRC(xF, yD)
        let exp = try SpecialFunctions.carlsonRC(Double(xF), yD)
        near(got, exp)
    }

    @Test("carlsonRF: Float + Double + Double → Double")
    func carlsonRF_Mixed() throws {
        let xF: Float = 0.2, yD: Double = 0.4, zD: Double = 0.6
        let got = try SpecialFunctions.carlsonRF(xF, yD, zD)
        let exp = try SpecialFunctions.carlsonRF(Double(xF), yD, zD)
        near(got, exp)
    }

    @Test("carlsonRD: Double + Float + Double → Double")
    func carlsonRD_Mixed() throws {
        let xD: Double = 0.2, yF: Float = 0.4, zD: Double = 0.6
        let got = try SpecialFunctions.carlsonRD(xD, yF, zD)
        let exp = try SpecialFunctions.carlsonRD(xD, Double(yF), zD)
        near(got, exp)
    }

    @Test("carlsonRJ: Double + Double + Float + Double → Double")
    func carlsonRJ_Mixed() throws {
        let xD: Double = 0.2, yD: Double = 0.3, zF: Float = 0.4, pD: Double = 0.5
        let got = try SpecialFunctions.carlsonRJ(xD, yD, zF, pD)
        let exp = try SpecialFunctions.carlsonRJ(xD, yD, Double(zF), pD)
        near(got, exp, 1e-11)
    }

    @Test("carlsonRG: Double + Double + Float → Double")
    func carlsonRG_Mixed() throws {
        let xD: Double = 0.2, yD: Double = 0.3, zF: Float = 0.4
        let got = try SpecialFunctions.carlsonRG(xD, yD, zF)
        let exp = try SpecialFunctions.carlsonRG(xD, yD, Double(zF))
        near(got, exp)
    }

    // MARK: - Splines (consistency across types)

    @Test("cardinal B-spline: Float vs Double (consistency)")
    func bSplineConsistency() throws {
        let n = 4
        let xF: Float = 0.25
        let xD: Double = Double(xF)
        let fF = try SpecialFunctions.cardinal_B_Spline(n, xF)
        let fD: Double = try SpecialFunctions.cardinal_B_Spline(n, xD)
        near(Double(fF), fD, 1e-6)
    }
}

#if arch(x86_64)
@Suite("Float80 Mixed Promotions (Beta, Elliptic, Splines)")
struct MixedPrecisionPromotionFloat80Tests {

    @inlinable
    func near80(_ a: Float80, _ b: Float80, _ tol: Float80 = 1e-16) {
        #expect(Swift.abs(a - b) <= tol)
    }

    // MARK: - Beta (Float80 mixes)

    @Test("beta: Float80 + Double → Float80")
    func betaL_Double() throws {
        let aL: Float80 = 2.25, bD: Double = 1.75
        let got = try SpecialFunctions.beta(aL, bD)
        let exp = try SpecialFunctions.beta(aL, Float80(bD))
        near80(got, exp)
    }

    @Test("regularizedIncompleteBeta: Float80 + Double → Float80")
    func regIncBetaL_Double() throws {
        let aL: Float80 = 2.5, bD: Double = 3.0, xD: Double = 0.4
        let got = try SpecialFunctions.regularizedIncompleteBeta(aL, bD, x: xD)
        let exp = try SpecialFunctions.regularizedIncompleteBeta(aL, Float80(bD), x: Float80(xD))
        near80(got, exp, 1e-16)
    }

    @Test("inverseRegularizedIncompleteBeta: Float80 + Float → Float80")
    func invRegIncBetaL_Float() throws {
        let aL: Float80 = 1.5, bF: Float = 2.0, pF: Float = 0.3
        let got = try SpecialFunctions.inverseRegularizedIncompleteBeta(aL, bF, p: pF)
        let exp = try SpecialFunctions.inverseRegularizedIncompleteBeta(aL, Float80(bF), p: Float80(pF))
        near80(got, exp, 1e-16)
    }

    @Test("inverseComplementaryRegularizedIncompleteBeta: Double + Float80 → Float80")
    func invCompRegIncBetaD_L() throws {
        let aD: Double = 1.25, bL: Float80 = 0.75, pD: Double = 0.4
        let got = try SpecialFunctions.inverseComplementaryRegularizedIncompleteBeta(aD, bL, p: pD)
        let exp = try SpecialFunctions.inverseComplementaryRegularizedIncompleteBeta(Float80(aD), bL, p: Float80(pD))
        near80(got, exp, 1e-16)
    }

    @Test("solveAForRegularizedIncompleteBeta: b Float80, x/p Double → Float80")
    func solveARegIncBetaL_D() {
        let bL: Float80 = 2, xD: Double = 0.4, pD: Double = 0.3
        let got = SpecialFunctions.solveAForRegularizedIncompleteBeta(b: bL, x: xD, p: pD)
        let exp = SpecialFunctions.solveAForRegularizedIncompleteBeta(b: bL, x: Float80(xD), p: Float80(pD))
        near80(got, exp, 1e-14)
    }

    // MARK: - Elliptic Legendre (Float80 mixes)

    @Test("elliptic F: Float80 + Double → Float80")
    func ellipticF_L_D() throws {
        let kL: Float80 = 0.6, phiD: Double = 0.9
        let got = try SpecialFunctions.incompleteEllipticIntegralF(kL, phi: phiD)
        let exp = try SpecialFunctions.incompleteEllipticIntegralF(kL, phi: Float80(phiD))
        near80(got, exp, 1e-16)
    }

    @Test("elliptic E: Double + Float80 → Float80")
    func ellipticE_D_L() throws {
        let kD: Double = 0.5, phiL: Float80 = 0.7
        let got = try SpecialFunctions.incompleteEllipticIntegralE(kD, phi: phiL)
        let exp = try SpecialFunctions.incompleteEllipticIntegralE(Float80(kD), phi: phiL)
        near80(got, exp, 1e-16)
    }

    @Test("elliptic Π (incomplete): Float + Float80 + Double → Float80")
    func ellipticPiIncomplete_F_L_D() {
        let kF: Float = 0.3, nuL: Float80 = 0.1, phiD: Double = 0.75
        let got = SpecialFunctions.incompleteEllipticIntegralPi(kF, nuL, phiD)
        let exp = SpecialFunctions.incompleteEllipticIntegralPi(Float80(kF), nuL, Float80(phiD))
        near80(got, exp, 1e-16)
    }

    @Test("elliptic Π (complete): Float80 + Float → Float80")
    func ellipticPiComplete_L_F() {
        let kL: Float80 = 0.4, nuF: Float = 0.2
        let got = SpecialFunctions.completeEllipticIntegralPi(kL, nuF)
        let exp = SpecialFunctions.completeEllipticIntegralPi(kL, Float80(nuF))
        near80(got, exp, 1e-16)
    }

    // MARK: - Elliptic Carlson (Float80 mixes)

    @Test("carlsonRC: Float80 + Double → Float80")
    func carlsonRC_L_D() throws {
        let xL: Float80 = 0.25, yD: Double = 0.5
        let got = try SpecialFunctions.carlsonRC(xL, yD)
        let exp = try SpecialFunctions.carlsonRC(xL, Float80(yD))
        near80(got, exp, 1e-16)
    }

    @Test("carlsonRF: Float80 + Double + Double → Float80")
    func carlsonRF_L_DD() throws {
        let xL: Float80 = 0.2, yD: Double = 0.4, zD: Double = 0.6
        let got = try SpecialFunctions.carlsonRF(xL, yD, zD)
        let exp = try SpecialFunctions.carlsonRF(xL, Float80(yD), Float80(zD))
        near80(got, exp, 1e-16)
    }

    @Test("carlsonRD: Double + Float80 + Double → Float80")
    func carlsonRD_D_L_D() throws {
        let xD: Double = 0.2, yL: Float80 = 0.4, zD: Double = 0.6
        let got = try SpecialFunctions.carlsonRD(xD, yL, zD)
        let exp = try SpecialFunctions.carlsonRD(Float80(xD), yL, Float80(zD))
        near80(got, exp, 1e-16)
    }

    @Test("carlsonRJ: Float + Float + Float80 + Float → Float80")
    func carlsonRJ_F_F_L_F() throws {
        let xF: Float = 0.2, yF: Float = 0.3, zL: Float80 = 0.4, pF: Float = 0.5
        let got = try SpecialFunctions.carlsonRJ(xF, yF, zL, pF)
        let exp = try SpecialFunctions.carlsonRJ(Float80(xF), Float80(yF), zL, Float80(pF))
        near80(got, exp, 1e-16)
    }

    @Test("carlsonRG: Float80 + Float + Float → Float80")
    func carlsonRG_L_FF() throws {
        let xL: Float80 = 0.2, yF: Float = 0.3, zF: Float = 0.4
        let got = try SpecialFunctions.carlsonRG(xL, yF, zF)
        let exp = try SpecialFunctions.carlsonRG(xL, Float80(yF), Float80(zF))
        near80(got, exp, 1e-16)
    }

    // MARK: - Splines (Float80 basic consistency)

    @Test("cardinal B-spline: Float80 vs Double (consistency)")
    func bSplineFloat80Consistency() throws {
        let n = 5
        let xL: Float80 = 0.125
        let fL = try SpecialFunctions.cardinal_B_Spline(n, xL)
        let fD: Double = try SpecialFunctions.cardinal_B_Spline(n, Double(xL))
        near80(fL, Float80(fD), 1e-12)
    }
}
#endif
