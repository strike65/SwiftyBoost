import Testing
@testable import SwiftyBoost

@Suite("Mixed-Precision Promotions (Chebyshev, Digamma/Polygamma, Error, E-Integrals)")
struct MixedPrecisionPromotionMorePermutationsTests {

    @inlinable
    func near(_ a: Double, _ b: Double, _ tol: Double = 1e-12) { #expect(Swift.abs(a - b) <= tol) }

    // MARK: - Chebyshev T/U basic cross-type consistency

    @Test("chebyshevT: Float vs Double consistency")
    func chebyshevT_Float_Double() throws {
        let n = 7
        let xF: Float = 0.3
        let tF = try SpecialFunctions.chebyshevT(n, xF)
        let tD = try SpecialFunctions.chebyshevT(n, Double(xF))
        near(Double(tF), tD, 1e-6)
    }
    
    @Test("chebyshevU: Float vs Double consistency")
    func chebyshevU_Float_Double() throws {
        let n = 6
        let xF: Float = 0.25
        let uF = try SpecialFunctions.chebyshevU(n, xF)
        let uD = try SpecialFunctions.chebyshevU(n, Double(xF))
        near(Double(uF), uD, 1e-6)
    }

    // MARK: - Chebyshev Clenshaw mixed promotions

    @Test("chebyshevClenshaw: [Float] + Double x → Double")
    func chebyshevClenshaw_FloatCoeffs_DoubleX() throws {
        let coeffsF: [Float] = [1.0, 0.5, -0.25, 0.125]
        let xD: Double = 0.4
        let got = try SpecialFunctions.chebyshevClenshawRecurrence(coeffsF, x: xD)
        let exp = try SpecialFunctions.chebyshevClenshawRecurrence(coeffsF.map(Double.init), x: xD)
        near(got, exp)
    }

    @Test("chebyshevClenshaw: [Double] + Float x → Double")
    func chebyshevClenshaw_DoubleCoeffs_FloatX() throws {
        let coeffsD: [Double] = [1.0, -0.3, 0.2]
        let xF: Float = 0.6
        let got = try SpecialFunctions.chebyshevClenshawRecurrence(coeffsD, x: xF)
        let exp = try SpecialFunctions.chebyshevClenshawRecurrence(coeffsD, x: Double(xF))
        near(got, exp)
    }

    // MARK: - Digamma / Trigamma / Polygamma / Zeta

    @Test("digamma: Float vs Double")
    func digamma_F_D() throws {
        let xF: Float = 3.5
        let dF = try SpecialFunctions.digamma(xF)
        let dD = try SpecialFunctions.digamma(Double(xF))
        near(Double(dF), dD, 1e-6)
    }

    @Test("trigamma: Float vs Double")
    func trigamma_F_D() throws {
        let xF: Float = 2.25
        let tF = try SpecialFunctions.trigamma(xF)
        let tD = try SpecialFunctions.trigamma(Double(xF))
        near(Double(tF), tD, 1e-6)
    }

    @Test("polygamma (n=2): Float vs Double")
    func polygamma2_F_D() throws {
        let xF: Float = 4.0
        let pF = try SpecialFunctions.polygamma(order: 2, xF)
        let pD = try SpecialFunctions.polygamma(order: 2, Double(xF))
        near(Double(pF), pD, 1e-5)
    }

    @Test("riemannZeta: Float vs Double")
    func zeta_F_D() throws {
        let xF: Float = 2.0
        let zF = try SpecialFunctions.riemannZeta(xF)
        let zD = try SpecialFunctions.riemannZeta(Double(xF))
        near(Double(zF), zD, 1e-6)
    }

    // MARK: - Error functions

    @Test("erf/erfc: Float vs Double")
    func errorFunctions_F_D() throws {
        let xF: Float = 0.7
        let erfF = try SpecialFunctions.errorFunction(xF)
        let erfD = try SpecialFunctions.errorFunction(Double(xF))
        near(Double(erfF), erfD, 1e-6)

        let erfcF = try SpecialFunctions.complementaryErrorFunction(xF)
        let erfcD = try SpecialFunctions.complementaryErrorFunction(Double(xF))
        near(Double(erfcF), erfcD, 1e-6)
    }

    @Test("inverse erf: Float vs Double")
    func inverseErf_F_D() throws {
        let zF: Float = 0.5
        let invF = try SpecialFunctions.inverseErrorFunction(zF)
        let invD = try SpecialFunctions.inverseErrorFunction(Double(zF))
        near(Double(invF), invD, 1e-6)
    }

    // MARK: - Exponential integrals

    @Test("Ei: Float vs Double")
    func ei_F_D() throws {
        let xF: Float = 1.0
        let eiF = try SpecialFunctions.exponentialIntegralEi(xF)
        let eiD = try SpecialFunctions.exponentialIntegralEi(Double(xF))
        near(Double(eiF), eiD, 1e-6)
    }

    @Test("En: Float vs Double")
    func en_F_D() throws {
        let n = 2
        let xF: Float = 1.2
        let enF = try SpecialFunctions.exponentialIntegralEn(n, xF)
        let enD = try SpecialFunctions.exponentialIntegralEn(n, Double(xF))
        near(Double(enF), enD, 1e-6)
    }

    // MARK: - Gegenbauer mix and Laguerre/Legendre consistency

    @Test("gegenbauer: Float lambda + Double x → Double")
    func gegenbauer_Flambda_Dx() throws {
        let n = 4
        let lamF: Float = 0.75
        let xD: Double = 0.3
        let got = try SpecialFunctions.gegenbauer(n: n, lambda: lamF, x: xD)
        let exp = try SpecialFunctions.gegenbauer(n: n, lambda: Double(lamF), x: xD)
        near(got, exp)
    }

    @Test("gegenbauerPrime: Double lambda + Float x → Double")
    func gegenbauerPrime_Dlambda_Fx() throws {
        let n = 5
        let lamD: Double = 1.25
        let xF: Float = 0.2
        let got = try SpecialFunctions.gegenbauerPrime(n: n, lambda: lamD, x: xF)
        let exp = try SpecialFunctions.gegenbauerPrime(n: n, lambda: lamD, x: Double(xF))
        near(got, exp)
    }

    @Test("gegenbauerDerivative: Float lambda + Double x → Double")
    func gegenbauerDeriv_Flambda_Dx() throws {
        let n = 3, k = 2
        let lamF: Float = 0.6
        let xD: Double = 0.4
        let got = try SpecialFunctions.gegenbauerDerivative(n: n, lambda: lamF, x: xD, k: k)
        let exp = try SpecialFunctions.gegenbauerDerivative(n: n, lambda: Double(lamF), x: xD, k: k)
        near(got, exp, 1e-11)
    }

    @Test("laguerre: Float vs Double consistency")
    func laguerre_F_D() throws {
        let n = 6
        let xF: Float = 0.8
        let lF = try SpecialFunctions.laguerre(n, xF)
        let lD = try SpecialFunctions.laguerre(n, Double(xF))
        near(Double(lF), lD, 1e-6)
    }

    @Test("assocLaguerre: Float vs Double consistency")
    func assocLaguerre_F_D() throws {
        let n = 4, m = 2
        let xF: Float = 0.9
        let lF = try SpecialFunctions.assocLaguerre(n, m, xF)
        let lD = try SpecialFunctions.assocLaguerre(n, m, Double(xF))
        near(Double(lF), lD, 1e-6)
    }

    @Test("legendre P and P': Float vs Double consistency")
    func legendreP_consistency() throws {
        let n = 5
        let xF: Float = 0.35
        near(Double(try SpecialFunctions.legendreP(n, xF)), try SpecialFunctions.legendreP(n, Double(xF)), 1e-6)
        near(Double(try SpecialFunctions.legendrePPrime(n, xF)), try SpecialFunctions.legendrePPrime(n, Double(xF)), 1e-6)
        // Associated
        let m = 2
        near(Double(try SpecialFunctions.associatedLegendreP(n, m, xF)), try SpecialFunctions.associatedLegendreP(n, m, Double(xF)), 1e-6)
    }

    // MARK: - Owen's T and Spherical Harmonics

    @Test("Owen's T: Float + Double → Double")
    func owensT_F_D() throws {
        let hF: Float = 0.7
        let aD: Double = 1.25
        let got = try SpecialFunctions.owensT(h: hF, a: aD)
        let exp = try SpecialFunctions.owensT(h: Double(hF), a: aD)
        near(got, exp, 1e-12)
    }

    @Test("SphericalHarmonics: mixed Float/Double angles")
    func sphHarmonics_mixed() {
        let n: UInt = 3
        let m = 1
        let thetaF: Float = 0.6
        let phiD: Double = 1.2
        let got = SpecialFunctions.sphericalHarmonic(n: n, m: m, theta: thetaF, phi: phiD)
        let exp = SpecialFunctions.sphericalHarmonic(n: n, m: m, theta: Double(thetaF), phi: phiD)
        near(Double(got.real), Double(exp.real), 1e-12)
        near(Double(got.imag), Double(exp.imag), 1e-12)
    }

    // MARK: - Hypergeometric promotions

    @Test("1F0: Float a + Double z → Double")
    func hyper1F0_F_D() throws {
        let aF: Float = 0.5
        let zD: Double = 0.3
        let got = try SpecialFunctions.hypergeometric1F0(a: aF, z: zD)
        let exp = try SpecialFunctions.hypergeometric1F0(a: Double(aF), z: zD)
        near(got, exp, 1e-12)
    }

    @Test("0F1: Double b + Float z → Double")
    func hyper0F1_D_F() throws {
        let bD: Double = 1.5
        let zF: Float = 0.2
        let got = try SpecialFunctions.hypergeometric0F1(b: bD, z: zF)
        let exp = try SpecialFunctions.hypergeometric0F1(b: bD, z: Double(zF))
        near(got, exp)
    }

    @Test("1F1: single Float param mix → Double")
    func hyper1F1_mix() throws {
        let aF: Float = 0.4, bD: Double = 1.3, zD: Double = 0.25
        let got = try SpecialFunctions.hypergeometric1F1(a: aF, b: bD, z: zD)
        let exp = try SpecialFunctions.hypergeometric1F1(a: Double(aF), b: bD, z: zD)
        near(got, exp, 1e-12)
    }

    @Test("2F0: single Float param mix → Double (z=0)")
    func hyper2F0_mix() throws {
        let aD: Double = 0.3, bF: Float = 0.2, zD: Double = 0.0
        let got = try SpecialFunctions.hypergeometric2F0(a: aD, b: bF, z: zD)
        let exp = try SpecialFunctions.hypergeometric2F0(a: aD, b: Double(bF), z: zD)
        near(got, exp)
    }

    @Test("pFq: [Float] + z Double → Double")
    func hyperPFQ_mix() {
        let aF: [Float] = [0.5]
        let bF: [Float] = [1.5]
        let zD: Double = 0.2
        let got = SpecialFunctions.hypergeometricPFQ(a: aF, b: bF, z: zD)
        let exp = SpecialFunctions.hypergeometricPFQ(a: aF.map(Double.init), b: bF.map(Double.init), z: zD)
        near(got, exp)
    }
}

#if arch(x86_64)
@Suite("Float80 Mix (Chebyshev, Digamma/Polygamma, Error, E-Integrals)")
struct MixedPrecisionPromotionMorePermutationsFloat80Tests {

    @inlinable
    func near80(_ a: Float80, _ b: Float80, _ tol: Float80 = 1e-16) { #expect(Swift.abs(a - b) <= tol) }

    @Test("chebyshevT: Float80 vs Double")
    func chebyshevT_L_D() throws {
        let n = 5
        let xL: Float80 = 0.45
        let tL = try SpecialFunctions.chebyshevT(n, xL)
        let tD = try SpecialFunctions.chebyshevT(n, Double(xL))
        near80(tL, Float80(tD), 1e-15)
    }

    @Test("chebyshevClenshaw: [Float80] + Double x → Float80")
    func chebyshevClenshaw_L_Dx() throws {
        let coeffsL: [Float80] = [1.0, 0.2, -0.1]
        let xD: Double = 0.3
        let got = try SpecialFunctions.chebyshevClenshawRecurrence(coeffsL, x: xD)
        let exp = try SpecialFunctions.chebyshevClenshawRecurrence(coeffsL, x: Float80(xD))
        near80(got, exp)
    }

    @Test("digamma/trigamma/polygamma: Float80 vs Double")
    func psi_L_D() throws {
        let xL: Float80 = 3.25
        near80(try SpecialFunctions.digamma(xL), Float80(try SpecialFunctions.digamma(Double(xL))), 1e-16)
        near80(try SpecialFunctions.trigamma(xL), Float80(try SpecialFunctions.trigamma(Double(xL))), 1e-16)
        near80(try SpecialFunctions.polygamma(order: 2, xL), Float80(try SpecialFunctions.polygamma(order: 2, Double(xL))), 1e-16)
    }

    @Test("riemannZeta: Float80 vs Double")
    func zeta_L_D() throws {
        let xL: Float80 = 2.0
        near80(try SpecialFunctions.riemannZeta(xL), Float80(try SpecialFunctions.riemannZeta(Double(xL))), 1e-16)
    }

    @Test("erf/erfc/inv: Float80 vs Double")
    func erf_L_D() throws {
        let xL: Float80 = 0.65
        near80(try SpecialFunctions.errorFunction(xL), Float80(try SpecialFunctions.errorFunction(Double(xL))), 1e-16)
        near80(try SpecialFunctions.complementaryErrorFunction(xL), Float80(try SpecialFunctions.complementaryErrorFunction(Double(xL))), 1e-16)
        let zL: Float80 = 0.3
        near80(try SpecialFunctions.inverseErrorFunction(zL), Float80(try SpecialFunctions.inverseErrorFunction(Double(zL))), 1e-16)
    }

    @Test("Ei/En: Float80 vs Double")
    func ei_en_L_D() throws {
        let xL: Float80 = 1.1
        near80(try SpecialFunctions.exponentialIntegralEi(xL), Float80(try SpecialFunctions.exponentialIntegralEi(Double(xL))), 1e-16)
        let n = 2
        near80(try SpecialFunctions.exponentialIntegralEn(n, xL), Float80(try SpecialFunctions.exponentialIntegralEn(n, Double(xL))), 1e-16)
    }

    // MARK: - Gegenbauer/Legendre/Laguerre Float80 mixes

    @Test("gegenbauer: Float80 lambda + Double x → Float80")
    func gegenbauer_L_D() throws {
        let n = 3
        let lamL: Float80 = 0.7
        let xD: Double = 0.2
        let got = try SpecialFunctions.gegenbauer(n: n, lambda: lamL, x: xD)
        let exp = try SpecialFunctions.gegenbauer(n: n, lambda: lamL, x: Float80(xD))
        near80(got, exp)
    }

    @Test("gegenbauerPrime: Double lambda + Float80 x → Float80")
    func gegenbauerPrime_D_L() throws {
        let n = 4
        let lamD: Double = 1.1
        let xL: Float80 = 0.25
        let got = try SpecialFunctions.gegenbauerPrime(n: n, lambda: lamD, x: xL)
        let exp = try SpecialFunctions.gegenbauerPrime(n: n, lambda: Float80(lamD), x: xL)
        near80(got, exp)
    }

    @Test("gegenbauerDerivative: Float + Float80 → Float80")
    func gegenbauerDeriv_F_L() throws {
        let n = 2, k = 1
        let lamF: Float = 0.6
        let xL: Float80 = 0.4
        let got = try SpecialFunctions.gegenbauerDerivative(n: n, lambda: lamF, x: xL, k: k)
        let exp = try SpecialFunctions.gegenbauerDerivative(n: n, lambda: Float80(lamF), x: xL, k: k)
        near80(got, exp)
    }

    @Test("legendre/laguerre: Float80 vs Double consistency")
    func legendreLaguerre_L_D() throws {
        let n = 5, m = 2
        let xL: Float80 = 0.45
        near80(try SpecialFunctions.legendreP(n, xL), Float80(try SpecialFunctions.legendreP(n, Double(xL))), 1e-16)
        near80(try SpecialFunctions.associatedLegendreP(n, m, xL), Float80(try SpecialFunctions.associatedLegendreP(n, m, Double(xL))), 1e-16)
        near80(try SpecialFunctions.laguerre(n, xL), Float80(try SpecialFunctions.laguerre(n, Double(xL))), 1e-16)
        near80(try SpecialFunctions.assocLaguerre(n, m, xL), Float80(try SpecialFunctions.assocLaguerre(n, m, Double(xL))), 1e-16)
    }

    @Test("Owen's T: Float80 mixes → Float80")
    func owensT_L_mix() throws {
        let hL: Float80 = 0.7
        let aD: Double = 1.1
        let got = try SpecialFunctions.owensT(h: hL, a: aD)
        let exp = try SpecialFunctions.owensT(h: hL, a: Float80(aD))
        near80(got, exp)
    }

    @Test("SphericalHarmonics: Float80 + Double angles")
    func sphHarmonics_L_D() {
        let n: UInt = 3, m = 1
        let thetaL: Float80 = 0.7
        let phiD: Double = 1.3
        let got = SpecialFunctions.Y(n: n, m: m, theta: thetaL, phi: phiD)
        let exp = SpecialFunctions.Y(n: n, m: m, theta: thetaL, phi: Float80(phiD))
        near80(got.real, exp.real)
        near80(got.imag, exp.imag)
    }

    @Test("1F1: Float80 mixes → Float80")
    func hyper1F1_L_mix() throws {
        let aL: Float80 = 0.4, bD: Double = 1.25, zD: Double = 0.2
        let got = try SpecialFunctions.hypergeometric1F1(a: aL, b: bD, z: zD)
        let exp = try SpecialFunctions.hypergeometric1F1(a: aL, b: Float80(bD), z: Float80(zD))
        near80(got, exp)
    }
}
#endif
