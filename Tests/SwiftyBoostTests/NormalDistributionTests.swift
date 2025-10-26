import Testing
@testable import SwiftyBoost
import Foundation

@Suite("Normal distribution tests")
struct NormalDistributionTests {

    @Test("PDF / CDF / quantile consistency")
    func pdfCdfQuantile() throws {
        let mean: Double = 0.75
        let sd: Double = 1.4
        let normal = try Distribution.Normal<Double>(mean: mean, sd: sd)

        let xs: [Double] = [-2.0, -0.25, 0.0, 0.5, 1.5, 3.0]
        for x in xs {
            let pdf = try normal.pdf(x)
            let expected = referencePdf(x, mean: mean, sd: sd)
            #expect(abs(pdf - expected) <= 1e-12)

            let cdf = try normal.cdf(x)
            let expectedCdf = 0.5 * (1 + erf((x - mean) / (sd * Foundation.sqrt(2.0))))
            #expect(abs(cdf - expectedCdf) <= 1e-12)
        }

        let ps: [Double] = [1e-9, 1e-4, 0.1, 0.25, 0.75, 0.9, 1 - 1e-9]
        for p in ps {
            let x = try normal.quantile(p)
            let c = try normal.cdf(x)
            #expect(abs(c - p) <= 1e-10, "round trip failed for p=\(p)")
        }
    }

    @Test("Moments reflect parameters")
    func momentsMatchParameters() throws {
        let mean: Double = -0.35
        let sd: Double = 2.1
        let n = try Distribution.Normal<Double>(mean: mean, sd: sd)
        #expect(n.mean == mean)
        #expect(abs((n.variance ?? .nan) - sd * sd) <= 1e-12)
        #expect(n.median == mean)
        #expect(n.mode == mean)
        #expect(n.skewness == 0)
        #expect(n.kurtosisExcess == 0)
    }

    @Test("KL divergence closed form (typed)")
    func klDivergenceMatchesClosedForm() throws {
        let p = try Distribution.Normal<Double>(mean: 0.2, sd: 0.9)
        let q = try Distribution.Normal<Double>(mean: -0.4, sd: 1.7)
        let expected = analyticKL(meanP: 0.2, sdP: 0.9, meanQ: -0.4, sdQ: 1.7)
        let maybeKL = try p.klDivergence(relativeTo: q)
        let kl = try #require(maybeKL)
        #expect(abs(kl - expected) <= 1e-13)
    }

    @Test("Dynamic normal matches typed KL divergence")
    func dynamicNormalKLMatchesTyped() throws {
        let meanP: Double = 0.8
        let sdP: Double = 1.1
        let meanQ: Double = -1.25
        let sdQ: Double = 0.75
        let typedP = try Distribution.Normal<Double>(mean: meanP, sd: sdP)
        let typedQ = try Distribution.Normal<Double>(mean: meanQ, sd: sdQ)
        let expectedTyped = try typedP.klDivergence(relativeTo: typedQ)
        let expected = try #require(expectedTyped)

        let dynamicP = try Distribution.Dynamic<Double>(
            distributionName: "normal",
            parameters: ["mean": meanP, "sd": sdP]
        )
        let dynamicQ = try Distribution.Dynamic<Double>(
            distributionName: "normal",
            parameters: ["mean": meanQ, "sd": sdQ]
        )
        let maybeDynamicKL = try dynamicP.klDivergence(relativeTo: dynamicQ, options: .automatic())
        let dynamicKL = try #require(maybeDynamicKL)
        #expect(abs(dynamicKL - expected) <= 5e-8)
    }

    private func referencePdf(_ x: Double, mean: Double, sd: Double) -> Double {
        let z = (x - mean) / sd
        let coeff = 1 / (sd * Foundation.sqrt(2 * Double.pi))
        return coeff * Foundation.exp(-0.5 * z * z)
    }

    private func analyticKL(meanP: Double, sdP: Double, meanQ: Double, sdQ: Double) -> Double {
        let logTerm = Foundation.log(sdQ / sdP)
        let meanDiff = meanP - meanQ
        let varianceTerm = (sdP * sdP + meanDiff * meanDiff) / (2 * sdQ * sdQ)
        return logTerm + varianceTerm - 0.5
    }
}
