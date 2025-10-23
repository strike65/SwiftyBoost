import Testing
@testable import SwiftyBoost
import Foundation

@Suite("DistributionProtocol – Geometric")
struct GeometricDistributionTests {

    @Test("PDF and CDF match closed form")
    func pdfAndCdfClosedForm() throws {
        let p: Double = 0.23
        let g = try Distribution.Geometric<Double>(probabibilityOfSuccess: p)
        let ks = [0, 1, 2, 5, 10]
        for k in ks {
            let x = Double(k)
            let expectedPdf = pow(1 - p, x) * p
            let expectedCdf = 1 - pow(1 - p, x + 1)
            #expect(abs((try g.pdf(x)) - expectedPdf) <= 1e-12)
            #expect(abs((try g.cdf(x)) - expectedCdf) <= 1e-12)
            #expect(abs((try g.sf(x)) - (1 - expectedCdf)) <= 1e-12)
        }
    }

    @Test("Quantile round-trip and integrality")
    func quantileRoundTrip() throws {
        let p: Double = 0.41
        let g = try Distribution.Geometric<Double>(probabibilityOfSuccess: p)
        let probs: [Double] = [1e-9, 1e-6, 0.1, 0.5, 0.9, 1 - 1e-9]
        for prob in probs {
            let q = try g.quantile(prob)
            let c = try g.cdf(q)
            #expect(c >= prob - 1e-12)

            let floorQ = Foundation.floor(q)
            if floorQ >= 1 {
                let prev = try g.cdf(floorQ - 1)
                #expect(prev <= prob + 1e-12)
            }
        }
    }

    @Test("Moments align with analytical expressions")
    func moments() throws {
        let p: Double = 0.65
        let g = try Distribution.Geometric<Double>(probabibilityOfSuccess: p)
        let expectedMean = (1 - p) / p
        let expectedVariance = (1 - p) / (p * p)
        let expectedSkewness = (2 - p) / (1 - p).squareRoot()
        let expectedExcess = 6 + (p * p) / (1 - p)
        let expectedEntropy = -((1 - p) / p) * Foundation.log(1 - p) - Foundation.log(p)

        #expect(abs((g.mean ?? .nan) - expectedMean) <= 1e-12)
        #expect(abs((g.variance ?? .nan) - expectedVariance) <= 1e-12)
        #expect(abs((g.mode ?? .nan) - 0) <= 1e-12)
        if let skew = g.skewness {
            #expect(abs(skew - expectedSkewness) <= 1e-10)
        }
        if let kurtosisExcess = g.kurtosisExcess {
            #expect(abs(kurtosisExcess - expectedExcess) <= 1e-10)
        }
        if let entropy = g.entropy {
            #expect(abs(entropy - expectedEntropy) <= 1e-12)
        } else {
            Issue.record("Expected entropy to be available for geometric distribution")
        }
    }
}
